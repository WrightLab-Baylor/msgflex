#!/usr/bin/env python3
"""
MSGFLEX: Conventional MS-GF+ Pipeline.
"""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
import psutil
import time
import uuid
import threading
import pandas as pd

from typing import Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

from msgflex.core.process_registry import managed_popen, managed_run
from msgflex.core import ProjectPaths, SparxError, get_logger
from msgflex.modules import masic_merger
from msgflex.modules.dbcurator import remove_redundancy
from msgflex.modules.inputfile_generator import generate_or_update_tsv
from msgflex.modules.MASIC_wrapper import run_masic_batch
from msgflex.modules.tic_plot import generate_tic_plot
from msgflex.tools import tool

log = get_logger("pipelines.conventional")

# ---helper---
def _normalize(p):
    return os.path.abspath(os.path.expanduser(p))

# Runtime hardware detection 
def _get_core_config(total_ram_gb):
    """Detect core count and NUMA topology; derive safe thread/task/job settings."""
    total_cores = os.cpu_count() or 4  # fallback if detection fails

    numa_nodes = 1
    try:
        numa_path = "/sys/devices/system/node"
        if os.path.isdir(numa_path):
            nodes = [d for d in os.listdir(numa_path) if d.startswith("node") and d[4:].isdigit()]
            numa_nodes = max(1, len(nodes))
    except Exception:
        pass

    cores_per_node = max(1, total_cores // numa_nodes)
    return total_cores, numa_nodes, cores_per_node

NUMACTL_AVAILABLE = shutil.which("numactl") is not None

#  Module-level helpers 
def run_msgfplus(mzML, db_name, i, datapath, database_path, conf_loc, java_mem="4G",
                  msgf_threads=None, msgf_tasks=None, numa_node=None):
    """
    Run MS-GF+ for a single mzML file against a given database.
    """
    sample_name = os.path.splitext(mzML)[0]
    mzid_file = sample_name + ".mzid"
    mzid_path = os.path.join(datapath, mzid_file)

    if os.path.exists(mzid_path):
        log.info(f"Skipping {mzML}: {mzid_file} already exists.")
        return True

    temp_fasta = os.path.join(
        database_path,
        f"temp_database_{uuid.uuid4().hex[:6]}.fasta"
    )

    try:
        shutil.copyfile(os.path.join(database_path, db_name), temp_fasta)

        cmd = [
            "java", f"-Xmx{java_mem}",
            "-jar", tool("MSGFPlus.jar"),
            "-s", os.path.join(datapath, mzML),
            "-d", temp_fasta,
            "-conf", str(conf_loc),
            "-o", mzid_path,
            "-addFeatures", "1",
            "-tda", "1",
        ]

        if msgf_threads is not None:
            cmd += ["-thread", str(msgf_threads)]
        if msgf_tasks is not None:
            cmd += ["-tasks", str(msgf_tasks)]

        if numa_node is not None and NUMACTL_AVAILABLE:
            cmd = ["numactl", f"--cpunodebind={numa_node}", f"--membind={numa_node}"] + cmd
        elif numa_node is not None:
            log.debug(f"numactl not found — running {mzML} without NUMA pinning.")

        log.info(f"Running MSGFPlus for {mzML} with {temp_fasta}")

        process = managed_popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

        def _stream(pipe, log_func):
            for line in iter(pipe.readline, ''):
                if line:
                    log_func("[%s] %s", mzML, line.rstrip())
            pipe.close()

        t_out = threading.Thread(target=_stream, args=(process.stdout, log.info))
        t_err = threading.Thread(target=_stream, args=(process.stderr, log.error))

        t_out.start()
        t_err.start()

        process.wait()
        t_out.join()
        t_err.join()

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)

        log.info(f"Success: {mzML}")
        return True

    except Exception:
        log.exception(f"Failed: {mzML}")
        return False

    finally:
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)

def run(paths, conf_loc, java_mem = "4G",
        msgf_tasks = None,
        stop_event = None):
    """
    Run the conventional MS-GF+ pipeline end-to-end.
    """
    conf_loc = str(conf_loc)
    log.info(f"Conventional pipeline -- base: {paths.base}")
    log.info(f"Params file: {conf_loc}")
    log.info(f"Java heap: {java_mem}")

    # Pre-flight checks
    paths.require_sparx_inputs()
    paths.ensure_dirs()

    # Generate inputfile.tsv if missing
    if not paths.inputfile.exists():
        log.info(f"Input file not found at {paths.inputfile}, generating it...")
        generate_or_update_tsv(str(paths.inputfile), str(paths.base))
        if not paths.inputfile.exists():
            raise SparxError(
                f"Failed to generate input file at {paths.inputfile}"
            )
        log.info(f"Generated input file: {paths.inputfile}")

    datapath = str(paths.data) + os.sep
    database_path = str(paths.database) + os.sep
    QCdir = str(paths.qcdir) + os.sep
    SICdir = str(paths.sicdir)
    resultpath = str(paths.results) + os.sep
    PHRPOut = str(paths.phrp_out)
    MasicOut = str(paths.masic_out)

    decoder = pd.read_csv(paths.inputfile, sep="\t", header=0)
    log.info(f"Loaded {len(decoder)} rows from inputfile.tsv")

    start_time = time.time()

    # Step 0: RAW → mzML conversion + QC plots 
    VALID_EXTENSIONS = (".mzML", ".mgf", ".mzXML", ".ms2")

    log.info("=== Step 0: File conversion + QC plots ===")
    processed_db: set[str] = set()
    for _, row in decoder.iterrows():
        msfile = str(row["msfilename"])
        msraw = msfile + ".raw"
        msmzml = msfile + ".mzML"
        db_name = str(row["database"])
        outname = str(row["QCplot"]) + ".png"

        log.info("-" * 35)
        log.info(f"{msraw}")
        spec_exists = any(
            os.path.exists(os.path.join(datapath, msfile + ext)) 
            for ext in VALID_EXTENSIONS
        )

        if db_name not in processed_db:
            input_fasta = os.path.join(database_path, db_name)
            try:
                log.info(f"Processing {db_name} for redundancy...")
                remove_redundancy(input_fasta, input_fasta)
                log.info(f"Finished removing redundancy for {db_name}")
                processed_db.add(db_name)
            except FileNotFoundError:
                log.warning(f"File {input_fasta} not found. Skipping redundancy removal.")
                continue

        # Only convert from RAW if no valid spectrum file is found
        if not spec_exists:
            log.info(f"No compatible spectrum files detected; converting {msraw}")
            managed_run(
                [
                    "mono",
                    tool("ThermoRawFileParser/ThermoRawFileParser.exe"),
                    "-i", os.path.join(datapath, msraw),
                    "-f", "2",
                    "-L", "1-",
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE
            )
        mzml_full = os.path.join(datapath, msmzml)
        if os.path.exists(mzml_full):
            try:
                spectrum_count = 0
                # errors='ignore' prevents crashes from weird vendor-encoded characters
                with open(mzml_full, 'r', errors='ignore') as f:
                    for line in f:
                        if "<spectrum " in line:  # Note the space to avoid matching <spectrumList>
                            spectrum_count += 1
                log.info(f"Conversion complete. Detected {spectrum_count} spectra in {msmzml}")
            except Exception as e:
                log.warning(f"Could not count spectra in {msmzml}: {e}")

        if os.path.exists(os.path.join(datapath, msmzml)):
            generate_tic_plot(f"{datapath}{msmzml}", f"{QCdir}{outname}")

    # Step 1: MS-GF+ database search
    log.info("=== Step 1: Running MS-GF+ ===")

    peak_files = [f for f in os.listdir(datapath) if f.endswith(VALID_EXTENSIONS)]
    num_files  = len(peak_files)
   
    if num_files == 0:
        raise SparxError(f"No .mzML files found in {datapath}")

    # --- Universal Resource Guard ---
    try:
        total_ram_gb = psutil.virtual_memory().total / (1024**3)
    except ImportError:
        log.warning("psutil not installed. Assuming 32GB RAM for safety limits.")
        total_ram_gb = 32.0
    mem_per_job_gb = int(java_mem.upper().replace("G", ""))
    safe_ram_gb = total_ram_gb * 0.90
    max_jobs_by_ram = max(1, int(safe_ram_gb / mem_per_job_gb))

    total_cores, numa_nodes, cores_per_node = _get_core_config(total_ram_gb)
   
    MIN_THREADS_PER_JOB = 4
    max_jobs_by_cpu = max(1, total_cores // MIN_THREADS_PER_JOB)

    num_jobs = max(1, min(max_jobs_by_ram, max_jobs_by_cpu, num_files))
    jobs_per_node = max(1, -(-num_jobs // numa_nodes))  # Ceiling division to distribute jobs across NUMA nodes
    threads_per_job = max(1, cores_per_node // jobs_per_node)

    TASK_MULTIPLIER = 2
    if msgf_tasks is not None:
        tasks_per_job = max(1, min(msgf_tasks, threads_per_job))
    else:
        tasks_per_job = threads_per_job * TASK_MULTIPLIER

    log.info(f"System Resources -> Total RAM: {total_ram_gb}GB | Java/job: {mem_per_job_gb}GB")
    log.info(f"Concurrency Limits -> By RAM: {max_jobs_by_ram} | By CPU Cap: {max_jobs_by_cpu} | Final Jobs: {num_jobs}")
    log.info(f"Spectrum files: {len(peak_files)}  |  parallel jobs: {len(num_jobs)}")
    log.info(f"Thread Allocation -> Total cores: {total_cores} | jobs: {num_jobs} | Threasds/job: {cores_per_node} | Tasks/job: {tasks_per_job}")
    log.info(f"Spectrum files: {len(peak_files)}  |  Parallel jobs: {num_jobs}")

    # Dynamically map whichever format actually exists to its database
    peak_to_db = {}
    for _, row in decoder.iterrows():
        base_name = str(row["msfilename"])
        db_name = str(row["database"])
        
        for ext in VALID_EXTENSIONS:
            if os.path.exists(os.path.join(datapath, base_name + ext)):
                peak_to_db[base_name + ext] = db_name
                break
        else:
            log.warning(f"No compatible spectrum file found on disk for entry: {base_name}")

    with ThreadPoolExecutor(max_workers=num_jobs) as executor:
        futures = []
        for i, peak_file in enumerate(peak_files):
            if peak_file not in peak_to_db:
                log.warning(f"Skipping {peak_file}: Could not find matching database in inputfile.tsv")
                continue

            numa_node = i % numa_nodes if numa_nodes > 1 else None

            futures.append(
                executor.submit(
                    run_msgfplus,
                    peak_file, peak_to_db[peak_file], i,
                    datapath, database_path, conf_loc, java_mem,
                    msgf_threads=threads_per_job,
                    msgf_tasks=tasks_per_job,
                    numa_node=numa_node,
                )
            )

        for future in as_completed(futures):
            future.result()

    # Clean up temp FASTA files
    for f in os.listdir(database_path):
        if "temp_database" in f:
            try:
                os.remove(os.path.join(database_path, f))
                log.debug("Deleted temporary database file: %s", f)
            except Exception as e:
                log.warning("Error deleting %s: %s", f, e)

    log.info("All MS-GF+ processes completed.")

    # Step 2: mzid to TSV
    log.info("=== Step 2: mzid → TSV conversion ===")
    mzid_files = glob.glob(os.path.join(datapath, "*.mzid"))
    for mzid in mzid_files:
        filename_no_ext = os.path.splitext(os.path.basename(mzid))[0]
        tsv_out = os.path.join(resultpath, f"{filename_no_ext}.tsv")
        managed_run(
            [
                "java", f"-Xmx{java_mem}",
                "-cp", tool("MSGFPlus.jar"),
                "edu.ucsd.msjava.ui.MzIDToTsv",
                "-i", mzid,
                "-o", tsv_out,
                "-showDecoy", "1",
                "-unroll", "1",
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE
        )
        log.info(f"Converted: {mzid} to {tsv_out}")

    # Step 3: PHRP
    # --------------------------------------------------------
    log.info("=== Step 3: PHRP processing ===")
    os.makedirs(PHRPOut, exist_ok=True)

    tsv_files = [f for f in os.listdir(resultpath) if f.endswith(".tsv")]
    for tsv_file in tsv_files:
        tsv_path = os.path.join(resultpath, tsv_file)
        base_name = os.path.splitext(tsv_file)[0]
        sample_output_dir = os.path.join(PHRPOut, base_name)
        os.makedirs(sample_output_dir, exist_ok=True)

        log.info(f"PHRP: {tsv_file}")

        db_name = None
        for _, row in decoder.iterrows():
            if str(row["msfilename"]) == base_name:
                db_name = str(row["database"])
                break
        if db_name is None:
            log.warning(f"No DB match for {tsv_file}, using first entry.")
            db_name = str(decoder.iloc[0]["database"])

        log.info(f"Using database: {db_name}")
        log.info(f"PHRP input path exists? {os.path.exists(tsv_path)}" )
        log.info(f"PHRP input path = {tsv_path}")

        managed_run(
            [
                "mono",
                _normalize(tool("PHRP/PeptideHitResultsProcRunner.exe")),
                f"/I:{_normalize(tsv_path)}",
                f"/M:{_normalize(tool('PHRP/MSGFDB_Mods.txt'))}",
                f"/N:{_normalize(conf_loc)}",
                f"/T:{_normalize(tool('PHRP/Mass_Correction_Tags.txt'))}",
                f"/F:{_normalize(os.path.join(database_path, db_name))}",
                f"/O:{_normalize(sample_output_dir)}",
            ],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        log.info(f"PHRP done for {tsv_file}.")

    for sample_dir in os.listdir(PHRPOut):
        sample_dir_path = os.path.join(PHRPOut, sample_dir)
        if os.path.isdir(sample_dir_path):
            for f in os.listdir(sample_dir_path):
                shutil.move(os.path.join(sample_dir_path, f), os.path.join(PHRPOut, f))
            os.rmdir(sample_dir_path)
    log.info("All PHRP outputs moved to PHRPOut.")

    # STEP 4 — MASIC SIC generation
    os.makedirs(MasicOut, exist_ok=True)
    msraw_files = [f for f in os.listdir(datapath) if f.endswith(".raw")]
    if not msraw_files:
        log.warning(f"No .raw files found in {datapath}. MASIC step will be skipped.")

    _step4_masic_done = (
        os.path.isdir(MasicOut) and bool(os.listdir(MasicOut))
        and all(
            any(f.startswith(os.path.splitext(raw)[0]) for f in os.listdir(MasicOut))
            for raw in msraw_files
        )
    )
    if _step4_masic_done:
        log.info("[SKIP] Step-4: MASIC outputs already exist.")
    else:
        log.info("Starting MASIC SIC generation ...")
        successful, failed, skipped = run_masic_batch(
            datapath, MasicOut,
            workers=4,
            masic_exe=tool("MASIC/MASIC_Console.exe"),
            masic_params=tool("MASICParameters.xml"),
            stop_event=stop_event,
        )
        log.info(f"MASIC: {successful} succeeded, {failed} failed, {skipped} skipped.")
    
    #  Step 5: MASIC merger 
    log.info("=== Step 5: MASIC merger ===")
    try:
        temp_merger = os.path.join(resultpath, "temp_merger")
        masic_merger.copy_files_to_temp(PHRPOut, MasicOut, temp_merger)
        masic_merger.merge_files(temp_merger, SICdir)
        log.info("All intensity merging completed.")
    except Exception as e:
        log.error(f"Error in MASIC merger: {e}")
        raise SparxError(f"MASIC merger failed: {e}") from e

    # Summary
    elapsed_hours = (time.time() - start_time) / 3600
    log.info("*" * 60)
    log.info(f"Binning pipeline completed! Time: {elapsed_hours:.2f} hours")
    log.info("*" * 60)