#!/usr/bin/env python3
"""
Standalone parallel MASIC processing script for proteomics raw files.
"""
import os
import logging
import subprocess
import time
import multiprocessing
import traceback
import psutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import threading

logger = logging.getLogger(__name__)

# Tuneable constants 
MAX_RETRIES = 2
RETRY_BASE_DELAY = 30
WORKER_STAGGER_DELAY = 5
LARGE_FILE_GB = 1.75
MIN_FREE_RAM_GB = 4.0
MEM_POLL_INTERVAL = 10
MASIC_TIMEOUT = 7_200
LARGE_FILE_TIMEOUT = 18_000
LARGE_FILE_MAX_WORKERS = 2
EST_RAM_PER_LARGE_FILE_GB = 16.0

_launch_lock = None
_mem_semaphore = None
_stop_flag = None

#  Output file helpers (single source of truth) 
def _masic_output_files(base_name, output_dir):
    """Return categorised output paths for one raw file."""
    stems = {
        "all": [
            f"{base_name}_SICs.xml",
            f"{base_name}_SICstats.txt",
            f"{base_name}_ScanStats.txt",
            f"{base_name}_ScanStatsEx.txt",
            f"{base_name}_ScanStatsConstant.txt",
            f"{base_name}_MSMS_scans.csv",
            f"{base_name}_MS_scans.csv",
        ],
        "required": [
            f"{base_name}_SICstats.txt",
            f"{base_name}_ScanStats.txt",
        ],
    }
    return {
        key: [os.path.join(output_dir, f) for f in names]
        for key, names in stems.items()
    }

def _files_ok(paths):
    return all(os.path.exists(p) and os.path.getsize(p) > 0 for p in paths)

# Memory gate
def _wait_for_ram(needed_gb: float, poll: int = MEM_POLL_INTERVAL) -> None:
    """
    Block until the system has at least `needed_gb` of free RAM.
    """
    while True:
        free_gb = psutil.virtual_memory().available / (1024 ** 3)
        if free_gb >= needed_gb:
            return
        logger.info(f"RAM gate: only {free_gb:.1f} GB free, need {needed_gb:.1f} GB; waiting {poll}s")
        time.sleep(poll)

def _pool_initializer(lock, mem_sem, stop_flag):
    global _launch_lock, _mem_semaphore, _stop_flag
    _launch_lock = lock
    _mem_semaphore = mem_sem
    _stop_flag = stop_flag

def _kill_tree(pid: int) -> None:
    """Kill a process and all its children."""
    try:
        parent = psutil.Process(pid)
        children = parent.children(recursive=True)
        for child in children:
            try:
                child.kill()
            except psutil.NoSuchProcess:
                pass
        parent.kill()
    except psutil.NoSuchProcess:
        pass


def masic_console(msraw, datapath, masic_out, masic_exe, masic_params,
                  stagger_delay=WORKER_STAGGER_DELAY):
    """
    Run MASIC for one raw file. Checks module-level _stop_flag every 2s.
    """
    msraw_basename = os.path.basename(msraw)
    msraw_path = os.path.join(datapath, msraw_basename)

    if not os.path.isfile(msraw_path):
        logger.error(f"Raw file not found: {msraw_path}")
        return False, 2

    file_size_gb = os.path.getsize(msraw_path) / (1024 ** 3)
    is_large = file_size_gb >= LARGE_FILE_GB
    logger.info(f"Starting {msraw_basename} ({file_size_gb:2f}.2f GB, large={is_large})")

    os.makedirs(masic_out, exist_ok=True)
    progress_file = os.path.join(masic_out, f"{msraw_basename}_progress.log")

    with open(progress_file, 'w') as pf:
        pf.write(f"Started : {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        pf.write(f"File : {msraw_path}\n")
        pf.write(f"Size : {file_size_gb:.2f} GB\n")
        pf.write(f"=" * 60 + "\nMASIC OUTPUT:\n" + "=" * 60 + "\n")

    masic_cmd = [
        'mono', masic_exe,
        f'/I:{msraw_path}',
        f'/O:{masic_out}',
        f'/P:{masic_params}'
    ]

    timeout = LARGE_FILE_TIMEOUT if is_large else MASIC_TIMEOUT

    acquired_sem = False
    if is_large and _mem_semaphore is not None:
        _mem_semaphore.acquire()
        acquired_sem = True

    try:
        if is_large:
            _wait_for_ram(MIN_FREE_RAM_GB)

        with open(progress_file, "a") as masic_log:
            if _launch_lock is not None:
                with _launch_lock:
                    process = subprocess.Popen(
                        masic_cmd, stdout=masic_log, stderr=masic_log,
                    )
                    time.sleep(stagger_delay)
            else:
                process = subprocess.Popen(
                    masic_cmd, stdout=masic_log, stderr=masic_log,
                )

            deadline = time.monotonic() + timeout
            while True:
                try:
                    process.wait(timeout=2.0)
                    break
                except subprocess.TimeoutExpired:
                    if _stop_flag is not None and _stop_flag.value:
                        logger.warning(f"Stop requested : killing MASIC process {msraw}")
                        _kill_tree(process.pid)
                        return False, -2
                    if time.monotonic() >= deadline:
                        _kill_tree(process.pid)
                        logger.error(f"{msraw} timed out : killed")
                        return False, -1
    finally:
        if acquired_sem:
            _mem_semaphore.release()

    rc = process.returncode
    with open(progress_file, "a") as pf:
        pf.write("=" * 60 + f"\nFinished: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        pf.write(f"Return code: {rc}\n")

    if rc != 0:
        reason = {
            32: "File access error — locked or I/O contention",
             1: "General MASIC error — check progress log",
             2: "File not found or invalid input",
            -1: "Killed by timeout guard",
            -2: "Stopped by user",
        }.get(rc, "Unknown error")
        logger.error(f"{msraw_basename}: MASIC failed rc={rc} ({reason})")
        with open(progress_file, "a") as pf:
            pf.write(f"FAILED: rc={rc} — {reason}\n")
        return False, rc

    logger.info(f"Completed {msraw_basename}")
    return True, 0

def process_masic(
    msraw,
    datapath,
    masic_out,
    masic_exe,
    masic_params,
    stagger_delay = WORKER_STAGGER_DELAY
    ):

    msraw_basename = os.path.basename(msraw)
    base_name = msraw_basename.rsplit(".", 1)[0]
    outputs = _masic_output_files(base_name, masic_out)

    try:
        if _files_ok(outputs["all"]):
            logger.info(f"All outputs exist for {msraw_basename}: skipping")
            return msraw_basename, True, "Skipped - all files exist"

        if _files_ok(outputs["required"]):
            logger.info(f"Required outputs exist for {msraw_basename}: skipping")
            return msraw_basename, True, "Skipped - required files exist"

        last_rc = None
        for attempt in range(1, MAX_RETRIES + 1):
            if _stop_flag is not None and _stop_flag.value:
                return msraw_basename, False, "Stopped by user"

            if attempt > 1:
                delay = RETRY_BASE_DELAY * (2 ** (attempt - 2))
                logger.warning(f"RETRY {msraw_basename}: attempt {attempt}/{MAX_RETRIES}: back-off {delay}s (prev rc={last_rc})")
                time.sleep(delay)

            success, rc = masic_console(
                msraw, datapath, masic_out,
                masic_exe, masic_params,
                stagger_delay=stagger_delay,
            )
            last_rc = rc

            if _files_ok(outputs["required"]):
                if rc in (-1, -2):
                    pass
                else:
                    if not success:
                        logger.info(f"{msraw_basename}: rc={rc} but required files OK; treating as success")
                    return msraw_basename, True, None

            if rc not in (32, -1):
                break

        missing = [
            os.path.basename(f) for f in outputs["required"]
            if not os.path.exists(f)
        ]
        return (
            msraw_basename, False,
            f"Failed after {attempt} attempt(s) rc={last_rc}; "
            f"missing: {', '.join(missing)}"
        )

    except Exception:
        logger.exception(f"Exception in {msraw_basename}")
        return msraw_basename, False, traceback.format_exc()

def _classify_files(file_list, input_dir):
    normal, large = [], []
    for f in file_list:
        size = os.path.getsize(os.path.join(input_dir, f))
        (large if size / (1024 ** 3) >= LARGE_FILE_GB else normal).append(f)
    return normal, large

def run_masic_batch(
    input_dir,
    output_dir,
    workers=4,
    masic_exe="MASIC/MASIC_Console.exe",
    masic_params="MASICParameters.xml",
    stop_event = None
    ):
    """
    Run MASIC on all .raw files in input_dir
    """
    for label, path, check in [
        ("Input directory", input_dir, os.path.isdir),
        ("MASIC executable", masic_exe, os.path.isfile),
        ("MASIC parameters", masic_params, os.path.isfile),
    ]:
        if not check(path):
            raise ValueError(f"{label} not found: {path}")

    os.makedirs(output_dir, exist_ok=True)

    msraw_files = [
        e.name for e in os.scandir(input_dir)
        if e.is_file() and e.name.lower().endswith(".raw")
    ]
    if not msraw_files:
        raise ValueError(f"No .raw files found in {input_dir}")

    logger.info("=" * 60)
    logger.info("MASIC Parallel Processing")
    logger.info(f"Input : {input_dir}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Files : {len(msraw_files)}")
    logger.info(f"Workers: {workers}")
    logger.info("=" * 60)

    files_to_process, files_already_done = [], []

    for msraw in msraw_files:
        base_name = msraw.rsplit(".", 1)[0]
        outputs = _masic_output_files(base_name, output_dir)

        if _files_ok(outputs["all"]):
            files_already_done.append(msraw)
            logger.debug(f"[SKIP-ALL] {msraw}")
        elif _files_ok(outputs["required"]):
            files_already_done.append(msraw)
            logger.debug(f"[SKIP-REQ] {msraw}")
        else:
            files_to_process.append(msraw)
            missing = [
                os.path.basename(f) for f in outputs["required"]
                if not os.path.exists(f)
            ]
            logger.debug(f"[QUEUE] {msraw} missing: {', '.join(missing)}")

    logger.info(f"Already done: {len(files_already_done)} | To process: {len(files_to_process)}")

    if not files_to_process:
        logger.info("Nothing to do.")
        return [], [], files_already_done

    normal_files, large_files = _classify_files(files_to_process, input_dir)

    if large_files:
        total_ram = psutil.virtual_memory().total / (1024 ** 3)
        ram_limited = max(1, int(total_ram // EST_RAM_PER_LARGE_FILE_GB))
        safe_workers = min(
            LARGE_FILE_MAX_WORKERS,
            ram_limited,
            len(large_files),
        )
        logger.info(f"Large files detected ({len(large_files)}). System RAM: {total_ram:.1f}GB => max {safe_workers} concurrent")
        for f in large_files:
            gb = os.path.getsize(os.path.join(input_dir, f)) / (1024 ** 3)
            logger.info(f" [LARGE] {f} ({gb:.2f} GB)", f, gb)
    else:
        safe_workers = 1

    manager = multiprocessing.Manager()
    launch_lock = manager.Lock()
    mem_semaphore = manager.Semaphore(safe_workers)
    stop_flag = multiprocessing.Value('b', False)

    if stop_event is not None:
        def _sync_stop():
            stop_event.wait()
            stop_flag.value = True
            logger.warning(
                "MASIC stop flag set — workers will exit after current poll",
            )
        threading.Thread(target=_sync_stop, daemon=True).start()

    def _run_batch(file_list, max_workers, stagger):
        if not file_list:
            return []
        pfunc = partial(
            process_masic,
            datapath=input_dir,
            masic_out=output_dir,
            masic_exe=masic_exe,
            masic_params=masic_params,
            stagger_delay=stagger,
        )
        executor = ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_pool_initializer,
            initargs=(launch_lock, mem_semaphore, stop_flag),
        )
        try:
            futures = {executor.submit(pfunc, f): f for f in file_list}
            results = []
            for fut in as_completed(futures):
                if stop_flag.value:
                    break
                try:
                    results.append(fut.result())
                except Exception as exc:
                    fname = futures[fut]
                    logger.error(f"[CRASH] {fname}: {exc}")
                    results.append((fname, False, f"Worker crashed: {exc}"))
            return results
        finally:
            executor.shutdown(wait=True, cancel_futures=True)

    start = time.time()

    all_results: list = []
    if normal_files:
        logger.info(f"Normal batch: {len(normal_files)} file(s), {workers} worker(s)")
        all_results += _run_batch(normal_files, workers, WORKER_STAGGER_DELAY)

    if large_files and not stop_flag.value:
        logger.info(f"Large batch : {len(large_files)} file(s), {safe_workers} worker(s)")
        all_results += _run_batch(
            large_files, safe_workers, WORKER_STAGGER_DELAY * 2,
        )

    manager.shutdown()

    elapsed = time.time() - start
    successful, failed = [], []
    skipped = list(files_already_done)

    logger.info("=" * 60 + "\nPROCESSING SUMMARY\n" + "=" * 60)
    for basename, success, message in all_results:
        if success:
            if message and "Skipped" in message:
                skipped.append(basename)
                logger.info(f"[SKIPPED] {basename}")
            else:
                successful.append(basename)
                logger.info(f"[SUCCESS] {basename}")
        else:
            failed.append(basename)
            logger.error(f"[FAILED] {basename} : {message}")

    logger.info(
        f"Total: {len(msraw_files)} | Done: {len(files_already_done)} | New: {len(successful)} | Failed: {len(failed)}| Elapsed:{elapsed / 60:.1f} min"
    )
    return successful, failed, skipped
