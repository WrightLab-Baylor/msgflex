#!/usr/bin/env python3
import os
import shutil
import pandas as pd

import logging
log = logging.getLogger(__name__)

def copy_files_to_temp(phrp_out, masic_out, temp_merger):
    os.makedirs(temp_merger, exist_ok=True)
    for file in os.listdir(phrp_out):
        if file.endswith("_fht.txt"):
            shutil.copy(os.path.join(phrp_out, file), os.path.join(temp_merger, file))
    for file in os.listdir(masic_out):
        if file.endswith("_ScanStats.txt") or file.endswith("_SICstats.txt"):
            shutil.copy(os.path.join(masic_out, file), os.path.join(temp_merger, file))

def compute_peak_width_minutes(scan_stats, fht_row):
    try:
        start_raw = fht_row["PeakScanStart"]
        end_raw   = fht_row["PeakScanEnd"]

        # NaN means MASIC found no peak for this scan, return 0
        if pd.isna(start_raw) or pd.isna(end_raw):
            return 0

        start_scan = int(start_raw)
        end_scan   = int(end_raw)
    except (ValueError, TypeError):
        log.debug(f"Could not parse PeakScanStart/End: {fht_row.get('PeakScanStart')}, {fht_row.get('PeakScanEnd')}")
        return 0

    # Get the corresponding scan stats for start and end scans from scan_stats dictionary
    start_scan_stats = scan_stats.get(start_scan)
    end_scan_stats = scan_stats.get(end_scan)

    # If no stats are found for the scans, return 0
    if not start_scan_stats or not end_scan_stats:
        log.debug(f"Scan stats not found for start_scan: {start_scan}, end_scan: {end_scan}")
        return 0

    try:
        # ElutionTime in scan_stats corresponds to ScanTime in ScanStats
        start_time_minutes = float(start_scan_stats["ElutionTime"])
        end_time_minutes = float(end_scan_stats["ElutionTime"])
    except (ValueError, KeyError) as e:
        log.error(f"Error parsing ElutionTime: {e}")
        return 0

    # Calculate the peak width in minutes
    return end_time_minutes - start_time_minutes

def merge_files(temp_merger, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    file_groups = {}
    for file in os.listdir(temp_merger):
        if file.endswith("_fht.txt"):
            base_name = file.replace("_fht.txt", "")
        elif file.endswith("_ScanStats.txt"):
            base_name = file.replace("_ScanStats.txt", "")
        elif file.endswith("_SICstats.txt"):
            base_name = file.replace("_SICstats.txt", "")
        else:
            continue
        if base_name not in file_groups:
            file_groups[base_name] = []
        file_groups[base_name].append(file)

    if not file_groups:
        log.error(f"No files found to process.")
        return

    for base_name, files in file_groups.items():
        fht_file = next((f for f in files if f.endswith("_fht.txt")), None)
        scanstats_file = next((f for f in files if f.endswith("_ScanStats.txt")), None)
        sicstats_file = next((f for f in files if f.endswith("_SICstats.txt")), None)

        if not fht_file or not scanstats_file or not sicstats_file:
            log.info(f"Skipping {base_name}: Missing required files.")
            continue

        try:
            fht_df = pd.read_csv(os.path.join(temp_merger, fht_file), sep='\t')
            scanstats_df = pd.read_csv(os.path.join(temp_merger, scanstats_file), sep='\t')
            sicstats_df = pd.read_csv(os.path.join(temp_merger, sicstats_file), sep='\t')
        except Exception as e:
            log.error(f"Error reading files for {base_name}: {e}")
            continue

        # Rename column in SICStats
        scanstats_df.rename(columns={"ScanTime":"ElutionTime"}, inplace=True)
        sicstats_df.rename(columns={"OptimalPeakApexScanNumber": "OptimalScanNumber"}, inplace=True)
        sicstats_df.rename(columns={"MZ": "ParentIonMZ"}, inplace=True)

        # Create dictionaries
        scan_stats_dict = scanstats_df.set_index("ScanNumber")[
            ["ElutionTime", "ScanType", "TotalIonIntensity", "BasePeakIntensity", "BasePeakMZ"]
        ].to_dict("index")

        sic_stats_dict = sicstats_df.set_index("FragScanNumber")[
            [
                "OptimalScanNumber", "PeakMaxIntensity", "PeakSignalToNoiseRatio", 
                "FWHMInScans", "PeakArea", "ParentIonIntensity", "ParentIonMZ", 
                "StatMomentsArea", "PeakScanStart", "PeakScanEnd"
            ]
        ].to_dict("index")

        for column in ["ElutionTime", "ScanType", "TotalIonIntensity", "BasePeakIntensity", "BasePeakMZ"]:
            fht_df[column] = fht_df["Scan"].apply(lambda x: scan_stats_dict.get(x, {}).get(column, None))

        # Create new columns from sic_stats_dict (OptimalScanNumber, PeakMaxIntensity, etc.)
        for column in [
            "OptimalScanNumber", "PeakMaxIntensity", "PeakSignalToNoiseRatio", 
            "FWHMInScans", "PeakArea", "ParentIonIntensity", "ParentIonMZ", 
            "StatMomentsArea", "PeakScanStart", "PeakScanEnd"
        ]:
            fht_df[column] = fht_df["Scan"].apply(lambda x: sic_stats_dict.get(x, {}).get(column, None))

        fht_df["PeakWidthMinutes"] = fht_df.apply(
            lambda row: round(compute_peak_width_minutes(scan_stats_dict, row), 4), axis=1
        )

        # Save the merged DataFrame
        output_file = os.path.join(output_dir, f"{base_name}_fht_PlusSICStats.txt")
        fht_df.to_csv(output_file, sep='\t', index=False)
        log.info(f"Merged file saved {output_file}:")

    # Delete temp_merger folder after processing
    shutil.rmtree(temp_merger)
    log.info(f"Deleted {temp_merger} directory.")