"""
QUANTIX module 1 - Prepare SICs (Selected Ion Chromatograms).

Reads PHRP+MASIC merged files from SICdir, selects key columns,
picks monoisotopic peaks, and enriches with:
  - is_decoy (0/1 integer, consistent with XPECTRA)
  - MSMSScore = -log10(SpecEValue)
  - absPPM = abs(DelM_PPM)

Author: Tulasi Rao Relangi, PhD
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

COLUMNS_LIST = [
    "Scan", "PrecursorMZ", "DelM_PPM", "Peptide", "Protein",
    "MSGFScore", "MSGFDB_SpecEValue", "EValue", "QValue", "PepQValue",
    "ElutionTime", "ParentIonIntensity", "PeakArea", "StatMomentsArea",
]

FILE_EXTENSION = "_PlusSICStats"
DECOY_PREFIX = "XXX_"

def prepare_SICs(
    input_dir: str,
    output_dir: str,
    columns_list: list[str] = COLUMNS_LIST,
    file_extension: str = FILE_EXTENSION,
    delimiter: str = "\t",
    decoy_prefix: str = DECOY_PREFIX,
) -> None:
    """Extract selected columns, pick monoisotopic peaks, add scores.

    Writes one *_PlusSICStats.tsv per input file into output_dir.
    """
    os.makedirs(output_dir, exist_ok=True)

    for file_name in os.listdir(input_dir):
        if file_extension not in file_name:
            continue

        # print(f"Processed: {file_name}")
        input_file_path = os.path.join(input_dir, file_name)
        output_file_path = os.path.join(
            output_dir,
            file_name.replace("_fht_PlusSICStats.txt", "_PlusSICStats.tsv"),
        )

        df = pd.read_csv(input_file_path, delimiter=delimiter)

        # Pick monoisotopic peak per scan (highest intensity)
        mono_idx = df.groupby("Scan")["ParentIonIntensity"].idxmax()
        selected = df.loc[mono_idx, columns_list].copy()

        # Rename for downstream consistency
        selected.rename(
            columns={"MSGFDB_SpecEValue": "SpecEValue", "Scan": "ScanNum"},
            inplace=True,
        )

        selected["is_decoy"] = (
            selected["Protein"]
            .astype(str)
            .str.startswith(decoy_prefix)
            .astype(int)
        )

        spec_e = pd.to_numeric(selected["SpecEValue"], errors="coerce")
        selected["MSMSScore"] = -np.log10(spec_e.replace(0, np.nan).fillna(1))

        # absPPM
        if "DelM_PPM" in selected.columns:
            selected["absPPM"] = pd.to_numeric(
                selected["DelM_PPM"], errors="coerce"
            ).abs()

        selected.to_csv(output_file_path, sep="\t", index=False)
        # print(f"Exported to: {output_file_path}")