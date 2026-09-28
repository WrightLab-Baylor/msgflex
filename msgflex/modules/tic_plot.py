#!/usr/bin/en python3

import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pyteomics import mzml

def generate_tic_plot(input_mzml, output_png, *, max_minutes: float=999.0, dpi=300):
    """
    Generate TIC plot using Pyteomics
    """
    input_mzml = Path(input_mzml)
    output_png = Path(output_png)

    if not input_mzml.exists():
        raise FileNotFoundError(f"mzML file not found: {input_mzml}")

    retention_times = []
    intensities = []

    # Streaming read (IMPORTANT for large files)
    with mzml.MzML(str(input_mzml)) as reader:
        for spectrum in reader:
            # Only MS1 spectra
            if spectrum.get("ms level") != 1:
                continue

            rt = spectrum["scanList"]["scan"][0]["scan start time"]

            # Convert to minutes if needed
            if spectrum["scanList"]["scan"][0].get("unitName") == "second":
                rt = rt / 60.0

            intensity_array = spectrum.get("intensity array")
            if intensity_array is None:
                continue

            tic = float(np.sum(intensity_array))

            retention_times.append(rt)
            intensities.append(tic)

    if not retention_times:
        raise ValueError(f"No MS1 spectra found in {input_mzml}")

    retention_times = np.array(retention_times)
    intensities = np.array(intensities)

    # Normalize RT to start at 0
    retention_times = retention_times - retention_times[0]

    # Filter by max_minutes
    mask = retention_times < max_minutes
    retention_times = retention_times[mask]
    intensities = intensities[mask]

    if len(retention_times) == 0:
        raise ValueError("No data left after RT filtering")

    # Plot
    plt.figure(figsize=(6, 4))
    plt.plot(retention_times, intensities, color="#00283F", linestyle='-')

    plt.xlabel("Time (min)")
    plt.ylabel("Intensity (cps)")
    plt.xlim([0, float(retention_times.max())])
    plt.ylim(float(intensities.min()), float(intensities.max()))

    plt.tight_layout()
    plt.savefig(output_png, dpi=dpi)
    plt.close()

def main():
    if len(sys.argv) != 3:
        print("Usage: tic_plot.py <input.mzML> <output.png>")
        sys.exit(1)

    generate_tic_plot(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()