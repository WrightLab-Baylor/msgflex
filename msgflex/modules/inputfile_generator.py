#!/usr/bin/env python3
import os

def generate_or_update_tsv(output_file, base_dir="."):
    data_dir = os.path.join(base_dir, "data")
    database_dir = os.path.join(base_dir, "database")

    # Get the .fasta or .faa file from the "database" directory
    database_files = [f for f in os.listdir(database_dir) if f.endswith(".fasta") or f.endswith(".faa")]
    if len(database_files) != 1:
        raise ValueError("There should be exactly one .fasta or .faa file in the 'database' directory.")
    database_file = database_files[0]

    # Try to regenerate from .raw files first
    if os.path.isdir(data_dir):
        ms_files = sorted([f for f in os.listdir(data_dir) if f.endswith(".raw")])
    else:
        ms_files = []

    if ms_files:
        # print(f"Generating new '{output_file}' from raw files in '{data_dir}'.")
        with open(output_file, "w") as file:
            file.write("msfilename\tprotein_input\tdatabase\tQCplot\n")
            for ms_file in ms_files:
                ms_file_name = os.path.splitext(ms_file)[0]
                file.write(f"{ms_file_name}\tgenome\t{database_file}\t{ms_file_name}\n")

    elif os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        # Fall back to updating existing file
        print(f"No raw files found. Updating database column in existing '{output_file}' to '{database_file}'.")
        updated_lines = []

        with open(output_file, "r") as file:
            header = file.readline().strip()
            updated_lines.append(header)

            for line in file:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                updated_lines.append(f"{parts[0]}\t{parts[1]}\t{database_file}\t{parts[3]}")

        with open(output_file, "w") as file:
            file.write("\n".join(updated_lines) + "\n")

    else:
        raise FileNotFoundError(
            f"No raw files found in '{data_dir}' and no existing '{output_file}' to update."
        )

