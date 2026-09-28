#!/usr/bin/env python3
from pathlib import Path
import re

# ---- Helpers ----
def get_gene_from_info(info):
    if not isinstance(info, str) or not info:
        return "Unknown"
    m = re.search(r'(?:^|\s|;|\|)(?:GN|Gene)=([^\s;|]+)', info, flags=re.IGNORECASE)
    return m.group(1) if m else "Unknown"

def get_protein_function(protein_id, proteins):
    protein_info = proteins.get(protein_id)
    if not protein_info:
        return {'function': 'Unknown', 'gene': 'Unknown'}
    function_string = protein_info.get('Function', '') or ''
    if "OS" in function_string:
        parts = function_string.split('OS', 1)
        function = parts[0].strip()
        info = "OS" + parts[1].strip() if len(parts) > 1 else ''
    else:
        function = function_string.strip()
        info = ''

    return {
        'function': function or 'Unknown',
        'gene': get_gene_from_info(info) or 'Unknown'
    }

def build_protein_info_map(fasta_paths):
    """
    Parse one or more FASTA files and return a combined map:
    """
    prot2info = {}

    # Ensure fasta_paths is a list
    if isinstance(fasta_paths, (str, Path)):
        fasta_paths = [fasta_paths]

    for fasta_path in fasta_paths:
        fasta_path = Path(fasta_path)
        if not fasta_path.exists():
            print(f"Warning: FASTA not found, skipping: {fasta_path}")
            continue

        # If directory, glob for fasta/faa files
        if fasta_path.is_dir():
            files = sorted(fasta_path.glob("*.fasta")) + sorted(fasta_path.glob("*.faa"))
            if not files:
                print(f"Warning: No FASTA files found in directory: {fasta_path}")
                continue
        else:
            files = [fasta_path]

        for f in files:
            with open(f, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line.startswith(">"):
                        continue
                    header = line[1:].rstrip("\n")
                    protein_id = header.split()[0]
                    rest = header[len(protein_id):].lstrip()
                    prot2info[protein_id] = {"Function": rest}

    return prot2info
