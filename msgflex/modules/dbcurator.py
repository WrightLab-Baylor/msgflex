#!/usr/bin/env python3
from Bio import SeqIO
import argparse

def remove_redundancy(input_file, output_file):
    # Dictionary to store unique sequences
    unique_sequences = {}

    # Open input FASTA file
    with open(input_file, 'r') as f:
        # Iterate over sequences in the input FASTA file
        for record in SeqIO.parse(f, 'fasta'):
            # Add to dictionary if sequence is unique
            if str(record.seq) not in unique_sequences:
                unique_sequences[str(record.seq)] = record

    # Get unique records
    unique_records = unique_sequences.values()

    # Open output FASTA file
    with open(output_file, 'w') as f:
        # Write unique records to the output FASTA file
        SeqIO.write(unique_records, f, 'fasta')