#!/usr/bin/env python3
"""
Select maximum 3 species per genus from FASTA file.
Usage: python select_max3_per_genus.py <input_fasta> <output_fasta>
"""

from Bio import SeqIO
from collections import defaultdict
import sys

def extract_genus(description):
    """Extract genus from FASTA description."""
    parts = description.split()
    if len(parts) >= 2:
        return parts[1]  # Second word is genus
    return "Unknown"

def select_max_per_genus(input_file, output_file, max_per_genus=3):
    """Select maximum N species per genus."""
    genus_records = defaultdict(list)

    # Group records by genus
    for record in SeqIO.parse(input_file, "fasta"):
        genus = extract_genus(record.description)
        genus_records[genus].append(record)

    # Select up to max_per_genus from each genus
    selected_records = []

    print("=" * 60)
    print(f"Selecting max {max_per_genus} species per genus")
    print("=" * 60)
    print(f"{'Genus':<20} {'Total':>8} {'Selected':>10}")
    print("-" * 60)

    for genus in sorted(genus_records.keys()):
        records = genus_records[genus]
        total = len(records)
        selected = records[:max_per_genus]
        selected_records.extend(selected)
        print(f"{genus:<20} {total:>8} {len(selected):>10}")

    print("-" * 60)
    print(f"{'TOTAL':<20} {sum(len(r) for r in genus_records.values()):>8} {len(selected_records):>10}")
    print()

    # Write selected records
    SeqIO.write(selected_records, output_file, "fasta")
    print(f"Written {len(selected_records)} sequences to {output_file}")

    # List selected sequences
    print()
    print("=" * 60)
    print("Selected sequences:")
    print("=" * 60)
    for record in selected_records:
        print(f"  {record.id} {record.description.split(None, 1)[1][:60]}...")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        input_file = "ref/Orobanchaceae_cp_genomes_representative.fasta"
        output_file = "ref/Orobanchaceae_cp_genomes_max3.fasta"
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]

    select_max_per_genus(input_file, output_file, max_per_genus=3)
