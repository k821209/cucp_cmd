#!/usr/bin/env python3
"""
Count species from FASTA file using metadata CSV for organism names.
Usage: python count_species.py <fasta_file> <metadata_csv>

Examples:
    python count_species.py ref/Cuscuta_rbcL_with_outgroups.fasta input_rbcL/Cuscuta_rbcL_all_metadata.csv
"""

import csv
import sys
from collections import Counter
from Bio import SeqIO

def load_metadata(csv_file):
    """Load accession -> organism mapping from metadata CSV."""
    accession_to_organism = {}

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            accession = row.get('Accession', '')
            organism = row.get('Organism', 'Unknown')
            accession_to_organism[accession] = organism

    return accession_to_organism

def count_species_from_fasta(fasta_file, accession_to_organism):
    """Count sequences per species in FASTA file using metadata mapping."""
    species_counter = Counter()
    unmapped = []
    total_seqs = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        total_seqs += 1
        accession = record.id

        if accession in accession_to_organism:
            organism = accession_to_organism[accession]
        else:
            # Try to extract from description for outgroups
            desc = record.description
            parts = desc.split(None, 1)
            if len(parts) > 1:
                # Extract genus species from description
                words = parts[1].split()
                if len(words) >= 2:
                    organism = f"{words[0]} {words[1]}"
                else:
                    organism = parts[1]
            else:
                organism = "Unknown"
            unmapped.append(accession)

        species_counter[organism] += 1

    return species_counter, total_seqs, unmapped

def main():
    if len(sys.argv) < 3:
        fasta_file = "ref/Cuscuta_rbcL_with_outgroups.fasta"
        csv_file = "input_rbcL/Cuscuta_rbcL_all_metadata.csv"
    else:
        fasta_file = sys.argv[1]
        csv_file = sys.argv[2]

    print(f"FASTA file: {fasta_file}")
    print(f"Metadata CSV: {csv_file}")
    print()

    # Load metadata
    accession_to_organism = load_metadata(csv_file)
    print(f"Loaded {len(accession_to_organism)} accessions from metadata")
    print()

    # Count species
    species_counter, total_seqs, unmapped = count_species_from_fasta(fasta_file, accession_to_organism)

    # Summary
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total sequences in FASTA: {total_seqs}")
    print(f"Unique species/taxa: {len(species_counter)}")
    if unmapped:
        print(f"Unmapped accessions (outgroups): {len(unmapped)}")
        for acc in unmapped:
            print(f"  - {acc}")
    print()

    # Sort by count (descending), then by name
    sorted_species = sorted(species_counter.items(), key=lambda x: (-x[1], x[0]))

    # Print table
    print("=" * 60)
    print("SPECIES COUNT (sorted by abundance)")
    print("=" * 60)
    print(f"{'Species':<50} {'Count':>6}")
    print("-" * 60)

    for species, count in sorted_species:
        print(f"{species:<50} {count:>6}")

    count_sum = sum(species_counter.values())
    print("-" * 60)
    print(f"{'TOTAL':<50} {count_sum:>6}")

if __name__ == "__main__":
    main()
