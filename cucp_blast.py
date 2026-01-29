#!/usr/bin/env python3
"""
CUCP BLAST Analysis Command-Line Tool

A standalone command-line tool for running BLAST analysis and phylogenetic tree
construction for Cuscuta (dodder) comparative phylogenomics.

Usage:
    python cucp_blast.py --reference ref.fa --query query.fa --output results/
    python cucp_blast.py -r ref.fa -q query.fa -o results/ --run-iqtree
"""

import argparse
import json
import subprocess
import os
from pathlib import Path
import pandas as pd
import numpy as np
import io
import base64
import sys
from datetime import datetime

# Optional imports for visualization
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sb
    from Bio import Phylo
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False
    print("Warning: matplotlib, seaborn, or biopython not available. Visualization will be skipped.")

# Get the directory of this script
SCRIPT_DIR = Path(__file__).parent.resolve()

def get_binary_path(binary_name):
    """
    Get the path to a binary from system PATH (conda environment).
    """
    try:
        result = subprocess.run(['which', binary_name], capture_output=True, text=True)
        if result.returncode == 0:
            return Path(result.stdout.strip())
    except:
        pass
    return None

# Binary paths from conda environment (system PATH)
BLAST_PATH = get_binary_path("blastn")
MAKEBLASTDB_PATH = get_binary_path("makeblastdb")
IQTREE_PATH = get_binary_path("iqtree")


def setup_permissions():
    """Set executable permissions for binary tools."""
    for binary in [BLAST_PATH, MAKEBLASTDB_PATH, IQTREE_PATH]:
        if binary and binary.exists():
            os.chmod(binary, 0o755)


def read_fasta(fasta_path):
    """Read a FASTA file and return a dictionary of {name: sequence}."""
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]  # Take first word after >
                current_seq = []
            else:
                current_seq.append(line)

        if current_name:
            sequences[current_name] = ''.join(current_seq)

    return sequences


def write_fasta(sequences, output_path):
    """Write sequences dictionary to a FASTA file."""
    with open(output_path, 'w') as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")


def compute_pairwise_distances(df_profile):
    """Compute pairwise distance matrix from genotype profile."""
    profile_values = df_profile.values
    n = profile_values.shape[0]
    dist_matrix = np.zeros((n, n))

    for i in range(n):
        comparisons = (profile_values[i] == profile_values).astype(int)
        matching_fractions = comparisons.mean(axis=1)
        dist_matrix[i] = matching_fractions
        dist_matrix[:, i] = matching_fractions

    df_dist = pd.DataFrame(dist_matrix, index=df_profile.index, columns=df_profile.index)
    return df_dist


def generate_fasta_from_dataframe(df_var):
    """Generate FASTA content from genotype matrix dataframe."""
    fasta_lines = []
    sample_columns = [col for col in df_var.columns if col != 'query']

    for sample_name in sample_columns:
        fasta_lines.append(f">{sample_name}")
        sequence = ''.join(df_var[sample_name].astype(str).values)
        fasta_lines.append(sequence)

    return '\n'.join(fasta_lines)


def run_blast(query_fasta, reference_fasta, work_dir, sample_name="query"):
    """
    Run BLAST analysis.

    Args:
        query_fasta: Path to query FASTA file
        reference_fasta: Path to reference FASTA file
        work_dir: Working directory for intermediate files
        sample_name: Name for the sample

    Returns:
        dict with status and results
    """
    try:
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        # Read and merge reference with query
        merged_path = work_dir / "merged.fa"

        # Copy reference sequences
        ref_seqs = read_fasta(reference_fasta)
        query_seqs = read_fasta(query_fasta)

        # Merge sequences
        merged_seqs = {**ref_seqs, **query_seqs}
        write_fasta(merged_seqs, merged_path)

        print(f"Merged {len(ref_seqs)} reference sequences with {len(query_seqs)} query sequences")

        # Create BLAST database
        db_path = work_dir / "blastdb"
        makeblastdb_cmd = [
            str(MAKEBLASTDB_PATH),
            "-in", str(merged_path),
            "-dbtype", "nucl",
            "-out", str(db_path)
        ]

        print("Creating BLAST database...")
        subprocess.run(makeblastdb_cmd, check=True, capture_output=True)

        # Need a base sequence for BLAST query
        # Use the first reference sequence as the base
        base_name = list(ref_seqs.keys())[0]
        base_path = work_dir / "base.fa"
        write_fasta({base_name: ref_seqs[base_name]}, base_path)

        # Run BLAST
        blast_output = work_dir / f"{sample_name}_blast_results.txt"
        blast_cmd = [
            str(BLAST_PATH),
            "-task", "blastn",
            "-query", str(base_path),
            "-db", str(db_path),
            "-outfmt", "13",  # JSON output
            "-out", str(blast_output),
            "-evalue", "1e-5"
        ]

        print("Running BLAST...")
        subprocess.run(blast_cmd, check=True, capture_output=True)

        # Read BLAST results
        json_output = str(blast_output) + "_1.json"
        if os.path.exists(json_output):
            with open(json_output, 'r') as f:
                blast_results = f.read()
        else:
            # Try without _1 suffix
            with open(str(blast_output), 'r') as f:
                blast_results = f.read()

        return {
            "status": "success",
            "blast_results": blast_results,
            "base_sequence": ref_seqs[base_name],
            "query_names": list(query_seqs.keys())
        }

    except Exception as e:
        return {
            "status": "error",
            "error": str(e)
        }


def process_blast_results(blast_json_str, base_sequence, query_names, run_iqtree=False, work_dir=None, output_dir=None):
    """
    Process BLAST results to generate variation matrix and optionally run IQ-TREE.

    Args:
        blast_json_str: JSON string of BLAST results
        base_sequence: Base sequence string
        query_names: List of query sequence names
        run_iqtree: Whether to run IQ-TREE for phylogenetic analysis
        work_dir: Working directory
        output_dir: Output directory for results

    Returns:
        dict with processed results
    """
    try:
        blast_data = json.loads(blast_json_str)

        # Initialize variations dictionary
        variations = {i+1: {'query': base.upper()}
                     for i, base in enumerate(base_sequence)}

        # Extract query info
        query_title = blast_data['BlastOutput2']['report']['results']['search']['query_title']
        query_length = blast_data['BlastOutput2']['report']['results']['search']['query_len']

        # Process hits
        for hit in blast_data['BlastOutput2']['report']['results']['search']['hits']:
            title_parts = hit['description'][0]['title'].split()

            # Extract scientific name
            if len(title_parts) >= 3:
                scientific_name = '_'.join(title_parts[0:3])
            else:
                scientific_name = '_'.join(title_parts)

            for hsp in hit['hsps']:
                alignment_length = hsp['align_len']
                if alignment_length < 100:
                    continue

                query_seq = hsp['qseq']
                hit_seq = hsp['hseq']
                query_start = hsp['query_from']

                query_pos = query_start
                for q, h in zip(query_seq, hit_seq):
                    if q != '-':
                        if query_pos not in variations:
                            variations[query_pos] = {'query': q.upper()}
                        else:
                            if scientific_name not in variations[query_pos]:
                                variations[query_pos][scientific_name] = h.upper() if h != '-' else '-'
                        query_pos += 1

        # Create DataFrame
        df = pd.DataFrame.from_dict(variations, orient='index')
        df.index.name = 'Query_Position'
        df = df.sort_index()

        # Move query column to first position
        cols = ['query'] + [col for col in df.columns if col != 'query']
        df = df[cols]

        # Filter for polymorphic sites
        variant_cols = df.columns[1:]
        if len(variant_cols) < 2:
            return {
                "status": "success",
                "data": {
                    "query_title": query_title,
                    "query_length": query_length,
                    "message": "Not enough variant data for analysis."
                }
            }

        multiple_var_mask = df[variant_cols].apply(lambda x: len(set(x)) > 1, axis=1)
        df_filtered = df[multiple_var_mask].dropna()

        if df_filtered.empty or len(df_filtered.columns) <= 1:
            return {
                "status": "success",
                "data": {
                    "query_title": query_title,
                    "query_length": query_length,
                    "message": "Not enough variant data after filtering."
                }
            }

        df_var = df_filtered.copy()

        # Generate FASTA content
        fasta_content = generate_fasta_from_dataframe(df_var)

        # Initialize results
        treefile_content = None
        tree_visualization = None

        # Run IQ-TREE if requested
        if run_iqtree and IQTREE_PATH.exists():
            print("Running IQ-TREE phylogenetic analysis...")

            # Use absolute paths for IQ-TREE
            work_dir_abs = Path(work_dir).resolve()
            fasta_file_path = work_dir_abs / "genotype_alignment.fasta"

            if 'query' in df_var.columns:
                genotype_df_var = df_var.drop(["query"], axis=1)
            else:
                genotype_df_var = df_var

            if not genotype_df_var.empty:
                # Write alignment FASTA
                with open(fasta_file_path, 'w') as f:
                    for c in genotype_df_var.columns:
                        seq = ''.join(genotype_df_var[c].astype(str).values)
                        f.write(f">{c}\n{seq}\n")

                treefile_path = str(fasta_file_path) + ".treefile"

                try:
                    iqtree_cmd = [
                        str(IQTREE_PATH), "-s", str(fasta_file_path),
                        "-B", "1000", "-redo", "-nt", "AUTO", "-safe",
                        "-m", "MFP", "-mrate", "G,I,I+G,R"
                    ]

                    result = subprocess.run(
                        iqtree_cmd, check=True, cwd=str(work_dir_abs),
                        capture_output=True, text=True, timeout=1800
                    )

                    if os.path.exists(treefile_path):
                        with open(treefile_path, 'r') as f:
                            treefile_content = f.read()
                        print("IQ-TREE analysis completed successfully")

                        # Generate tree visualization if available
                        if VISUALIZATION_AVAILABLE:
                            tree_visualization = visualize_tree(
                                treefile_content, query_names, output_dir
                            )
                    else:
                        print(f"Warning: Treefile not found at {treefile_path}")

                except subprocess.TimeoutExpired:
                    print("IQ-TREE timed out after 30 minutes")
                    treefile_content = "IQ-TREE analysis timed out"
                except subprocess.CalledProcessError as e:
                    print(f"IQ-TREE error: {e.stderr}")
                    treefile_content = f"IQ-TREE error: {e.stderr}"

        # Generate clustermap if visualization available
        clustermap_path = None
        if VISUALIZATION_AVAILABLE:
            df_dist_input = df_var.drop(["query"], axis=1) if 'query' in df_var.columns else df_var

            if not df_dist_input.empty and df_dist_input.shape[1] >= 2:
                df_dist = compute_pairwise_distances(df_dist_input.T)

                num_samples = len(df_dist.columns)
                clustermap_height = max(10, num_samples * 0.4)
                clustermap_width = max(10, num_samples * 0.4)

                try:
                    clustermap = sb.clustermap(
                        df_dist, method='complete',
                        figsize=(clustermap_width, clustermap_height),
                        tree_kws=dict(linewidths=1.5, colors=(0.2, 0.2, 0.4))
                    )

                    # Highlight query samples
                    for label in clustermap.ax_heatmap.get_xticklabels():
                        if label.get_text() in query_names:
                            label.set_color('red')
                            label.set_weight('bold')

                    for label in clustermap.ax_heatmap.get_yticklabels():
                        if label.get_text() in query_names:
                            label.set_color('red')
                            label.set_weight('bold')

                    clustermap_path = output_dir / "clustermap.png"
                    clustermap.savefig(clustermap_path, format='png', bbox_inches='tight', dpi=150)
                    plt.close()
                    print(f"Clustermap saved to {clustermap_path}")

                except Exception as e:
                    print(f"Error generating clustermap: {e}")

        # Save results
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            # Save variant matrix
            variant_csv_path = output_dir / "variant_matrix.csv"
            df_var.to_csv(variant_csv_path)
            print(f"Variant matrix saved to {variant_csv_path}")

            # Save alignment FASTA
            alignment_fasta_path = output_dir / "alignment.fasta"
            with open(alignment_fasta_path, 'w') as f:
                f.write(fasta_content)
            print(f"Alignment FASTA saved to {alignment_fasta_path}")

            # Save tree file
            if treefile_content and not treefile_content.startswith("IQ-TREE"):
                tree_path = output_dir / "phylogenetic_tree.nwk"
                with open(tree_path, 'w') as f:
                    f.write(treefile_content)
                print(f"Phylogenetic tree saved to {tree_path}")

        return {
            "status": "success",
            "data": {
                "query_title": query_title,
                "query_length": query_length,
                "num_variants": len(df_var),
                "num_samples": len(df_var.columns) - 1,
                "treefile": treefile_content,
                "tree_visualization": tree_visualization,
                "fasta_content": fasta_content
            }
        }

    except Exception as e:
        import traceback
        print(f"Error processing BLAST results: {e}")
        print(traceback.format_exc())
        return {
            "status": "error",
            "error": str(e)
        }


def visualize_tree(newick_string, query_names, output_dir):
    """Generate tree visualization using BioPython."""
    if not VISUALIZATION_AVAILABLE:
        return None

    if not newick_string or not newick_string.strip():
        return None

    try:
        handle = io.StringIO(newick_string)
        tree = Phylo.read(handle, "newick")

        num_samples = len(tree.get_terminals())
        height = max(15, num_samples * 0.4)

        fig = plt.figure(figsize=(12, height), dpi=150)
        ax = fig.add_subplot(1, 1, 1)

        # Try to root with outgroup
        outgroup_keywords = ['Cassytha', 'Persea', 'Ipomoea', 'Coffea', 'Nicotiana']

        for keyword in outgroup_keywords:
            matching_terminals = [l for l in tree.get_terminals()
                                if l.name and keyword in l.name]
            if matching_terminals:
                try:
                    if len(matching_terminals) == 1:
                        tree.root_with_outgroup(matching_terminals[0])
                    else:
                        outgroup_mrca = tree.common_ancestor(matching_terminals)
                        tree.root_with_outgroup(outgroup_mrca)
                    print(f"Rooted tree with outgroup: {keyword}")
                    break
                except:
                    continue

        # Define label function
        def label_func(clade):
            if clade.confidence is not None:
                bootstrap_value = int(round(float(clade.confidence)))
                if bootstrap_value >= 70:
                    return str(bootstrap_value)
                return ""
            elif clade.is_terminal():
                return clade.name if clade.name else "Unnamed"
            return ""

        # Draw tree
        Phylo.draw(tree, axes=ax, label_func=label_func, do_show=False)

        # Highlight query sequences
        for label in ax.texts:
            label_text = label.get_text().strip()
            if label_text in query_names:
                label.set_color('red')
                label.set_weight('bold')

        plt.xlabel("Branch length")
        plt.ylabel("")
        plt.tight_layout()

        # Save figure
        tree_img_path = output_dir / "phylogenetic_tree.png"
        plt.savefig(tree_img_path, format='png', bbox_inches='tight', dpi=150)
        plt.close()
        print(f"Tree visualization saved to {tree_img_path}")

        return str(tree_img_path)

    except Exception as e:
        print(f"Error visualizing tree: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="CUCP BLAST Analysis Tool - Cuscuta Comparative Phylogenomics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic BLAST analysis
  python cucp_blast.py -r reference.fa -q query.fa -o results/

  # With phylogenetic tree construction
  python cucp_blast.py -r reference.fa -q query.fa -o results/ --run-iqtree

  # Specify sample name
  python cucp_blast.py -r reference.fa -q query.fa -o results/ -n "Cuscuta_sample1"
        """
    )

    parser.add_argument(
        '-r', '--reference',
        required=True,
        help='Path to reference FASTA file'
    )

    parser.add_argument(
        '-q', '--query',
        required=True,
        help='Path to query FASTA file'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory for results'
    )

    parser.add_argument(
        '-n', '--name',
        default='query_sample',
        help='Sample name (default: query_sample)'
    )

    parser.add_argument(
        '--run-iqtree',
        action='store_true',
        help='Run IQ-TREE for phylogenetic tree construction'
    )

    parser.add_argument(
        '-w', '--work-dir',
        default=None,
        help='Working directory for intermediate files (default: output/work)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.reference):
        print(f"Error: Reference file not found: {args.reference}")
        sys.exit(1)

    if not os.path.exists(args.query):
        print(f"Error: Query file not found: {args.query}")
        sys.exit(1)

    # Setup paths
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    work_dir = Path(args.work_dir) if args.work_dir else output_dir / "work"
    work_dir.mkdir(parents=True, exist_ok=True)

    # Check for binaries
    setup_permissions()

    if not BLAST_PATH or not BLAST_PATH.exists():
        print("Error: blastn not found in PATH")
        print("\nInstall via conda:")
        print("  conda install -c bioconda blast")
        sys.exit(1)

    if not MAKEBLASTDB_PATH or not MAKEBLASTDB_PATH.exists():
        print("Error: makeblastdb not found in PATH")
        print("\nInstall via conda:")
        print("  conda install -c bioconda blast")
        sys.exit(1)

    if args.run_iqtree and (not IQTREE_PATH or not IQTREE_PATH.exists()):
        print("Warning: iqtree not found in PATH")
        print("Install via conda: conda install -c bioconda iqtree")
        print("Tree construction will be skipped")
        args.run_iqtree = False

    print("=" * 60)
    print("CUCP BLAST Analysis Tool")
    print("=" * 60)
    print(f"Reference: {args.reference}")
    print(f"Query: {args.query}")
    print(f"Output: {output_dir}")
    print(f"Sample name: {args.name}")
    print(f"Run IQ-TREE: {args.run_iqtree}")
    print("=" * 60)

    start_time = datetime.now()

    # Run BLAST
    print("\n[Step 1/2] Running BLAST analysis...")
    blast_results = run_blast(
        args.query,
        args.reference,
        work_dir,
        args.name
    )

    if blast_results["status"] != "success":
        print(f"Error: BLAST failed - {blast_results.get('error', 'Unknown error')}")
        sys.exit(1)

    # Process results
    print("\n[Step 2/2] Processing BLAST results...")
    processed = process_blast_results(
        blast_results["blast_results"],
        blast_results["base_sequence"],
        blast_results["query_names"],
        run_iqtree=args.run_iqtree,
        work_dir=work_dir,
        output_dir=output_dir
    )

    if processed["status"] != "success":
        print(f"Error: Processing failed - {processed.get('error', 'Unknown error')}")
        sys.exit(1)

    # Summary
    elapsed = datetime.now() - start_time
    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    print(f"Time elapsed: {elapsed}")

    if "data" in processed:
        data = processed["data"]
        print(f"Query length: {data.get('query_length', 'N/A')}")
        print(f"Variant positions: {data.get('num_variants', 'N/A')}")
        print(f"Samples analyzed: {data.get('num_samples', 'N/A')}")

    print(f"\nOutput files saved to: {output_dir}")
    print("  - variant_matrix.csv")
    print("  - alignment.fasta")
    if args.run_iqtree:
        print("  - phylogenetic_tree.nwk")
        if VISUALIZATION_AVAILABLE:
            print("  - phylogenetic_tree.png")
            print("  - clustermap.png")


if __name__ == "__main__":
    main()
