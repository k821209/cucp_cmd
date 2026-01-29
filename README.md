# CUCP BLAST Command-Line Tool

A standalone command-line tool for running BLAST analysis and phylogenetic tree construction for Cuscuta (dodder) comparative phylogenomics.

## Installation

```bash
# Create conda environment with all dependencies
conda create -n cucp python=3.11 -y
conda activate cucp

# Install BLAST+ and IQ-TREE from bioconda
conda install -c bioconda blast iqtree -y

# Install Python dependencies
pip install -r requirements.txt
```

## How It Works

1. **Merge sequences**: Reference and sample sequences are merged into a single FASTA file
2. **Create BLAST database**: A BLAST database is created from the merged sequences
3. **Run BLAST**: Align all sequences against the first sequence in the reference file
4. **Extract variants**: Genotype positions are mapped based on the first sequence of reference.fa. Polymorphic positions are extracted and saved as a variant matrix
5. **Build tree** (optional): IQ-TREE constructs a phylogenetic tree from the variant alignment

## Usage

### Basic BLAST Analysis

```bash
python cucp_blast.py -r reference.fa -q query.fa -o results/
```

### With Phylogenetic Tree Construction

```bash
python cucp_blast.py -r reference.fa -q query.fa -o results/ --run-iqtree
```

### Full Options

```bash
python cucp_blast.py \
  --reference reference.fa \
  --query query.fa \
  --output results/ \
  --name "Cuscuta_sample1" \
  --run-iqtree \
  --work-dir /tmp/work
```

## Arguments

| Argument | Short | Required | Description |
|----------|-------|----------|-------------|
| `--reference` | `-r` | Yes | Path to reference FASTA file (first sequence used as base) |
| `--query` | `-q` | Yes | Path to sample FASTA file (your sequences to analyze) |
| `--output` | `-o` | Yes | Output directory for results |
| `--name` | `-n` | No | Sample name (default: query_sample) |
| `--run-iqtree` | | No | Run IQ-TREE for phylogenetic tree construction |
| `--work-dir` | `-w` | No | Working directory for intermediate files |

## Output Files

The tool generates the following output files:

**Always generated:**
- `variant_matrix.csv` - Matrix of variant positions across samples
- `alignment.fasta` - Multiple sequence alignment in FASTA format
- `clustermap.png` - Distance-based hierarchical clustering heatmap

**With `--run-iqtree`:**
- `phylogenetic_tree.nwk` - Newick format tree file
- `phylogenetic_tree.png` - Tree visualization

## Example

Using reference sequences from the CUCP database:

```bash
# Using CP genome reference
python cucp_blast.py \
  -r ref_250905_nocassytha_noam711639.fa \
  -q my_sample.fa \
  -o output/ \
  --run-iqtree

# Using trnL-F reference
python cucp_blast.py \
  -r ref_trnlf_250904_costea_2025.fa \
  -q my_trnlf_sample.fa \
  -o output_trnlf/
```

## Reference Databases

Available reference databases:

| File | Description |
|------|-------------|
| `ref_250905_noam711639.fa` | CP genome reference (with Cassytha) |
| `ref_250905_nocassytha_noam711639.fa` | CP genome reference (without Cassytha) |
| `ref_trnlf_250904_costea_2025.fa` | trnL-F reference database |

## Preparing a Custom Reference File

The **first sequence** in your reference file is used as the anchoring sequence for genotype position mapping. All other sequences will be aligned against this first sequence.

### Reference File Format

```
>Anchoring_Sequence_Name description
ATGCATGCATGC...
>Reference_Species_1 description
ATGCATGCATGC...
>Reference_Species_2 description
ATGCATGCATGC...
```

### Guidelines

1. **First sequence = Anchoring sequence**: Place your best quality, complete sequence first. This sequence defines the coordinate system for all variant positions.

2. **Sequence requirements**:
   - All sequences should be the same genomic region (e.g., complete chloroplast genome or trnL-F region)
   - The anchoring sequence should be complete without gaps
   - Use consistent naming format: `>Accession Species_name additional_info`

3. **Example**: Creating a custom reference from NCBI sequences

```bash
# Download sequences and concatenate (anchoring sequence first)
cat anchoring_sequence.fa other_references.fa > my_reference.fa

# Verify the first sequence
head -2 my_reference.fa
```

4. **Recommended anchoring sequences**:
   - For CP genome analysis: Use a well-annotated complete chloroplast genome (e.g., *Cuscuta australis* MN866891.1)
   - For trnL-F analysis: Use a high-quality trnL-F sequence from a closely related species
