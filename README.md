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

### Single Sample Analysis

```bash
# With explicit sample name
python cucp_blast.py -r reference.fa -q query.fa -o results/ -n "Cuscuta_sample1"

# Name derived from filename (will use "query" as sample name)
python cucp_blast.py -r reference.fa -q query.fa -o results/
```

### Multiple Samples Analysis

```bash
# Multiple samples with explicit names
python cucp_blast.py -r reference.fa \
  -q sample1.fa sample2.fa sample3.fa \
  -o results/ \
  -n "Sample_A" "Sample_B" "Sample_C"

# Multiple samples with names derived from filenames
python cucp_blast.py -r reference.fa \
  -q sample1.fa sample2.fa sample3.fa \
  -o results/
# Names will be: sample1, sample2, sample3
```

### With Phylogenetic Tree Construction

```bash
python cucp_blast.py -r reference.fa -q query.fa -o results/ --run-iqtree
```

### Full Options

```bash
python cucp_blast.py \
  --reference reference.fa \
  --query sample1.fa sample2.fa \
  --output results/ \
  --name "Sample_A" "Sample_B" \
  --run-iqtree \
  --work-dir /tmp/work
```

## Arguments

| Argument | Short | Required | Description |
|----------|-------|----------|-------------|
| `--reference` | `-r` | Yes | Path to reference FASTA file (first sequence used as anchoring base) |
| `--query` | `-q` | Yes | Path to query FASTA file(s). Multiple files = multiple samples |
| `--output` | `-o` | Yes | Output directory for results |
| `--name` | `-n` | No | Sample name(s). Must match number of query files. If not provided, names are derived from filenames |
| `--run-iqtree` | | No | Run IQ-TREE for phylogenetic tree construction |
| `--work-dir` | `-w` | No | Working directory for intermediate files |

### Multiple Contigs/Scaffolds

Each query file can contain multiple sequences (contigs/scaffolds). All sequences within a single file are treated as belonging to the **same sample** and will be merged into one column in the output matrix.

If your sample has multiple scaffold files, concatenate them first:
```bash
cat scaffold1.fa scaffold2.fa scaffold3.fa > my_sample.fa
python cucp_blast.py -r reference.fa -q my_sample.fa -o results/ -n "My_Sample"
```

## Output Files

The tool generates the following output files:

**Always generated:**
- `variant_matrix.csv` - Matrix of variant positions across samples
- `alignment.fasta` - Multiple sequence alignment in FASTA format
- `clustermap.png` - Distance-based hierarchical clustering heatmap

**With `--run-iqtree`:**
- `phylogenetic_tree.nwk` - Newick format tree file
- `phylogenetic_tree.png` - Tree visualization

## Examples

### Single Sample with CP Genome Reference

```bash
python cucp_blast.py \
  -r ref/ref_250905_nocassytha_noam711639.fa \
  -q my_sample.fa \
  -o output/ \
  -n "Cuscuta_japonica_Korea" \
  --run-iqtree
```

### Multiple Samples

```bash
python cucp_blast.py \
  -r ref/ref_250905_nocassytha_noam711639.fa \
  -q sample_korea.fa sample_japan.fa sample_china.fa \
  -o output/ \
  -n "C_japonica_Korea" "C_japonica_Japan" "C_chinensis_China" \
  --run-iqtree
```

### Batch Analysis (names from filenames)

```bash
# Analyze all .fa files in a directory
python cucp_blast.py \
  -r ref/ref_250905_nocassytha_noam711639.fa \
  -q samples/*.fa \
  -o batch_results/
```

## Reference Databases

Reference databases are located in the `ref/` directory.

### Included Reference

| File | Description |
|------|-------------|
| `ref/ref_250905_nocassytha_noam711639.fa` | CP genome reference (without Cassytha, recommended for most analyses) |

### Additional References (from CUCP web application)

Additional reference databases can be obtained from the main CUCP project (`blast/data/`):

| File | Description |
|------|-------------|
| `ref_250905_noam711639.fa` | CP genome reference (with Cassytha) |
| `ref_trnlf_250904_costea_2025.fa` | trnL-F reference database (Costea 2025 expanded) |

## Preparing a Custom Reference File

The **first sequence** in your reference file is used as the anchoring sequence for genotype position mapping. All other sequences will be aligned against this first sequence.

### Reference File Format

```
>Accession Genus species description
ATGCATGCATGC...
```

Example:
```
>AM711640.1 Cuscuta reflexa complete chloroplast genome
ATGCATGCATGC...
>MN866891.1 Cuscuta australis chloroplast complete genome
ATGCATGCATGC...
```

The tool parses headers and outputs sample names as `Genus_species_(Accession)`:
- `>AM711640.1 Cuscuta reflexa ...` → `Cuscuta_reflexa_(AM711640.1)`

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
