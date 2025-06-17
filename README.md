Protein Vocabulary Reduction for Gene2Vec Embeddings

Implementation of vocabulary reduction strategy for Gene2Vec protein language models, achieving 5.6-fold compression while maintaining 28.68% sequence diversity coverage.

## Overview

This repository contains the analysis pipeline for reducing Gene2Vec vocabulary from ~560,000 to 99,571 gene families through multi-stage clustering and abundance-based selection. The method identifies the most prevalent protein families that capture substantial genomic diversity with minimal computational overhead.

## Key Results

- **Vocabulary compression**: 5.6-fold (560,000 → 99,571 gene families)
- **Diversity coverage**: 28.68% with only 0.38% of clusters
- **Efficiency**: 75-fold enrichment in sequences per selected cluster
- **Statistical validation**: Spearman ρ = -0.804 between cluster size and identity

## Installation

```bash
# Clone repository
git clone https://github.com/yourusername/protein-vocabulary-reduction.git
cd protein-vocabulary-reduction

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
External Dependencies


CD-HIT v4.8.1
MMseqs2 v13.45111
DIAMOND v2.0.15

Data Availability
Due to file size constraints, the following large files are not included in this repository:
Required Input Files - those files are on the server:

CD-HIT output: *.clstr file (~GB size)
MMseqs2 results: mmseqs_s1_c0.5_cluster.tsv (~GB size)
Protein FASTA: *.faa file (~GB size)

Small Data Files (included):

data/advanced_cluster_stats.csv - Final cluster statistics (11.2 MB)
data/summary_statistics.txt - Summary metrics
data/sample_input/ - Small example files for testing

# Run vocabulary selection
python scripts/vocabulary_selection.py \
    --cdhit-clusters $CDHIT_CLUSTERS \
    --fasta $PROTEIN_FASTA \
    --mmseqs-results $MMSEQS_RESULTS \
    --output-dir results/

3. Reproduce Figures from Provided Statistics
bash# Generate all figures using the included cluster statistics
python scripts/generate_figures.py \
    --stats-file data/advanced_cluster_stats.csv \
    --output-dir figures/

Repository Structure
├── scripts/
│   ├── vocabulary_selection.py      # Main pipeline
│   ├── cluster_statistics.py        # Statistical analysis
│   ├── process_diamond_output.py    # DIAMOND processing
│   └── generate_figures.py          # Create publication figures
│
├── data/
│   ├── advanced_cluster_stats.csv   # Final statistics (included)
│   ├── summary_statistics.txt       # Summary metrics (included)
│   ├── sample_input/               # Small example files
│   └── README.md                   # Instructions for obtaining full data
│
├── figures/                         # Publication figures

Expected Input Formats
CD-HIT Cluster File (.clstr)
>Cluster 0
0       2799aa, >protein1... *
1       2554aa, >protein2... at 80.05%
MMseqs2 Output (TSV)
cluster_1    protein1
cluster_1    protein2
Outputs
The pipeline generates:

selected_vocabulary.txt: 99,571 selected gene families
vocabulary_statistics.csv: Detailed metrics
diversity_coverage_report.txt: Coverage analysis


Author: Saar Gerassi
Email: saargerassi@mail.tau.ac.il
Lab: Dr. David Burstein's Lab, Tel Aviv University
