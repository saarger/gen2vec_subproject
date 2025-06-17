# 🧬 Protein Vocabulary Reduction for Gene2Vec Embeddings

Implementation of a vocabulary reduction strategy for Gene2Vec protein language models, achieving **5.6× compression** while maintaining **28.68% sequence diversity coverage**.

---

## 📌 Overview

This repository contains the analysis pipeline for reducing the Gene2Vec vocabulary from ~560,000 to **99,571 gene families** through multi-stage clustering and abundance-based selection.  
The method prioritizes highly prevalent protein families to capture genomic diversity with minimal computational overhead.

---

## ✅ Key Results

- **Vocabulary compression**: 5.6× (560,000 → 99,571 gene families)
- **Diversity coverage**: 28.68% with only 0.38% of clusters
- **Efficiency**: 75-fold enrichment in sequences per selected cluster
- **Validation**: Spearman ρ = **-0.804** between cluster size and identity

---

## ⚙️ Installation

```bash
# Clone repository
git clone https://github.com/yourusername/protein-vocabulary-reduction.git
cd protein-vocabulary-reduction

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

---

## 🔗 External Dependencies

- **[CD-HIT v4.8.1](https://github.com/weizhongli/cdhit)**
- **[MMseqs2 v13.45111](https://github.com/soedinglab/MMseqs2)**
- **[DIAMOND v2.0.15](https://github.com/bbuchfink/diamond)**

---

## 📂 Data Availability

### Not Included (stored externally due to size):
- `*.clstr` file (CD-HIT clusters)
- `mmseqs_s1_c0.5_cluster.tsv` (MMseqs2 results)
- `*.faa` file (protein FASTA)

### Included in repo:
- `data/advanced_cluster_stats.csv` – Final cluster statistics (11.2 MB)
- `data/summary_statistics.txt` – Summary metrics
- `data/sample_input/` – Example input files

---

## 🚀 Usage

### 1. Run Vocabulary Selection
```bash
python scripts/vocabulary_selection.py \
    --cdhit-clusters $CDHIT_CLUSTERS \
    --fasta $PROTEIN_FASTA \
    --mmseqs-results $MMSEQS_RESULTS \
    --output-dir results/
```

### 2. Generate Figures
```bash
python scripts/generate_figures.py \
    --stats-file data/advanced_cluster_stats.csv \
    --output-dir figures/
```

---

## 🧪 MMseqs2 Parameter Exploration

| Sensitivity (s) | Coverage (c) | Total Clusters |
|------------------|--------------|----------------|
| 1.0              | 0.50         | 26,247,182     |
| 1.0              | 0.65         | 28,399,675     |
| 1.0              | 0.80         | 32,015,052     |
| 4.0              | 0.50         | 21,093,901     |
| 4.0              | 0.65         | 23,407,319     |
| 4.0              | 0.80         | 27,546,184     |
| 7.5              | 0.50         | 19,660,851     |
| 7.5              | 0.65         | 22,080,886     |

---

## 📁 Repository Structure

```
├── scripts/
│   ├── vocabulary_selection.py      # Main pipeline
│   ├── cluster_statistics.py        # Statistical analysis
│   ├── process_diamond_output.py    # DIAMOND processing
│   └── generate_figures.py          # Create publication figures
│
├── data/
│   ├── advanced_cluster_stats.csv   # Final statistics (included)
│   ├── summary_statistics.txt       # Summary metrics (included)
│   ├── sample_input/                # Small example files
│   └── README.md                    # Instructions for obtaining full data
│
├── figures/                         # Publication figures
```

---

## 📥 Expected Input Formats

### 📄 CD-HIT `.clstr` Example:
```
>Cluster 0
0       2799aa, >protein1... *
1       2554aa, >protein2... at 80.05%
```

### 📄 MMseqs2 `.tsv` Example:
```
cluster_1    protein1
cluster_1    protein2
```

---

## 📤 Outputs

- `selected_vocabulary.txt`: List of 99,571 selected gene families
- `vocabulary_statistics.csv`: Cluster-level metrics
- `diversity_coverage_report.txt`: Global diversity summary

---

## 👨‍🔬 Author

**Saar Gerassi**  
Email: [saargerassi@mail.tau.ac.il](mailto:saargerassi@mail.tau.ac.il)  
Lab: Dr. David Burstein's Lab, Tel Aviv University

---