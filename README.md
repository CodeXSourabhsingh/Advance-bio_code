# 🧬 Advance-bio_code

> A connected bioinformatics portfolio for sequence analysis, protein structure assessment, clinical cohort filtering, and viral variant tracking.

[![Live Demo](https://img.shields.io/badge/🚀%20Live%20Demo-CRISPR%20Gene%20Editing%20Simulator-2ea44f?style=for-the-badge)](https://advance-biocode-whdkzngbtm3ojrz4du4kgc.streamlit.app/)
[![PRISM Demo](https://img.shields.io/badge/📊%20PRISM%20Demo-Advanced%20Clinical%20Pipeline-2ea44f?style=for-the-badge)](https://advance-biocode-95t2yzvng6j7hqordph6bc.streamlit.app/)
[![COVID-19 Demo](https://img.shields.io/badge/🦠%20COVID--19%20Demo-Variant%20Mutation%20Tracker-2ea44f?style=for-the-badge)](https://advance-biocode-2dvc5hdwriquhcezve3wyz.streamlit.app/)
[![Protein Mutation Demo](https://img.shields.io/badge/🧫%20Protein%20Mutation%20Demo-Structure%20Analyzer-2ea44f?style=for-the-badge)](https://advance-biocode-dmogu9dlhr6l5kjm8peqfz.streamlit.app/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.13](https://img.shields.io/badge/Python-3.13-blue.svg)](https://www.python.org/)

Advance-bio_code is a Python-based bioinformatics toolkit designed for research, learning, and exploration across multiple layers of molecular and clinical analysis. The project brings together five integrated tools spanning genetic engineering, clinical trial analysis, viral surveillance, protein structure assessment, and drug-target interaction prediction.

---

## 📋 Table of Contents

- [Project at a Glance](#project-at-a-glance)
- [The Ecosystem](#-the-ecosystem)
- [Tools & Features](#tools--features)
- [Live Demo](#live-demo)
- [Tech Stack](#tech-stack)
- [Prerequisites](#prerequisites)
- [Quick Start](#quick-start)
- [Project Structure](#project-structure)
- [Performance Benchmarks](#performance-benchmarks)
- [Scaling & Deployment](#scaling--deployment)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Notes & Disclaimer](#notes--disclaimer)
- [Author](#author)
- [License](#license)

---

## Project at a Glance

- **5 integrated interactive applications**
- Sequence analysis using public genomic data
- Mutation and cleavage simulation for CRISPR workflows
- Clinical cohort filtering and trial-fit evaluation
- Viral variant tracking for SARS-CoV-2
- Protein structure impact assessment
- Drug-target binding prediction and scoring
- MySQL-powered audit logging and data persistence
- Built with Python, Streamlit, Pandas, NumPy, and Biopython
- Optimized for standard consumer hardware (tested on mid-range laptop)
- Benchmark: Process 10 million patient records in under 2 minutes

---

## 🌍 The Ecosystem

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    ADVANCE-BIO_CODE BIOINFORMATICS ECOSYSTEM                 │
└─────────────────────────────────────────────────────────────────────────────┘

                           ┌─────────────────────┐
                           │   Data Collection   │
                           │     & Analysis      │
                           └──────────┬──────────┘
                                      │
                    ┌─────────────────┼─────────────────┐
                    │                 │                 │
                    ▼                 ▼                 ▼
            ┌──────────────┐   ┌──────────────┐  ┌──────────────┐
            │    NCBI      │   │  PubChem/    │  │     PDB      │
            │  GenBank     │   │  ChEMBL      │  │  Database    │
            └──────┬───────┘   └──────┬───────┘  └──────┬───────┘
                   │                  │                 │
        ┌──────────┴──────────┬───────┴────────┬────────┴─────────┐
        │                     │                │                  │
        ▼                     ▼                ▼                  ▼
   ┌─────────────┐   ┌──────────────┐  ┌─────────────┐  ┌──────────────┐
   │   CRISPR    │   │    PRISM     │  │  COVID-19   │  │   Protein    │
   │ Gene Editor │   │  Clinical    │  │   Variant   │  │  Structure   │
   │ Simulator   │   │  Pipeline    │  │  Tracker    │  │  Analyzer    │
   └──────┬──────┘   └──────┬───────┘  └──────┬──────┘  └──────┬───────┘
          │                 │                  │               │
          └─────────────────┼──────────────────┼───────────────┘
                            │
                   ┌────────┴────────┐
                   │                 │
                   ▼                 ▼
            ┌────────────────┐  ┌──────────────┐
            │   Drug-Target  │  │    MySQL     │
            │  Interaction   │  │    Logging   │
            │   Explorer     │  │  & Audit     │
            └────────────────┘  └──────────────┘
                    │
                    ▼
        ┌──────────────────────────────┐
        │  Research & Clinical Outputs │
        │  - Mutation Classifications  │
        │  - Trial Fit Predictions     │
        │  - Variant Risk Scoring      │
        │  - Structure Impact Analysis │
        │  - Drug Binding Predictions  │
        └──────────────────────────────┘
```

### **Workflow Integration**

The five tools are designed as a connected research pipeline:

```
1. CRISPR Simulator
   ↓ (Identifies sequence mutations)
   
2. Protein Structure Analyzer
   ↓ (Evaluates structural impact)
   
3. PRISM Clinical Pipeline
   ↓ (Finds relevant patient cohorts)
   
4. Drug-Target Explorer
   ↓ (Predicts drug interactions)
   
5. COVID-19 Tracker
   └─→ (Extends analysis to viral evolution)

   ALL RESULTS → MySQL Logging & Audit Trail
```

This creates a comprehensive research workflow spanning:
- **Molecular Biology** — Sequence editing and analysis
- **Structural Biology** — Protein impact assessment
- **Clinical Relevance** — Patient cohort identification
- **Pharmacology** — Drug-target interactions
- **Viral Surveillance** — Variant monitoring and mutation scoring

---

## Tools & Features

### 🧬 CRISPR Gene Editing Simulator

**Purpose:** Design and simulate CRISPR-Cas9 gene editing workflows

- Fetches real gene sequences from NCBI GenBank
- Scans both DNA strands for PAM motifs (NGG, NGG variants)
- Calculates mismatch patterns and specificity scores
- Simulates Cas9 cleavage and off-target potential
- Classifies likely DNA repair outcomes (NHEJ, HDR)
- Exports guide RNA designs and cleavage predictions

**Entry Point:** `CRISPR.py`

---

### 🧪 PRISM — Advanced Clinical Pipeline

**Purpose:** Filter and analyze large patient cohorts for trial-fit and biomarker discovery

- Fetches public clinical trial metadata and patient cohorts
- Filters by disease, stage, age, genetic mutations, biomarkers
- Identifies statistical outliers and potential trial-fit candidates
- Supports analysis of 500K–10M patient records
- Benchmarked for speed: 10M patients in <2 minutes
- Generates cohort reports with demographic and genetic breakdowns
- MySQL integration for full audit and result persistence

**Entry Point:** `PRISM.py`

---

### 🦠 COVID-19 Variant Mutation Tracker

**Purpose:** Monitor SARS-CoV-2 evolution and assess variant risk

- Fetches latest SARS-CoV-2 sequence data from GISAID/NCBI
- Aligns Spike protein region against Wuhan-Hu-1 reference genome
- Detects substitutions, insertions, and deletions
- Scores mutations for immune escape potential (ACE2 binding, epitope impact)
- Ranks transmission risk using phylogenetic distance and prevalence
- Tracks emerging variants and mutation patterns over time

**Entry Point:** `COVID_19_Tracker.py`

---

### 🧫 Protein Structure Analyzer

**Purpose:** Predict how mutations affect protein structure and function

- Fetches protein structures from the Protein Data Bank (PDB)
- Locates mutation sites on protein chains in 3D space
- Calculates physicochemical property changes (charge, polarity, hydrophobicity)
- Estimates structural destabilization risk (RSA, DSSP properties)
- Classifies mutations: benign, disease-associated, or deleterious
- Interactive 3D visualization of mutation sites
- Comparison against known disease databases

**Entry Point:** `Protein_mutation_analyzer.py`

---

### 💊 Drug-Target Interaction Explorer

**Purpose:** Predict drug-target binding affinity and mechanism

- Fetches drug-target binding data from ChEMBL database
- Classifies interactions by binding strength (IC50, Ki, Kd values)
- Visualizes primary protein targets in 3D using PDB structures
- Supports target filtering by gene, protein family, or pathway
- Scoring framework for therapeutic potential and off-target risk
- Full audit trail logging to MySQL for regulatory compliance
- Enables drug repurposing and polypharmacology discovery

**Entry Point:** `Drug_Target_Explorer.py`

---

## Live Demo

Try the deployed applications here:

| Tool | Live Demo |
|------|-----------|
| 🧬 CRISPR Gene Editing Simulator | [Launch Demo](https://advance-biocode-whdkzngbtm3ojrz4du4kgc.streamlit.app/) |
| 🧪 PRISM Clinical Pipeline | [Launch Demo](https://advance-biocode-95t2yzvng6j7hqordph6bc.streamlit.app/) |
| 🦠 COVID-19 Variant Tracker | [Launch Demo](https://advance-biocode-2dvc5hdwriquhcezve3wyz.streamlit.app/) |
| 🧫 Protein Mutation Analyzer | [Launch Demo](https://advance-biocode-dmogu9dlhr6l5kjm8peqfz.streamlit.app/) |

> **⚠️ Cloud Deployment Note:** The Streamlit Cloud deployment runs without MySQL logging due to networking limitations. For the complete pipeline with persistent logging, run locally (see [Quick Start](#quick-start)).

---

## Tech Stack

| Technology | Version | Role |
|---|---|---|
| **Python** | 3.13+ | Core runtime and computational engine |
| **Streamlit** | Latest | Interactive web dashboards and UI |
| **Pandas** | Latest | Tabular data handling and analysis |
| **NumPy** | Latest | Numerical arrays and vectorized operations |
| **Matplotlib** | Latest | Static plotting and data visualization |
| **Biopython** | 1.80+ | Sequence analysis, alignment, structure parsing |
| **Requests** | Latest | HTTP client for API access |
| **MySQL Connector** | 8.0+ | Database connectivity and logging |
| **Py3Dmol** | Latest | Interactive 3D protein structure rendering |

**External APIs:**
- NCBI Entrez (GenBank, GeneID, Protein)
- PubChem/ChEMBL (Drug-target data)
- Protein Data Bank (PDB structures)
- GISAID/NCBI (SARS-CoV-2 sequences)

---

## Prerequisites

### System Requirements

- **Python:** 3.13 or higher
- **OS:** Windows, macOS, or Linux
- **RAM:** 4 GB minimum (16 GB recommended for PRISM scaling)
- **Storage:** 500 MB for code and dependencies
- **GPU:** Not required (CPU-based computation)

### Software Dependencies

- Git (for cloning the repository)
- pip (Python package manager)
- MySQL Server 8.0+ (optional, for logging)

### API Requirements

- **NCBI Entrez Email:** Free NCBI account (register at [NCBI](https://www.ncbi.nlm.nih.gov/))
- **MySQL Credentials:** If using database logging locally

---

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/CodeXSourabhsingh/Advance-bio_code.git
cd Advance-bio_code
```

### 2. Create Virtual Environment (Recommended)

```bash
# On Windows
python -m venv venv
venv\Scripts\activate

# On macOS/Linux
python3 -m venv venv
source venv/bin/activate
```

### 3. Install Dependencies

```bash
pip install -r requirements.txt
```

Or install manually:

```bash
pip install biopython pandas numpy matplotlib streamlit mysql-connector-python requests py3dmol
```

### 4. Configure the Application

Create a `config.py` file in the project root:

```python
# config.py
MYSQL_HOST = "localhost"
MYSQL_USER = "your_mysql_user"
MYSQL_PASSWORD = "your_mysql_password"
MYSQL_DATABASE = "bioinformatics_db"
MYSQL_PORT = 3306

ENTREZ_EMAIL = "your.email@example.com"
ENTREZ_API_KEY = "optional_ncbi_api_key"

# Streamlit Cloud deployment (set to False locally)
USE_CLOUD_DEPLOYMENT = False
```

**⚠️ Security:** Never commit `config.py` with real credentials. Add it to `.gitignore`:

```bash
echo "config.py" >> .gitignore
```

### 5. Verify Installation

Test that all dependencies are installed:

```bash
python -c "import streamlit, biopython, pandas, numpy; print('✓ All dependencies loaded')"
```

### 6. Run an Application

From the project root:

```bash
# CRISPR Gene Editing Simulator
streamlit run CRISPR.py

# PRISM Clinical Pipeline
streamlit run PRISM.py

# COVID-19 Variant Tracker
streamlit run COVID_19_Tracker.py

# Protein Structure Analyzer
streamlit run Protein_mutation_analyzer.py

# Drug-Target Interaction Explorer
streamlit run Drug_Target_Explorer.py
```

Each app will launch at `http://localhost:8501`.

---

## Project Structure

```
Advance-bio_code/
├── README.md                          # This file
├── LICENSE                            # MIT License
├── requirements.txt                   # Python dependencies
├── config.py                          # Configuration (NOT committed)
├── .gitignore                         # Git ignore rules
│
├── CRISPR.py                          # CRISPR Gene Editing Simulator
├── PRISM.py                           # Clinical Pipeline Tool
├── COVID_19_Tracker.py                # Viral Variant Tracker
├── Protein_mutation_analyzer.py       # Protein Structure Tool
├── Drug_Target_Explorer.py            # Drug-Target Interaction Tool
│
├── utils/
│   ├── ncbi_fetcher.py               # NCBI GenBank queries
│   ├── sequence_analysis.py          # Alignment and mutation detection
│   ├── protein_structure.py          # PDB parsing and analysis
│   ├── pam_scanner.py                # PAM motif scanning
│   ├── clinical_filters.py           # Cohort filtering logic
│   ├── drug_target_scorer.py         # Binding prediction
│   └── mysql_logger.py               # Database logging
│
├── data/
│   ├── reference_genomes/            # Reference sequences (FASTA)
│   ├── mutation_databases/           # Known disease mutations
│   ├── pdb_cache/                    # Cached protein structures
│   └── sample_data/                  # Example datasets
│
└── docs/
    ├── ARCHITECTURE.md               # System design
    ├── API_GUIDE.md                  # API usage examples
    └── CONTRIBUTING.md               # Contribution guidelines
```

---

## Performance Benchmarks

### Test Environment

- **Device:** HP Laptop 15-fr0xxx
- **Processor:** 13th Gen Intel Core i5-13420H (2.10 GHz, 8 cores)
- **RAM:** 16 GB DDR4
- **Storage:** 477 GB SSD
- **GPU:** Intel UHD Graphics (integrated, unused)
- **OS:** Windows 11 (64-bit)

### PRISM Cohort Filtering Performance

| Patients | Time | Memory |
|----------|------|--------|
| 500,000 | 7–8 sec | ~500 MB |
| 1,000,000 | 10–15 sec | ~800 MB |
| 2,000,000 | 20 sec | ~1.2 GB |
| 5,000,000 | 50 sec | ~1.5 GB |
| 10,000,000 | <2 min | ~1.9 GB |

**Key Insight:** The bottleneck is Python object generation, not computational hardware. The pipeline is efficient enough for standard consumer hardware without requiring GPU acceleration or distributed computing.

---

## Scaling & Deployment

### Local Development

For development and testing with full features (including MySQL logging):

```bash
streamlit run PRISM.py
```

To increase dataset limits locally, edit the relevant `max_value` parameter in each tool's configuration.

### Cloud Deployment (Streamlit Cloud)

The project is deployed on Streamlit Cloud with the following limitations:

- **Patient Limit:** 5,000 (vs. 10M locally) due to free-tier memory constraints
- **Logging:** MySQL disabled (networking restrictions)
- **Updates:** Re-deploy after git push

Deploy your own fork:

1. Fork this repository
2. Connect to Streamlit Cloud: https://share.streamlit.io/
3. Deploy from your fork

### Production Deployment

For large-scale production use:

- **Local Execution:** Recommended for processing >1M records
- **Database:** Set up MySQL server and update `config.py`
- **Scaling:** Increase `max_value` parameters as needed
- **Monitoring:** Check MySQL logs for query performance
- **Caching:** Implement Redis for API response caching

---

## Troubleshooting

### Issue: "ModuleNotFoundError: No module named 'streamlit'"

**Solution:**
```bash
pip install --upgrade streamlit
```

### Issue: NCBI API Rate Limiting

**Solution:** Add API key to `config.py`:
```python
ENTREZ_API_KEY = "your_api_key_here"
```

Register free at: https://www.ncbi.nlm.nih.gov/account/

### Issue: MySQL Connection Refused

**Solution:** Verify MySQL is running:
```bash
# Windows
net start MySQL80

# macOS
brew services start mysql

# Linux
sudo systemctl start mysql
```

Or disable MySQL in `config.py`:
```python
USE_MYSQL = False
```

### Issue: Streamlit Port Already in Use

**Solution:** Specify a different port:
```bash
streamlit run CRISPR.py --server.port 8502
```

### Issue: "Sequence Too Long" or Memory Errors

**Solution:** Reduce dataset size in the app UI or locally in code:
```python
# In PRISM.py, change:
max_value=10000000  # to:
max_value=1000000   # for 1 million
```

---

## Contributing

We welcome contributions! Please follow these steps:

1. **Fork the repository**
2. **Create a feature branch:** `git checkout -b feature/your-feature`
3. **Make changes and test locally**
4. **Commit:** `git commit -m "Add: your feature description"`
5. **Push:** `git push origin feature/your-feature`
6. **Open a Pull Request** with a clear description

### Code Standards

- Follow PEP 8 style guide
- Add docstrings to functions
- Test all changes locally before submitting
- Update documentation as needed

---

## Notes & Disclaimer

⚠️ **Important:** Results and predictions from Advance-bio_code are intended for **research and educational exploration only**. They should **NOT** be used as a substitute for:

- Professional clinical diagnosis
- Medical advice from qualified healthcare providers
- Regulatory or therapeutic decision-making
- Patient treatment planning

Always validate computational predictions with experimental data and consult appropriate domain experts.

---

## Author

**Sourabh Singh**

- **GitHub:** [@CodeXSourabhsingh](https://github.com/CodeXSourabhsingh)
- **LinkedIn:** [Sourabh Singh](https://www.linkedin.com/in/sourabh-singh-7b1249434/)
- **Email:** sourabh@example.com (update as needed)

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

You are free to:
- ✅ Use commercially
- ✅ Modify and distribute
- ✅ Use privately

You must:
- 📋 Include license and copyright notice
- 📋 State changes

---

## Acknowledgments

- NCBI GenBank and Entrez API
- Protein Data Bank (PDB) and structure tools
- ChEMBL for drug-target data
- GISAID for SARS-CoV-2 sequences
- Streamlit for rapid web app development
- Biopython community

---

**Last Updated:** September 2024  
**Status:** Active Development  
**Version:** 1.0.0
