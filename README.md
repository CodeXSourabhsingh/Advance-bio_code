# 🧬 Advance-bio_code

A Python-based bioinformatics toolkit for sequence analysis, protein structure assessment, clinical cohort filtering, and viral mutation tracking.

This repository brings together multiple interactive Streamlit applications for research, education, and exploratory analysis across molecular biology, structural biology, clinical data analysis, and drug discovery.

[![Live Demo](https://img.shields.io/badge/🚀%20Live%20Demo-CRISPR%20Gene%20Editing%20Simulator-2ea44f?style=for-the-badge)](https://advance-biocode-whdkzngbtm3ojrz4du4kgc.streamlit.app/)
[![PRISM Demo](https://img.shields.io/badge/📊%20PRISM%20Demo-Advanced%20Clinical%20Pipeline-2ea44f?style=for-the-badge)](https://advance-biocode-95t2yzvng6j7hqordph6bc.streamlit.app/)
[![COVID-19 Demo](https://img.shields.io/badge/🦠%20COVID--19%20Demo-Variant%20Mutation%20Tracker-2ea44f?style=for-the-badge)](https://advance-biocode-2dvc5hdwriquhcezve3wyz.streamlit.app/)
[![Protein Mutation Demo](https://img.shields.io/badge/🧫%20Protein%20Mutation%20Demo-Structure%20Analyzer-2ea44f?style=for-the-badge)](https://advance-biocode-dmogu9dlhr6l5kjm8peqfz.streamlit.app/)
[![Drug-Target Demo](https://img.shields.io/badge/💊%20Drug--Target%20Demo-Interaction%20Explorer-2ea44f?style=for-the-badge)](https://advance-biocode-ewjbtfap7ensteshckapp8.streamlit.app/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.13](https://img.shields.io/badge/Python-3.13-blue.svg)](https://www.python.org/)

---

## Overview

Advance-bio_code is designed to support several connected bioinformatics workflows:

- CRISPR guide design and cleavage simulation
- Clinical cohort filtering and trial-fit analysis
- SARS-CoV-2 variant tracking and mutation scoring
- Protein mutation impact assessment
- Drug-target interaction exploration and prioritization

The project is intended for research and educational exploration, not for clinical decision-making.

---

## Applications

### 1. CRISPR Gene Editing Simulator
Purpose: design and evaluate CRISPR-Cas9 editing workflows.

Features:
- Fetches gene sequences from public databases
- Scans DNA strands for PAM motifs
- Simulates cleavage positions and mismatch patterns
- Estimates editing outcomes and specificity
- Exports guide design summaries

Entry point: `CRISPR.py`

### 2. PRISM Clinical Pipeline
Purpose: filter and analyze large clinical or cohort datasets.

Features:
- Filters records by disease, mutation, age, stage, and biomarker conditions
- Identifies candidate patient cohorts
- Supports large-scale data exploration
- Produces summary reports for trial-fit and biomarker investigations

Entry point: `PRISM.py`

### 3. COVID-19 Variant Mutation Tracker
Purpose: monitor SARS-CoV-2 mutations and variant patterns.

Features:
- Retrieves viral sequence data
- Detects substitutions, insertions, and deletions
- Scores mutation risk and potential immune-escape impact
- Tracks variant patterns over time

Entry point: `COVID_19_Tracker.py`

### 4. Protein Structure Analyzer
Purpose: evaluate how sequence variants may affect protein function or structure.

Features:
- Retrieves protein structures from PDB
- Maps mutation sites in 3D space
- Estimates physicochemical changes
- Highlights potential structural destabilization
- Visualizes mutation hotspots

Entry point: `Protein_mutation_analyzer.py`

### 5. Drug-Target Interaction Explorer
Purpose: explore candidate drug-target relationships.

Features:
- Retrieves compound and target metadata
- Scores interaction potential
- Helps prioritize binding hypotheses
- Supports exploratory therapeutic analysis

Entry point: `Drug_Target_Explorer.py`

---

## Workflow

The applications are structured as an end-to-end research workflow:

```text
CRISPR Simulator
   ↓
Protein Structure Analyzer
   ↓
PRISM Clinical Pipeline
   ↓
Drug-Target Explorer
   ↓
COVID-19 Variant Tracker
   ↓
Research insights and audit logging
```

This workflow connects:
- sequence editing and mutation analysis
- protein structure interpretation
- cohort-level patient filtering
- therapeutic prioritization
- pathogen evolution tracking

---

## Live Demos

| Tool | Demo |
|------|------|
| 🧬 CRISPR Gene Editing Simulator | [Launch Demo](https://advance-biocode-whdkzngbtm3ojrz4du4kgc.streamlit.app/) |
| 🧪 PRISM Clinical Pipeline | [Launch Demo](https://advance-biocode-95t2yzvng6j7hqordph6bc.streamlit.app/) |
| 🦠 COVID-19 Variant Tracker | [Launch Demo](https://advance-biocode-2dvc5hdwriquhcezve3wyz.streamlit.app/) |
| 🧫 Protein Mutation Analyzer | [Launch Demo](https://advance-biocode-dmogu9dlhr6l5kjm8peqfz.streamlit.app/) |
| 💊 Drug-Target Interaction Explorer | [Launch Demo](https://advance-biocode-ewjbtfap7ensteshckapp8.streamlit.app/) |

> Note: cloud-hosted Streamlit apps may run without full local features such as persistent MySQL logging or larger-scale dataset processing.

---

## Tech Stack

| Technology | Purpose |
|-----------|---------|
| Python 3.13+ | Core application runtime |
| Streamlit | Interactive web dashboards |
| Pandas | Data processing and tabular analysis |
| NumPy | Numerical computation |
| Matplotlib | Visualization |
| Biopython | Sequence and structure analysis |
| Requests | API access |
| MySQL Connector | Optional database logging |
| Py3Dmol | 3D molecular visualization |

External data sources:
- NCBI Entrez
- Protein Data Bank (PDB)
- ChEMBL / PubChem
- GISAID / NCBI viral data

---

## Prerequisites

### System requirements
- Python 3.13 or newer
- Windows, macOS, or Linux
- 4 GB RAM minimum; 8–16 GB recommended
- 500 MB+ free storage
- No GPU required

### Software dependencies
- Git
- pip
- MySQL Server 8.0+ (optional, for local logging)

### API access
- NCBI Entrez account/email for public sequence queries
- MySQL credentials if local database logging is enabled

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/CodeXSourabhsingh/Advance-bio_code.git
cd Advance-bio_code
```

### 2. Create a virtual environment

```bash
# Windows
python -m venv venv
venv\Scripts\activate

# macOS / Linux
python3 -m venv venv
source venv/bin/activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

If you want to install manually:

```bash
pip install biopython pandas numpy matplotlib streamlit mysql-connector-python requests py3dmol
```

### 4. Create configuration file

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

USE_CLOUD_DEPLOYMENT = False
```

Important:
- Do not commit real credentials
- Add `config.py` to `.gitignore`

Example:

```bash
echo "config.py" >> .gitignore
```

### 5. Verify installation

```bash
python -c "import streamlit, biopython, pandas, numpy; print('Dependencies loaded successfully')"
```

### 6. Run the applications

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

Each app launches locally on:

```text
http://localhost:8501
```

---

## Project Structure

```text
Advance-bio_code/
├── README.md
├── LICENSE
├── requirements.txt
├── config.py
├── .gitignore
│
├── CRISPR.py
├── PRISM.py
├── COVID_19_Tracker.py
├── Protein_mutation_analyzer.py
├── Drug_Target_Explorer.py
│
├── utils/
│   ├── ncbi_fetcher.py
│   ├── sequence_analysis.py
│   ├── protein_structure.py
│   ├── pam_scanner.py
│   ├── clinical_filters.py
│   ├── drug_target_scorer.py
│   └── mysql_logger.py
│
├── data/
│   ├── reference_genomes/
│   ├── mutation_databases/
│   ├── pdb_cache/
│   └── sample_data/
│
└── docs/
    ├── ARCHITECTURE.md
    ├── API_GUIDE.md
    └── CONTRIBUTING.md
```

---

## Typical Use Cases

- Investigate CRISPR guide candidates for a target gene
- Compare mutation patterns across a protein sequence
- Filter clinical records by disease and biomarker criteria
- Rank potential drug-target candidates
- Monitor and assess emerging SARS-CoV-2 variants

---

## Benchmarks and Performance Notes

This project is intended for exploratory and research-scale workloads. Local performance depends on machine configuration, dataset size, and API latency.

Example local-scale benchmark areas:
- large cohort filtering
- sequence scanning for motif detection
- mutation scoring
- protein feature evaluation

These use standard CPU-based processing, so they are suitable for laptops and moderate workstation environments without requiring GPU acceleration.

> Benchmark numbers should be interpreted as environment-specific estimates rather than guaranteed production metrics.

---

## Deployment Notes

### Local development
For full functionality, run locally with:
- local data access
- MySQL logging enabled
- larger dataset processing

### Streamlit Cloud
The live cloud deployment is useful for demos and lightweight usage, but may have limitations such as:
- limited memory compared to local machines
- disabled or restricted MySQL connectivity
- lower dataset capacity

### Production use
For larger-scale or production workflows:
- run locally or in a managed environment
- configure MySQL properly
- monitor API rate limits
- validate results against reliable domain datasets before use in decision-making

---

## Troubleshooting

### ModuleNotFoundError: No module named 'streamlit'
```bash
pip install --upgrade streamlit
```

### NCBI API rate limiting
Add your API key in `config.py`:

```python
ENTREZ_API_KEY = "your_api_key_here"
```

### MySQL connection refused
Verify the service is running:

```bash
# Windows
net start MySQL80

# macOS
brew services start mysql

# Linux
sudo systemctl start mysql
```

Or disable database logging in config if you do not need it.

### Streamlit port already in use
```bash
streamlit run CRISPR.py --server.port 8502
```

### Memory or sequence length issues
Reduce dataset size in the app UI or by adjusting local limits in the code.

---

## Contributing

Contributions are welcome.

### Suggested workflow
1. Fork the repository
2. Create a new branch
3. Make your changes
4. Test locally
5. Commit with a clear message
6. Open a pull request

### Code standards
- Follow PEP 8
- Add docstrings where appropriate
- Keep functions modular and readable
- Update documentation when behavior changes

---

## Notes & Disclaimer

⚠️ Advance-bio_code is for research and educational exploration only.

It should not be used as a substitute for:
- professional clinical diagnosis
- medical advice
- therapeutic decision-making
- regulatory or patient-care planning

All computational outputs should be validated with experimental data and domain expertise.

---

## Author

**Sourabh Singh**

- GitHub: [@CodeXSourabhsingh](https://github.com/CodeXSourabhsingh)
- LinkedIn: [Sourabh Singh](https://www.linkedin.com/in/sourabh-singh-7b1249434/)

---

## License

This project is licensed under the [MIT License](LICENSE).

You are free to:
- use commercially
- modify and distribute
- use privately

You must:
- include the license and copyright notice
- state changes made

---

## Acknowledgments

- NCBI GenBank and Entrez API
- Protein Data Bank (PDB)
- ChEMBL / PubChem
- GISAID and public viral sequence resources
- Streamlit community
- Biopython contributors

---
