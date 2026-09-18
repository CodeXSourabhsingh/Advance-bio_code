# 🧬 Advance-bio_code

> A connected bioinformatics portfolio for sequence analysis, protein structure assessment, clinical cohort filtering, and viral variant tracking.

[![Live Demo](https://img.shields.io/badge/🚀%20Live%20Demo-CRISPR%20Gene%20Editing%20Simulator-2ea44f?style=for-the-badge)](https://advance-biocode-whdkzngbtm3ojrz4du4kgc.streamlit.app/)
[![PRISM Demo](https://img.shields.io/badge/📊%20PRISM%20Demo-Advanced%20Clinical%20Pipeline-2ea44f?style=for-the-badge)](https://advance-biocode-95t2yzvng6j7hqordph6bc.streamlit.app/)

**Live demo:** [CRISPR Gene Editing Simulator](https://advance-biocode-whdkzngbtm3ojrz4du4kgc.streamlit.app/) and [PRISM — Advanced Clinical Pipeline](https://advance-biocode-95t2yzvng6j7hqordph6bc.streamlit.app/)

"Cloud demo is capped at 5,000 patients due to Streamlit Cloud free-tier resource limits. Local installation verified up to 10,000,000 patients in under 2 minutes on an i5-13420H + 16 GB RAM setup."

## Overview

Advance-bio_code brings together four interactive tools that connect genomic, structural, clinical, and viral mutation analysis into one bioinformatics ecosystem. The applications use publicly available data sources and computational methods to support research and educational exploration.

> **Note:** Results are intended for research and educational exploration. They should not be used as a substitute for clinical diagnosis or medical advice.

## Tools

### 🧬 CRISPR Gene Editing Simulator

- Fetches real gene sequences from NCBI.
- Scans both DNA strands for PAM sites.
- Calculates mismatch patterns and simulates Cas9 cleavage.
- Classifies likely DNA repair outcomes.

### 🧪 PRISM — Advanced Clinical Pipeline

- Fetches real clinical trial data.
- Filters patient cohorts by disease, stage, age, and mutation.
- Identifies biomarker outliers and potential trial-fit candidates.

### 🦠 COVID-19 Variant Mutation Tracker

- Fetches SARS-CoV-2 sequences.
- Aligns the Spike region against the Wuhan-Hu-1 reference.
- Detects substitutions and deletions.
- Scores potential immune escape and transmission-related effects.

### 🧫 Protein Structure Analyzer

- Fetches protein structures from the Protein Data Bank.
- Locates mutation sites on protein chains.
- Quantifies physicochemical changes.
- Classifies mutations according to their predicted structural impact.

## Live Demo

Try the deployed applications here:

- **[Open the CRISPR Gene Editing Simulator](https://advance-biocode-whdkzngbtm3ojrz4du4kgc.streamlit.app/)**
- **[Open PRISM — Advanced Clinical Pipeline](https://advance-biocode-95t2yzvng6j7hqordph6bc.streamlit.app/)**

The cloud demo runs without MySQL logging because of cloud networking limitations. To test the complete pipeline, including MySQL logging, clone the repository and run it locally with your own MySQL installation.

## Tech Stack

| Technology | Role |
| --- | --- |
| Python 3.13 | Core application runtime |
| Streamlit | Interactive dashboards and visualization interfaces |
| Pandas | Data handling and tabular analysis |
| NumPy | Numerical calculations and array operations |
| Matplotlib | Static plotting for analysis results |
| MySQL | Persistent logging of pipeline outputs |
| Biopython | Sequence analysis, alignment, and structure parsing |
| Requests | API access to public data sources |
| Py3Dmol | 3D protein structure visualization |

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/CodeXSourabhsingh/Advance-bio_code.git
cd Advance-bio_code
```

### 2. Install dependencies

```bash
pip install biopython pandas numpy matplotlib streamlit mysql-connector-python requests py3dmol
```

### 3. Configure the application

Create a `config.py` file in the project root with your MySQL and NCBI Entrez settings:

```python
MYSQL_HOST = "localhost"
MYSQL_USER = "your_user"
MYSQL_PASSWORD = "your_password"
MYSQL_DATABASE = "your_database"
ENTREZ_EMAIL = "your.email@example.com"
```

Keep credentials private and do not commit `config.py` if it contains secrets.

### 4. Run an application

Run the relevant Streamlit app from the repository root:

```bash
streamlit run <app_file>.py
```

Replace `<app_file>.py` with the Streamlit entry point for the tool you want to use.

## The Ecosystem

The four tools are designed as a connected bioinformatics pipeline:

1. **CRISPR Gene Editing Simulator** identifies mutations at the sequence level.
2. **Protein Structure Analyzer** evaluates whether a mutation may affect protein structure or function.
3. **PRISM** helps identify relevant clinical cohorts and potential trial-fit candidates.
4. **COVID-19 Variant Mutation Tracker** extends mutation analysis to viral evolution and variant monitoring.

## Author

**Sourabh Singh**  

- GitHub: [CodeXSourabhsingh](https://github.com/CodeXSourabhsingh)
- 
- LinkedIn:https://www.linkedin.com/in/sourabh-singh-7b1249434/

## License

This project is licensed under the MIT License.
