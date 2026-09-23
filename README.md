# 🧬 Advance-bio_code

> A connected bioinformatics portfolio for sequence analysis, protein structure assessment, clinical cohort filtering, and viral variant tracking.

[![Live Demo](https://img.shields.io/badge/🚀%20Live%20Demo-CRISPR%20Gene%20Editing%20Simulator-2ea44f?style=for-the-badge)](https://advance-biocode-whdkzngbtm3ojrz4du4kgc.streamlit.app/)
[![PRISM Demo](https://img.shields.io/badge/📊%20PRISM%20Demo-Advanced%20Clinical%20Pipeline-2ea44f?style=for-the-badge)](https://advance-biocode-95t2yzvng6j7hqordph6bc.streamlit.app/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Advance-bio_code is a Python-based bioinformatics toolkit designed for research, learning, and exploration across multiple layers of molecular and clinical analysis. The project brings together four interactive applications that connect genomic sequencing, structural biology, clinical cohort analysis, and viral mutation tracking into a single ecosystem.

## Project at a Glance

- 4 interactive applications
- Sequence analysis using public genomic data
- Mutation and cleavage simulation for CRISPR workflows
- Clinical cohort filtering and trial-fit evaluation
- Viral variant tracking for SARS-CoV-2
- Protein structure impact assessment
- Built with Python, Streamlit, Pandas, NumPy, and Biopython

---

## Tools

### 🧬 CRISPR Gene Editing Simulator

- Fetches real gene sequences from NCBI
- Scans both DNA strands for PAM motifs
- Calculates mismatch patterns and simulates Cas9 cleavage
- Classifies likely DNA repair outcomes

### 🧪 PRISM — Advanced Clinical Pipeline

- Fetches public clinical trial data
- Filters patient cohorts by disease, stage, age, and mutation
- Identifies biomarker outliers and potential trial-fit candidates
- Supports large-scale cohort screening for research exploration

### 🦠 COVID-19 Variant Mutation Tracker

- Fetches SARS-CoV-2 sequence data
- Aligns the Spike region against the Wuhan-Hu-1 reference
- Detects substitutions and deletions
- Scores mutations related to immune escape and transmission potential

### 🧫 Protein Structure Analyzer

- Fetches protein structures from the Protein Data Bank
- Locates mutation sites on protein chains
- Quantifies physicochemical changes caused by substitutions
- Classifies mutations based on predicted structural impact

---

## Live Demo

Try the deployed applications here:

- [Open the CRISPR Gene Editing Simulator](https://advance-biocode-whdkzngbtm3ojrz4du4kgc.streamlit.app/)
- [Open PRISM — Advanced Clinical Pipeline](https://advance-biocode-95t2yzvng6j7hqordph6bc.streamlit.app/)

> Note: The cloud-hosted demo runs without MySQL logging because of Streamlit Cloud networking limitations. For the complete pipeline, including MySQL-powered logging, run the project locally with your own MySQL instance.

---

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

---

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

### 4. Run an app

From the repository root, run the relevant Streamlit app:

```bash
streamlit run <app_file>.py
```

Replace `<app_file>.py` with the entry point for the tool you want to use, such as:

```bash
streamlit run CRISPR.py
streamlit run PRISM.py
```

---

## The Ecosystem

The four tools are designed as a connected bioinformatics workflow:

1. **CRISPR Gene Editing Simulator** identifies sequence-level mutations and likely editing outcomes.
2. **Protein Structure Analyzer** evaluates whether a mutation may affect protein structure or function.
3. **PRISM** helps identify relevant clinical cohorts and potential trial-fit candidates.
4. **COVID-19 Variant Mutation Tracker** extends mutation analysis to viral evolution and variant monitoring.

This creates a research-oriented pipeline spanning molecular biology, structural interpretation, clinical relevance, and viral surveillance.

---

## Local Benchmark Environment

The 10,000,000-patient benchmark was tested on:

- Device: HP Laptop 15-fr0xxx
- Processor: 13th Gen Intel Core i5-13420H (2.10 GHz, 8 cores, 12 threads)
- RAM: 16 GB
- Storage: 477 GB SSD
- GPU: Intel UHD Graphics (integrated; not used for computation)
- OS: Windows 11, 64-bit

No dedicated GPU is required. No cloud compute is required. No distributed infrastructure is required. This workflow ran efficiently on a mid-range consumer laptop using pure Python, Pandas, and NumPy.

### Performance

- 500,000 patients → 7–8 seconds
- 1,000,000 patients → 10–15 seconds
- 2,000,000 patients → 20 seconds
- 5,000,000 patients → 50 seconds
- 10,000,000 patients → under 2 minutes

Peak memory: approximately 1.9 GB at 10 million patients.

This indicates that the pipeline is efficient enough for standard hardware. The main bottleneck is the Python object-generation loop rather than computational hardware limitations.

---

## Demo and Scaling Notes

The Streamlit Cloud deployment is suitable for testing and demonstration. For larger-scale analysis, especially with very large datasets:

- Local execution is recommended
- Increase dataset limits locally if needed
- MySQL logging is supported when running on your own machine

Example local scaling note:

> Cloud demo is capped at 5,000 patients due to Streamlit Cloud free-tier limits. To process up to 10,000,000 patients locally, change `max_value=5000` to `max_value=10000000` in `PRISM.py` and run it locally.

---

## Notes

> Results are intended for research and educational exploration. They should not be used as a substitute for clinical diagnosis or medical advice.

---

## Author

**Sourabh Singh**

- GitHub: [CodeXSourabhsingh](https://github.com/CodeXSourabhsingh)
- LinkedIn: [Sourabh Singh](https://www.linkedin.com/in/sourabh-singh-7b1249434/)

---

## License

This project is licensed under the MIT License.
