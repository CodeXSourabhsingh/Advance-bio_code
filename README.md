# 🧬 Advance-bio_code

Advance-bio_code is a connected bioinformatics portfolio that brings together CRISPR editing analysis, structural mutation assessment, clinical cohort filtering, and viral variant tracking in a single research workflow. The four Python applications use public data sources and sequence inputs to connect molecular observations with patient impact, making the repository a compact bioinformatics toolkit for studying how variants emerge, damage proteins, and influence disease outcomes.

## Tools

- CRISPR Gene Editing Simulator — Fetches real gene sequences from NCBI, scans both DNA strands for PAM sites, calculates mismatch patterns, simulates Cas9 cleavage, and classifies likely repair outcomes.
- PRISM - Advanced Clinical Pipeline — Fetches real clinical trial data, filters patient cohorts by disease, stage, age, and mutation, and identifies biomarker outliers and trial-fit candidates.
- COVID-19 Variant Mutation Tracker — Fetches SARS-CoV-2 sequences, aligns the Spike region to the Wuhan-Hu-1 reference, detects substitutions and deletions, and scores immune escape and transmissibility.
- Protein Structure Analyzer — Fetches PDB structures, locates mutation sites on protein chains, quantifies physicochemical changes, and classifies mutations based on structural impact.

## Live Demo

CRISPR Gene Editing Simulator: https://codexsourabhsingh-advance-bio-code.streamlit.app

The live cloud demo runs without MySQL logging due to cloud networking limitations. To test the full pipeline including MySQL logging, clone the repository and run locally with your own MySQL instance.

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
| Py3Dmol | 3D protein structure viewing |

## Quick Start

1. git clone the repo
2. pip install biopython pandas numpy matplotlib streamlit mysql-connector-python requests py3dmol
3. create config.py with MYSQL_HOST, MYSQL_USER, MYSQL_PASSWORD, MYSQL_DATABASE, ENTREZ_EMAIL
4. run each app with streamlit run

Example configuration:

```python
MYSQL_HOST = "localhost"
MYSQL_USER = "your_user"
MYSQL_PASSWORD = "your_password"
MYSQL_DATABASE = "your_database"
ENTREZ_EMAIL = "your.email@example.com"
```

## The Ecosystem

The four tools work as a connected bioinformatics pipeline. CRISPR finds the mutation at the sequence level, Protein Analyzer predicts whether that mutation is structurally damaging, PRISM identifies the patient populations most likely affected by the mutation, and COVID tracker monitors how those mutations evolve over time in circulating viral sequences.

## Roadmap

- V1 shipped for all four tools
- V2 planned for full Pango lineage database in COVID
- V3 planned for conservation scoring in Protein
- V4 planned for ML-based binding affinity prediction in Drug-Target Explorer (next project)

## Author

Sourabh Singh is a self-taught bioinformatics developer building applied tools at the intersection of genomics, clinical data, and computational biology.

- GitHub: https://github.com/CodeXSourabhsingh
- LinkedIn: https://www.linkedin.com/in/sourabh-singh

## License

MIT

