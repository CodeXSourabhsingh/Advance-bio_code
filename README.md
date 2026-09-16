# Advance-bio_code

# CRISPR Gene Editing Simulator

A bioinformatics tool that fetches real gene sequences from NCBI, scans both DNA strands for PAM sites, calculates gRNA mismatches, simulates Cas9 cuts, and predicts the repair pathway (NHEJ/HDR). Results are logged to MySQL and visualized in a Streamlit dashboard.

## What It Does

1. **Fetches** real gene sequences from NCBI using Biopython's Entrez API.
2. **Scans** both the forward and reverse complement strands for PAM sites (NGG, NAG, NGA, NGC).
3. **Matches** a user-provided gRNA against each PAM site using Hamming distance.
4. **Filters** candidate sites by GC content (default 30-70%).
5. **Simulates** the Cas9 cut (3 bp upstream of PAM) and the repair pathway (NHEJ or HDR).
6. **Logs** every run to a MySQL vault (`crispr_vault`).
7. **Visualizes** the results in a Streamlit dashboard with 4 charts.

## The Biology

- **PAM Recognition:** Cas9 requires a PAM sequence (NGG for SpCas9) adjacent to the target.
- **gRNA Matching:** The 20-base protospacer upstream of the PAM is compared to the gRNA. 0 mismatches = on-target, 1-2 = off-target, 3+ = ignored.
- **Cut Simulation:** SpCas9 cleaves 3 bp upstream of the PAM.
- **Repair Pathways:** NHEJ (random indels, gene knockout) or HDR (precise edit with a donor template).

## Test Results

| Gene | Sequence Length | PAMs Scanned | On-Target Cuts | Efficiency |
|------|----------------|--------------|----------------|------------|
| HBB  | 10,106 bp      | 857          | 1              | 0.12%      |
| VEGFA| 23,272 bp      | 4,763        | 1              | 0.02%      |

*Tested on HBB gRNA: `CTTGCCCCACAGGGCAGTAA` (20 bases, 60% GC content).*

## Tech Stack

- **Python 3.13**
- **Biopython** — NCBI Entrez API integration
- **Pandas / NumPy** — data processing and mismatch math
- **Matplotlib** — 4 diverse visualizations
- **Streamlit** — interactive dashboard with `st.session_state`
- **MySQL** — run logging and persistence

## Setup

1. Clone the repo:
   ```bash
   git clone https://github.com/CodeXSourabhsingh/bio_code.git
   cd bio_code
