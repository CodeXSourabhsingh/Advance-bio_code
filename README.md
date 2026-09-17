# 🧬 Advance-bio_code

## ⚡ CRISPR Gene Editing Simulator

A cutting-edge bioinformatics tool that fetches real gene sequences from NCBI, scans both DNA strands for PAM sites, calculates gRNA mismatches, simulates Cas9 cuts, and predicts the repair pathway (NHEJ or HDR) in an interactive Streamlit dashboard with integrated MySQL logging.

---

## 🎯 What It Does

- ✅ **Fetches** real gene sequences from NCBI using Biopython's Entrez API
- 🔍 **Scans** both the forward and reverse complement strands for PAM sites (NGG, NAG, NGA, NGC)
- 🎯 **Matches** a user-provided gRNA against each PAM site using Hamming distance
- 📊 **Filters** candidate sites by GC content (default 30–70%)
- ⚔️ **Simulates** the Cas9 cut (3 bp upstream of PAM) and the repair pathway (NHEJ or HDR)
- 💾 **Logs** every run to a MySQL vault (`crispr_vault`)
- 📈 **Visualizes** the results in a Streamlit dashboard with 4 interactive charts

---

## 🔬 The Biology

- **PAM Recognition:** Cas9 requires a PAM sequence (NGG for SpCas9) adjacent to the target
- **gRNA Matching:** The 20-base protospacer upstream of the PAM is compared to the gRNA
  - 0 mismatches = on-target ✅
  - 1–2 mismatches = off-target ⚠️
  - 3+ mismatches = ignored ❌
- **Cut Simulation:** SpCas9 cleaves 3 bp upstream of the PAM
- **Repair Pathways:** NHEJ (random indels, gene knockout) or HDR (precise edit with donor template)

---

## 📊 Test Results

| Gene  | Sequence Length | PAMs Scanned | On-Target Cuts | Efficiency |
|-------|-----------------|--------------|----------------|-----------|
| HBB   | 10,106 bp       | 857          | 1              | 0.12%     |
| VEGFA | 23,272 bp       | 4,763        | 1              | 0.02%     |
| EGFR  | 1,575 bp        | 257          | 0              | 0%        |
| TP53  | 32,772 bp       | 4,826        | 0              | 0%        |

*EGFR, TP53, and KRAS returned 0 cuts because the HBB-specific gRNA does not perfectly match their sequence. This is correct behavior, not a bug.*

*Tested with HBB gRNA: `CTTGCCCCACAGGGCAGTAA` (20 bases, 60% GC content)*

---

## 🛠️ Tech Stack

| Technology | Purpose |
|-----------|---------|
| **Python 3.13** | Core language |
| **Biopython** | NCBI Entrez API integration |
| **Pandas / NumPy** | Data processing and mismatch calculations |
| **Matplotlib** | Visualizations |
| **Streamlit** | Dashboard with session state |
| **MySQL** | Run logging and data persistence |

---

## 🚀 Quick Start
## 🚀 Live Demo

🔴 **Live Demo:** [advance-biocode-whdikzngbtm3ojrz4du4kgc.streamlit.app](https://advance-biocode-whdikzngbtm3ojrz4du4kgc.streamlit.app)

**Note on Database Logging:** The live cloud demo runs without MySQL logging due to cloud networking limitations. To test the full pipeline including MySQL logging, follow the local installation instructions below.
### Prerequisites
- Python 3.13+
- MySQL Server (for data logging)

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/CodeXSourabhsingh/Advance-bio_code.git
   cd Advance-bio_code
   ```

2. **Install dependencies:**
   ```bash
   pip install biopython pandas numpy matplotlib streamlit mysql-connector-python
   ```

3. **Create a `config.py` file:**
   ```python
   MYSQL_HOST = "localhost"
   MYSQL_USER = "your_username"
   MYSQL_PASSWORD = "your_password"
   MYSQL_DATABASE = "crispr_vault"
   ```

4. **Run the Streamlit app:**
   ```bash
   streamlit run CRISPR_simulator.py
   ```

5. **Open in your browser:**
   Navigate to `http://localhost:8501`

---

## 📖 How to Use

The Streamlit application presents an interactive interface:

1. **Enter gene name(s):** Input a single gene (e.g., `HBB`) or comma-separated list (e.g., `HBB, VEGFA, EGFR`)
2. **Select PAM sequence:** Choose from NGG, NAG, NGA, or NGC
3. **Enter gRNA:** Provide a 20-base guide RNA sequence (e.g., `CTTGCCCCACAGGGCAGTAA`)
4. **Adjust GC content filter:** Use the slider to set the min/max GC content range (default 30–70%)
5. **Click "Run simulation":** The app fetches sequences, scans for PAM sites, calculates mismatches, and logs results to MySQL
6. **View results:**
   - Metrics: Total PAMs scanned and on-target cuts
   - 4 interactive Matplotlib charts
   - Complete results dataframe
   - CSV download button

---

## 🎓 Understanding the Output

### Results DataFrame Columns:

| Column | Description |
|--------|-------------|
| `gene_name` | Gene symbol |
| `strand` | Forward (+) or reverse (−) |
| `pam_position` | Position of PAM in sequence (bp) |
| `target` | 20-base protospacer |
| `gc_content` | GC percentage of protospacer |
| `mismatches` | Hamming distance to gRNA (0 = perfect match) |
| `match_type` | On-target, Off-target, or No Cuts |
| `cut_position` | Position 3 bp upstream of PAM (on-target only) |
| `repair_type` | NHEJ or HDR (randomly assigned on-target) |

### Visualizations Generated:

1. **PAMs vs Cuts** — Bar chart showing total PAM sites scanned vs on-target cuts
2. **Mismatch Distribution** — Histogram of mismatch counts across all PAM sites
3. **Cut Position Map** — Scatter plot of on-target cut positions vs mismatches
4. **Strand Distribution** — Pie chart showing cut distribution between forward and reverse strands

---

## ⚙️ Configuration

Create a `config.py` file in the repository root with these variables:

```python
MYSQL_HOST = "localhost"
MYSQL_USER = "your_username"
MYSQL_PASSWORD = "your_password"
MYSQL_DATABASE = "crispr_vault"
```

The Streamlit interface provides inputs for:
- Gene name(s)
- PAM sequence selection (NGG, NAG, NGA, NGC)
- 20-base gRNA sequence
- GC content filter (slider: 0–100%)

---

## 📊 Database Schema

The `crispr_vault` table automatically stores:

| Column | Type | Description |
|--------|------|-------------|
| `id` | INT (AUTO_INCREMENT PRIMARY KEY) | Unique run identifier |
| `gene_name` | VARCHAR(255) | Gene symbol |
| `gene_id` | VARCHAR(255) | NCBI gene ID |
| `sequence_length` | INT | Length of fetched sequence (bp) |
| `total_pams` | INT | Total PAM sites found on both strands |
| `total_cuts` | INT | On-target cuts found |
| `efficiency` | FLOAT | Efficiency percentage (cuts / PAMs) |

**Query your runs:**
```sql
SELECT * FROM crispr_vault WHERE gene_name = 'HBB' ORDER BY id DESC;
```

---

## ⚠️ Known Issues

**Windows Smart App Control:** Some Windows 11 users may experience an `ImportError` for `_codonaligner` due to Smart App Control blocking Biopython's compiled files.

**Solution:** Temporarily disable Smart App Control in Windows settings.

---

## 🗺️ Roadmap

- **V1 (Shipped):** NCBI fetch, dual-strand PAM scan, Hamming distance mismatch calculation, Cas9 cut simulation, NHEJ/HDR repair assignment, MySQL logging, 4 visualizations
- **V2 (October):** Indel simulation for NHEJ with premature stop codon detection, HDR donor template input for precise gene correction
- **V3 (2027):** Genome-wide off-target risk assessment

---

## 👨‍💻 Author

**Sourabh Singh** — Self-taught bioinformatics developer

- **GitHub:** https://github.com/CodeXSourabhsingh
- **LinkedIn:** https://www.linkedin.com/in/sourabh-singh-7b124934/

---

## 📝 License & Citation

Built with ❤️ for bioinformatics research

If you use this tool in your research, please cite:
```bibtex
@software{advance_bio_code,
  title={Advance-bio_code: CRISPR Gene Editing Simulator},
  author={Sourabh Singh},
  url={https://github.com/CodeXSourabhsingh/Advance-bio_code},
  year={2024}
}
```

---

## 📧 Support

For issues, questions, or feature requests, please open an [Issue](https://github.com/CodeXSourabhsingh/Advance-bio_code/issues) on GitHub.
