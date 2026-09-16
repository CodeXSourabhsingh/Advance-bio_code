# 🧬 Advance-bio_code

# ⚡ CRISPR Gene Editing Simulator

A cutting-edge bioinformatics tool that fetches real gene sequences from NCBI, scans both DNA strands for PAM sites, calculates gRNA mismatches, simulates Cas9 cuts, and predicts the repair pathway (NHEJ/HDR).

---

## 🎯 What It Does
  
- ✅ **Fetches** real gene sequences from NCBI using Biopython's Entrez API
- 🔍 **Scans** both the forward and reverse complement strands for PAM sites (NGG, NAG, NGA, NGC)
- 🎯 **Matches** a user-provided gRNA against each PAM site using Hamming distance
- 📊 **Filters** candidate sites by GC content (default 30-70%)
- ⚔️ **Simulates** the Cas9 cut (3 bp upstream of PAM) and the repair pathway (NHEJ or HDR)
- 💾 **Logs** every run to a MySQL vault (`crispr_vault`)
- 📈 **Visualizes** the results in a Streamlit dashboard with 4 interactive charts

---

## 🔬 The Biology

- **PAM Recognition:** Cas9 requires a PAM sequence (NGG for SpCas9) adjacent to the target
- **gRNA Matching:** The 20-base protospacer upstream of the PAM is compared to the gRNA
  - 0 mismatches = on-target ✅
  - 1-2 mismatches = off-target ⚠️
  - 3+ mismatches = ignored ❌
- **Cut Simulation:** SpCas9 cleaves 3 bp upstream of the PAM
- **Repair Pathways:** NHEJ (random indels, gene knockout) or HDR (precise edit with donor template)

---

## 📊 Test Results

| Gene  | Sequence Length | PAMs Scanned | On-Target Cuts | Efficiency |
|-------|-----------------|--------------|----------------|-----------|
| HBB   | 10,106 bp       | 857          | 1              | 0.12%     |
| VEGFA | 23,272 bp       | 4,763        | 1              | 0.02%     |

*Tested on HBB gRNA: `CTTGCCCCACAGGGCAGTAA` (20 bases, 60% GC content)*

---

## 🛠️ Tech Stack

| Technology | Purpose |
|-----------|---------|
| **Python 3.13** | Core language |
| **Biopython** | NCBI Entrez API integration |
| **Pandas / NumPy** | Data processing & mismatch calculations |
| **Matplotlib** | Interactive visualizations |
| **Streamlit** | Dynamic dashboard with session state |
| **MySQL** | Run logging & data persistence |

---

## 🚀 Quick Start

### Prerequisites
- Python 3.13+
- MySQL Server (for data logging)
- NCBI API key (free, from [NCBI](https://www.ncbi.nlm.nih.gov/account/settings/))

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/CodeXSourabhsingh/Advance-bio_code.git
   cd Advance-bio_code
   ```

2. **Create a virtual environment:**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

4. **Configure MySQL connection:**
   - Update `config.py` or `.env` with your MySQL credentials:
     ```python
     DB_HOST = "localhost"
     DB_USER = "your_username"
     DB_PASSWORD = "your_password"
     DB_NAME = "crispr_vault"
     ```

5. **Set your NCBI API key:**
   ```bash
   export NCBI_API_KEY="your_ncbi_api_key"
   ```

---

## 📖 How to Use

### Option 1: Command Line Interface

```bash
python crispr_simulator.py \
  --gene HBB \
  --grna "CTTGCCCCACAGGGCAGTAA" \
  --pam_types NGG \
  --gc_min 30 \
  --gc_max 70 \
  --repair NHEJ
```

**Parameters:**
- `--gene`: Gene accession ID (e.g., `HBB`, `VEGFA`)
- `--grna`: 20-base guide RNA sequence
- `--pam_types`: PAM sites to scan (NGG, NAG, NGA, NGC)
- `--gc_min`, `--gc_max`: GC content filter range (%)
- `--repair`: Repair pathway (NHEJ or HDR)

### Option 2: Streamlit Dashboard (Interactive)

```bash
streamlit run app.py
```

**Features:**
- 🎨 Real-time sequence visualization
- 📊 Interactive PAM site location plots
- 📈 Mismatch distribution charts
- 💾 Export results as CSV/JSON
- 🔍 Filter and sort results dynamically

---

## 💡 Examples

### Example 1: Target HBB Gene (Sickle Cell Disease)

**Goal:** Find optimal CRISPR targets in the HBB gene

```python
from crispr_simulator import CRISPRSimulator

# Initialize simulator
simulator = CRISPRSimulator(
    gene_id="NM_000518",      # HBB gene
    grna="CTTGCCCCACAGGGCAGTAA",
    pam="NGG",
    gc_range=(30, 70)
)

# Run simulation
results = simulator.run()

# View results
print(f"Total PAM sites found: {len(results)}")
print(f"On-target matches: {results['on_target'].sum()}")
print(results.head())

# Visualize
simulator.visualize()
simulator.save_to_mysql()
```

**Output:**
```
Total PAM sites found: 857
On-target matches: 1

   Position  PAM_Sequence  Strand  Mismatches  GC_Content  Repair_Type
0     5421  NGG           +      0           60%         NHEJ
1     7834  NGG           -      2           55%         HDR
2     9102  NGG           +      4           65%         Ignored
```

---

### Example 2: Multiple Gene Comparison

**Goal:** Compare CRISPR efficiency across multiple genes

```python
genes = ["NM_000518", "NM_003376"]  # HBB, VEGFA
grna = "CTTGCCCCACAGGGCAGTAA"

results_summary = []

for gene_id in genes:
    simulator = CRISPRSimulator(gene_id=gene_id, grna=grna)
    results = simulator.run()
    
    summary = {
        "Gene": gene_id,
        "Total_PAMs": len(results),
        "On_Target": results['on_target'].sum(),
        "Efficiency": (results['on_target'].sum() / len(results)) * 100
    }
    results_summary.append(summary)

import pandas as pd
comparison_df = pd.DataFrame(results_summary)
print(comparison_df)
```

**Output:**
```
         Gene  Total_PAMs  On_Target  Efficiency
0  NM_000518        857         1         0.12%
1  NM_003376       4763         1         0.02%
```

---

### Example 3: HDR-Based Precise Editing

**Goal:** Design a precise edit using Homology-Directed Repair (HDR)

```python
# For sickle cell correction
simulator = CRISPRSimulator(
    gene_id="NM_000518",
    grna="CTTGCCCCACAGGGCAGTAA",
    repair_pathway="HDR",
    donor_template="CTTGCCCCACAGGGCAGTAA"  # Corrected sequence
)

results = simulator.run()
on_target_sites = results[results['mismatches'] == 0]

print(f"Precise edit targets: {len(on_target_sites)}")
for idx, site in on_target_sites.iterrows():
    print(f"Position {site['position']}: {site['pam_sequence']} ({site['gc_content']}% GC)")

# Generate repair prediction
repair_report = simulator.predict_repair_efficiency()
print(repair_report)
```

---

### Example 4: Export & Share Results

```python
# Run simulation
simulator = CRISPRSimulator(gene_id="NM_000518", grna="CTTGCCCCACAGGGCAGTAA")
results = simulator.run()

# Export formats
results.to_csv("crispr_results.csv", index=False)
results.to_json("crispr_results.json")

# Generate visual report
simulator.generate_report(output_file="crispr_report.pdf")

# Save to MySQL vault
simulator.save_to_mysql(run_name="HBB_Experiment_v1")
print("✅ Results saved to database!")
```

---

## 🎓 Understanding the Output

### Results DataFrame Columns:

| Column | Description |
|--------|-------------|
| `Position` | Location in the gene sequence |
| `PAM_Sequence` | Actual PAM found (NGG, NAG, etc.) |
| `Strand` | Forward (+) or Reverse (-) |
| `Mismatches` | Hamming distance to gRNA (0 = perfect match) |
| `GC_Content` | Percentage of G/C bases (30-70% optimal) |
| `Repair_Type` | Predicted repair pathway (NHEJ or HDR) |
| `Off_Target_Risk` | Estimated off-target activity score |

### Visualizations Generated:

1. **PAM Site Distribution** — Histogram of PAM locations
2. **Mismatch Heatmap** — gRNA alignment quality
3. **GC Content Analysis** — Distribution across binding sites
4. **Efficiency Metrics** — On-target vs off-target ratios

---

## ⚙️ Configuration

Edit `config.py` to customize:

```python
# PAM Site Definitions
PAM_SITES = {
    'NGG': 'SpCas9 (Streptococcus pyogenes)',
    'NAG': 'Off-target NGG',
    'NGA': 'Alternative PAM',
    'NGC': 'Rare PAM'
}

# GC Content Range (%)
GC_MIN = 30
GC_MAX = 70

# Mismatch Tolerance
MAX_MISMATCHES = 4  # 0 = on-target, 3+ = ignore

# MySQL Settings
DB_HOST = "localhost"
DB_PORT = 3306
DB_USER = "root"
```

---

## 📊 Database Schema

The MySQL `crispr_vault` automatically stores:
- Simulation metadata (timestamp, user, parameters)
- Complete results table (PAM sites, mismatches, predictions)
- Visualization snapshots
- Run history for comparison

Query your runs:
```sql
SELECT * FROM crispr_vault WHERE gene = 'HBB' ORDER BY created_at DESC;
```

---

## ⚠️ Known Issues

**Windows Smart App Control:** Some Windows 11 users may experience an `ImportError` for `_codonaligner` due to Smart App Control blocking Biopython's compiled files. 

**Solution:** Temporarily disable Smart App Control in Windows settings.

---

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit your changes (`git commit -m "Add your feature"`)
4. Push to the branch (`git push origin feature/your-feature`)
5. Open a Pull Request

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
