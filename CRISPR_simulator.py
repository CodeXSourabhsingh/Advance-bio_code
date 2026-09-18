import pandas as pd
from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
import streamlit as st
import time, random
import mysql.connector
from config import ENTREZ_EMAIL
try:
    MYSQL_HOST = st.secrets["MYSQL_HOST"]
    MYSQL_PORT = st.secrets["MYSQL_PORT"]
    MYSQL_USER = st.secrets["MYSQL_USER"]
    MYSQL_PASSWORD = st.secrets["MYSQL_PASSWORD"]
    MYSQL_DATABASE = st.secrets["MYSQL_DATABASE"]
except Exception:
   
    from config import MYSQL_HOST, MYSQL_USER, MYSQL_PASSWORD, MYSQL_DATABASE
    MYSQL_PORT = 3306

Entrez.email = "ENTREZ_EMAIL"

def get_db_connection():
    return mysql.connector.connect(
        host=MYSQL_HOST,
        port=MYSQL_PORT,
        user=MYSQL_USER,
        password=MYSQL_PASSWORD,
        database=MYSQL_DATABASE,
    )

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement.get(base, base) for base in reversed(seq))

gene_input = st.text_input("Enter gene name(s) — single or comma-separated:", "HBB")
B = st.selectbox("PAM sequence: ", ["NGG", "NAG", "NGA", "NGC"])
C = st.text_input("Enter gRNA (20 base):", "CTTGCCCCACAGGGCAGTAA")
gc_min, gc_max = st.slider("GC content filter:", 0, 100, (30, 70))
run_button = st.button("Run simulation")

if "crispr_results" not in st.session_state:
    st.session_state.crispr_results = None
if "total_pams" not in st.session_state:
    st.session_state.total_pams = 0

if run_button:
    gene_list = [g.strip() for g in gene_input.split(",") if g.strip()]
    all_results = []
    total_pams_found_across_genes = 0

    for gene_name in gene_list:
        try:
            handle = Entrez.esearch(db="nucleotide", term=f"{gene_name}[Gene Name] AND Homo sapiens[Organism] AND RefSeq[Filter]", retmax=1)
            record = Entrez.read(handle)
            handle.close()

            if not record.get("IdList"):
                # Fallback to general search if RefSeq fails
                handle = Entrez.esearch(db="nucleotide", term=f"{gene_name}[All Fields] AND Homo sapiens[Organism]", retmax=1)
                record = Entrez.read(handle)
                handle.close()
                if not record.get("IdList"):
                    st.warning(f"No gene ID found for {gene_name}.")
                    continue

            gene_id = record["IdList"][0]
            time.sleep(0.5)

            handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
            sequence = SeqIO.read(handle, "fasta")
            handle.close()
            dna = str(sequence.seq).upper()

            st.success(f"Fetched {gene_name}! Gene ID: {gene_id} | Length: {len(dna)} bp")

            pam_pairs = {"NGG": "GG", "NAG": "AG", "NGA": "GA", "NGC": "GC"}
            target_pair = pam_pairs[B]
            
            sequences_to_scan = [("forward", dna), ("reverse", reverse_complement(dna))]
            gene_pams = 0
            
            for strand_name, seq in sequences_to_scan:
                pam_positions = [pos for pos in range(len(seq)-2) if seq[pos+1:pos+3] == target_pair]
                gene_pams += len(pam_positions)
                total_pams_found_across_genes += len(pam_positions)

                for i in pam_positions:
                    target = seq[i-20:i]
                    if len(target) < 20:
                        continue
                    
                    
                    gc_count = target.count('G') + target.count('C')
                    gc_content = (gc_count / 20) * 100
                    if gc_content < gc_min or gc_content > gc_max:
                        continue

                    mismatches = sum(1 for a, b in zip(target, C) if a != b)

                    if mismatches == 0:
                        match_type = "On-target"
                        cut_position = i - 3
                        repair_type = random.choice(["NHEJ", "HDR"])
                    elif mismatches <= 2:
                        match_type = "Off-target"
                        cut_position = None
                        repair_type = None
                    else:
                        continue

                    all_results.append({
                        "gene_name": gene_name,
                        "strand": strand_name,
                        "pam_position": i,
                        "target": target,
                        "gc_content": round(gc_content, 2),
                        "mismatches": mismatches,
                        "match_type": match_type,
                        "cut_position": cut_position,
                        "repair_type": repair_type
                    })

            gene_cuts = [r for r in all_results if r["gene_name"] == gene_name]
            total_cuts = len(gene_cuts)
            efficiency = (total_cuts / gene_pams) * 100 if gene_pams > 0 else 0

               
            mydb = get_db_connection()
            cursor = mydb.cursor()
            insert_query = """
                INSERT INTO crispr_vault (gene_name, gene_id, sequence_length, total_pams, total_cuts, efficiency) 
                VALUES (%s, %s, %s, %s, %s, %s)
            """
            values = (str(gene_name), str(gene_id), int(len(dna)), int(gene_pams), int(total_cuts), float(efficiency))
            cursor.execute(insert_query, values)
            mydb.commit()
            cursor.close()
            mydb.close()

        except Exception as e:
            st.error(f"Error on {gene_name}: {e}")
            continue

    st.session_state.total_pams = total_pams_found_across_genes
    if all_results:
        st.session_state.crispr_results = pd.DataFrame(all_results)
    else:
        st.session_state.crispr_results = pd.DataFrame([{
            "gene_name": gene_input, "strand": "N/A", "pam_position": 0, "target": "None",
            "gc_content": 0, "mismatches": 0, "match_type": "No Cuts", "cut_position": None, "repair_type": None
        }])


if st.session_state.crispr_results is not None:
    df = st.session_state.crispr_results

    total_scanned = st.session_state.total_pams
    total_on_target = len(df[df["match_type"] == "On-target"]) if not df.empty else 0

    if total_on_target == 0:
        st.warning("0 on-target cuts found. The gRNA does not perfectly match this gene on either strand.")
    else:
        st.success(f"Found {total_on_target} On-target cut(s).")

    st.metric("Total PAMs Scanned", total_scanned)
    st.metric("On-Target Cuts", total_on_target)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(["PAMs", "Cuts"], [max(1, total_scanned), total_on_target], color=["blue", "green"])
    st.pyplot(fig)
    plt.close(fig)

    fig2, ax2 = plt.subplots(figsize=(6, 4))
    df['mismatches'].value_counts().sort_index().plot(kind='bar', color='purple', ax=ax2)
    ax2.set_title("Mismatch Distribution")
    ax2.set_xlabel("Number of Mismatches")
    ax2.set_ylabel("PAM Site Count")
    st.pyplot(fig2)
    plt.close(fig2)

    fig3, ax3 = plt.subplots(figsize=(8, 4))
    df_on = df[df['match_type'] == 'On-target']
    if not df_on.empty:
     ax3.scatter(df_on['pam_position'], df_on['mismatches'], color='red', alpha=0.7)
     ax3.set_title("Cut Position Map (On-Target)")
     ax3.set_xlabel("Position in Gene (bp)")
     ax3.set_ylabel("Mismatches")
     st.pyplot(fig3)
    plt.close(fig3)

    if 'strand' in df.columns and not df.empty:
        fig4, ax4 = plt.subplots(figsize=(6, 6))
        df['strand'].value_counts().plot(kind='pie', autopct='%1.1f%%', ax=ax4)
        ax4.set_title("Cuts by DNA Strand")
        ax4.set_ylabel("")
        st.pyplot(fig4)
        plt.close(fig4)

    st.dataframe(df)
    st.download_button("Download CSV", df.to_csv(index=False), "crispr_results.csv")