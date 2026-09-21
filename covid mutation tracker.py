import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mysql.connector
import time
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner

try:
    MYSQL_HOST = st.secrets["MYSQL_HOST"]
    MYSQL_PORT = st.secrets["MYSQL_PORT"]
    MYSQL_USER = st.secrets["MYSQL_USER"]
    MYSQL_PASSWORD = st.secrets["MYSQL_PASSWORD"]
    MYSQL_DATABASE = st.secrets["MYSQL_DATABASE"]
    ENTREZ_EMAIL = st.secrets["ENTREZ_EMAIL"]
except Exception:
    from config import MYSQL_HOST, MYSQL_USER, MYSQL_PASSWORD, MYSQL_DATABASE, ENTREZ_EMAIL
    MYSQL_PORT = 3306

Entrez.email = ENTREZ_EMAIL

SPIKE_START = 21562
SPIKE_END = 25383

ESCAPE_MUTATIONS = {"E484K", "L452R", "K417N", "F486V", "R493Q"}
TRANSMISSIBILITY_MUTATIONS = {"N501Y", "N481K", "D614G"}
LINEAGE_SIGNATURES = {
    "BA.3.2.2": {"R493Q", "N529"},
    "PQ.16.1.1": {"D253G", "N417T", "D420N", "I478T"},
    "BA.3.2": {"N501Y", "D614G", "E484K"},
}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

@st.cache_data
def fetch_covid_sequence(accession_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq).upper()
    except Exception as e:
        st.error(f"Error fetching {accession_id}: {e}")
        return None

@st.cache_data
def search_variant(variant_name):
    try:
        search_term = f"SARS-CoV-2[Organism] AND {variant_name}[All Fields] AND complete genome[Title]"
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1, sort="date")
        record = Entrez.read(handle)
        handle.close()
        if not record.get("IdList"):
            return None
        return record["IdList"][0]
    except Exception as e:
        st.error(f"Search error for {variant_name}: {e}")
        return None

@st.cache_data
def get_reference():
    time.sleep(1)
    return fetch_covid_sequence("NC_045512.2")

def find_aligned_index(ref_aligned, ref_pos):
    count = 0
    for j in range(len(ref_aligned)):
        if ref_aligned[j] != '-':
            if count == ref_pos:
                return j
            count += 1
    return 0

def align_spike(ref_spike, sample_spike):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    alignment = aligner.align(ref_spike, sample_spike)[0]
    ref_aligned = str(alignment[0])
    sample_aligned = str(alignment[1])

    mutations = []
    ref_pos = 0

    for i in range(len(ref_aligned)):
        r = ref_aligned[i]
        s = sample_aligned[i]

        if r == '-':
            continue

        absolute_pos = SPIKE_START + ref_pos
        aa_pos = (ref_pos // 3) + 1

        if s == '-':
            mutations.append({
                "position": absolute_pos,
                "ref_base": r,
                "alt_base": "del",
                "mutation_type": "Deletion",
                "gene": "Spike",
                "aa_change": f"del{aa_pos}"
            })
            ref_pos += 1
            continue

        if r != s:
            codon_start_ref = (ref_pos // 3) * 3
            codon_start_align = find_aligned_index(ref_aligned, codon_start_ref)
            ref_codon = ref_aligned[codon_start_align:codon_start_align + 3].replace('-', '')
            alt_codon = sample_aligned[codon_start_align:codon_start_align + 3].replace('-', '')

            if len(ref_codon) == 3 and len(alt_codon) == 3:
                ref_aa = CODON_TABLE.get(ref_codon, "X")
                alt_aa = CODON_TABLE.get(alt_codon, "X")
                aa_change = f"{ref_aa}{aa_pos}{alt_aa}"
            else:
                aa_change = "N/A"

            mutations.append({
                "position": absolute_pos,
                "ref_base": r,
                "alt_base": s,
                "mutation_type": "Substitution",
                "gene": "Spike",
                "aa_change": aa_change
            })
        ref_pos += 1

    return mutations

def match_lineage(mutation_set):
    for lineage, signature in LINEAGE_SIGNATURES.items():
        if signature.issubset(mutation_set):
            return lineage
    return "Unassigned"

def score_variant(mutations):
    aa_changes = [m.get("aa_change", "") for m in mutations]
    escape_score = len(set(aa_changes) & ESCAPE_MUTATIONS)
    transmissibility_score = len(set(aa_changes) & TRANSMISSIBILITY_MUTATIONS)
    return {
        "escape_score": escape_score,
        "transmissibility_score": transmissibility_score
    }

@st.cache_resource
def get_db_connection():
    return mysql.connector.connect(
        host=MYSQL_HOST,
        port=MYSQL_PORT,
        user=MYSQL_USER,
        password=MYSQL_PASSWORD,
        database=MYSQL_DATABASE
    )

def log_to_mysql(sample_data, mutation_list):
    try:
        mydb = get_db_connection()
        cursor = mydb.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS variant_tracker (
                id INT AUTO_INCREMENT PRIMARY KEY,
                sample_id VARCHAR(255),
                lineage VARCHAR(255),
                total_mutations INT,
                immune_escape_score INT,
                transmissibility_score INT,
                logged_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS mutation_details (
                id INT AUTO_INCREMENT PRIMARY KEY,
                sample_id VARCHAR(255),
                position INT,
                ref_base VARCHAR(10),
                alt_base VARCHAR(10),
                aa_change VARCHAR(50),
                logged_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        cursor.execute("""
            INSERT INTO variant_tracker
            (sample_id, lineage, total_mutations, immune_escape_score, transmissibility_score)
            VALUES (%s, %s, %s, %s, %s)
        """, (
            sample_data["sample_id"],
            sample_data["lineage"],
            sample_data["total_mutations"],
            sample_data["escape_score"],
            sample_data["transmissibility_score"]
        ))
        for m in mutation_list:
            cursor.execute("""
                INSERT INTO mutation_details
                (sample_id, position, ref_base, alt_base, aa_change)
                VALUES (%s, %s, %s, %s, %s)
            """, (
                sample_data["sample_id"],
                m["position"],
                m["ref_base"],
                m["alt_base"],
                m.get("aa_change", "N/A")
            ))
        mydb.commit()
        cursor.close()
    except Exception as e:
        st.error(f"Database Error: {e}")

st.title("COVID-19 Variant Mutation Tracker")

if "results" not in st.session_state:
    st.session_state.results = None

st.subheader("Sample Input")

input_mode = st.radio("Search by:", ["Variant Name", "Accession ID"], horizontal=True)

if input_mode == "Variant Name":
    user_input = st.text_input("Enter variant name (e.g., Omicron, Delta):", "Omicron")
else:
    user_input = st.text_input("Enter NCBI Accession ID:", "OL672836.1")

run_button = st.button("Fetch and Analyze")

if run_button:
    with st.spinner("Fetching sequences from NCBI..."):
        ref_seq = get_reference()

        if input_mode == "Variant Name":
            target_id = search_variant(user_input)
            if target_id is None:
                st.error(f"No sequence found for {user_input}.")
                st.stop()
            st.info(f"Found accession ID: {target_id}")
            time.sleep(1)
            current_seq = fetch_covid_sequence(target_id)
        else:
            target_id = user_input
            time.sleep(1)
            current_seq = fetch_covid_sequence(user_input)

    if ref_seq is None or current_seq is None:
        st.error("Failed to fetch. Try again in 30 seconds.")
        st.stop()

    st.write(f"Reference: {len(ref_seq)} bp | Sample: {len(current_seq)} bp")

    if len(current_seq) < SPIKE_END:
        st.error(f"Sample is too short. Needs at least {SPIKE_END} bp.")
        st.stop()

    ref_spike = ref_seq[SPIKE_START:SPIKE_END]
    sample_spike = current_seq[SPIKE_START:SPIKE_END]

    with st.spinner("Aligning Spike region..."):
        mutations = align_spike(ref_spike, sample_spike)

    mutation_set = {m["aa_change"] for m in mutations}
    lineage = match_lineage(mutation_set)
    scores = score_variant(mutations)

    sample_id = f"SAMPLE-{pd.Timestamp.now().strftime('%Y%m%d%H%M%S')}"

    sample_data = {
        "sample_id": sample_id,
        "lineage": lineage,
        "total_mutations": len(mutations),
        "escape_score": scores["escape_score"],
        "transmissibility_score": scores["transmissibility_score"]
    }

    log_to_mysql(sample_data, mutations)

    st.session_state.results = {
        "sample_id": sample_id,
        "lineage": lineage,
        "mutations": mutations,
        "scores": scores,
        "sample_data": sample_data
    }

if st.session_state.results is not None:
    res = st.session_state.results

    st.divider()
    st.subheader("Analysis Results")

    tab1, tab2, tab3 = st.tabs(["Mutation Table", "Variant Summary", "Charts"])

    with tab1:
        df = pd.DataFrame(res["mutations"])
        st.dataframe(df)

    with tab2:
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Lineage", res["lineage"])
        with col2:
            st.metric("Escape Score", res["scores"]["escape_score"])
        with col3:
            st.metric("Transmissibility Score", res["scores"]["transmissibility_score"])

    with tab3:
        df = pd.DataFrame(res["mutations"])

        if not df.empty:
            type_counts = df["mutation_type"].value_counts()

            fig1, ax1 = plt.subplots(figsize=(6, 4))
            ax1.bar(type_counts.index, type_counts.values, color="steelblue")
            ax1.set_title("Mutation Types")
            ax1.set_xlabel("Type")
            ax1.set_ylabel("Count")
            st.pyplot(fig1)
            plt.close(fig1)

            if "gene" in df.columns:
                gene_counts = df["gene"].value_counts()
                fig2, ax2 = plt.subplots(figsize=(6, 6))
                ax2.pie(gene_counts.values, labels=gene_counts.index, autopct='%1.1f%%')
                ax2.set_title("Mutations by Gene")
                st.pyplot(fig2)
                plt.close(fig2)
        else:
            st.info("No mutations detected.")

    st.dataframe(df)
    st.download_button("Download CSV", df.to_csv(index=False), "covid_results.csv")