import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mysql.connector
from sklearn.decomposition import PCA

try:
    MYSQL_HOST = st.secrets["MYSQL_HOST"]
    MYSQL_PORT = st.secrets["MYSQL_PORT"]
    MYSQL_USER = st.secrets["MYSQL_USER"]
    MYSQL_PASSWORD = st.secrets["MYSQL_PASSWORD"]
    MYSQL_DATABASE = st.secrets["MYSQL_DATABASE"]
except Exception:
    from config import MYSQL_HOST, MYSQL_USER, MYSQL_PASSWORD, MYSQL_DATABASE
    MYSQL_PORT = 3306

TAXONOMY_MAP = {
    "OTU_1": "Bacteroides",
    "OTU_2": "Firmicutes",
    "OTU_3": "Proteobacteria",
    "OTU_4": "Actinobacteria",
    "OTU_5": "Fusobacteria",
    "OTU_6": "Verrucomicrobia",
    "OTU_7": "Euryarchaeota",
    "OTU_8": "Cyanobacteria",
    "OTU_9": "Spirochaetes",
    "OTU_10": "Chloroflexi",
    "OTU_11": "Acidobacteria",
    "OTU_12": "Planctomycetes"
}

FUNCTION_MAP = {
    "Bacteroides": "Polysaccharide Degradation",
    "Firmicutes": "Butyrate Production",
    "Proteobacteria": "Nitrogen Fixation",
    "Actinobacteria": "Antibiotic Production",
    "Fusobacteria": "Inflammatory Signaling",
    "Verrucomicrobia": "Mucin Degradation",
    "Euryarchaeota": "Methanogenesis",
    "Cyanobacteria": "Photosynthesis",
    "Spirochaetes": "Motility",
    "Chloroflexi": "Carbon Fixation",
    "Acidobacteria": "Organic Acid Metabolism",
    "Planctomycetes": "Ammonia Oxidation"
}

@st.cache_data
def generate_otu_table(num_samples, num_taxa):
    np.random.seed(42)
    num_taxa = min(num_taxa, 12)
    data = np.random.gamma(shape=2, scale=100, size=(num_samples, num_taxa)).astype(int)
    data[np.random.rand(num_samples, num_taxa) < 0.4] = 0
    for i in range(num_samples):
        if data[i].sum() == 0:
            data[i, np.random.randint(0, num_taxa)] = np.random.randint(50, 200)
    taxa_names = [f"OTU_{i+1}" for i in range(num_taxa)]
    sample_names = [f"Sample_{i+1}" for i in range(num_samples)]
    return pd.DataFrame(data, index=sample_names, columns=taxa_names)

def clean_otu_table(df):
    df = df.loc[:, (df != 0).any(axis=0)]
    totals = df.sum(axis=1)
    df = df[totals > 0]
    return df

def to_relative_abundance(df):
    totals = df.sum(axis=1)
    rel = df.div(totals, axis=0)
    rel = rel.div(rel.sum(axis=1), axis=0)
    return rel.fillna(0)

def alpha_diversity(raw_df, rel_df):
    results = []
    for sample in rel_df.index:
        row = rel_df.loc[sample]
        row_nonzero = row[row > 0]
        observed = len(row_nonzero)
        if observed == 0:
            continue
        shannon = -np.sum(row_nonzero * np.log(row_nonzero))
        simpson = 1 - np.sum(row_nonzero ** 2)
        counts = raw_df.loc[sample]
        f1 = int(np.sum(counts == 1))
        f2 = int(np.sum(counts == 2))
        if f2 > 0:
            chao1 = observed + (f1 ** 2) / (2 * f2)
        else:
            chao1 = observed + (f1 * (f1 - 1)) / 2
        evenness = shannon / np.log(observed) if observed > 1 else 0
        results.append({
            "sample_id": sample,
            "observed_otus": observed,
            "shannon": round(shannon, 4),
            "simpson": round(simpson, 4),
            "chao1": round(chao1, 2),
            "evenness": round(evenness, 4)
        })
    return pd.DataFrame(results)

def bray_curtis_matrix(rel_df):
    samples = rel_df.index.tolist()
    matrix = pd.DataFrame(0.0, index=samples, columns=samples)
    for i, a in enumerate(samples):
        for j, b in enumerate(samples):
            if i == j:
                continue
            row_a = rel_df.loc[a].values
            row_b = rel_df.loc[b].values
            numerator = np.sum(np.abs(row_a - row_b))
            denominator = np.sum(row_a + row_b)
            if denominator > 0:
                value = numerator / denominator
                matrix.loc[a, b] = min(max(value, 0.0), 1.0)
    return matrix

def jaccard_matrix(df):
    presence = (df > 0).astype(int)
    samples = presence.index.tolist()
    matrix = pd.DataFrame(0.0, index=samples, columns=samples)
    for i, a in enumerate(samples):
        for j, b in enumerate(samples):
            if i == j:
                continue
            a_set = set(presence.loc[a][presence.loc[a] == 1].index)
            b_set = set(presence.loc[b][presence.loc[b] == 1].index)
            intersection = len(a_set & b_set)
            union = len(a_set | b_set)
            if union > 0:
                value = 1 - (intersection / union)
                matrix.loc[a, b] = min(max(value, 0.0), 1.0)
    return matrix

def pcoa_projection(bray_matrix):
    pca = PCA(n_components=2)
    coords = pca.fit_transform(bray_matrix.values)
    return pd.DataFrame(coords, columns=["PC1", "PC2"], index=bray_matrix.index)

def map_taxonomy(otu_names):
    return [TAXONOMY_MAP.get(name, "Unknown") for name in otu_names]

def taxonomy_abundance(rel_df, taxonomy_list, otu_names):
    tax_abundance = {}
    for otu_name, tax_name in zip(otu_names, taxonomy_list):
        mean_abundance = rel_df[otu_name].mean()
        tax_abundance[tax_name] = tax_abundance.get(tax_name, 0) + mean_abundance
    return pd.Series(tax_abundance).sort_values(ascending=False)

def predict_function(taxonomy_names):
    unique_tax = set(taxonomy_names)
    functions = []
    for name in unique_tax:
        fn = FUNCTION_MAP.get(name, "Unknown Function")
        if fn not in functions:
            functions.append(fn)
    return functions

@st.cache_resource
def get_db_connection():
    return mysql.connector.connect(
        host=MYSQL_HOST,
        port=MYSQL_PORT,
        user=MYSQL_USER,
        password=MYSQL_PASSWORD,
        database=MYSQL_DATABASE
    )

def log_to_mysql(alpha_df, beta_matrix):
    try:
        mydb = get_db_connection()
        cursor = mydb.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS microbiome_samples (
                id INT AUTO_INCREMENT PRIMARY KEY,
                sample_id VARCHAR(100),
                observed_otus INT,
                shannon FLOAT,
                simpson FLOAT,
                chao1 FLOAT,
                evenness FLOAT,
                logged_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS beta_diversity (
                id INT AUTO_INCREMENT PRIMARY KEY,
                sample_a VARCHAR(100),
                sample_b VARCHAR(100),
                bray_curtis FLOAT,
                logged_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        for _, row in alpha_df.iterrows():
            cursor.execute("""
                INSERT INTO microbiome_samples
                (sample_id, observed_otus, shannon, simpson, chao1, evenness)
                VALUES (%s, %s, %s, %s, %s, %s)
            """, (
                row["sample_id"], int(row["observed_otus"]),
                float(row["shannon"]), float(row["simpson"]),
                float(row["chao1"]), float(row["evenness"])
            ))
        samples = beta_matrix.index.tolist()
        for i, a in enumerate(samples):
            for j, b in enumerate(samples):
                if i < j:
                    cursor.execute("""
                        INSERT INTO beta_diversity
                        (sample_a, sample_b, bray_curtis)
                        VALUES (%s, %s, %s)
                    """, (a, b, float(beta_matrix.loc[a, b])))
        mydb.commit()
        cursor.close()
    except Exception as e:
        st.error(f"Database Error: {e}")

st.title("Microbiome Diversity Dashboard")

if "results" not in st.session_state:
    st.session_state.results = None

st.subheader("Input")

input_mode = st.radio("Data source:", ["Use Synthetic Data", "Upload CSV"], horizontal=True)

uploaded_file = None
num_samples = 20
num_taxa = 12

if input_mode == "Upload CSV":
    uploaded_file = st.file_uploader("Upload OTU table (CSV)", type=["csv"])
    st.caption("Format: rows = samples, columns = taxa, values = read counts.")
else:
    col1, col2 = st.columns(2)
    with col1:
        num_samples = st.number_input("Number of Samples:", min_value=5, max_value=100, value=20)
    with col2:
        num_taxa = st.number_input("Number of Taxa:", min_value=5, max_value=12, value=12)

run_button = st.button("Run Analysis")

if run_button:
    if input_mode == "Upload CSV":
        if uploaded_file is None:
            st.error("Please upload a CSV file.")
            st.stop()
        raw_df = pd.read_csv(uploaded_file, index_col=0)
    else:
        raw_df = generate_otu_table(num_samples, num_taxa)

    if raw_df.empty:
        st.error("Empty OTU table.")
        st.stop()

    cleaned = clean_otu_table(raw_df)

    if cleaned.empty:
        st.error("No valid data after cleaning.")
        st.stop()

    rel_df = to_relative_abundance(cleaned)

    with st.spinner("Computing alpha diversity..."):
        alpha_df = alpha_diversity(cleaned, rel_df)

    with st.spinner("Computing beta diversity..."):
        bray = bray_curtis_matrix(rel_df)
        jaccard = jaccard_matrix(cleaned)

    with st.spinner("Running PCoA..."):
        pcoa_df = pcoa_projection(bray)

    otu_names = cleaned.columns.tolist()
    taxonomy_names = map_taxonomy(otu_names)
    tax_series = taxonomy_abundance(rel_df, taxonomy_names, otu_names)
    functions = predict_function(taxonomy_names)

    log_to_mysql(alpha_df, bray)

    st.session_state.results = {
        "raw_df": raw_df,
        "cleaned": cleaned,
        "rel_df": rel_df,
        "alpha_df": alpha_df,
        "bray": bray,
        "jaccard": jaccard,
        "pcoa_df": pcoa_df,
        "taxonomy": taxonomy_names,
        "tax_series": tax_series,
        "functions": functions
    }

if st.session_state.results is not None:
    res = st.session_state.results

    st.divider()
    st.subheader("Analysis Results")

    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        ["Alpha Diversity", "Beta Diversity", "PCoA", "Taxonomy", "Function"]
    )

    with tab1:
        st.dataframe(res["alpha_df"])

        fig1, ax1 = plt.subplots(figsize=(10, 5))
        ax1.bar(res["alpha_df"]["sample_id"], res["alpha_df"]["shannon"], color="steelblue")
        ax1.set_title("Shannon Index per Sample")
        ax1.set_xlabel("Sample")
        ax1.set_ylabel("Shannon Index")
        plt.xticks(rotation=45, ha="right")
        st.pyplot(fig1)
        plt.close(fig1)

        fig2, ax2 = plt.subplots(figsize=(10, 5))
        ax2.bar(res["alpha_df"]["sample_id"], res["alpha_df"]["observed_otus"], color="darkorange")
        ax2.set_title("Observed OTUs per Sample")
        ax2.set_xlabel("Sample")
        ax2.set_ylabel("Richness")
        plt.xticks(rotation=45, ha="right")
        st.pyplot(fig2)
        plt.close(fig2)

    with tab2:
        st.write("Bray-Curtis Dissimilarity Matrix")
        st.dataframe(res["bray"].round(3))
        st.write("Jaccard Dissimilarity Matrix")
        st.dataframe(res["jaccard"].round(3))

        fig3, ax3 = plt.subplots(figsize=(8, 6))
        im = ax3.imshow(res["bray"].values, cmap="viridis", aspect="auto", vmin=0, vmax=1)
        ax3.set_xticks(range(len(res["bray"].index)))
        ax3.set_yticks(range(len(res["bray"].index)))
        ax3.set_xticklabels(res["bray"].index, rotation=45, ha="right")
        ax3.set_yticklabels(res["bray"].index)
        ax3.set_title("Bray-Curtis Heatmap")
        plt.colorbar(im, ax=ax3)
        st.pyplot(fig3)
        plt.close(fig3)

    with tab3:
        pcoa = res["pcoa_df"]
        fig4, ax4 = plt.subplots(figsize=(8, 6))
        ax4.scatter(pcoa["PC1"], pcoa["PC2"], s=80, c="teal", alpha=0.7)
        for sample in pcoa.index:
            ax4.annotate(sample, (pcoa.loc[sample, "PC1"], pcoa.loc[sample, "PC2"]),
                         fontsize=8, alpha=0.7)
        ax4.set_title("PCoA Projection (Bray-Curtis)")
        ax4.set_xlabel("PC1")
        ax4.set_ylabel("PC2")
        ax4.grid(True, alpha=0.3)
        st.pyplot(fig4)
        plt.close(fig4)

    with tab4:
        tax_df = pd.DataFrame({
            "OTU": res["cleaned"].columns.tolist(),
            "Taxonomy": res["taxonomy"]
        })
        st.dataframe(tax_df)

        fig5, ax5 = plt.subplots(figsize=(8, 6))
        tax_series = res["tax_series"]
        ax5.pie(tax_series.values, labels=tax_series.index, autopct='%1.1f%%')
        ax5.set_title("Mean Relative Abundance by Taxon")
        st.pyplot(fig5)
        plt.close(fig5)

    with tab5:
        st.write("Predicted Metabolic Functions")
        for fn in res["functions"]:
            st.write(f"• {fn}")