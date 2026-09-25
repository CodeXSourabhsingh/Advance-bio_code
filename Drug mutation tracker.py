import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mysql.connector
import requests
import urllib.request
from Bio.PDB import PDBParser
from io import StringIO
import py3Dmol
import streamlit.components.v1 as components

try:
    MYSQL_HOST = st.secrets["MYSQL_HOST"]
    MYSQL_PORT = st.secrets["MYSQL_PORT"]
    MYSQL_USER = st.secrets["MYSQL_USER"]
    MYSQL_PASSWORD = st.secrets["MYSQL_PASSWORD"]
    MYSQL_DATABASE = st.secrets["MYSQL_DATABASE"]
except Exception:
    from config import MYSQL_HOST, MYSQL_USER, MYSQL_PASSWORD, MYSQL_DATABASE
    MYSQL_PORT = 3306

ACTIVITY_TYPE_MAP = {
    "IC50": "Inhibition",
    "Ki": "Inhibition",
    "Kd": "Inhibition",
    "EC50": "Activation",
    "AC50": "Activation",
    "Potency": "Modulation",
    "Log IC50": "Inhibition",
    "pIC50": "Inhibition",
    "Ratio": "Modulation",
    "Inhibition": "Inhibition",
    "% Inhibition": "Inhibition",
    "Kb": "Inhibition"
}

@st.cache_data
def search_drug(drug_name):
    try:
        url = "https://www.ebi.ac.uk/chembl/api/data/molecule/search"
        params = {"q": drug_name, "format": "json"}
        response = requests.get(url, params=params, timeout=30)
        if response.status_code != 200:
            raise Exception(f"API error: {response.status_code}")
        data = response.json()
        molecules = data.get("molecules", [])
        if not molecules:
            return None
        return molecules[0].get("molecule_chembl_id")
    except requests.exceptions.Timeout:
        st.error("ChEMBL search timed out. Try again.")
        return None
    except Exception as e:
        st.error(f"Search error: {e}")
        return None

@st.cache_data
def fetch_activities(chembl_id):
    try:
        url = "https://www.ebi.ac.uk/chembl/api/data/activity"
        params = {
            "molecule_chembl_id": chembl_id,
            "format": "json",
            "limit": 100
        }
        response = requests.get(url, params=params, timeout=30)
        if response.status_code != 200:
            raise Exception(f"API error: {response.status_code}")
        data = response.json()
        activities = data.get("activities", [])
        parsed = []
        for act in activities:
            parsed.append({
                "target_name": act.get("target_pref_name", "Unknown"),
                "target_chembl_id": act.get("target_chembl_id"),
                "standard_type": act.get("standard_type", "Unknown"),
                "standard_value": act.get("standard_value"),
                "standard_units": act.get("standard_units"),
                "assay_description": act.get("assay_description", "")
            })
        return parsed
    except requests.exceptions.Timeout:
        st.error("ChEMBL activity fetch timed out. Try again.")
        return []
    except Exception as e:
        st.error(f"Activity fetch error: {e}")
        return []

@st.cache_data
def get_uniprot_from_target(target_chembl_id):
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/target/{target_chembl_id}?format=json"
        response = requests.get(url, timeout=30)
        if response.status_code != 200:
            return None
        data = response.json()
        components_list = data.get("target_components", [])
        for comp in components_list:
            accession = comp.get("accession")
            if accession:
                return accession
        return None
    except Exception:
        return None

@st.cache_data
def uniprot_to_pdb(uniprot_id):
    try:
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": uniprot_id
                }
            },
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": 1}}
        }
        response = requests.post(url, json=query, timeout=30)
        if response.status_code != 200:
            return None
        data = response.json()
        results = data.get("result_set", [])
        if results:
            return results[0].get("identifier")
        return None
    except Exception:
        return None

@st.cache_data
def fetch_pdb(pdb_id):
    try:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = urllib.request.urlopen(url, timeout=30)
        raw_text = response.read().decode("utf-8")
        return raw_text
    except Exception as e:
        st.error(f"PDB fetch error for {pdb_id}: {e}")
        return None

def convert_to_nM(value, unit):
    if value is None or unit is None:
        return None
    try:
        v = float(value)
    except (ValueError, TypeError):
        return None

    unit_lower = str(unit).lower().strip()

    if "micromolar" in unit_lower or unit_lower == "um" or unit_lower == "µm":
        return v * 1000
    elif "millimolar" in unit_lower or unit_lower == "mm":
        return v * 1000000
    elif "nanomolar" in unit_lower or unit_lower == "nm":
        return v
    elif "picomolar" in unit_lower or unit_lower == "pm":
        return v / 1000
    elif "molar" in unit_lower and "micro" not in unit_lower and "milli" not in unit_lower and "nano" not in unit_lower and "pico" not in unit_lower:
        return v * 1000000000
    else:
        return None

def classify_strength(ic50_values):
    conditions = [
        ic50_values < 10,
        ic50_values < 100,
        ic50_values < 1000,
        ic50_values < 10000
    ]
    choices = ["Very Strong", "Strong", "Moderate", "Weak"]
    return np.select(conditions, choices, default="Inactive")

def map_activity_type(standard_type):
    if standard_type is None:
        return "Unknown"
    return ACTIVITY_TYPE_MAP.get(str(standard_type).strip(), "Unknown")

@st.cache_resource
def get_db_connection():
    return mysql.connector.connect(
        host=MYSQL_HOST,
        port=MYSQL_PORT,
        user=MYSQL_USER,
        password=MYSQL_PASSWORD,
        database=MYSQL_DATABASE
    )

def log_to_mysql(drug_name, chembl_id, df):
    try:
        mydb = get_db_connection()
        cursor = mydb.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS drug_targets (
                id INT AUTO_INCREMENT PRIMARY KEY,
                drug_name VARCHAR(255),
                chembl_id VARCHAR(50),
                target_name VARCHAR(500),
                ic50_nm FLOAT,
                activity_type VARCHAR(50),
                binding_strength VARCHAR(50),
                logged_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        insert_query = """
            INSERT INTO drug_targets
            (drug_name, chembl_id, target_name, ic50_nm, activity_type, binding_strength)
            VALUES (%s, %s, %s, %s, %s, %s)
        """
        for _, row in df.iterrows():
            cursor.execute(insert_query, (
                drug_name,
                chembl_id,
                str(row.get("target_name", "Unknown"))[:500],
                float(row["ic50_nm"]) if pd.notna(row["ic50_nm"]) else None,
                str(row.get("activity_type", "Unknown")),
                str(row.get("binding_strength", "Unknown"))
            ))
        mydb.commit()
        cursor.close()
    except Exception as e:
        st.error(f"Database Error: {e}")

st.title("Drug-Target Interaction Explorer")

if "results" not in st.session_state:
    st.session_state.results = None

st.subheader("Input")
drug_name = st.text_input("Enter drug name:", "Imatinib")
run_button = st.button("Fetch and Analyze")

if run_button:
    with st.spinner(f"Searching ChEMBL for {drug_name}..."):
        chembl_id = search_drug(drug_name)

    if chembl_id is None:
        st.error(f"Drug '{drug_name}' not found in ChEMBL.")
        st.stop()

    st.success(f"Found ChEMBL ID: {chembl_id}")

    with st.spinner("Fetching activity records..."):
        activities = fetch_activities(chembl_id)

    if not activities:
        st.error("No activity data found for this drug.")
        st.stop()

    df = pd.DataFrame(activities)
    df = df.dropna(subset=["standard_value", "standard_units"])

    if df.empty:
        st.error("No valid activity records with values and units.")
        st.stop()

    df["ic50_nm"] = df.apply(
        lambda row: convert_to_nM(row["standard_value"], row["standard_units"]),
        axis=1
    )
    df = df.dropna(subset=["ic50_nm"])

    if df.empty:
        st.error("No convertible IC50 values found.")
        st.stop()

    df["binding_strength"] = classify_strength(df["ic50_nm"].values)
    df["activity_type"] = df["standard_type"].apply(map_activity_type)

    log_to_mysql(drug_name, chembl_id, df)

    pdb_id = None
    st.info("Searching for 3D structure across targets...")

    for tid in df["target_chembl_id"].dropna().unique()[:5]:
        uniprot_id = get_uniprot_from_target(tid)
        if uniprot_id:
            found_pdb = uniprot_to_pdb(uniprot_id)
            if found_pdb:
                pdb_id = found_pdb
                st.info(f"Found PDB ID {pdb_id}")
                break

    if not pdb_id:
        st.warning("No PDB structure found for any target of this drug.")

    st.session_state.results = {
        "drug_name": drug_name,
        "chembl_id": chembl_id,
        "df": df,
        "pdb_id": pdb_id
    }

if st.session_state.results is not None:
    res = st.session_state.results
    df = res["df"]

    st.divider()
    st.subheader("Analysis Results")

    tab1, tab2, tab3 = st.tabs(["Interactions", "Summary", "Charts"])

    with tab1:
        display_cols = ["target_name", "ic50_nm", "activity_type", "binding_strength"]
        st.dataframe(df[display_cols])

    with tab2:
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Total Targets", len(df))
        with col2:
            very_strong = len(df[df["binding_strength"] == "Very Strong"])
            st.metric("Very Strong Binders", very_strong)
        with col3:
            inhibition_count = len(df[df["activity_type"] == "Inhibition"])
            st.metric("Inhibition Count", inhibition_count)

    with tab3:
        if not df.empty:
            top_targets = df.nlargest(10, "ic50_nm").sort_values("ic50_nm")

            fig1, ax1 = plt.subplots(figsize=(8, 5))
            ax1.barh(top_targets["target_name"].astype(str), top_targets["ic50_nm"], color="steelblue")
            ax1.set_xscale("log")
            ax1.set_xlabel("IC50 (nM) - log scale")
            ax1.set_title("Top 10 Targets by IC50")
            st.pyplot(fig1)
            plt.close(fig1)

            strength_counts = df["binding_strength"].value_counts()

            fig2, ax2 = plt.subplots(figsize=(6, 6))
            ax2.pie(strength_counts.values, labels=strength_counts.index, autopct='%1.1f%%')
            ax2.set_title("Binding Strength Distribution")
            st.pyplot(fig2)
            plt.close(fig2)

            activity_counts = df["activity_type"].value_counts()

            fig3, ax3 = plt.subplots(figsize=(6, 4))
            ax3.bar(activity_counts.index, activity_counts.values, color="darkorange")
            ax3.set_title("Activity Type Distribution")
            ax3.set_xlabel("Activity Type")
            ax3.set_ylabel("Count")
            st.pyplot(fig3)
            plt.close(fig3)

        if res["pdb_id"]:
            st.subheader("3D Structure View")
            with st.spinner(f"Loading {res['pdb_id']}..."):
                pdb_text = fetch_pdb(res["pdb_id"])

            if pdb_text:
                try:
                    view = py3Dmol.view(width=700, height=500)
                    view.addModel(pdb_text, "pdb")
                    view.setStyle({"cartoon": {"color": "spectrum"}})
                    view.addLabel(
                        f"{res['drug_name']} → {res['pdb_id']}",
                        {
                            "fontSize": 14,
                            "fontColor": "black",
                            "backgroundColor": "white",
                            "backgroundOpacity": 0.7
                        }
                    )
                    view.zoomTo()
                    components.html(view._make_html(), height=520)
                except Exception as e:
                    st.error(f"3D render error: {e}")
            else:
                st.warning("PDB structure could not be loaded.")
        else:
            st.info("No 3D structure available for this drug's primary target.")