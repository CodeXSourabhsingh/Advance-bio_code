import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mysql.connector
import urllib.request
import py3Dmol
import streamlit.components.v1 as components
from builtins import ValueError
from io import StringIO
from Bio.PDB import PDBParser
from config import MYSQL_HOST, MYSQL_USER, MYSQL_PASSWORD, MYSQL_DATABASE

HYDROPHOBICITY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

CHARGE = {
    "A": "neutral", "R": "positive", "N": "neutral", "D": "negative", "C": "neutral",
    "Q": "neutral", "E": "negative", "G": "neutral", "H": "positive", "I": "neutral",
    "L": "neutral", "K": "positive", "M": "neutral", "F": "neutral", "P": "neutral",
    "S": "neutral", "T": "neutral", "W": "neutral", "Y": "neutral", "V": "neutral",
}

MOLECULAR_WEIGHT = {
    "A": 89.1, "R": 174.2, "N": 132.1, "D": 133.1, "C": 121.2,
    "Q": 146.2, "E": 147.1, "G": 75.1, "H": 155.2, "I": 131.2,
    "L": 131.2, "K": 146.2, "M": 149.2, "F": 165.2, "P": 115.1,
    "S": 105.1, "T": 119.1, "W": 204.2, "Y": 181.2, "V": 117.1,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

BLOSUM62 = {
    "A": {"A": 4, "R": -1, "N": -2, "D": -2, "C": 0, "Q": -1, "E": -1, "G": 0, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1, "T": 0, "W": -3, "Y": -2, "V": 0},
    "R": {"A": -1, "R": 5, "N": 0, "D": -2, "C": -3, "Q": 1, "E": 0, "G": -2, "H": 0, "I": -3, "L": -2, "K": 2, "M": -1, "F": -3, "P": -2, "S": -1, "T": -1, "W": -3, "Y": -2, "V": -3},
    "N": {"A": -2, "R": 0, "N": 6, "D": 1, "C": -3, "Q": 0, "E": 0, "G": 0, "H": 1, "I": -3, "L": -3, "K": 0, "M": -2, "F": -3, "P": -2, "S": 1, "T": 0, "W": -4, "Y": -2, "V": -3},
    "D": {"A": -2, "R": -2, "N": 1, "D": 6, "C": -3, "Q": 0, "E": 2, "G": -1, "H": -1, "I": -3, "L": -4, "K": -1, "M": -3, "F": -3, "P": -1, "S": 0, "T": -1, "W": -4, "Y": -3, "V": -3},
    "C": {"A": 0, "R": -3, "N": -3, "D": -3, "C": 9, "Q": -3, "E": -4, "G": -3, "H": -3, "I": -1, "L": -1, "K": -3, "M": -1, "F": -2, "P": -3, "S": -1, "T": -1, "W": -2, "Y": -2, "V": -1},
    "Q": {"A": -1, "R": 1, "N": 0, "D": 0, "C": -3, "Q": 5, "E": 2, "G": -2, "H": 0, "I": -3, "L": -2, "K": 1, "M": 0, "F": -3, "P": -1, "S": 0, "T": -1, "W": -2, "Y": -1, "V": -2},
    "E": {"A": -1, "R": 0, "N": 0, "D": 2, "C": -4, "Q": 2, "E": 5, "G": -2, "H": 0, "I": -3, "L": -3, "K": 1, "M": -2, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "G": {"A": 0, "R": -2, "N": 0, "D": -1, "C": -3, "Q": -2, "E": -2, "G": 6, "H": -2, "I": -4, "L": -4, "K": -2, "M": -3, "F": -3, "P": -2, "S": 0, "T": -2, "W": -2, "Y": -3, "V": -3},
    "H": {"A": -2, "R": 0, "N": 1, "D": -1, "C": -3, "Q": 0, "E": 0, "G": -2, "H": 8, "I": -3, "L": -3, "K": -1, "M": -2, "F": -1, "P": -2, "S": -1, "T": -2, "W": -2, "Y": 2, "V": -3},
    "I": {"A": -1, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -3, "E": -3, "G": -4, "H": -3, "I": 4, "L": 2, "K": -3, "M": 1, "F": 0, "P": -3, "S": -2, "T": -1, "W": -3, "Y": -1, "V": 3},
    "L": {"A": -1, "R": -2, "N": -3, "D": -4, "C": -1, "Q": -2, "E": -3, "G": -4, "H": -3, "I": 2, "L": 4, "K": -2, "M": 2, "F": 0, "P": -3, "S": -2, "T": -1, "W": -2, "Y": -1, "V": 1},
    "K": {"A": -1, "R": 2, "N": 0, "D": -1, "C": -3, "Q": 1, "E": 1, "G": -2, "H": -1, "I": -3, "L": -2, "K": 5, "M": -1, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "M": {"A": -1, "R": -1, "N": -2, "D": -3, "C": -1, "Q": 0, "E": -2, "G": -3, "H": -2, "I": 1, "L": 2, "K": -1, "M": 5, "F": 0, "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V": 1},
    "F": {"A": -2, "R": -3, "N": -3, "D": -3, "C": -2, "Q": -3, "E": -3, "G": -3, "H": -1, "I": 0, "L": 0, "K": -3, "M": 0, "F": 6, "P": -4, "S": -2, "T": -2, "W": 1, "Y": 3, "V": -1},
    "P": {"A": -1, "R": -2, "N": -2, "D": -1, "C": -3, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -3, "L": -3, "K": -1, "M": -2, "F": -4, "P": 7, "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2},
    "S": {"A": 1, "R": -1, "N": 1, "D": 0, "C": -1, "Q": 0, "E": 0, "G": 0, "H": -1, "I": -2, "L": -2, "K": 0, "M": -1, "F": -2, "P": -1, "S": 4, "T": 1, "W": -3, "Y": -2, "V": -2},
    "T": {"A": 0, "R": -1, "N": 0, "D": -1, "C": -1, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1, "T": 5, "W": -2, "Y": -2, "V": 0},
    "W": {"A": -3, "R": -3, "N": -4, "D": -4, "C": -2, "Q": -2, "E": -3, "G": -2, "H": -2, "I": -3, "L": -2, "K": -3, "M": -1, "F": 1, "P": -4, "S": -3, "T": -2, "W": 11, "Y": 2, "V": -3},
    "Y": {"A": -2, "R": -2, "N": -2, "D": -3, "C": -2, "Q": -1, "E": -2, "G": -3, "H": 2, "I": -1, "L": -1, "K": -2, "M": -1, "F": 3, "P": -3, "S": -2, "T": -2, "W": 2, "Y": 7, "V": -1},
    "V": {"A": 0, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -2, "E": -2, "G": -3, "H": -3, "I": 3, "L": 1, "K": -2, "M": 1, "F": -1, "P": -2, "S": -2, "T": 0, "W": -3, "Y": -1, "V": 4},
}

@st.cache_data
def open_pdb_text(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = urllib.request.urlopen(url)
    return response.read().decode("utf-8")

@st.cache_data
def fetch_structure(pdb_id):
    try:
        raw_text = open_pdb_text(pdb_id)
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, StringIO(raw_text))
        return structure, raw_text
    except Exception as e:
        st.error(f"Error fetching {pdb_id}: {e}")
        return None, None

def parse_secondary_structure(raw_text):
    helices = []
    sheets = []
    for line in raw_text.split("\n"):
        if line.startswith("HELIX"):
            try:
                chain = line[19]
                start = int(line[21:25])
                end = int(line[33:37])
                helices.append((chain, start, end))
            except Exception:
                continue
        elif line.startswith("SHEET"):
            try:
                chain = line[21]
                start = int(line[22:26])
                end = int(line[33:37])
                sheets.append((chain, start, end))
            except Exception:
                continue
    return helices, sheets

def get_secondary_structure(position, chain_id, helices, sheets):
    for ch, start, end in helices:
        if ch == chain_id and start <= position <= end:
            return "Helix"
    for ch, start, end in sheets:
        if ch == chain_id and start <= position <= end:
            return "Sheet"
    return "Loop"

def locate_residue(structure, position, chain_id="A"):
    try:
        chain = structure[0][chain_id]
    except KeyError:
        return None
    for residue in chain:
        if residue.id[0] == " " and residue.id[1] == position:
            return residue
    return None

def get_ca_coordinates(residue):
    if "CA" in residue:
        return residue["CA"].get_coord()
    return None

def find_neighbors(structure, target_coord, chain_id="A", cutoff=5.0):
    neighbors = []
    try:
        chain = structure[0][chain_id]
    except KeyError:
        return neighbors
    for residue in chain:
        if residue.id[0] != " ":
            continue
        if "CA" not in residue:
            continue
        coord = residue["CA"].get_coord()
        distance = float(np.linalg.norm(coord - target_coord))
        if 0.1 < distance < cutoff:
            neighbors.append({
                "residue": THREE_TO_ONE.get(residue.resname, "X"),
                "position": residue.id[1],
                "distance": round(distance, 2)
            })
    return neighbors

def classify_exposure(neighbor_count, cutoff_buried=20):
    if neighbor_count >= cutoff_buried:
        return "Buried"
    return "Exposed"

def score_mutation(ref_aa, alt_aa):
    hydrophobicity_delta = HYDROPHOBICITY[alt_aa] - HYDROPHOBICITY[ref_aa]
    charge_change = CHARGE[ref_aa] != CHARGE[alt_aa]
    size_delta = MOLECULAR_WEIGHT[alt_aa] - MOLECULAR_WEIGHT[ref_aa]
    return {
        "hydrophobicity_delta": round(hydrophobicity_delta, 2),
        "charge_change": charge_change,
        "size_delta": round(size_delta, 2)
    }

def conservation_score(ref_aa, alt_aa):
    score = BLOSUM62.get(ref_aa, {}).get(alt_aa, 0)
    if score <= -3:
        return "Highly Conserved", score
    elif score <= 0:
        return "Conserved", score
    else:
        return "Variable", score

def classify_mutation(scores, exposure, conservation_label):
    h = abs(scores["hydrophobicity_delta"])
    c = scores["charge_change"]
    s = abs(scores["size_delta"])

    if conservation_label == "Highly Conserved" and (h > 2 or c or s > 30):
        return "Destabilizing"

    if exposure == "Buried":
        if h > 3.0:
            return "Destabilizing"
        if c:
            return "Destabilizing"
        if s > 50:
            return "Destabilizing"

    if h > 5.0:
        return "Destabilizing"
    if s > 80:
        return "Destabilizing"

    return "Neutral"

@st.cache_resource
def get_db_connection():
    return mysql.connector.connect(
        host=MYSQL_HOST,
        user=MYSQL_USER,
        password=MYSQL_PASSWORD,
        database=MYSQL_DATABASE
    )

def log_to_mysql(data):
    try:
        mydb = get_db_connection()
        cursor = mydb.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS protein_mutations (
                id INT AUTO_INCREMENT PRIMARY KEY,
                pdb_id VARCHAR(20),
                mutation VARCHAR(20),
                position INT,
                ref_aa VARCHAR(5),
                alt_aa VARCHAR(5),
                secondary_structure VARCHAR(20),
                exposure VARCHAR(20),
                conservation VARCHAR(30),
                blosum_score INT,
                hydrophobicity_delta FLOAT,
                charge_change BOOLEAN,
                size_delta FLOAT,
                classification VARCHAR(30),
                logged_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        cursor.execute("""
            INSERT INTO protein_mutations
            (pdb_id, mutation, position, ref_aa, alt_aa, secondary_structure,
             exposure, conservation, blosum_score, hydrophobicity_delta,
             charge_change, size_delta, classification)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """, (
            data["pdb_id"],
            data["mutation"],
            data["position"],
            data["ref_aa"],
            data["alt_aa"],
            data["secondary_structure"],
            data["exposure"],
            data["conservation"],
            data["blosum_score"],
            data["hydrophobicity_delta"],
            data["charge_change"],
            data["size_delta"],
            data["classification"]
        ))
        mydb.commit()
        cursor.close()
    except Exception as e:
        st.error(f"Database Error: {e}")

st.title("Protein Structure Analyzer")
st.caption("Analyze the structural impact of point mutations.")

if "results" not in st.session_state:
    st.session_state.results = None

st.subheader("Input")
col1, col2 = st.columns(2)
with col1:
    pdb_id = st.text_input("PDB ID (e.g., 1HHO):", "1HHO")
with col2:
    mutation = st.text_input("Mutation (e.g., E6V):", "E6V")

chain_id = st.text_input("Chain ID:", "B")

run_button = st.button("Analyze Mutation")

if run_button:
    if len(mutation) < 3:
        st.error("Mutation must be in format like E6V.")
        st.stop()

    ref_aa = mutation[0].upper()
    alt_aa = mutation[-1].upper()

    try:
        position = int(mutation[1:-1])
    except ValueError:
        st.error("Position must be a number.")
        st.stop()

    if ref_aa not in HYDROPHOBICITY or alt_aa not in HYDROPHOBICITY:
        st.error("Invalid amino acid codes.")
        st.stop()

    with st.spinner(f"Fetching {pdb_id} from RCSB PDB..."):
        structure, raw_text = fetch_structure(pdb_id)

    if structure is None:
        st.stop()

    residue = locate_residue(structure, position, chain_id)
    if residue is None:
        st.error(f"Residue {position} not found in chain {chain_id}.")
        st.stop()

    actual_ref_aa = THREE_TO_ONE.get(residue.resname, "X")
    if actual_ref_aa != ref_aa:
        st.error(f"Mismatch: Position {position} is {actual_ref_aa}, not {ref_aa}.")
        st.stop()

    ca_coord = get_ca_coordinates(residue)
    if ca_coord is None:
        st.error("No CA atom found for this residue.")
        st.stop()

    neighbors = find_neighbors(structure, ca_coord, chain_id, cutoff=5.0)
    exposure = classify_exposure(len(neighbors))

    helices, sheets = parse_secondary_structure(raw_text)
    secondary_structure = get_secondary_structure(position, chain_id, helices, sheets)

    scores = score_mutation(ref_aa, alt_aa)
    conservation_label, blosum_score = conservation_score(ref_aa, alt_aa)
    classification = classify_mutation(scores, exposure, conservation_label)

    result_data = {
        "pdb_id": pdb_id,
        "mutation": mutation.upper(),
        "position": position,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "secondary_structure": secondary_structure,
        "exposure": exposure,
        "conservation": conservation_label,
        "blosum_score": blosum_score,
        "hydrophobicity_delta": scores["hydrophobicity_delta"],
        "charge_change": scores["charge_change"],
        "size_delta": scores["size_delta"],
        "classification": classification
    }

    log_to_mysql(result_data)

    st.session_state.results = {
        "data": result_data,
        "neighbors": neighbors,
        "structure": structure,
        "chain_id": chain_id,
        "raw_text": raw_text
    }

if st.session_state.results is not None:
    res = st.session_state.results
    data = res["data"]
    chain_id = res["chain_id"]
    position = data["position"]

    st.divider()
    st.subheader("Analysis Results")

    tab1, tab2, tab3 = st.tabs(["Summary", "Neighbors", "Charts"])

    with tab1:
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Position", data["position"])
        with col2:
            st.metric("Mutation", f"{data['ref_aa']} → {data['alt_aa']}")
        with col3:
            st.metric("Classification", data["classification"])

        st.write(f"**Secondary Structure:** {data['secondary_structure']}")
        st.write(f"**Exposure:** {data['exposure']}")
        st.write(f"**Conservation:** {data['conservation']} (BLOSUM62: {data['blosum_score']})")
        st.write(f"**Hydrophobicity Delta:** {data['hydrophobicity_delta']}")
        st.write(f"**Charge Change:** {data['charge_change']}")
        st.write(f"**Size Delta:** {data['size_delta']} Da")

    with tab2:
        if res["neighbors"]:
            neighbor_df = pd.DataFrame(res["neighbors"])
            st.dataframe(neighbor_df)
        else:
            st.info("No neighbors within 5 Angstroms.")

    with tab3:
        fig1, ax1 = plt.subplots(figsize=(6, 4))
        labels = ["Hydrophobicity Δ", "Charge Δ", "Size Δ"]
        values = [
            data["hydrophobicity_delta"],
            1 if data["charge_change"] else 0,
            data["size_delta"] / 10
        ]
        colors = ["red" if v < 0 else "green" for v in values]
        ax1.bar(labels, values, color=colors)
        ax1.axhline(0, color="black", linewidth=0.8)
        ax1.set_title("Mutation Impact Scores")
        ax1.set_ylabel("Magnitude")
        st.pyplot(fig1)
        plt.close(fig1)

        structure = res["structure"]
        window = 20
        positions = []
        hydro_values = []

        try:
            chain = structure[0][chain_id]
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                pos = residue.id[1]
                if position - window <= pos <= position + window:
                    one_letter = THREE_TO_ONE.get(residue.resname, "X")
                    if one_letter in HYDROPHOBICITY:
                        positions.append(pos)
                        hydro_values.append(HYDROPHOBICITY[one_letter])
        except Exception:
            pass

        if positions:
            fig2, ax2 = plt.subplots(figsize=(8, 4))
            ax2.plot(positions, hydro_values, marker="o", color="steelblue")
            ax2.axvline(position, color="red", linestyle="--", label=f"Mutation at {position}")
            ax2.set_title("Hydrophobicity Profile (±20 residues)")
            ax2.set_xlabel("Residue Position")
            ax2.set_ylabel("Kyte-Doolittle Value")
            ax2.legend()
            st.pyplot(fig2)
            plt.close(fig2)

        fig3, ax3 = plt.subplots(figsize=(6, 6))
        classification_counts = pd.Series([data["classification"]]).value_counts()
        color_map = {"Destabilizing": "red", "Neutral": "gray", "Stabilizing": "green"}
        chart_colors = [color_map.get(c, "gray") for c in classification_counts.index]
        ax3.pie(classification_counts.values, labels=classification_counts.index,
                autopct='%1.0f%%', colors=chart_colors)
        ax3.set_title("Classification")
        st.pyplot(fig3)
        plt.close(fig3)

        st.subheader("3D Structure View")
        try:
            pdb_text = open_pdb_text(data["pdb_id"])
            view = py3Dmol.view(width=700, height=500)
            view.addModel(pdb_text, "pdb")
            view.setStyle({"cartoon": {"color": "spectrum"}})
            view.addStyle(
                {"chain": chain_id, "resi": position},
                {"stick": {"color": "red", "radius": 0.3}}
            )
            view.addStyle(
                {"chain": chain_id, "resi": position},
                {"sphere": {"color": "red", "radius": 0.8}}
            )
            view.zoomTo({"chain": chain_id, "resi": position})
            view.zoom(0.9)
            components.html(view._make_html(), height=520)
        except Exception as e:
            st.error(f"3D render error: {e}")