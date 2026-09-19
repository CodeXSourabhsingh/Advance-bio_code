import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
import requests
import random
import numpy as np
import mysql.connector
try:
    MYSQL_HOST = st.secrets["MYSQL_HOST"]
    MYSQL_PORT = st.secrets["MYSQL_PORT"]
    MYSQL_USER = st.secrets["MYSQL_USER"]
    MYSQL_PASSWORD = st.secrets["MYSQL_PASSWORD"]
    MYSQL_DATABASE = st.secrets["MYSQL_DATABASE"]
except Exception:
    from config import MYSQL_HOST, MYSQL_USER, MYSQL_PASSWORD, MYSQL_DATABASE
    MYSQL_PORT = 3306

@st.cache_resource
def get_db_connection():
    return mysql.connector.connect(
        host=MYSQL_HOST,
        port=MYSQL_PORT,
        user=MYSQL_USER,
        password=MYSQL_PASSWORD,
        database=MYSQL_DATABASE,
        connection_timeout=3
    )
class Patient:
    def __init__(self, id, age, sex, disease_state, cancer_stage, crispr_mutation, biomarker_value, is_outliers=False):
        self.id = id
        self.age = age
        self.sex = sex
        self.disease_state = disease_state
        self.cancer_stage = cancer_stage
        self.crispr_mutation = crispr_mutation
        self.biomarker_value = biomarker_value
        self.is_outliers = is_outliers

@st.cache_data
def generate_synthetic_data(num_patients):
    patient_list = []
    diseases = ['NSCLC','Breast Cancer','Melanoma','Diabetes',"Alzheimer's",
    "Parkinson's",
    'Cardiovascular Disease',
    'Chronic Kidney Disease',
    'Leukemia',
    'Rheumatoid Arthritis'
    ]
   
    stages = ['1', '2', '3', '4']
    mutations = ['EGFR Exon 19 Del', 'ALK Fusion', 'KRAS G12C', 'None']
    sexes = ['Male', 'Female']

    for i in range(num_patients):
        p = Patient(
            id=i + 1,
            age=random.randint(30, 85),
            sex=random.choice(sexes),
            disease_state=random.choice(diseases),
            cancer_stage=random.choice(stages),
            crispr_mutation=random.choice(mutations),
                        biomarker_value=(
                                round(random.gauss(150, 20), 2)
                                if random.random() < 0.95
                                else round(random.choice([40, 300]), 2)
                        ),
        )
        patient_list.append(p)

    return pd.DataFrame([vars(p) for p in patient_list])

@st.cache_data(ttl=3600)
def fetch_clinical_trials(condition, page_size=10):
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.cond": condition,
        "pageSize": page_size,
        "fields": "NCTId,BriefTitle,Phase,OverallStatus,EligibilityCriteria"
    }
    try:
        response = requests.get(url, params=params, timeout= 20)
        if response.status_code != 200:
            raise Exception(f"API returned status {response.status_code}")
        data = response.json()
        studies = data.get("studies", [])
        if not studies:
            raise Exception(f"No trials found for condition: {condition}")

        trials = []
        for study in studies:
            protocol = study.get("protocolSection", {})
            identification = protocol.get("identificationModule", {})
            status = protocol.get("statusModule", {})
            eligibility = protocol.get("eligibilityModule", {})
            design = protocol.get("designModule", {})
            phases = design.get("phases", ["Unknown"])

            trials.append({
                "nct_id": identification.get("nctId", "Unknown"),
                "title": identification.get("briefTitle", "Unknown"),
                "phase": phases[0] if phases else "Unknown",
                "status": status.get("overallStatus", "Unknown"),
                "eligibility": eligibility.get("eligibilityCriteria", "No criteria provided")
            })
        return pd.DataFrame(trials)

    except requests.exceptions.Timeout:
        st.error("API request timed out. Try again.")
        return None
    except requests.exceptions.RequestException as e:
        st.error(f"API connection error: {e}")
        return None
    except Exception as e:
        st.error(f"Error fetching trials: {e}")
        return None

st.title("PRISM: Advanced Clinical Pipeline")

if "filtered_cohort" not in st.session_state:
    st.session_state.filtered_cohort = None
if "summary_stats" not in st.session_state:
    st.session_state.summary_stats = None
if "outliers" not in st.session_state:
    st.session_state.outliers = None
if "trials_df" not in st.session_state:
    st.session_state.trials_df = None

st.subheader("Clinical Trial Fetcher (ClinicalTrials.gov)")
condition = st.text_input("Enter Condition (e.g., NSCLC)", "NSCLC")
fetch_button = st.button("Fetch Trials")

if fetch_button:
    with st.spinner(f"Fetching trials for {condition}..."):
        trials_df = fetch_clinical_trials(condition)
        if trials_df is not None:
            st.session_state.trials_df = trials_df
            st.success(f"Fetched {len(trials_df)} trials for {condition}.")

if st.session_state.trials_df is not None:
    selected_nct = st.selectbox(
        "Select a Trial",
        st.session_state.trials_df["nct_id"].tolist()
    )
    selected_trial = st.session_state.trials_df[
        st.session_state.trials_df["nct_id"] == selected_nct
    ].iloc[0]
    st.write(f"**{selected_trial['title']}**")
    st.write(f"Phase: {selected_trial['phase']} | Status: {selected_trial['status']}")
    with st.expander("View Eligibility Criteria"):
        st.write(selected_trial["eligibility"])

st.subheader("Cohort Filters")
selected_diseases = st.multiselect("Select Disease", ['NSCLC', 'Breast Cancer', 'Melanoma', 'Diabetes', "Alzheimer's",
    "Parkinson's",
    'Cardiovascular Disease',
    'Chronic Kidney Disease',
    'Leukemia',
    'Rheumatoid Arthritis'], default=['NSCLC']) 

selected_stages = st.multiselect("Select Stage", ['1', '2', '3', '4'], default=['1', '2'])
min_age, max_age = st.slider("Age Range", 18, 90, (40, 70))
num_patients = st.number_input("Number of Patients to Simulate", min_value=100, max_value=5000, value=2000, step=100)

st.subheader("CRISPR Integration (Optional)")
use_crispr = st.checkbox("Filter by CRISPR Mutation")
crispr_input = ""
if use_crispr:
    crispr_input = st.text_input("Enter Mutation (from CRISPR Simulator)", "EGFR Exon 19 Del")

run_button = st.button("Run Pipeline")

if run_button:
    df = generate_synthetic_data(num_patients)
    df = df.reset_index(drop=True)
    cohort = df.loc[
        (df['disease_state'].isin(selected_diseases)) &
        (df['cancer_stage'].isin(selected_stages)) &
        (df['age'].between(min_age, max_age))
    ].copy()
    if use_crispr and crispr_input:
        cohort = cohort[cohort['crispr_mutation'] == crispr_input]

    if not cohort.empty:
        mean_val = cohort['biomarker_value'].mean()
        std_val = cohort['biomarker_value'].std()
        lower_bound = mean_val - (2 * std_val)
        upper_bound = mean_val + (2 * std_val)
        outliers = cohort[(cohort['biomarker_value'] < lower_bound) | (cohort['biomarker_value'] > upper_bound)]
        st.session_state.filtered_cohort = cohort
        st.session_state.summary_stats = {
            'mean': mean_val,
            'std': std_val,
            'lower': lower_bound,
            'upper': upper_bound
        }
        st.session_state.outliers = outliers

    if not cohort.empty:
        try:
            conn = get_db_connection()
            cursor = conn.cursor()
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS clinical_pipeline_vault (
                    id INT AUTO_INCREMENT PRIMARY KEY,
                    disease VARCHAR(255),
                    cohort_size INT,
                    mean_biomarker FLOAT,
                    std_dev FLOAT,
                    outlier_count INT,
                    logged_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            insert_query = """
                INSERT INTO clinical_pipeline_vault
                (disease, cohort_size, mean_biomarker, std_dev, outlier_count)
                VALUES (%s, %s, %s, %s, %s)
            """
            disease_str = ", ".join(selected_diseases) if selected_diseases else "All"
            values = (disease_str, len(cohort), mean_val, std_val, len(outliers))
            cursor.execute(insert_query, values)
            conn.commit()
            cursor.close()
        except Exception as e:
            st.error(f"Database Error: {e}")
    else:
        st.warning("No patients match the selected criteria. Please broaden your filter.")


if st.session_state.filtered_cohort is not None:
    st.divider()
    st.subheader("PRISM Dashboard")

    tab1, tab2, tab3 = st.tabs(["Cohort Explorer", "Baseline Summary", "Safety Alerts"])

    with tab1:
        st.dataframe(st.session_state.filtered_cohort)
        fig4, ax4 = plt.subplots(figsize=(6, 6))
        st.session_state.filtered_cohort.groupby('cancer_stage').size().plot(kind='pie', autopct='%1.1f%%', ax=ax4)
        ax4.set_title("Cohort by Cancer Stage")
        ax4.set_ylabel("")
        st.pyplot(fig4)
        plt.close(fig4)

    with tab2:
        stats = st.session_state.summary_stats
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Cohort Size", len(st.session_state.filtered_cohort))
        with col2:
            st.metric("Mean Biomarker", f"{stats['mean']:.2f}")
        with col3:
            st.metric("Normal Range", f"{stats['lower']:.2f} - {stats['upper']:.2f}")

        fig2, ax2 = plt.subplots(figsize=(8, 4))
        ax2.hist(st.session_state.filtered_cohort['biomarker_value'], bins=20, color='steelblue', edgecolor='black')
        ax2.axvline(stats['mean'], color='red', linestyle='--', label=f"Mean: {stats['mean']:.2f}")
        ax2.axvline(stats['lower'], color='orange', linestyle=':', label=f"Lower: {stats['lower']:.2f}")
        ax2.axvline(stats['upper'], color='orange', linestyle=':', label=f"Upper: {stats['upper']:.2f}")
        ax2.set_title("Biomarker Distribution")
        ax2.set_xlabel("Biomarker Value")
        ax2.set_ylabel("Patient Count")
        ax2.legend()
        st.pyplot(fig2)
        plt.close(fig2)

        fig3, ax3 = plt.subplots(figsize=(8, 4))
        ax3.scatter(st.session_state.filtered_cohort['age'], st.session_state.filtered_cohort['biomarker_value'], alpha=0.5, color='teal')
        ax3.set_title("Age vs Biomarker Value")
        ax3.set_xlabel("Age")
        ax3.set_ylabel("Biomarker Value")
        st.pyplot(fig3)
        plt.close(fig3)

    with tab3:
        outliers = st.session_state.outliers
        if outliers is None:
            st.warning("No outlier data is available. Please run the pipeline first.")
        else:
            st.dataframe(outliers)

            fig, ax = plt.subplots(figsize=(6, 3))
            normal_count = len(st.session_state.filtered_cohort) - len(outliers)
            outlier_count = len(outliers)
            ax.bar(["Normal", "Outlier"], [normal_count, outlier_count], color=["blue", "red"])
            ax.set_title("Cohort Safety Profile")
            ax.set_ylabel("Patient Count")
            st.pyplot(fig)
            plt.close(fig)

