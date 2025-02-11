from admet_ai import ADMETModel
import math
import matplotlib.pyplot as plt
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
from io import StringIO

model = ADMETModel()
# Function to predict ADMET properties
def admet_from_smiles(smiles):
    preds = model.predict(smiles=smiles)
    return preds

admet_from_smiles("C1CCCC1")

upper_limits = {
    "molecular_weight": 600,
    "logP": 3,
    "hydrogen_bond_acceptors": 12,
    "hydrogen_bond_donors": 7,
    "Lipinski": 5,
    "QED": 1,
    "stereo_centers": 10,
    "tpsa": 140
}

st.set_page_config(layout="wide")
admet_tab, virutal_screening, de_novo = st.tabs(["ADMET", "Virtual Screening", "De Novo"])
data = pd.read_csv("consolidated_output.csv")
data = data.dropna()

def add_admet_properties_to_data(ligand_data):
    admet_results = []

    for smiles in ligand_data["Ligand SMILES"]:
        try:
            admet_data = admet_from_smiles(smiles)  
            if admet_data and isinstance(admet_data, dict):
                admet_results.append(admet_data)
            else:
                admet_results.append(None)
        except:
            continue
        
    admet_df = pd.DataFrame(admet_results)

    ligand_data = ligand_data.reset_index(drop=True)

    ligand_data = pd.concat([ligand_data[["Ligand SMILES", "Protein", "pIC50"]], admet_df], axis=1)
    
    return ligand_data

with admet_tab:
    col1, col2 = st.columns(2)

    with col1:
        sequence = st.text_input("SMILES Sequence", "CCCC")
        st.write("ADMET Results: ")
        results = admet_from_smiles(sequence)

        st.dataframe(results, use_container_width=True)

    with col2:
        m = Chem.MolFromSmiles(sequence)
        st.image(Draw.MolToImage(m))

        st.write("Radar Chart:")
        # fig = radar_chart(results, upper_limits)
        # st.pyplot(fig)

with de_novo:
    st.write("De Novo")    

with virutal_screening:
    st.write("## Virtual Screening")

    proteins = data["Protein"].unique()
    
    selected_protein = st.selectbox("Select Protein", proteins)

    # Initialize session state for filtered_ligands if not already set
    if "filtered_ligands" not in st.session_state or st.session_state.get("last_selected_protein") != selected_protein:
        st.session_state.filtered_ligands = data[data["Protein"] == selected_protein].drop(
            columns=["ProteinID"] + ["rpssm" + str(x) for x in range(110)] + ["IC50 (nM)"]
        )
        st.session_state.last_selected_protein = selected_protein

    st.write("May take up to " + str(math.ceil((len(st.session_state.filtered_ligands) * 0.1) / 60)) + " minute(s) to complete.")

    if st.button("Run ADMET Predictions"):
        st.session_state.filtered_ligands = add_admet_properties_to_data(st.session_state.filtered_ligands)
        
        st.write("Ligands with ADMET Properties")
        st.dataframe(st.session_state.filtered_ligands, use_container_width=True)

        st.download_button(
            label="Download Updated Data with ADMET",
            data=st.session_state.filtered_ligands.to_csv(index=False),
            file_name="updated_ligands_with_admet.csv",
            mime="text/csv"
        )

    # Ensure filtering works persistently
    filtered_ligands = st.session_state.filtered_ligands.copy()
    
    numeric_columns = filtered_ligands.select_dtypes(include=['number']).columns
    filters = {}

    for column in numeric_columns:
        min_val, max_val = filtered_ligands[column].min(), filtered_ligands[column].max()
        if min_val == max_val:
            max_val = min_val + 0.1

        # Store slider values in session state to persist across reruns
        if f"{column}_slider" not in st.session_state:
            st.session_state[f"{column}_slider"] = (float(min_val), float(max_val))

        filters[column] = st.slider(
            f"Select range for {column}", 
            min_value=float(min_val), 
            max_value=float(max_val), 
            value=st.session_state[f"{column}_slider"], 
            step=0.1
        )

        st.session_state[f"{column}_slider"] = filters[column]  # Store updated values

    # Apply filters based on stored values
    for column, range_values in filters.items():
        min_value, max_value = range_values
        filtered_ligands = filtered_ligands[(filtered_ligands[column] >= min_value) & (filtered_ligands[column] <= max_value)]

    # Display the filtered ligands
    st.write(f"### Ligands for {selected_protein}")
    st.dataframe(filtered_ligands[["Ligand SMILES", "pIC50"]])

    # Allow user to select a ligand
    selected_ligand_smiles = st.selectbox("Select Ligand (SMILES)", filtered_ligands["Ligand SMILES"])
    
    # Display the molecule structure
    selected_ligand = Chem.MolFromSmiles(selected_ligand_smiles)
    if selected_ligand:
        st.image(Draw.MolToImage(selected_ligand), caption="Selected Ligand Structure")