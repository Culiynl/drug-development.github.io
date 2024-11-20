from admet_ai import ADMETModel

def admet_from_smiles(smiles):
    model = ADMETModel()
    preds = model.predict(smiles=smiles)

    categories = {
        "Physicochemical": [
            "molecular_weight",
            "logP",
            "hydrogen_bond_acceptors",
            "hydrogen_bond_donors",
            "Lipinski",
            "QED",
            "stereo_centers",
            "tpsa",
        ],
        "Absorption": [
            "HIA_Hou",
            "Bioavailability_Ma",
            "Solubility_AqSolDB",
            "Lipophilicity_AstraZeneca",
            "HydrationFreeEnergy_FreeSolv",
            "Caco2_Wang",
            "PAMPA_NCATS",
            "Pgp_Broccatelli",
        ],
        "Distribution": [
            "BBB_Martins",
            "PPBR_AZ",
            "VDss_Lombardo",
        ],
        "Metabolism": [
            "CYP1A2_Veith",
            "CYP2C19_Veith",
            "CYP2C9_Veith",
            "CYP2D6_Veith",
            "CYP3A4_Veith",
            "CYP2C9_Substrate_CarbonMangels",
            "CYP2D6_Substrate_CarbonMangels",
            "CYP3A4_Substrate_CarbonMangels",
        ],
        "Excretion": [
            "Half_Life_Obach",
            "Clearance_Hepatocyte_AZ",
            "Clearance_Microsome_AZ",
        ],
        "Toxicity": [
            "hERG",
            "ClinTox",
            "AMES",
            "DILI",
            "Carcinogens_Lagunin",
            "LD50_Zhu",
            "Skin_Reaction",
            "NR-AR",
            "NR-AR-LBD",
            "NR-AhR",
            "NR-Aromatase",
            "NR-ER",
            "NR-ER-LBD",
            "NR-PPAR-gamma",
            "SR-ARE",
            "SR-ATAD5",
            "SR-HSE",
            "SR-MMP",
            "SR-p53",
        ],
    }

    organized_predictions = {}
    for category, props in categories.items():
        organized_predictions[category] = {prop: preds[prop] for prop in props if prop in preds}

    return organized_predictions

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

st.set_page_config(layout="wide")
admet_tab, de_novo = st.tabs(["ADMET", "De Novo"])

with admet_tab:
    col1, col2 = st.columns(2)

    with col1:        
        sequence = st.text_input("SMILES Sequence", "CCCC")
        st.write("ADMET Results: ")
        results = admet_from_smiles(sequence)
        
        # Create a unified DataFrame
        combined_data = []
        
        # Combine Physicochemical data
        for key, value in results["Physicochemical"].items():
            combined_data.append({"Property": key, "Category": "Physicochemical", "Value": value})
        
        # Combine Absorption data
        for key, value in results["Absorption"].items():
            combined_data.append({"Property": key, "Category": "Absorption", "Value": value})
        
        # Combine Distribution data
        for key, value in results["Distribution"].items():
            combined_data.append({"Property": key, "Category": "Distribution", "Value": value})
        
        # Combine Metabolism data
        for key, value in results["Metabolism"].items():
            combined_data.append({"Property": key, "Category": "Metabolism", "Value": value})
        
        # Combine Excretion data
        for key, value in results["Excretion"].items():
            combined_data.append({"Property": key, "Category": "Excretion", "Value": value})
        
        # Combine Toxicity data
        for key, value in results["Toxicity"].items():
            combined_data.append({"Property": key, "Category": "Toxicity", "Value": value})
        
        # Create DataFrame from combined data
        combined_df = pd.DataFrame(combined_data)
        
        # Display the table with a fixed width (use_container_width)
        st.dataframe(combined_df, use_container_width=True)

    with col2:
        m = Chem.MolFromSmiles(sequence)
        st.image(Draw.MolToImage(m))


with de_novo:
    st.write("De Novo")

    