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

title = st.text_input("SMILES Sequence", "CCCC")
st.write("ADMET Results: ")
results = admet_from_smiles(title)
st.table(results["Physicochemical"])
st.table(results["Absorption"])
st.table(results["Distribution"])
st.table(results["Metabolism"])
st.table(results["Excretion"])