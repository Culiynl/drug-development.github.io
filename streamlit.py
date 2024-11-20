from admet_ai import ADMETModel
import math
import matplotlib.pyplot as plt
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

# Function to predict ADMET properties
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

    # Organize predictions
    organized_predictions = {}
    for category, props in categories.items():
        organized_predictions[category] = {prop: preds[prop] for prop in props if prop in preds}

    return organized_predictions

# Function to create a radar chart
def radar_chart(data, upper_limits):
    pi = math.pi
    categories = list(data.keys())
    values = [data[key] for key in categories]

    # Number of variables
    N = len(categories)
    
    # Compute angle of each axis
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]  # To close the circle

    # Plotting
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
    ax.set_theta_offset(pi / 2)  # Rotate the chart so that the first category is at the top
    ax.set_theta_direction(-1)  # Clockwise

    # Draw one axis per variable and add labels
    ax.set_rlabel_position(0)
    ax.set_xticks(angles[:-1])  # Remove the last angle
    ax.set_xticklabels(categories, color='black', size=10)

    # Adjust yticks based on upper limits
    max_value = max(upper_limits.values())  # Find the maximum upper limit for scaling
    ax.set_yticks([i * max_value / 5 for i in range(1, 6)])  # Adjust as necessary
    ax.set_yticklabels([str(int(i * max_value / 5)) for i in range(1, 6)], color='gray', size=8)

    # Scale the values according to their respective upper limits
    scaled_values = [value / upper_limits[categories[i]] * max_value for i, value in enumerate(values)]
    scaled_values += scaled_values[:1]  # Close the circle by repeating the first value

    # Plot the data
    ax.plot(angles, scaled_values, color='b', linewidth=2, linestyle='solid')
    ax.fill(angles, scaled_values, color='b', alpha=0.3)

    return fig

# Define upper limits for scaling the values
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

# Streamlit app setup
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

        # Combine data from all categories
        for category in results:
            for key, value in results[category].items():
                combined_data.append({"Property": key, "Category": category, "Value": value})

        # Create DataFrame from combined data
        combined_df = pd.DataFrame(combined_data)

        # Display the table with a fixed width
        st.dataframe(combined_df, use_container_width=True)

    with col2:
        m = Chem.MolFromSmiles(sequence)
        st.image(Draw.MolToImage(m))

        st.write("Radar Chart:")
        fig = radar_chart(results["Physicochemical"], upper_limits)
        st.pyplot(fig)

with de_novo:
    st.write("De Novo")
