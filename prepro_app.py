import streamlit as st
import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client
import base64

# Define a Streamlit app function
def main():
    st.title("Data Preprocessing App")
    st.write("Welcome to the Data Preprocessing App for QSAR analysis.")
    st.write("This app allows you to preprocess data for QSAR analysis based on a ChemBL ID.")
    st.write("Please enter your ChemBL ID below to get started.")

    # User input for ChemBL ID
    chembl_id = st.text_input("Enter ChemBL ID:")

    if st.button("Preprocess Data"):
        # Preprocess data based on user input
        preprocess_data(chembl_id)

# Define the preprocess_data function
def preprocess_data(chembl_id):
    try:
        # Target search
        target = new_client.target
        target_query = target.search(chembl_id)
        targets = pd.DataFrame.from_dict(target_query)

        selected_target = targets.target_chembl_id[0]

        # Activity search and filtering
        activity = new_client.activity
        res = activity.filter(target_chembl_id=selected_target, standard_type="IC50", standard_relation="=")

        # Step 1: Save the raw data to 'Raw_Data.csv'
        raw_data = pd.DataFrame.from_dict(res)
        raw_data.to_csv("Raw_Data.csv", index=False)

        st.write("1. **Raw_Data.csv**: Contains the original data retrieved from ChEMBL.")

        # Step 2: Filter data and save to 'Filtered_Data.csv'
        filtered_data = raw_data.dropna(subset=['standard_value'])
        filtered_data = filtered_data.drop_duplicates(subset=['canonical_smiles', 'molecule_chembl_id'])
        filtered_data.to_csv("Filtered_Data.csv", index=False)

        st.write("2. **Filtered_Data.csv**: Contains data after applying filters for standard type as IC50, standard relation as '=', and removing empty standard values and duplicates.")

        # Step 3: Calculate pIC50 and save to 'Preprocessed_Data.csv'
        preprocessed_data = filtered_data[['molecule_chembl_id', 'canonical_smiles', 'standard_value']]
        preprocessed_data['pIC50'] = -np.log10(filtered_data['standard_value'].astype(float) * 1e-9)
        preprocessed_data.to_csv("Preprocessed_Data.csv", index=False)

        st.write("3. **Preprocessed_Data.csv**: Contains the final preprocessed data with selected columns including calculated pIC50 values.")

        # Provide download links to individual CSV files
        st.markdown(get_table_download_link(raw_data, "Raw_Data.csv"), unsafe_allow_html=True)
        st.markdown(get_table_download_link(filtered_data, "Filtered_Data.csv"), unsafe_allow_html=True)
        st.markdown(get_table_download_link(preprocessed_data, "Preprocessed_Data.csv"), unsafe_allow_html=True)

        st.write("\n---\n")
        st.write("Powered by Parth Sanghavi")
        st.write('Preprocessed_Data.csv can be uploaded to the QSAR webapp (link coming soon) to generate a robust 2D-QSAR model')

        # For feedback
        st.write("For feedback and inquiries, please email: unprofessor.edu@gmail.com")

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")

# Function to create a download link for a DataFrame
def get_table_download_link(df, filename):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">Download {filename}</a>'
    return href

# Run the Streamlit app
if __name__ == "__main__":
    main()
