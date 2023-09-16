import streamlit as st
import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client

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

        # Create a directory to store CSV files
        output_dir = f"preprocessed_data_{chembl_id}"
        
        # Step 1: Save the raw data to 'Raw_Data.csv'
        raw_data = pd.DataFrame.from_dict(res)
        raw_data.to_csv(f"{output_dir}/Raw_Data.csv", index=False)

        st.write("1. **Raw_Data.csv**: Contains the original data retrieved from ChEMBL.")
        
        # Step 2: Filter data and save to 'Filtered_Data.csv'
        filtered_data = raw_data[['molecule_chembl_id', 'canonical_smiles', 'standard_value']]
        filtered_data = filtered_data.dropna(subset=['standard_value'])
        filtered_data = filtered_data.drop_duplicates(subset=['canonical_smiles', 'molecule_chembl_id'])
        filtered_data.to_csv(f"{output_dir}/Filtered_Data.csv", index=False)

        st.write("2. **Filtered_Data.csv**: Contains data after applying filters for standard type as IC50, standard relation as '=', and removing empty standard values and duplicates.")
        
        # Step 3: Calculate pIC50 and save to 'Preprocessed_Data.csv'
        filtered_data['pIC50'] = -np.log10(filtered_data['standard_value'].astype(float) * 1e-9)
        filtered_data.to_csv(f"{output_dir}/{chembl_id}_Preprocessed_Data.csv", index=False)

        st.write("3. **Preprocessed_Data.csv**: Contains the final preprocessed data with selected columns including calculated pIC50 values.")

        # Provide download links to individual CSV files
        st.markdown(f"**Download Raw Data:** [{chembl_id}_Raw_Data.csv](/{output_dir}/Raw_Data.csv)")
        st.markdown(f"**Download Filtered Data:** [{chembl_id}_Filtered_Data.csv](/{output_dir}/Filtered_Data.csv)")
        st.markdown(f"**Download Preprocessed Data:** [{chembl_id}_Preprocessed_Data.csv](/{output_dir}/{chembl_id}_Preprocessed_Data.csv)")

        st.write("\n---\n")
        st.write("Powered by Parth Sanghavi")

        # For feedback
        st.write("For feedback and inquiries, please email: unprofessor.edu@gmail.com")

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")

# Run the Streamlit app
if __name__ == "__main__":
    main()
