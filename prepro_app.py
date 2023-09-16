import streamlit as st
import pandas as pd
import numpy as np
import os
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
        res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50").filter(unit="nM")

        # Create a directory to store CSV files
        output_dir = f"preprocessed_data_{chembl_id}"
        os.makedirs(output_dir, exist_ok=True)

        # Step 1: Save the raw data to 'Raw_Data.csv'
        raw_data = pd.DataFrame.from_dict(res)
        raw_data.to_csv(f"{output_dir}/Raw_Data.csv", index=False)

        st.write("1. **Raw_Data.csv**: Contains the original data retrieved from ChEMBL.")
        st.write("2. **Filtered_Data.csv**: Contains data after applying filters for standard type as IC50,  standard units as nM, and '=' as the standard relation.")
        st.write("3. **Preprocessed_Data.csv**: Contains the final preprocessed data with selected columns.")

        # Perform IC50 to pIC50 conversion (you can add this part if needed)

        st.success("Data preprocessing completed.")

        # Provide download links to individual CSV files
        download_links = {
            "Raw Data": f"{output_dir}/Raw_Data.csv",
            "Filtered Data": f"{output_dir}/Filtered_Data.csv",
            "Preprocessed Data": f"{output_dir}/Preprocessed_Data.csv"
        }

        for file_name, file_path in download_links.items():
            st.markdown(get_download_link(file_path, file_name), unsafe_allow_html=True)

        st.write("\n---\n")
        st.write("Powered by Parth Sanghavi")

        # For feedback
        st.write("For feedback and inquiries, please email: unprofessor.edu@gmail.com")

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")

# Function to generate download links
def get_download_link(file_path, file_name):
    with open(file_path, "rb") as file:
        file_bytes = file.read()
    return f'<a href="data:file/csv;base64,{file_bytes.decode()}" download="{file_name}">Download {file_name}</a>'

# Run the Streamlit app
if __name__ == "__main__":
    main()
