import streamlit as st
import pandas as pd
import base64  # Import base64 module
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

        # Create a directory to store CSV files (if needed)
        output_dir = f"preprocessed_data_{chembl_id}"
        # os.makedirs(output_dir, exist_ok=True)  # Commented out to provide direct download links

        # Step 1: Save the raw data to 'Raw_Data.csv'
        raw_data = pd.DataFrame.from_dict(res)
        st.write("1. **Raw_Data.csv**: Contains the original data retrieved from ChEMBL.")
        st.markdown(get_table_download_link(raw_data, f"{chembl_id}_Raw_Data.csv"), unsafe_allow_html=True)

        # Step 2: Apply filters and save to 'Filtered_Data.csv' (add your code here)
        filtered_data = raw_data  # Replace with your filtered data
        st.write("2. **Filtered_Data.csv**: Contains data after applying filters.")
        st.markdown(get_table_download_link(filtered_data, f"{chembl_id}_Filtered_Data.csv"), unsafe_allow_html=True)

        # Step 3: Perform additional preprocessing and save to 'Preprocessed_Data.csv' (add your code here)
        preprocessed_data = raw_data  # Replace with your preprocessed data
        st.write("3. **Preprocessed_Data.csv**: Contains the final preprocessed data with selected columns.")
        st.markdown(get_table_download_link(preprocessed_data, f"{chembl_id}_Preprocessed_Data.csv"), unsafe_allow_html=True)

        # Perform IC50 to pIC50 conversion (you can add this part if needed)

        st.success("Data preprocessing completed.")

        st.write("\n---\n")
        st.write("Powered by Parth Sanghavi")

        # For feedback
        st.write("For feedback and inquiries, please email: unprofessor.edu@gmail.com")

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")

# Function to generate a download link for a DataFrame
def get_table_download_link(df, filename):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">Download {filename}</a>'
    return href

# Run the Streamlit app
if __name__ == "__main__":
    main()
