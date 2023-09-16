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
        st.markdown(get_table_download_link(raw_data, f"Raw_Data.csv"), unsafe_allow_html=True)

        # Step 2: Apply filters and save to 'Filtered_Data.csv'
        filtered_data = raw_data.copy()
        filtered_data = filtered_data[filtered_data['standard_relation'] == '=']
        filtered_data = filtered_data.dropna(subset=['standard_value'])
        filtered_data = filtered_data.drop_duplicates(subset=['canonical_smiles', 'molecule_chembl_id'])

        filtered_data.to_csv(f"{output_dir}/Filtered_Data.csv", index=False)
        st.write("2. **Filtered_Data.csv**: Contains data after applying filters for standard type as IC50, standard relation as '=', and removing duplicates.")
        st.markdown(get_table_download_link(filtered_data, f"Filtered_Data.csv"), unsafe_allow_html=True)

        # Step 3: Get only 3 columns from the filtered_data
        filtered_data = filtered_data[['molecule_chembl_id', 'canonical_smiles', 'standard_value']].copy()

        # Add a column to convert standard_value to pIC50 values
        filtered_data['pIC50'] = -1 * (filtered_data['standard_value'].apply(lambda x: pd.np.log10(x / 1e9)))

        # Save the preprocessed data to 'Preprocessed_Data.csv'
        st.write("3. **Preprocessed_Data.csv**: Contains the final preprocessed data with selected columns.")
        st.markdown(get_table_download_link(filtered_data, f"Preprocessed_Data.csv"), unsafe_allow_html=True)

        st.success("Data preprocessing completed.")

        st.write("\n---\n")
        st.write("Powered by Parth Sanghavi")

        # For feedback
        st.write("For feedback and inquiries, please email: unprofessor.edu@gmail.com")

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
