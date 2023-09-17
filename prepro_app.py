import streamlit as st
import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client
import base64
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split, cross_val_predict
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from sklearn.metrics import mean_squared_error
from math import sqrt

# Define a Streamlit app function
def main():
    st.title("Data Preprocessing and QSAR Modeling App")
    st.write("Welcome to the Data Preprocessing and QSAR Modeling App.")
    st.write("This app allows you to preprocess data and build a QSAR regression model based on a ChemBL ID.")
    st.write("Please enter your ChemBL ID below to get started.")

    # User input for ChemBL ID
    chembl_id = st.text_input("Enter ChemBL ID:")


    if st.button("Preprocess Data and Build Model"):
        # Measure the start time
        start_time = time.time()

        # Preprocess data based on user input
        success, num_initial_molecules, num_filtered_molecules, num_omitted_molecules_fp, preprocessed_data, morgan_data, raw_data, filtered_data = preprocess_data(chembl_id)

        # Measure the end time
        end_time = time.time()

        # Calculate the processing time
        processing_time = end_time - start_time

                
        if success:
            # Display processing time and molecule information
            st.write(f"Data preprocessing completed in {processing_time:.2f} seconds.")
            st.write(f"Initial number of molecules: {num_initial_molecules}")
            st.write(f"Number of molecules after filtering: {num_filtered_molecules}")
            st.write(f"Molecules failed to convert to fingerprint: {num_omitted_molecules_fp}")

            # Measure the start time for model development
            model_start_time = time.time()

            # Create qsar_raw DataFrame
            qsar_raw = preprocessed_data[['pIC50']].copy()
            qsar_raw = pd.concat([qsar_raw, morgan_data.iloc[:, 1:]], axis=1)

            # Save the qsar_raw data to 'qsar_raw.csv'
            qsar_raw.to_csv("qsar_raw.csv", index=False)

            st.write("Data for QSAR analysis saved as **qsar_raw.csv**.")

            # Split the data into internal (80%) and external dataset (20%)
            internal_data, external_data = train_test_split(qsar_raw, test_size=0.2, random_state=42)

            # Split the internal dataset into internal train (80%) and internal test (20%)
            internal_train, internal_test = train_test_split(internal_data, test_size=0.2, random_state=42)

            # Build a random forest regression model
            rf_model = RandomForestRegressor(n_estimators=100, random_state=42)

            # Internal train and test sets
            X_internal_train = internal_train.drop(columns=['pIC50'])
            y_internal_train = internal_train['pIC50']
            X_internal_test = internal_test.drop(columns=['pIC50'])
            y_internal_test = internal_test['pIC50']


            # Perform 10 trials for model building and evaluation
            num_trials = 10
            internal_train_rmses = []
            qsquared_cv_values = []
            rmse_cv_values = []

            for _ in range(num_trials):
                # Fit the model on the internal train data
                rf_model.fit(X_internal_train, y_internal_train)

                # Calculate RMSE for internal train data
                internal_train_rmse = np.sqrt(mean_squared_error(y_internal_train, rf_model.predict(X_internal_train)))
                internal_train_rmses.append(internal_train_rmse)

                # Perform 10-fold cross-validation
                cv_predictions = cross_val_predict(rf_model, X_internal_train, y_internal_train, cv=10)

                # Calculate Q-squared (CV)
                r2_internal_train = rf_model.score(X_internal_train, y_internal_train)
                qsquared_cv = 1 - (np.var(y_internal_train - cv_predictions) / np.var(y_internal_train))
                qsquared_cv_values.append(qsquared_cv)

                # Calculate RMSE for cross-validation
                rmse_cv = np.sqrt(mean_squared_error(y_internal_train, cv_predictions))
                rmse_cv_values.append(rmse_cv)

            # Calculate the metrics
            internal_train_rmse_mean = np.mean(internal_train_rmses)
            qsquared_cv_mean = np.mean(qsquared_cv_values)
            qsquared_cv_std_dev = np.std(qsquared_cv_values)

            # Calculate RMSE for internal training set with std deviation
            internal_train_rmse_std_dev = np.std(internal_train_rmses)

            # Prepare data for external validation
            X_external = external_data.drop(columns=['pIC50'])
            y_external = external_data['pIC50']

            # Predict pIC50 for external data
            external_predictions = rf_model.predict(X_external)

            # Calculate R-squared for external validation
            r2_external = r2_score(y_external, external_predictions)

            # Calculate RMSE for external validation
            rmse_external = np.sqrt(mean_squared_error(y_external, external_predictions))

            # Calculate RMSE standard deviation for external validation
            rmse_external_std_dev = np.std([np.sqrt(mean_squared_error(y_external, rf_model.predict(X_external))) for _ in range(num_trials)])

            # Calculate R-Squared - Q-Squared (EV)
            r2_qsquared_ev = r2_internal_train - r2_external

            # Calculate R-Squared - Q-Squared (CV)
            r2_qsquared_cv = r2_internal_train - qsquared_cv_mean

            # Create a dataframe for results
            results_df = pd.DataFrame({
                'Descriptor class': ['Morgan'],
                'N': [2048],
                'Training set (R-Squared)': [f'{r2_internal_train:.2f} ± {np.std(cv_predictions):.2f}'],
                'RMSE_Tr': [f'{internal_train_rmse_mean:.2f} ± {internal_train_rmse_std_dev:.2f}'],
                '10-Fold CV (Q-squared)': [f'{qsquared_cv:.2f} ± {qsquared_cv_std_dev:.2f}'],
                'CV_RMSE': [f'{np.mean(rmse_cv_values):.2f} ± {np.std(rmse_cv_values):.2f}'],
                'External set (Q-Squared)': [f'{r2_external:.2f} ± {np.std(external_predictions):.2f}'],
                'EV-RMSE': [f'{rmse_external:.2f} ± {rmse_external_std_dev:.2f}'],
                'R2−Q2CV': [f'{r2_qsquared_cv:.2f}'],
                'R2−Q2Ext': [f' {r2_qsquared_ev:.2f}']
            })
            # Create a dataframe for experimental vs predicted values
            train_results = pd.DataFrame({
                'Experimental pIC50': y_internal_train,
                'Predicted pIC50 (Train)': rf_model.predict(X_internal_train)
            })

            cv_results = pd.DataFrame({
                'Experimental pIC50': y_internal_train,
                'Predicted pIC50 (CV)': cv_predictions
            })

            ev_results = pd.DataFrame({
                'Experimental pIC50': y_external,
                'Predicted pIC50 (EV)': external_predictions
            })
            # Save the results to 'qsar_results.csv'
            results_df.to_csv("qsar_results.csv", index=False)


            # Plot 1: Experimental vs Predicted for Train
# Replace these lines:
# st.pyplot()

            # With these lines:
            fig, ax = plt.subplots()
            sns.scatterplot(x='Experimental pIC50', y='Predicted pIC50 (Train)', data=train_results, ax=ax)
            plt.title('Experimental vs Predicted pIC50 for Train')
            plt.xlabel('Experimental pIC50')
            plt.ylabel('Predicted pIC50 (Train)')
            # Add a central line (y=x)
            plt.plot([train_results.min().min(), train_results.max().max()], [train_results.min().min(), train_results.max().max()], color='red', linestyle='--')
            st.pyplot(fig)


            # Plot 2: Experimental vs Predicted for CV
            fig2, ax2 = plt.subplots()
            sns.scatterplot(x='Experimental pIC50', y='Predicted pIC50 (CV)', data=cv_results, ax=ax2)
            plt.title('Experimental vs Predicted pIC50 for CV')
            plt.xlabel('Experimental pIC50')
            plt.ylabel('Predicted pIC50 (CV)')
            # Add a central line (y=x)
            plt.plot([cv_results.min().min(), cv_results.max().max()], [cv_results.min().min(), cv_results.max().max()], color='red', linestyle='--')
            st.pyplot(fig2)

            # Plot 3: Experimental vs Predicted for EV
            fig3, ax3 = plt.subplots()
            sns.scatterplot(x='Experimental pIC50', y='Predicted pIC50 (EV)', data=ev_results, ax=ax3)
            plt.title('Experimental vs Predicted pIC50 for EV')
            plt.xlabel('Experimental pIC50')
            plt.ylabel('Predicted pIC50 (EV)')
            # Add a central line (y=x)
            plt.plot([ev_results.min().min(), ev_results.max().max()], [ev_results.min().min(), ev_results.max().max()], color='red', linestyle='--')
            st.pyplot(fig3)

            
            # Measure the end time for model development
            model_development_time = time.time() - model_start_time

            st.write(f"Model development completed in {model_development_time:.2f} seconds.")

            st.write("\n---\n")
            st.write("QSAR Model Results:")
            st.write(results_df)

            # Generate a downloadable link for the model
            model_filename = "model.pkl"
            with open(model_filename, 'wb') as model_file:
                pickle.dump(rf_model, model_file)

            # Provide download links to individual CSV files
            st.markdown(get_table_download_link(raw_data, "Raw_Data.csv"), unsafe_allow_html=True)
            st.markdown(get_table_download_link(filtered_data, "Filtered_Data.csv"), unsafe_allow_html=True)
            st.markdown(get_table_download_link(preprocessed_data, "Preprocessed_Data.csv"), unsafe_allow_html=True)
            st.markdown(get_table_download_link(morgan_data, "Morgan_Fingerprints.csv"), unsafe_allow_html=True)
            st.markdown(get_table_download_link(qsar_raw, "qsar_raw.csv"), unsafe_allow_html=True)
            st.markdown(get_table_download_link(results_df, "QSAR_Results.csv"), unsafe_allow_html=True)
            st.markdown(get_model_download_link(model_filename), unsafe_allow_html=True)

            st.write("\n---\n")
            st.write("Powered by Parth Sanghavi")
            st.write('Preprocessed_Data.csv and Morgan_Fingerprints.csv can be uploaded to the QSAR webapp (link coming soon) to generate a robust 2D-QSAR model')

        else:
            # Display an error message for invalid ChemBL ID
            st.error("Please enter a valid ChemBL ID.")

# Define the preprocess_data function
def preprocess_data(chembl_id):
    try:
        # Remove spaces from the user input
        chembl_id = chembl_id.strip()

        # Target search
        target = new_client.target
        target_query = target.search(chembl_id)
        targets = pd.DataFrame.from_dict(target_query)

        if targets.empty:
            # Display an error message if no data is found
            st.error("No data found for the provided ChemBL ID.")
            return False, 0, 0, 0, None, None, None, None  # Return None for dataframes

        selected_target = targets.target_chembl_id[0]

        # Activity search and filtering
        activity = new_client.activity
        res = activity.filter(target_chembl_id=selected_target, standard_type="IC50", standard_relation="=")

        # Step 1: Save the raw data to 'Raw_Data.csv'
        raw_data = pd.DataFrame.from_dict(res)
        raw_data.to_csv("Raw_Data.csv", index=False)

        st.write("1. **Raw_Data.csv**: Contains the original data retrieved from ChEMBL.")

        num_initial_molecules = len(raw_data)

        # Step 2: Filter data and save to 'Filtered_Data.csv'
        filtered_data = raw_data.dropna(subset=['standard_value'])
        filtered_data = filtered_data.drop_duplicates(subset=['canonical_smiles', 'molecule_chembl_id'])
        filtered_data.to_csv("Filtered_Data.csv", index=False)

        st.write("2. **Filtered_Data.csv**: Contains data after applying filters for standard type as IC50, standard relation as '=', and removing empty standard values and duplicates.")

        num_filtered_molecules = len(filtered_data)

        # Step 3: Calculate pIC50 and save to 'Preprocessed_Data.csv'
        preprocessed_data = filtered_data[['molecule_chembl_id', 'canonical_smiles', 'standard_value']]
        preprocessed_data['pIC50'] = -np.log10(filtered_data['standard_value'].astype(float) * 1e-9)
        preprocessed_data.to_csv("Preprocessed_Data.csv", index=False)

        st.write("3. **Preprocessed_Data.csv**: Contains the final preprocessed data with selected columns including calculated pIC50 values.")

        # Calculate Morgan fingerprints
        morgan_data, num_omitted_molecules = calculate_morgan_fingerprints(preprocessed_data)

        return True, num_initial_molecules, num_filtered_molecules, num_omitted_molecules, preprocessed_data, morgan_data, raw_data, filtered_data

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
        return False, 0, 0, 0, None, None, None, None

# Function to calculate Morgan fingerprints
def calculate_morgan_fingerprints(data):
    morgan_data = []
    num_omitted_molecules_fp = 0

    for index, row in data.iterrows():
        smiles = row['canonical_smiles']
        molecule = Chem.MolFromSmiles(smiles)

        if molecule is not None:
            fingerprints = AllChem.GetMorganFingerprintAsBitVect(molecule, 3, nBits=2048)
            fingerprint_values = list(fingerprints.ToBitString())
            morgan_data.append([row['molecule_chembl_id']] + fingerprint_values)
        else:
            num_omitted_molecules_fp += 1

    morgan_data = pd.DataFrame(morgan_data, columns=['molecule_chembl_id'] + [f'morgan_{i}' for i in range(2048)])
    return morgan_data, num_omitted_molecules_fp


# Function to calculate internal train RMSE
def calculate_internal_train_rmse(model, X_train, y_train):
    y_train_pred = model.predict(X_train)
    rmse = sqrt(mean_squared_error(y_train, y_train_pred))
    return rmse

# Function to create a download link for the model
def get_model_download_link(model_filename):
    with open(model_filename, 'rb') as file:
        model_bytes = file.read()
    b64 = base64.b64encode(model_bytes).decode()
    href = f'<a href="data:application/octet-stream;base64,{b64}" download="{model_filename}">Download Model</a>'
    return href

# Function to create a download link for a DataFrame
def get_table_download_link(df, filename):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">Download {filename}</a>'
    return href

# Run the Streamlit app
if __name__ == "__main__":
    main()
