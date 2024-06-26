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
from sklearn.metrics import precision_score, accuracy_score, recall_score, confusion_matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve


# Define a Streamlit app function
def main():
    st.title("QSAR Modeling App")
    st.write("Welcome to the Data Preprocessing and QSAR Modeling App.")
    st.write("This app allows you to preprocess data and build a QSAR regression model based on the ChemBL ID.")
    st.write("If you encounter issues with the application, please drop an email to parthsanghavi017@gmail.com. Kindly include the CHEMBL ID, a screenshot of the issue, and a description if applicable.")
    st.write("Please enter your ChEMBL ID below to get started.")

    # User input for ChemBL ID
    chembl_id = st.text_input("Enter ChEMBL ID:")

    if st.button("Preprocess Data"):
        # Call the preprocess_data function to perform dataset preprocessing
        success, num_initial_molecules, num_filtered_molecules, num_omitted_molecules_fp, preprocessed_data, morgan_data, raw_data, filtered_data = preprocess_data(chembl_id)

        if success:
            # Display preprocessing results and ask the user if they want to continue with modeling
            st.write(f"Data preprocessing completed.")
            st.write(f"Initial number of molecules: {num_initial_molecules}")
            st.write(f"Number of molecules after filtering: {num_filtered_molecules}")
            st.write(f"Molecules failed to convert to fingerprint: {num_omitted_molecules_fp}")

            st.write("You can download the preprocessed data files now.")
            st.write("1. **Raw_Data.csv**: Contains the original data retrieved from ChEMBL.")
            st.markdown(get_table_download_link(raw_data, "Raw_Data.csv"), unsafe_allow_html=True)

            st.write("2. **Filtered_Data.csv**: Contains data after applying filters for standard type as IC50, standard relation as '=', and removing empty standard values and duplicates.")
            st.markdown(get_table_download_link(filtered_data, "Filtered_Data.csv"), unsafe_allow_html=True)

            st.write("3. **Preprocessed_Data.csv**: Contains the final preprocessed data with selected columns including calculated pIC50 values.")
            st.markdown(get_table_download_link(preprocessed_data, "Preprocessed_Data.csv"), unsafe_allow_html=True)

            st.write("4. **Morgan_Fingerprints.csv**: Contains ChEMBL ID, Morgan Fingerprints with nBits = 2048 and radius = 3")
            st.markdown(get_table_download_link(morgan_data, "Morgan_Fingerprints.csv"), unsafe_allow_html=True)

            # Ask the user if they want to continue with modeling
            continue_modeling = st.button("Continue with QSAR Modeling")

            if continue_modeling:
                # Call the function to perform QSAR modeling
                perform_qsar_modelling(chembl_id, preprocessed_data, morgan_data)
        else:
            # Display an error message for an invalid ChemBL ID
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
        num_initial_molecules = len(raw_data)

        # Step 2: Filter data and save to 'Filtered_Data.csv'
        filtered_data = raw_data.dropna(subset=['standard_value'])
        filtered_data = filtered_data.drop_duplicates(subset=['canonical_smiles', 'molecule_chembl_id'])
        filtered_data.to_csv("Filtered_Data.csv", index=False)
        num_filtered_molecules = len(filtered_data)

        # Step 3: Calculate pIC50 and save to 'Preprocessed_Data.csv'
        preprocessed_data = filtered_data[['molecule_chembl_id', 'canonical_smiles', 'standard_value']]

        # Filter out rows with empty canonical SMILES
        preprocessed_data = preprocessed_data.dropna(subset=['canonical_smiles'])

        # Calculate pIC50
        preprocessed_data['pIC50'] = -np.log10(preprocessed_data['standard_value'].astype(float) * 1e-9)

        # Drop rows with "inf" values in the 'pIC50' column
        preprocessed_data = preprocessed_data.replace([np.inf, -np.inf], np.nan).dropna(subset=['pIC50'])
        preprocessed_data.to_csv("Preprocessed_Data.csv", index=False)

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

    # Lists to store Morgan fingerprints and molecule IDs
    morgan_fps = []
    molecule_ids = []

    for index, row in data.iterrows():
        smiles = row['canonical_smiles']
        molecule = Chem.MolFromSmiles(smiles)

        if molecule is not None:
            try:
                fingerprints = AllChem.GetMorganFingerprintAsBitVect(molecule, 3, nBits=2048)
                if fingerprints:
                    fingerprint_values = list(fingerprints.ToBitString())
                    morgan_fps.append([row['molecule_chembl_id']] + fingerprint_values)
                    molecule_ids.append(row['molecule_chembl_id'])
                else:
                    num_omitted_molecules_fp += 1
            except Exception as e:
                num_omitted_molecules_fp += 1
        else:
            num_omitted_molecules_fp += 1

    # Create a DataFrame with 'molecule_chembl_id' and individual Morgan fingerprint columns
    morgan_data = pd.DataFrame(morgan_fps, columns=['molecule_chembl_id'] + [f'morgan_{i}' for i in range(2048)])

    return morgan_data, num_omitted_molecules_fp


def perform_qsar_modelling(chembl_id, preprocessed_data, morgan_data, random_seed=42):
    # Rest of the function remains the same

        try:
            # Measure the start time for model development
            model_start_time = time.time()

            # Create qsar_raw DataFrame by joining preprocessed_data and morgan_data on 'molecule_chembl_id'
            qsar_raw = preprocessed_data[['molecule_chembl_id', 'pIC50']].copy()
            qsar_raw = qsar_raw.set_index('molecule_chembl_id').join(morgan_data.set_index('molecule_chembl_id'))

            # Reset the index to bring back 'molecule_chembl_id' as a column
            qsar_raw = qsar_raw.reset_index()

            # Drop the 'molecule_chembl_id' column from the qsar_raw DataFrame
            qsar_raw = qsar_raw.drop(columns=['molecule_chembl_id'])

            # After Step 5: Save the qsar_raw data to 'QSAR_Raw.csv'
            qsar_raw.to_csv("QSAR_Raw.csv", index=False)
            st.write("5. **QSAR_Raw.csv**: Contains the pIC50 and Morgan Fingerprint Data")
            st.markdown(get_table_download_link(qsar_raw, "QSAR_Raw.csv"), unsafe_allow_html=True)

            # Split the data into internal (80%) and external dataset (20%)
            internal_data, external_data = train_test_split(qsar_raw, test_size=0.2, random_state=random_seed)

            # Split the internal dataset into internal train (80%) and internal test (20%)
            internal_train, internal_test = train_test_split(internal_data, test_size=0.2, random_state=random_seed)

            # Build a random forest regression model
            rf_model = RandomForestRegressor(n_estimators=100, random_state=random_seed)

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

            # After Step 6: Save the results_df data to 'QSAR_Results.csv'
            results_df.to_csv("QSAR_Results.csv", index=False)
            st.write("6. **QSAR_Results.csv**: Contains the results of the QSAR model")
            st.markdown(get_table_download_link(results_df, "QSAR_Results.csv"), unsafe_allow_html=True)

            # Calculate precision, accuracy, sensitivity, and specificity
            tn, fp, fn, tp = confusion_matrix(y_external > 5, external_predictions > 5).ravel()
            precision = tp / (tp + fp)
            accuracy = (tp + tn) / (tp + tn + fp + fn)
            sensitivity = tp / (tp + fn)
            specificity = tn / (tn + fp)

            # Create a dataframe for additional metrics
            additional_metrics_df = pd.DataFrame({
                'Metric': ['Precision', 'Accuracy', 'Sensitivity', 'Specificity'],
                'Value': [precision, accuracy, sensitivity, specificity]
            })

            # After Step 7: Save the additional metrics data to 'Additional_Metrics.csv'
            additional_metrics_df.to_csv("Additional_Metrics.csv", index=False)
            st.write("7. **Additional_Metrics.csv**: Contains additional metrics (Precision, Accuracy, Sensitivity, Specificity)")
            st.markdown(get_table_download_link(additional_metrics_df, "Additional_Metrics.csv"), unsafe_allow_html=True)


            # Display the additional metrics
            st.write("Additional Metrics:")
            st.write(additional_metrics_df)

            # Plot 4: Precision-Recall Curve
            fig4, ax4 = plt.subplots()
            precision, recall, _ = precision_recall_curve(y_external > 5, external_predictions)
            plt.plot(recall, precision, marker='.')
            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.title('Precision-Recall Curve')
            plot_filename4 = "plot_precision_recall.png"
            plt.savefig(plot_filename4)
            st.image(plot_filename4)

            # Plot 5: ROC Curve
            fig5, ax5 = plt.subplots()
            fpr, tpr, _ = roc_curve(y_external > 5, external_predictions)
            plt.plot(fpr, tpr, marker='.')
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('ROC Curve')
            plot_filename5 = "plot_roc.png"
            plt.savefig(plot_filename5)
            st.image(plot_filename5)

            # Measure the end time for model development
            model_development_time = time.time() - model_start_time

            # Calculate the model development time
            model_development_time_minutes = int(model_development_time // 60)
            model_development_time_seconds = int(model_development_time % 60)

            st.write(f"Model development completed in {model_development_time_minutes} minutes and {model_development_time_seconds} seconds.")
            st.write("\n---\n")
            st.write("QSAR Model Results:")
            st.write(results_df)

            # Generate a downloadable link for the model
            model_filename = "RF_Model.pkl"
            with open(model_filename, 'wb') as model_file:
                pickle.dump(rf_model, model_file)

            # Provide download links to the QSAR model
            st.write("7. **RF_Model.pkl**: Contains the trained Random Forest QSAR model")
            st.markdown(get_model_download_link(model_filename), unsafe_allow_html=True)


            st.write("\n---\n")
            st.write("Powered by Parth Sanghavi")

        except Exception as e:
            st.error(f"An error occurred during QSAR modeling: {str(e)}")




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
