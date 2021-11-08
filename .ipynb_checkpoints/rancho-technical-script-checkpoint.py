import pandas as pd
import numpy as np

def import_data(excel_file):
    
    ## ******NEED to get ask for input file and read to 
    
    ''' Read excel file and output pandas dataframes for each sheet '''
    
    patient_clinical_data = pd.read_excel('Technical Test - Data Wrangling.xlsx', sheet_name = 'Patient_clinical_data')
    tissue_sample_metadata = pd.read_excel('Technical Test - Data Wrangling.xlsx', sheet_name = 'Tissue Sample Metadata')
    serum_protein_data = pd.read_excel('Technical Test - Data Wrangling.xlsx', sheet_name = 'Serum Protein data')
    rna_seq_data = pd.read_excel('Technical Test - Data Wrangling.xlsx', sheet_name = 'RNA-seq (RPKM)')    
    return patient_clinical_data, tissue_sample_metadata, serum_protein_data, rna_seq_data

def munge_patient_clinical_data(patient_clinical_data):
    
    ''' Add and alter columns of the patient clinical data for final database specifications '''
    
    # Add unique_patient_id column to patient_clinical_data df
    patient_clinical_data['Unique_Patient_ID'] = patient_clinical_data['Study_ID'] + '_' + patient_clinical_data['Patient  Number'].astype(str)

    # Round age to nearest number and convert to integer datatype
    patient_clinical_data.Age = patient_clinical_data.Age.round(decimals = 0).astype(int)

    # Convert sex to say male or female instead of M or F 
    patient_clinical_data.Sex = patient_clinical_data.Sex.replace({'M': 'MALE', 'F': 'FEMALE'})

    #Rename Patient Number column to Patient_ID and convert to integer
    patient_clinical_data.rename({'Patient  Number': 'Patient_ID'}, axis = 1, inplace = True)
    patient_clinical_data.Patient_ID.astype(int)
    return patient_clinical_data


def munge_tissue_sample_metadata(tissue_sample_metadata):
    
    ''' Add and alter columns of the tissue sample metadata for final databse specifications''' 
    
    # Add results_units column to tissue_sample_metadata
    tissue_sample_metadata['Result_Units'] = 'RPKM'

    # Drop RUN and Total Read(millions) columnb
    tissue_sample_metadata.drop(labels = ['RIN', 'Total Reads(millions)'], axis = 1, inplace = True)

    #Rename Sample type column
    tissue_sample_metadata.rename({'Sample type': 'Sample_General_Pathology', 'Sample': 'Sample_ID', 'Material': 'Material_type', 'Patient  Number': 'Patient_ID'}, axis = 1, inplace = True)
    return tissue_sample_metadata

def munge_rna_seq_data(rna_seq_data):
    
    ''' Alter rna_seq_data dataframe to isolate ICAM1 results to meet final database specifications '''

    #Transpose rna_seq_data 
    rna_seq_data_t = rna_seq_data.transpose()

    # Name columns based on first row information (currently the name of the genes)
    rna_seq_data_t = rna_seq_data_t.rename(columns=rna_seq_data_t.iloc[0]).drop(rna_seq_data_t.index[0])

    # Get name of all genes 
    gene_list = rna_seq_data_t.columns
    
    # Set up empty list 
    df_list = []
    
    #Isolate each gene and get gene name into Gene_Symbol column and the results in the Result column
    for gene in gene_list:
        rna_seq_gene = rna_seq_data_t.loc[:, [gene]]

        #Create Gene_Symbol column
        rna_seq_gene['Gene_Symbol'] = gene

        #Rename ICAM1 column to Results
        rna_seq_gene.rename({gene: 'Result'}, axis = 1, inplace = True)

        df_list.append(rna_seq_gene)
        
    # Concatenate all the individual gene dataframes
    rna_seq_all = pd.concat(df_list)
    
    return rna_seq_all


def munge_serum_protein_data(serum_protein_data):
    
    ''' Add and alter columns from serum protein dataframe to meet specifications of final database '''
    
    #Add Material_type column
    serum_protein_data['Material_type'] = 'Serum'

    #Convert data type of the Serum IL-6 Receptor and Serum IL-6 columns from string to float, using error coercion
    serum_protein_data['Serum IL-6 Receptor (g/L)'] = pd.to_numeric(serum_protein_data['Serum IL-6 Receptor (mg/L)'], errors = 'coerce')
    serum_protein_data['Serum IL-6 (g/L)'] = pd.to_numeric(serum_protein_data['Serum IL-6 (g/L)'], errors = 'coerce')

    #Convert mg/L unit to g/L for the Serum IL-6 receptor column in serum_protein_data dataframe
    serum_protein_data['Serum IL-6 Receptor (g/L)'] = serum_protein_data['Serum IL-6 Receptor (g/L)'].apply(lambda x: x * 0.001)

    #Remove Serum IL-6 Receptor Column with mg/L units, since it has now been calculated as g/L units
    serum_protein_data.drop(labels = 'Serum IL-6 Receptor (mg/L)', axis = 1, inplace = True)

    #Add Result_Units column
    serum_protein_data['Result_Units'] = 'g/L'

    #Rename Serum IL-6 and Serum IL-6 Receptor columns to IL6 and IL6R
    serum_protein_data.rename({'Serum IL-6 Receptor (g/L)': 'IL6R', 'Serum IL-6 (g/L)': 'IL6', 'Patient': 'Patient_ID', 'Sample': 'Sample_ID'}, axis = 1, inplace = True)

    # Melt the IL6 and IL6R columns into rows
    serum_protein_data = serum_protein_data.melt(id_vars = ['Patient_ID', 'Sample_ID', 'Result_Units', 'Material_type'], var_name = 'Gene_Symbol', value_name = 'Result')
    
    return serum_protein_data

def merge_dfs(patient_clinical_data, tissue_sample_metadata, rna_seq_icam1, serum_protein_data):
    
    ''' Merge and concatenate munged dataframes to creat final database. Create new column to identify if Result is available or not. Return final dataframe.'''
    
    # Merge tissue_sample df with patient_clinical_data df
    clinical_tissue = patient_clinical_data.merge(tissue_sample_metadata)

    # Merge rna_seq df with clinical_tissue df
    clinical_tissue_rnaseq = clinical_tissue.merge(rna_seq_icam1, left_on = 'Sample_ID', right_index = True)

    # Merge patient_clinical_data df with serum_protein_data df b
    clinical_serum = patient_clinical_data.merge(serum_protein_data)

    # Concatenate clinical_serum and clinical_tissue_rnaseq dataframes
    all_clinical_info = pd.concat([clinical_serum, clinical_tissue_rnaseq])

    # Create new column for status based on if results are NaN or numeric
    all_clinical_info['Status'] = np.where(all_clinical_info.Result.isnull(), 'Not Done', 'NaN')

    #Sort the dataframe
    all_clinical_info.sort_values(['Patient_ID', 'Material_type', 'Sample_ID', 'Gene_Symbol'], inplace = True)
    
    # Rearrange columns
    col_list = all_clinical_info.columns.tolist()
    col_order = ['Study_ID', 'Patient_ID','Sex','Age','Unique_Patient_ID','Sample_ID','Sample_General_Pathology','Material_type','Gene_Symbol','Result', 'Result_Units', 'Status']
    all_clinical_info = all_clinical_info.reindex(columns = col_order)
    
    return all_clinical_info

def save_final_database(final_data):
    file_name = 'Technical Test Submission - Data Wrangling.xlsx'
    final_data.to_excel(file_name, 'Final Database',  na_rep = 'NaN', index = False)
    print('Data wrangling complete - please open ' + file_name)
    
# Read in data  
patient_clinical_data, tissue_sample_metadata, serum_protein_data, rna_seq_data = import_data('Technical Test - Data Wrangling.xlsx')

# Munge patient clinical data
patient_clinical_data_munged = munge_patient_clinical_data(patient_clinical_data)

# Munge tissue sample data
tissue_sample_metadata_munged  = munge_tissue_sample_metadata(tissue_sample_metadata)

# Munge RNA_seq data
rna_seq = munge_rna_seq_data(rna_seq_data)

# Munge serum_protein data
serum_protein_data_munged = munge_serum_protein_data(serum_protein_data)

# Combine munged dataframes
final_database = merge_dfs(patient_clinical_data_munged, tissue_sample_metadata_munged, rna_seq, serum_protein_data_munged)

# Save database
save_final_database(final_database)

### Finally save merged dataframes as new sheet in origingal xlsx
## **** NEED TO LOOK UP SYNTAX

## Look. up the. __main__ initialize syntax