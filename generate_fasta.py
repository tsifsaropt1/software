import pandas as pd
import os
import re # Make sure re is imported for cleaning sequences

def generate_fasta_from_csv(input_csv_path: str, output_fasta_path: str, id_column: str, cbm_family_column: str, sequence_column: str):
    """
    Reads a CSV file containing protein IDs (now UniProt Accessions), CBM families,
    and CBM domain sequences, then generates a single FASTA file.

    Args:
        input_csv_path (str): Path to the input CSV file (e.g., 'enriched_cbm_data.csv').
        output_fasta_path (str): Path for the output FASTA file (e.g., 'cbm_domains_for_alphafold.fasta').
        id_column (str): The name of the column in the CSV containing unique protein IDs (UniProt Accessions).
        cbm_family_column (str): The name of the column containing CBM family information.
        sequence_column (str): The name of the column containing the CBM domain amino acid sequences.
    """
    print(f"Loading data from {input_csv_path}...", flush=True)
    try:
        df = pd.read_csv(input_csv_path)
    except FileNotFoundError:
        print(f"Error: Input CSV file not found at '{input_csv_path}'. Please check the path and filename.", flush=True)
        return
    except Exception as e:
        print(f"Error reading CSV file: {e}", flush=True)
        return

    # Check if necessary columns exist
    for col in [id_column, cbm_family_column, sequence_column]:
        if col not in df.columns:
            print(f"Error: Required column '{col}' not found in the CSV. Available columns: {df.columns.tolist()}", flush=True)
            return

    fasta_entries = []
    skipped_count = 0

    print(f"Generating FASTA entries for {len(df)} proteins...", flush=True)
    for index, row in df.iterrows():
        # Get UniProt Accession for the FASTA header
        protein_id = row[id_column] 
        cbm_family = row[cbm_family_column]
        sequence = row[sequence_column]

        # Ensure UniProt Accession, CBM family, and sequence are valid
        if pd.isna(protein_id) or not str(protein_id).strip():
            print(f"  Skipping row {index + 1}: Missing or invalid UniProt Accession. (Original ID might be {row.get('APR_ID', 'N/A')}).", flush=True)
            skipped_count += 1
            continue

        if pd.isna(cbm_family) or not str(cbm_family).strip():
             cbm_family_str = "UnknownCBM" # Use a placeholder if family is missing
        else:
            cbm_family_str = str(cbm_family).strip().replace(" ", "") # Remove spaces if any

        if pd.isna(sequence) or not str(sequence).strip():
            print(f"  Skipping row {index + 1}: Empty or invalid sequence for UniProt Accession '{protein_id}'.", flush=True)
            skipped_count += 1
            continue
        
        # Ensure sequence is a string and clean it
        sequence = str(sequence).strip()
        cleaned_sequence = re.sub(r'[^A-Z]', '', sequence.upper()) 
        
        if not cleaned_sequence:
            print(f"  Skipping row {index + 1}: Sequence for UniProt Accession '{protein_id}' became empty after cleaning.", flush=True)
            skipped_count += 1
            continue

        # Create FASTA header: >UniProtAccession_CBMFamily (e.g., >P54583_CBM2)
        fasta_header = f">{protein_id}_{cbm_family_str}"
        fasta_entries.append(fasta_header)
        fasta_entries.append(cleaned_sequence) # Add the cleaned sequence

    if not fasta_entries:
        print("No valid FASTA entries were generated. Please check your input CSV data.", flush=True)
        return

    # Write all entries to the output FASTA file
    with open(output_fasta_path, 'w') as f:
        f.write('\n'.join(fasta_entries))

    print(f"\nFASTA file '{output_fasta_path}' created successfully with {len(fasta_entries) // 2} entries.", flush=True)
    if skipped_count > 0:
        print(f"Note: {skipped_count} entries were skipped due to missing or invalid data.", flush=True)
    print("You can now use this FASTA file directly with ColabFold or other protein structure prediction tools.", flush=True)

if __name__ == "__main__":
    # --- Configuration ---
    # IMPORTANT: Now using 'enriched_cbm_data.csv' as input,
    # as it contains the UniProt Accessions needed for FASTA headers.
    input_csv_file = "enriched_cbm_data.csv" 
    
    # The name of the FASTA file to be created - now more indicative.
    output_fasta_file = "cbm_domains_for_alphafold.fasta" 

    # IMPORTANT: Ensure these column names exactly match your 'enriched_cbm_data.csv' header
    protein_id_col = "UniProt_Accession" # This column must contain the UniProt Accession
    cbm_family_col = "CBM_Family" # This column must contain the CBM Family (e.g., 'CBM2', 'CBM64')
    sequence_col = "Amino_Acid_Sequence" # This column must contain the CBM domain amino acid sequence

    # --- Run the FASTA generation ---
    generate_fasta_from_csv(input_csv_file, output_fasta_file, protein_id_col, cbm_family_col, sequence_col)

