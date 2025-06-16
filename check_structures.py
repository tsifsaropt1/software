import pandas as pd
import requests
import time
import os
import json
from botocore.exceptions import ClientError # Import ClientError for more specific error handling

def download_pdb_file(protein_id: str, structure_id: str, source_type: str, base_download_dir: str) -> str | None:
    """
    Downloads a PDB or AlphaFoldDB structure file into a specific subdirectory.
    
    Args:
        protein_id (str): The UniProt Accession of the protein (used for naming the file).
        structure_id (str): The PDB ID (e.g., "1ECE") or AlphaFoldDB ID (which is the UniProt Accession).
        source_type (str): "PDB Available" or "AlphaFoldDB Available".
        base_download_dir (str): The *base* local directory (e.g., 'CBM Structures').
                                  Subdirectories ('PDB', 'AlphaFold DB') will be created/used within it.
        
    Returns:
        str: The local path to the downloaded file, or None if download fails.
    """
    specific_download_dir = ""
    file_name = ""
    download_url = ""

    if source_type == "PDB Available":
        specific_download_dir = os.path.join(base_download_dir, "PDB") # Changed to 'PDB'
        # PDB IDs can sometimes be comma-separated (e.g., "1ECE, 1VRX").
        # We'll download the first one listed, or iterate if needed later.
        first_pdb_id = structure_id.split(',')[0].strip()
        file_name = f"{protein_id}_{first_pdb_id}.pdb"
        download_url = f"https://files.rcsb.org/download/{first_pdb_id}.pdb"
        print(f"  Attempting to download PDB model for {protein_id} ({first_pdb_id})... ", end='', flush=True)
    elif source_type == "AlphaFoldDB Available":
        specific_download_dir = os.path.join(base_download_dir, "AlphaFold DB") # Changed to 'AlphaFold DB'
        # For AlphaFoldDB, the 'AlphaFoldDB_ID' from the CSV is actually the UniProt Accession.
        file_name = f"{protein_id}_AF.pdb"
        # Standard AlphaFoldDB URL format, typically model_v4 is the latest and best
        download_url = f"https://alphafold.ebi.ac.uk/files/AF-{structure_id}-F1-model_v4.pdb"
        print(f"  Attempting to download AlphaFoldDB model for {protein_id} ({structure_id})... ", end='', flush=True)
    else:
        print(f"  Unknown source type for download: {source_type}", flush=True)
        return None

    # Ensure the specific download directory exists
    if not os.path.exists(specific_download_dir):
        os.makedirs(specific_download_dir)
        print(f"  Created download subdirectory: {specific_download_dir}", flush=True)

    local_file_path = os.path.join(specific_download_dir, file_name)

    try:
        response = requests.get(download_url, stream=True)
        response.raise_for_status() # Raise an HTTPError for bad responses (4xx or 5xx)

        with open(local_file_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(" Success.", flush=True)
        return local_file_path
    except requests.exceptions.HTTPError as e:
        print(f" Failed (HTTP Error {e.response.status_code}).", flush=True)
        return None
    except requests.exceptions.ConnectionError as e:
        print(f" Failed (Connection Error).", flush=True)
        return None
    except Exception as e:
        print(f" Failed (Unexpected Error: {e}).", flush=True)
        return None


def check_existing_structure(uniprot_accession: str):
    """
    Checks UniProt for existing PDB or AlphaFoldDB structures for a given UniProt Accession.
    
    Args:
        uniprot_accession (str): The UniProt Accession ID of the protein (e.g., P54583, Q60029).
                                 This function expects a valid UniProt Accession, not an internal ID like 'APR1'.
                                 
    Returns:
        dict: A dictionary containing status and available structure IDs.
    """
    # Use the specific endpoint for detailed UniProt entry retrieval
    detail_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_accession}"
    
    # IMPORTANT CHANGE: Removed the 'fields' parameter entirely for this endpoint.
    # The /uniprotkb/{accession} endpoint is designed to return all details by default.
    # Explicitly requesting fields seems to cause issues for some entries.
    params = {
        "format": "json" 
    }

    try:
        response = requests.get(detail_url, params=params)
        response.raise_for_status() # Raise an HTTPError for bad responses (4xx or 5xx)
        data = response.json()

        pdb_ids = []
        # 'uniProtKBCrossReferences' is the correct field for cross-references like PDB
        if 'uniProtKBCrossReferences' in data:
            for xref in data['uniProtKBCrossReferences']:
                if xref.get('database') == 'PDB':
                    pdb_entry_id = xref.get('id')
                    if pdb_entry_id:
                        pdb_ids.append(pdb_entry_id)

        alphafold_id = None
        if 'uniProtKBCrossReferences' in data:
            for xref in data['uniProtKBCrossReferences']:
                if xref.get('database') == 'AlphaFoldDB':
                    alphafold_id = xref.get('id') # This should be the AF-ID if available
                    break

        if pdb_ids:
            return {
                "Status": "PDB Available",
                "PDB_ID": ", ".join(pdb_ids),
                "AlphaFoldDB_ID": alphafold_id if alphafold_id else "N/A"
            }
        elif alphafold_id:
            return {
                "Status": "AlphaFoldDB Available",
                "PDB_ID": "N/A",
                "AlphaFoldDB_ID": alphafold_id
            }
        else:
            return {
                "Status": "Needs Prediction",
                "PDB_ID": "N/A",
                "AlphaFoldDB_ID": "N/A"
            }

    except requests.exceptions.HTTPError as e:
        # Check if 404 Not Found, means accession itself wasn't found
        if e.response.status_code == 404:
            print(f"  DEBUG: UniProt Accession '{uniprot_accession}' Not Found (HTTP 404).", flush=True)
            return {
                "Status": "UniProt Accession Not Found",
                "PDB_ID": "N/A",
                "AlphaFoldDB_ID": "N/A"
            }
        print(f"  DEBUG: HTTP Error for {uniprot_accession}: {e.response.status_code} - {e.response.text[:200]}...", flush=True)
        return {
            "Status": f"HTTP Error ({e.response.status_code})",
            "PDB_ID": "N/A",
            "AlphaFoldDB_ID": "N/A"
        }
    except requests.exceptions.ConnectionError as e:
        print(f"  DEBUG: Connection Error for {uniprot_accession}: {e}", flush=True)
        return {
            "Status": "Connection Error",
            "PDB_ID": "N/A",
            "AlphaFoldDB_ID": "N/A"
        }
    except Exception as e:
        print(f"  DEBUG: Unexpected Error for {uniprot_accession}: {e}", flush=True)
        return {
            "Status": f"Error: {str(e)}",
            "PDB_ID": "N/A",
            "AlphaFoldDB_ID": "N/A"
        }

def process_cbm_for_structure_check(input_csv_path: str, output_csv_path: str, protein_id_col: str, 
                                    download_structures_flag: bool = False, base_download_directory: str = "downloaded_cbm_structures"):
    """
    Reads a CSV with protein IDs (which must be UniProt Accessions for this script to work),
    checks for existing structures, and saves the status to a new CSV.
    Optionally downloads the structure files into categorized subfolders.

    Args:
        input_csv_path (str): Path to the input CSV file.
        output_csv_path (str): Path where the output CSV with structure status will be saved.
        protein_id_col (str): The name of the column in the input CSV that contains
                              the UniProt Accessions.
        download_structures_flag (bool): If True, attempts to download PDB/AlphaFoldDB files.
        base_download_directory (str): The *base* folder where subfolders ('PDB_Structures',
                                       'AlphaFoldDB_Structures') for downloaded .pdb files will be created.
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

    if protein_id_col not in df.columns:
        print(f"Error: Column '{protein_id_col}' (containing UniProt Accessions) not found in the input CSV.", flush=True)
        print(f"Available columns: {df.columns.tolist()}", flush=True)
        return

    results = []
    total_proteins = len(df)
    print(f"Checking structure availability for {total_proteins} proteins...", flush=True)

    for index, row in df.iterrows():
        # Ensure that current_uniprot_accession is indeed a UniProt Accession (e.g., P54583),
        # not an internal ID like 'APR1'. The script relies on this.
        current_uniprot_accession = str(row[protein_id_col]).strip()
        cbm_family_val = row["CBM_Family"] if "CBM_Family" in row else "N/A" # Get CBM Family

        if pd.isna(current_uniprot_accession) or not current_uniprot_accession:
            print(f"({index + 1}/{total_proteins}) Skipping entry: Missing or invalid UniProt Accession. (Original ID might be {row.get('APR_ID', 'N/A')}).", flush=True)
            results.append({
                "Original_Protein_ID": row[protein_id_col] if protein_id_col in row else "N/A",
                "CBM_Family": cbm_family_val,
                "Structure_Status": "Skipped (Missing/Invalid Accession)",
                "PDB_ID": "N/A",
                "AlphaFoldDB_ID": "N/A",
                "Downloaded_File_Path": "N/A"
            })
            continue

        print(f"({index + 1}/{total_proteins}) Checking '{current_uniprot_accession}'... ", end='', flush=True)
        
        status_info = check_existing_structure(current_uniprot_accession)

        downloaded_file_path = "N/A"
        if download_structures_flag:
            if status_info["Status"] == "PDB Available":
                downloaded_file_path = download_pdb_file(current_uniprot_accession, status_info["PDB_ID"], "PDB Available", base_download_directory)
            elif status_info["Status"] == "AlphaFoldDB Available":
                downloaded_file_path = download_pdb_file(current_uniprot_accession, status_info["AlphaFoldDB_ID"], "AlphaFoldDB Available", base_download_directory)
            # For "Needs Prediction" or errors, downloaded_file_path remains "N/A"

        result_row = {
            "Original_Protein_ID": current_uniprot_accession,
            "CBM_Family": cbm_family_val,
            "Structure_Status": status_info["Status"],
            "PDB_ID": status_info["PDB_ID"],
            "AlphaFoldDB_ID": status_info["AlphaFoldDB_ID"],
            "Downloaded_File_Path": downloaded_file_path
        }
        results.append(result_row)
        print(f"Status: {status_info['Status']}", flush=True)
        
        time.sleep(0.1) # Be polite to the UniProt API

    output_df = pd.DataFrame(results)
    print(f"\nSaving structure availability data to {output_csv_path}...", flush=True)
    output_df.to_csv(output_csv_path, index=False)
    print("Structure availability check and download complete!", flush=True)
    if download_structures_flag:
        print(f"Downloaded structures can be found in subfolders within '{base_download_directory}' folder (e.g., '{base_download_directory}/PDB' and '{base_download_directory}/AlphaFold DB').")

if __name__ == "__main__":
    # --- Configuration ---
    # IMPORTANT: You must ensure the 'input_csv_file' and 'protein_id_for_check_col'
    # correctly point to the CSV file and the column within it that contains
    # the ACTUAL UniProt Accessions (e.g., P54583, Q60029), NOT your internal APR IDs (e.g., APR1).

    # Option A (Recommended): Use your 'enriched_cbm_data.csv' as input
    # if it contains a 'UniProt_Accession' column with valid UniProt Accessions.
    input_csv_file = "enriched_cbm_data.csv"
    protein_id_for_check_col = "UniProt_Accession"

    # Option B (Only if you are certain that your 'cbm_sequences_for_fasta.csv'
    # now has ACTUAL UniProt Accessions in its 'Protein ID' column):
    # input_csv_file = "cbm_sequences_for_fasta.csv"
    # protein_id_for_check_col = "Protein ID"

    # The output CSV file that will tell you which structures are available.
    output_csv_file = "cbm_structure_availability.csv" 

    # --- Download Settings ---
    # Set to True if you want to automatically download the PDB/AlphaFoldDB files.
    download_structures_enabled = True 
    # The BASE directory where downloaded .pdb files will be saved.
    # Subfolders 'PDB' and 'AlphaFold DB' will be created/used inside this.
    # This should be relative to where you run the script, or an absolute path.
    base_download_folder = "CBM Structures" # Matches the user's existing folder name
    
    # --- Run the structure availability check and optional download ---
    process_cbm_for_structure_check(
        input_csv_file,
        output_csv_file,
        protein_id_for_check_col,
        download_structures_enabled,
        base_download_folder
    )

