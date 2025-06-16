import pandas as pd
import requests
import time
import os
import re
import json # Import json for pretty printing

def get_uniprot_data(identifier: str):
    """
    Retrieves UniProt Accession, UniProt ID (Entry Name), Sequence, CBM Family,
    and Gene Names for a given identifier.
    It performs a two-stage search:
    1. Tries to find a UniProt Accession using various query strategies, requesting default fields.
    2. If a UniProt Accession is found, it uses that accession to retrieve
       all detailed information, including sequence and features/domains.
    
    Args:
        identifier (str): The identifier to search for in UniProt.
                          This can be a GenBank Accession (like BAA78554),
                          a UniProt Accession (like P0A7B8), or a protein name.
                          
    Returns:
        dict: A dictionary containing the retrieved protein data and a status message.
    """
    base_search_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # --- Stage 1: Find UniProt Accession ---
    # These queries will try to find a UniProt Accession.
    # We DO NOT request specific fields in this stage to avoid "Invalid fields" errors.
    # UniProt will return default fields which often include 'primaryAccession' or 'accession'.
    accession_finding_queries = [
        f"xref:GenBank-{identifier}",       # Primary attempt for GenBank IDs
        f"database:GenBank {identifier}",   # Alternative GenBank query
        identifier,                         # General keyword search (most likely to find TrEMBL primaryAccession)
        f"accession:{identifier}",          # Direct UniProt accession lookup (if it happens to be one)
        f"protein_name:{identifier}",       # Search by protein name
    ]
    
    found_uniprot_accession = None
    
    for query_string in accession_finding_queries:
        if not query_string.strip():
            continue

        params_stage1 = {
            "query": query_string,
            # IMPORTANT CHANGE: Removed the 'fields' parameter for Stage 1.
            # This allows UniProt to return its default set of fields,
            # which we then parse for 'accession' or 'primaryAccession'.
            "format": "json"
        }
        
        try:
            response = requests.get(base_search_url, params=params_stage1)
            response.raise_for_status() # Will raise HTTPError for 4xx/5xx responses
            data = response.json()

            if data and data.get('results'):
                # Try 'accession' first, then 'primaryAccession' from the default results
                result_entry = data['results'][0]
                found_uniprot_accession = result_entry.get('accession')
                if not found_uniprot_accession: 
                    found_uniprot_accession = result_entry.get('primaryAccession')

                if found_uniprot_accession:
                    print(f"  DEBUG: Stage 1 Success! Found UniProt Accession '{found_uniprot_accession}' with query '{query_string}'.", flush=True)
                    break # Found an accession, no need to try more Stage 1 queries
                else:
                    # This case should ideally not happen if a result was found
                    print(f"  DEBUG (Stage 1 - No Accession Extracted): Query '{query_string}' returned results but no expected accession was extracted. Raw data (first 200 chars): {json.dumps(data, indent=2)[:200]}...", flush=True)

            else:
                # Debug if results are empty for the current query string
                raw_data_str = json.dumps(data, indent=2) if data else "No data received"
                print(f"  DEBUG (Stage 1 - No results): Query '{query_string}' returned empty results. Raw data (first 200 chars): {raw_data_str[:200]}...", flush=True)
                
        except requests.exceptions.HTTPError as e:
            # Print HTTP error details and continue to next query in Stage 1
            print(f"  DEBUG (Stage 1 - HTTP Error): Query '{query_string}' for '{identifier}' failed with HTTP {e.response.status_code}. Response text (first 200 chars): {e.response.text[:200]}...", flush=True)
        except requests.exceptions.ConnectionError as e:
            # If a connection error occurs, it's likely a broader network issue, so return immediately.
            print(f"  DEBUG (Stage 1 - Connection Error): Connection error for '{identifier}' during query '{query_string}': {e}. Returning Connection Error.", flush=True)
            return {
                "UniProt_Accession": None, "UniProt_ID": None, "Amino_Acid_Sequence": None,
                "CBM_Family": None, "Gene_Names": None, "Status": "Connection Error"
            }
        except Exception as e:
            # Catch any other unexpected errors in Stage 1.
            print(f"  DEBUG (Stage 1 - General Error): An unexpected error occurred for '{identifier}' during query '{query_string}': {e}. Returning General Error.", flush=True)
            return {
                "UniProt_Accession": None, "UniProt_ID": None, "Amino_Acid_Sequence": None,
                "CBM_Family": None, "Gene_Names": None, "Status": f"Error: {str(e)}"
            }

    # If after all Stage 1 attempts, no UniProt Accession was found, return 'Not Found'.
    if not found_uniprot_accession:
        print(f"  DEBUG: Stage 1 failed for '{identifier}'. No UniProt Accession found after all attempts.", flush=True)
        return {
            "UniProt_Accession": None, "UniProt_ID": None, "Amino_Acid_Sequence": None,
            "CBM_Family": None, "Gene_Names": None, "Status": "Not Found after multiple attempts"
        }

    # --- Stage 2: Retrieve Full Details using Found UniProt Accession ---
    # IMPORTANT CHANGE: Use the /uniprotkb/{accession} endpoint for full details
    detail_url = f"https://rest.uniprot.org/uniprotkb/{found_uniprot_accession}"
    
    # IMPORTANT CHANGE: Removed the 'fields' parameter for Stage 2 as well.
    # The /uniprotkb/{accession} endpoint is designed to return all details by default.
    # Explicitly requesting fields seems to cause issues.
    
    params_stage2 = {
        "format": "json" # Only specify format for this endpoint
    }

    try:
        response = requests.get(detail_url, params=params_stage2) # Use detail_url and simplified params
        response.raise_for_status()
        data = response.json() # This response is the entry itself, not nested in a 'results' list

        # Parse the data directly from the response
        uniprot_accession = data.get('primaryAccession') or data.get('accession') # Robustly get accession
        uniprot_id = data.get('uniProtkbId') or data.get('id') # UniProtKB ID (entry name)
        sequence = data.get('sequence', {}).get('value') if data.get('sequence') else None # Sequence is often nested under 'sequence.value'

        # Extract gene names
        gene_names = ""
        if 'genes' in data and data['genes']:
            gene_values = []
            for gene_entry in data['genes']:
                if 'geneName' in gene_entry:
                    gene_values.append(gene_entry['geneName']['value'])
                if 'synonyms' in gene_entry:
                    for syn in gene_entry['synonyms']:
                        gene_values.append(syn['value'])
            gene_names = ", ".join(list(set(gene_values))) # Use set to avoid duplicates

        # Extract CBM Family from features
        cbm_family = ""
        if 'features' in data:
            for feature in data['features']:
                if feature.get('type') == 'Domain' and 'Carbohydrate-binding module' in feature.get('description', ''):
                    desc = feature.get('description')
                    match = re.search(r'family (\d+)', desc)
                    if match:
                        cbm_family = f"CBM{match.group(1)}"
                        break # Found a CBM family, typically one per protein

            print(f"  DEBUG: Stage 2 Success! Details retrieved for '{uniprot_accession}'.", flush=True)
            return {
                "UniProt_Accession": uniprot_accession,
                "UniProt_ID": uniprot_id,
                "Amino_Acid_Sequence": sequence,
                "CBM_Family": cbm_family,
                "Gene_Names": gene_names,
                "Status": f"Success (Found via '{found_uniprot_accession}')"
            }
        else:
            # This indicates an issue even after finding the accession, or UniProt no longer returns features for it.
            print(f"  DEBUG (Stage 2 - No results): Query to '{detail_url}' returned empty results or no valid data. Raw data (first 200 chars): {json.dumps(data, indent=2)[:200]}...", flush=True)
            return {
                "UniProt_Accession": found_uniprot_accession, "UniProt_ID": None, "Amino_Acid_Sequence": None,
                "CBM_Family": None, "Gene_Names": None, "Status": f"Found Accession but No Details (Accession: {found_uniprot_accession})"
            }

    except requests.exceptions.HTTPError as e:
        print(f"  DEBUG (Stage 2 - HTTP Error): Query to '{detail_url}' for '{identifier}' failed with HTTP {e.response.status_code}. Response text (first 200 chars): {e.response.text[:200]}...", flush=True)
        return {
            "UniProt_Accession": found_uniprot_accession, "UniProt_ID": None, "Amino_Acid_Sequence": None,
            "CBM_Family": None, "Gene_Names": None, "Status": f"HTTP Error in Stage 2 (Accession: {found_uniprot_accession})"
        }
    except requests.exceptions.ConnectionError as e:
        print(f"  DEBUG (Stage 2 - Connection Error): Connection error for '{identifier}' during query '{detail_url}': {e}.", flush=True)
        return {
            "UniProt_Accession": found_uniprot_accession, "UniProt_ID": None, "Amino_Acid_Sequence": None,
            "CBM_Family": None, "Gene_Names": None, "Status": "Connection Error in Stage 2"
        }
    except Exception as e:
        print(f"  DEBUG (Stage 2 - General Error): An unexpected error occurred for '{identifier}' during query '{detail_url}': {e}.", flush=True)
        return {
            "UniProt_Accession": None, "UniProt_ID": None, "Amino_Acid_Sequence": None,
            "CBM_Family": None, "Gene_Names": None, "Status": f"Error in Stage 2: {str(e)}"
        }

def process_cbm_list(input_csv_path: str, output_csv_path: str, id_column_name: str):
    """
    Reads a list of CBM IDs from an input CSV, retrieves UniProt data,
    and writes the combined data to an output CSV.
    
    Args:
        input_csv_path (str): Path to the input CSV file containing APR_IDs
                              and a column with identifiers for UniProt search.
        output_csv_path (str): Path where the enriched CSV will be saved.
        id_column_name (str): The name of the column in the input CSV that
                              contains the identifiers (e.g., 'GenBank_Accession',
                              'Protein_ID', or if directly APR IDs can be searched).
    """
    print(f"Loading data from {input_csv_path}...", flush=True)
    try:
        df = pd.read_csv(input_csv_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_csv_path}", flush=True)
        return
    except Exception as e:
        print(f"Error reading input CSV: {e}", flush=True)
        return

    if id_column_name not in df.columns:
        print(f"Error: Column '{id_column_name}' not found in the input CSV.", flush=True)
        print(f"Available columns: {df.columns.tolist()}", flush=True)
        return

    results = []
    total_proteins = len(df)
    print(f"Processing {total_proteins} proteins...", flush=True)

    for index, row in df.iterrows():
        # Ensure 'APR_ID' column exists in your input CSV.
        # If your primary ID column is different (e.g., just 'Protein_ID'),
        # adjust this line to match your actual column name.
        apr_id = row['APR_ID'] if 'APR_ID' in df.columns else f"Protein_{index+1}" 
        
        # Extract the actual identifier from the 'Header' column,
        # assuming it's the part before the first '|'
        raw_identifier = str(row[id_column_name]).strip() if pd.notna(row[id_column_name]) else ""
        
        # Split by '|' and take the first part
        identifier_for_uniprot = raw_identifier.split('|')[0].strip() if raw_identifier else ""
        
        if not identifier_for_uniprot:
            print(f"({index + 1}/{total_proteins}) Skipping {apr_id}: Missing or unparseable identifier (Raw Header: '{raw_identifier}').", flush=True)
            combined_row = row.to_dict()
            combined_row.update({
                "UniProt_Accession": None, "UniProt_ID": None, "Amino_Acid_Sequence": None,
                "CBM_Family": None, "Gene_Names": None, "Status": "Skipped (Empty/Invalid Identifier)"
            })
            results.append(combined_row)
            continue

        # Print the extracted ID and processing message without a newline
        print(f"({index + 1}/{total_proteins}) Processing {apr_id} (Extracted ID: '{identifier_for_uniprot}')... ", end='', flush=True) 
        
        # Call the UniProt retrieval function with the extracted identifier
        uniprot_info = get_uniprot_data(identifier_for_uniprot)

        # Print the status immediately after getting info, completing the line
        print(f"Status: {uniprot_info['Status']}", flush=True)

        # Combine original row data with new UniProt data
        combined_row = row.to_dict()
        combined_row.update(uniprot_info)
        results.append(combined_row)
        
        # Be polite to the UniProt server and avoid being blocked
        time.sleep(0.1) 

    # Create a new DataFrame from the results
    output_df = pd.DataFrame(results)

    # Save the enriched data to a new CSV file
    print(f"Saving enriched data to {output_csv_path}...", flush=True)
    output_df.to_csv(output_csv_path, index=False)
    print("Processing complete!", flush=True)

if __name__ == "__main__":
    input_file = "my_cbm_data.csv"
    output_file = "enriched_cbm_data.csv"
    column_for_uniprot_search = "Header" 
    
    # Dummy data creation block is commented out
    # if not os.path.exists(input_file):
    #     print(f"Creating a dummy input file '{input_file}' for demonstration.")
    #     print("Please replace this dummy file with your actual 'my_cbm_data.csv' before proceeding.")
    #     dummy_data = {
    #         'APR_ID': ['APR1', 'APR2', 'APR3', 'APR4', 'APR5', 'APR6', 'APR7', 'APR8', 'APR9', 'APR10'],
    #         'Header': [
    #             'BAA78554|Paramecium bursaria Chlorella virus CVK2|1-110', 
    #             'P0A7B8',     
    #             'A0A024RBG1', 
    #             'CELA_BACSU', 
    #             'NONEXISTENT_ID|Some Other Info', 
    #             'P07297|Example Protein A', 
    #             'Q9I6L1|Some Protein B', 
    #             'A7A786|Protein C from Organism D', 
    #             'Q83T69|Another Protein', 
    #             'P54530|Final Protein' 
    #         ]
    #     }
    #     pd.DataFrame(dummy_data).to_csv(input_file, index=False)
    #     print("Dummy input file created. Please ensure it matches your expected format.")
    #     # import sys; sys.exit(1)

    process_cbm_list(input_file, output_file, column_for_uniprot_search)
