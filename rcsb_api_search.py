import os
import json
import requests
import pandas as pd


def get_cofactor_ids(path_to_jsons):
    cofactors_query = []
    for filename in os.listdir(path_to_jsons):
        if filename.endswith(".json"):
            if filename == 'cofactors.json':  ## need to parse cofactor separately
                cofactors_query.append(filename)

    cofactor_query_path = path_to_jsons + cofactors_query[0]
    cofactor = pd.read_json(cofactor_query_path, orient='index').rename(columns={0: 'dict'})
    cofactor = cofactor['dict'].apply(lambda x: pd.Series([x['EC'], x['cofactors']])).rename(
        columns={0: 'EC', 1: 'cofactors'}).reset_index().drop('EC', axis=1) \
        .rename(columns={'index': 'names'}).reset_index(drop=True)

    names_cofacts_dict = {}
    for index, row in cofactor.iterrows():
        names = row[0]
        names_cofacts_dict[names] = row[1]

    return names_cofacts_dict

def call_api(end_url, query):
    # ok_status = 200  # HTTP Status 200 (OK) status ; https://search.rcsb.org/redoc/index.html
    response = requests.post(end_url, json=query, )
    status = response.status_code
    no_entry_found = 204
    if status == no_entry_found:
        return None, status
    else:
        response_json = response.json()
        return response_json, status


def query_ligands(ligand_name, outpath, output_name, write_to_file=False):
    # print('here')
    true = True
    false = False

    ligand_query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_instance_feature_summary.comp_id",
                        "operator": "exact_match",
                        "value": f"{ligand_name}",
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_instance_feature_summary.type",
                        "operator": "exact_match",
                        "value": "HAS_COVALENT_LINKAGE"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_instance_feature_summary.count",
                        "operator": "equals",
                        "value": 0
                    }
                },

                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                              "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                              "operator": "exact_match",
                              "value": "Protein (only)"}
                },

                {  # resolution from 1-3.25A
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "range",
                        "value": {
                            "from": 1,
                            "include_lower": true,
                            "to": 3.25,
                            "include_upper": true
                        }
                    }
                }
            ],
        },
        "request_options": {
            # "paginate": {
            #     "start": 0,
            #     "rows": 100
            # }
            "return_all_hits": true,  # return all matches
        },
        "return_type": "non_polymer_entity"
        # Returns a list of PDB IDs appended with entity IDs in the format of a [pdb_id]_[entity_id], corresponding to non-polymeric entities (or ligands).
    }

    query_end_url = "https://search.rcsb.org/rcsbsearch/v2/query"

    rcsb_api_call = call_api(query_end_url, ligand_query)
    status = rcsb_api_call[1]
    # print(rcsb_api_call[0])

    ok_status = 200
    no_entry_found = 204
    no_entry_found_IDs = []
    # print(status)

    if status != ok_status and status != no_entry_found:
        print("query schema did not match")
    if status == no_entry_found:
        no_entry_found_IDs.append(ligand_name)

    if rcsb_api_call[0] is not None:
        response_json = rcsb_api_call[0]

        if write_to_file:
            if output_name == ligand_name:
                filename = outpath + ligand_name + "_query.json"
                print('writing to file:: ', filename)
            else:
                filename = outpath + output_name + "_" + ligand_name + "_query.json"
                print('writing cofactor to file:: ', filename)
            with open(filename, 'w') as f:
                json.dump(response_json, f)

    return no_entry_found_IDs



if __name__ == '__main__':

    write_to_file = True

    ## PArse through co-factors
    path_to_jsons = "/Users/Anushriya/Documents/ligands_cofactors/"
    cofactors_dict = get_cofactor_ids(path_to_jsons)

    number_of_ids_no_entry_dict = {}
    number_of_ids_no_entry = []
    for key, value in cofactors_dict.items():
        # print(type(key), key)
        cofact_call_list = value
        out_file_name = key.replace(' ', '')
        for i in cofact_call_list:
            cofactor_query = query_ligands(i, path_to_jsons, out_file_name, write_to_file=write_to_file)
            number_of_ids_no_entry.append(cofactor_query)
            number_of_ids_no_entry_dict[key] = number_of_ids_no_entry


