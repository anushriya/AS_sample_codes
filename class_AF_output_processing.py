import os
import sys
import json
import math
import warnings
from pathlib import Path
import pandas as pd
import numpy as np
from pymol import cmd
import matplotlib.pyplot as plt
from UtilityFunctions import *
from postprocess_pdbs import write_filter_pdb

desired_width = 320
pd.set_option('display.width', desired_width)
pd.set_option('display.max_columns', 10)

warnings.filterwarnings("ignore", category=RuntimeWarning)


class ParseAFoutputs:
    def __init__(self, json_file_path, pkl_file_path, complex_name, complex_seq_len_dict):
        self.json_file_path = json_file_path
        self.pkl_file_path = pkl_file_path
        self.complex_name = complex_name
        self.complex_seq_len_dict = complex_seq_len_dict

    def get_top_rankings(self):
        """
        Parse the ranking_debug.json outputs from AlphaFold.
        Args:
            json_file_path:Path to json file
            complex_name: name of the complex

        Returns: 1) ranked_list = ordered list of ranked models,
                 2) ordered_iptms_dict = the ordered dictionary of all models+iPTM scores
                 3) top_rank_model_score = a list of top ranking model with score,
                 4) top_model_name[0] = the name of top ranking model
        """
        ranked_list = []
        json_file = self.json_file_path / 'ranking_debug.json'
        if json_file.exists():
            try:
                with open(json_file, 'r') as file:
                    # print('current dir:', sub_dir, 'json_found:', json_file.exists())
                    data = json.load(file)
                    iptms_list = data.get('iptm+ptm',
                                          {})  # Use get method to provide a default value if the key is not present
                    ranked_list = [v for k, v in data.items() if k != 'iptm+ptm']
            except Exception as e:
                print(f"Error loading JSON file {json_file}: {e}")
        else:
            print(f"File not found: {json_file}")

        top_rank_model_score = ()
        if ranked_list:
            top_ranked_model_name = ranked_list[0][0]
            if top_ranked_model_name in iptms_list.keys():
                top_rank_model_score = (self.complex_name, top_ranked_model_name, iptms_list[top_ranked_model_name])
            ordered_iptms_dict = {key: iptms_list[key] for key in ranked_list[0]}
            top_model_name = next(iter(ranked_list))
        else:
            ordered_iptms_dict = {}
            top_model_name = [None]

        return ranked_list, ordered_iptms_dict, top_rank_model_score, top_model_name[0]

    def all_models_pae(self):
        """
        Extract the various metrics from all models
        Args:
            pkl_file_path: path to al *.pkl file

        Returns: A dictionary for each model, key: model name, values: (plddt:number, pae: number,ptm:number,ranking_confidence: number
                                                                        max_predicted_aligned_error:number)
        """
        model_names = []
        try:
            for filename in os.listdir(self.pkl_file_path):
                # if filename.startswith("result_") and filename.endswith(".pkl"):
                if filename.startswith("clean_result_") and filename.endswith(".pkl"):
                    full_path = os.path.join(self.pkl_file_path, filename)
                    if os.path.exists(full_path):
                        model_names.append(filename)
                    else:
                        print(f"File not found: {full_path}")
        except OSError as e:
            print(f"Error while listing files: {e}")

        out = {}
        for i, name in enumerate(model_names):
            name = os.path.join(self.pkl_file_path, name)
            d = U.read_pkl(name)

            basename = os.path.basename(name)
            basename = basename[basename.index('model'):]

            for k, v in d.items():
                out[f'{basename}'] = {'plddt': v['plddt'],
                                      'pae': v['pae'],
                                      'ptm': v['ptm'],
                                      'iptm': v['iptm'],
                                      'ranking_confidence': v['ranking_confidence'],
                                      # 'max_predicted_aligned_error': d['max_predicted_aligned_error'],
                                      }

        return out

    def clean_comp_lenght_dict(self, negative_interaction=False):
        """
        Args:
            complex_name: name of the current complex
            complex_seq_len_dict: Dictionary of complexes and their sequence
            negative_interaction: If True - processing negative non-interacting complex ,
                                  If False= Processing positive interaction complex

        Returns: a nested list of list, each containing length of sequence for chain A and B for each complex.
        [[chainA end res number, chainB end res number]]

        """
        if negative_interaction:
            temp_protein_lenght_dict = {}
            split_name = self.complex_name.split("_")
            for key_tuple, value_list in complex_seq_len_dict.items():
                if len(key_tuple) == len(value_list):
                    for i in range(len(key_tuple)):
                        # print('key_tuple[i]', key_tuple[i])
                        temp_protein_lenght_dict[key_tuple[i]] = value_list[i]
                if len(value_list) == 1:
                    for i in key_tuple:
                        # print('i', i)
                        temp_protein_lenght_dict[i] = value_list[0]

            temp_tup = [None] * len(split_name)
            for i in range(len(split_name)):
                temp_tup[i] = temp_protein_lenght_dict.get(split_name[i])
            flat_comp_lengths = temp_tup

        else:
            curr_comp_lenghts = []
            split_name = self.complex_name.split("_")
            for k, v in self.complex_seq_len_dict.items():
                if k[0] == split_name[0] and k[1] == split_name[1] and split_name[0] == split_name[1]:  # homodimers
                    curr_comp_lenghts.append((v[0], v[0]))
                elif k[0] == split_name[0] and k[1] == split_name[1]:  # match the order of names, heterdimers
                    curr_comp_lenghts.append((v[0], v[1]))
            flat_comp_lengths = [item for tpl in curr_comp_lenghts for item in tpl]

        return flat_comp_lengths

    def get_residues(self, pae_per_model, ranked_list, threshold=90):
        """
        Filter the plddt for entries with score above given threshold
        Args:
            pae_per_model: The output from unpacking pkl (a dictionary with 'pae' and 'plddt' score
            ranked_list: List of ranked_models

            threshold: The min score cutoff value of plddt

        Returns: 1) a list (indexes) of residues at or above the threshold score
                2) the same list but with added residues based on res_cutoff
        """
        top_ranked_model_name = next(iter(ranked_list))

        indexes = []
        for model_name, value in pae_per_model.items():
            if model_name == top_ranked_model_name + ".pkl":
                greater_threshold_plddts = value['plddt'].copy()
                indexes = np.where(greater_threshold_plddts >= threshold)  # index in list=resiude number

        return indexes[0]


class GetAvgs(ParseAFoutputs):
    def __init__(self, index_overthreshold):
        super().__init__(json_file, pkl_file, complex_name, complex_seq_len_dict)
        self.index_overthreshold = index_overthreshold

    def get_avg_pae_plddt(self, pae_plddt_per_model,
                          comp_dict_clean,
                          top_ranked_model_name,
                          # index_overthreshold
                          ):
        """
        Get average Pae and plddt for top ranked model,
        Args:
            pae_plddt_per_model:
            complex_name:
            comp_dict_clean:
        Returns: 1) concat_avg_pae: list with 3 elements (1: chain A overall AVG,
                                                2) chain B overall AVG,
                                                3) off diagonal overall average
                 2) concat_avg_plddt: list with nested tuple [(chain A avg, pass/fail), (chain B avg, pass/fail)]
                 3) num_res_above_90: Tuple with (Chain A: number of residues with plddt>=90, Chain B: number of residues with plddt>=90,)
                 4) percent_res_above_90: Tuple with (Chain A: percent of residues with plddt>=90, Chain B: percent of residues with plddt>=90,)
        """

        chain_A = comp_dict_clean[0]
        chain_A_filtered_pae_avg = []
        chain_B_filtered_pae_avg = []
        chain_AB_top_filtered_pae_avg = []
        chain_AB_bottom_filtered_pae_avg = []

        ## pLDDT SCORES for selected res
        chain_A_plddt_scores_selected = []
        chain_B_plddt_scores_selected = []

        for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
            if model_name == top_ranked_model_name + ".pkl":
                current_pae = value['pae']
                # Parsing pae array for average pae per chains/off-diagonal
                chainA_only_pae = current_pae[0:chain_A, 0:chain_A]
                chainB_only_pae = current_pae[chain_A + 1:, chain_A + 1:]
                chain_AB_topdiagonal = current_pae[0:chain_A:, chain_A + 1:]
                chainAB_bottomdiagonal = current_pae[chain_A + 1:, 0:chain_A]

                # print("self.index_overthreshold[0]", self.index_overthreshold)
                for res in self.index_overthreshold:
                    # "Chain A pae average"
                    if res <= chain_A:
                        sub_array = chainA_only_pae[res - 1]  # counting offset (index =0, res=1..)
                        subarray_filter = [x for x in sub_array if x <= 15]
                        filtered_avg = np.nanmean(subarray_filter)
                        chain_A_filtered_pae_avg.append(filtered_avg)

                        # "Chain A/B TOP Diagonal pae average"
                        sub_array = chain_AB_topdiagonal[res - 1]
                        subarray_filter = [x for x in sub_array if x <= 15]
                        drop_nan_subarray = [x for x in subarray_filter if not np.isnan(x)]
                        AB_top_filtered_avg = np.nanmean(drop_nan_subarray)
                        chain_AB_top_filtered_pae_avg.append(AB_top_filtered_avg)

                        # "Avg plddt for res  in chainA"
                        temp_chainA_plddt = value['plddt'][res]
                        chain_A_plddt_scores_selected.append(temp_chainA_plddt)
                        # print("temp_chainA_plddt", temp_chainA_plddt)

                    # "Chain B pae average"
                    if res > chain_A:
                        chainB_res = res - (chain_A + 1)
                        sub_array = chainB_only_pae[chainB_res - 1]
                        subarray_filter = [x for x in sub_array if x <= 15]
                        filtered_avg = np.nanmean(subarray_filter)
                        chain_B_filtered_pae_avg.append(filtered_avg)

                        # "Chain A/B Bottom Diagonal pae average"
                        sub_array = chainAB_bottomdiagonal[chainB_res - 1]
                        subarray_filter = [x for x in sub_array if x <= 15]
                        AB_bottom_filtered_avg = np.nanmean(subarray_filter)
                        chain_AB_bottom_filtered_pae_avg.append(AB_bottom_filtered_avg)

                        # "Avg plddt for res  in chainB"
                        temp_chainB_plddt = value['plddt'][res]
                        chain_B_plddt_scores_selected.append(temp_chainB_plddt)

        """creates a tuple of (chain_avg, pass/fail flag) if residue is over plddt_threshold=pass, else fail"""
        """plddt lists checks if at least one value  is >=90, if so adds 'pass' flag, else 'fail' flag"""

        chain_A_plddt_avg = np.mean(chain_A_plddt_scores_selected)
        chain_B_plddt_avg = np.mean(chain_B_plddt_scores_selected)

        avg_chain_A_plddt = [
            (chain_A_plddt_avg, 'pass' if any(x >= 90 for x in chain_A_plddt_scores_selected) else 'fail')]
        avg_chain_B_plddt = [
            (chain_B_plddt_avg, 'pass' if any(x >= 90 for x in chain_B_plddt_scores_selected) else 'fail')]
        concat_avg_plddt = avg_chain_A_plddt + avg_chain_B_plddt  # list with nested tuple [(chain a avg, pass/fail), (chain b avg, pass/fail)]

        # "Number of residues with plddt score>=90"
        num_res_above_90 = (
            sum(1 for x in chain_A_plddt_scores_selected if x >= 90),
            sum(1 for x in chain_B_plddt_scores_selected if x >= 90))

        chainA_percent_res_plddtabove90 = (num_res_above_90[0] / chain_A) * 100
        chainB_percent_res_plddtabove90 = (num_res_above_90[1] / comp_dict_clean[
            1]) * 100  # comp_dict_clean[1]=chain B sequence length

        # "Percentage of residues with plddt score>=90"
        percent_res_above_90 = (chainA_percent_res_plddtabove90, chainB_percent_res_plddtabove90)

        clean_chain_A_filtered_pae_avg = np.mean([x for x in chain_A_filtered_pae_avg if not math.isnan(x)])
        clean_chain_B_filtered_pae_avg = np.mean([x for x in chain_B_filtered_pae_avg if not math.isnan(x)])

        clean_chain_AB_top_pae_avg = [x for x in chain_AB_top_filtered_pae_avg if not math.isnan(x)]
        clean_chain_AB_bottom_pae_avg = [x for x in chain_AB_bottom_filtered_pae_avg if not math.isnan(x)]

        chain_AB_BA_summed = clean_chain_AB_top_pae_avg + clean_chain_AB_bottom_pae_avg
        clean_chain_AB_BA_summed_avg = np.mean(chain_AB_BA_summed)

        concat_avg_pae = [clean_chain_A_filtered_pae_avg] + [clean_chain_B_filtered_pae_avg] + [
            clean_chain_AB_BA_summed_avg]

        return concat_avg_pae, concat_avg_plddt, num_res_above_90, percent_res_above_90

    def all_metrics(self, avg_pae,
                    avg_plddt,
                    numres_plddtover90,
                    percent_plddtover90,
                    pae_plddt_per_model,
                    top_model_name, save_path,
                    threshold=50,
                    aim1=False,
                    aim2=True
                    ):
        """
        Gather all AF metrics and make a dataframe and write to txt file.
        Args:
            complex_name: name of complex
            avg_pae: average pae score of complex (chain A, chain B, Diagonal pae)
            avg_plddt: average plddt score (chain A, chain B)
            numres_plddtover90: Number of residues with plddt over 90% (chain A, chain B)
            percent_plddtover90: Percents of residues with plddt over 90% (chain A, chain B)
            pae_plddt_per_model: Dictionary of all other metrics (iptm, ranking confidence)
            top_model_name: Top ranked model name
            save_path: directory path to save text file with all the AF metrics

        Returns:none

        """

        iptm_score = 0
        ptm_score = 0
        ranking_confidence = 0
        max_predicted_aligned_error = 0
        for model_name, value in pae_plddt_per_model.items():
            if top_model_name == model_name[:-4]:
                iptm_score = value['iptm']
                ptm_score = value['ptm']
                ranking_confidence = value['ranking_confidence']  # .80iptm+.20ptm score
                # max_predicted_aligned_error = value['max_predicted_aligned_error']

        chainA_avg_pae = avg_pae[0]
        # print("chainA_avg_pae", chainA_avg_pae)
        chainB_avg_pae = avg_pae[1]
        diag_AB_avg_pae = avg_pae[2]

        chainA_plddt_avg = avg_plddt[0][0]
        chainA_plddt_qual = avg_plddt[0][1]

        chainB_plddt_avg = avg_plddt[1][0]
        chainB_plddt_qual = avg_plddt[1][1]

        chainA_numres_over90 = numres_plddtover90[0]
        chainB_numres_over90 = numres_plddtover90[1]

        chainA_percentres_over90 = percent_plddtover90[0]
        chainB_percentres_over90 = percent_plddtover90[1]

        # percent_plddtover90 = numres_plddtover90

        metrics_df = pd.DataFrame()
        # metrics_df['complex_name'] = [complex_name]
        metrics_df['complex_name'] = [self.complex_name]

        metrics_df['iptm'] = [iptm_score]
        metrics_df['ptm'] = [ptm_score]
        metrics_df['ranking_confidence'] = [ranking_confidence]
        # metrics_df['max_predicted_aligned_error'] = [max_predicted_aligned_error]

        metrics_df['chainA_avg_pae'] = [chainA_avg_pae]
        metrics_df['chainB_avg_pae'] = [chainB_avg_pae]
        metrics_df['diagonal_AB_avg_pae'] = [diag_AB_avg_pae]

        metrics_df['chainA_plddt_avg'] = [chainA_plddt_avg]
        metrics_df['chainB_plddt_avg'] = [chainB_plddt_avg]

        metrics_df['chainA_plddt_qual'] = [chainA_plddt_qual]
        metrics_df['chainB_plddt_qual'] = [chainB_plddt_qual]

        metrics_df['chainA_numres_plddtover90'] = [chainA_numres_over90]
        metrics_df['chainB_numres_plddtover90'] = [chainB_numres_over90]

        metrics_df['chainA_percentres_plddtover90'] = [chainA_percentres_over90]
        metrics_df['chainB_percentres_plddtover90'] = [chainB_percentres_over90]

        # save_path = f"{save_path}/{complex_name}_metrics_plddt{threshold}_CORRECTPAE.txt"
        if aim1:
            save_path = f"{save_path}/{complex_name}_metrics_plddt{threshold}_CORRECTPAE.txt"
        if aim2:
            save_path = f"{save_path}/{complex_name}_metrics_plddt{threshold}_CORRECTPAE_aim2.txt"
        # print('metrics save path', save_path)
        metrics_df.to_csv(save_path, sep='\t', index=False)


class FindPLDDTdomain(ParseAFoutputs):
    def __init__(self, internal_plddt_threshold):
        super().__init__(json_file, pkl_file, complex_name, complex_seq_len_dict)

        self.internal_plddt_threshold = internal_plddt_threshold

    def find_plddt_domains(self, pae_per_model, ranked_list,
                           comp_dict_clean):

        top_ranked_model_name = next(iter(ranked_list))
        plddt_domains_chainA = []
        plddt_domains_chainB = []
        chain_A = comp_dict_clean[0]
        chain_B = comp_dict_clean[0] + comp_dict_clean[1]

        for model_name, value in pae_per_model.items():
            if model_name == top_ranked_model_name + ".pkl":
                plddt_dict = value['plddt'].copy()

                chainA_plddt_sub = plddt_dict[:chain_A]
                chainB_plddt_sub = plddt_dict[chain_A + 1:chain_B]

                chainA_index = np.where(chainA_plddt_sub >= self.internal_plddt_threshold)
                chainB_index = np.where(chainB_plddt_sub >= self.internal_plddt_threshold)

                temp_chainA_domains = []
                temp_chainB_domains = []

                for cur_res in chainA_index[0]:
                    if temp_chainA_domains and cur_res != temp_chainA_domains[-1] + 1:
                        plddt_domains_chainA.append(temp_chainA_domains)
                        temp_chainA_domains = []
                    temp_chainA_domains.append(cur_res)
                if temp_chainA_domains:
                    plddt_domains_chainA.append(temp_chainA_domains)

                for cur_res in chainB_index[0]:
                    if temp_chainB_domains and cur_res != temp_chainB_domains[-1] + 1:
                        plddt_domains_chainB.append(temp_chainB_domains)
                        temp_chainB_domains = []
                    temp_chainB_domains.append(cur_res)
                    # print("temp_chainB_domains", temp_chainB_domains)
                if temp_chainB_domains:
                    plddt_domains_chainB.append(temp_chainB_domains)

        temp_plddt_chainA_dom = split_sublists(plddt_domains_chainA, 100)
        temp_plddt_chainB_dom = split_sublists(plddt_domains_chainB, 100)

        return temp_plddt_chainA_dom, temp_plddt_chainB_dom

    def domains_avgPAE(self, chainA_plddtdom,
                       chainB_plddtdom,
                       comp_dict_clean,
                       pae_plddt_per_model,
                       top_ranked_model_name,
                       debug=False,
                       ):

        chain_A = comp_dict_clean[0]
        chain_B = comp_dict_clean[0] + comp_dict_clean[1]

        chainA_domain_avg = []
        chainB_domain_avg = []
        chainABtop_domain_avg = []

        for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
            if model_name == top_ranked_model_name + ".pkl":
                current_pae = value['pae']

                chainA_only_pae = current_pae[0:chain_A, 0:chain_A]
                chainB_only_pae = current_pae[chain_A + 1:, chain_A + 1:]
                chain_AB_topdiagonal = current_pae[0:chain_A:, chain_A + 1:]
                chainAB_bottomdiagonal = current_pae[chain_A + 1:, 0:chain_A]

                ## Debug
                # plt.imshow(chain_AB_topdiagonal, label=model_name, cmap="bwr", vmin=0, vmax=30)
                # plt.title('chain_AB_topdiagonal ')
                # plt.figure()
                # plt.imshow(chainA_only_pae, label=model_name, cmap="bwr", vmin=0, vmax=30)
                # plt.title('chainA_only_pae ')
                # plt.show()

                for i in range(len(chainA_plddtdom)):
                    current_res_list = chainA_plddtdom[i]
                    temp_subarray = []
                    temp_subarrayAB = []
                    for res in current_res_list:
                        sub_array = chainA_only_pae[res]

                        sub_arrayAB = chain_AB_topdiagonal[res]
                        # print("res",res, "subarray",sub_array)
                        temp_subarray.append(sub_array)

                        temp_subarrayAB.append(sub_arrayAB)

                    chainA_domain_avg.append(np.nanmean(temp_subarray))
                    chainABtop_domain_avg.append(np.nanmean(temp_subarrayAB))

                for i in range(len(chainB_plddtdom)):
                    current_res_list = chainB_plddtdom[i]
                    temp_subarray = []
                    for res in current_res_list:
                        abdiagbottom_array = chainAB_bottomdiagonal[res]

                        temp_subarray.append(np.nanmean(abdiagbottom_array))
                    chainB_domain_avg.append(np.nanmean(temp_subarray))

        print('chainA_domain_avg')
        print(chainA_domain_avg)

        print('chainB_domain_avg')
        print(chainB_domain_avg)

        ## "Debug plotting to check correct parsing"
        if debug:
            for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
                if model_name == top_ranked_model_name + ".pkl":
                    current_pae = value['pae']
                    chainA_whole_pae = current_pae[0:chain_A]
                    chainB_whole_pae = current_pae[chain_A + 1:chain_B]

                    chainA_only = current_pae[0:chain_A, 0:chain_A]
                    chainB_only = current_pae[chain_A + 1:, chain_A + 1:]

                    chainAB_topdiagonal = current_pae[0:chain_A:, chain_A + 1:]
                    chainAB_bottomdiagonal = current_pae[chain_A + 1:, 0:chain_A]

                    plt.imshow(chainA_whole_pae, label=model_name, cmap="bwr", vmin=0, vmax=30)
                    plt.title('chainA_whole_pae only')

                    plt.figure()

                    plt.imshow(chainA_only, label=model_name, cmap="bwr", vmin=0, vmax=30)
                    plt.title('chainA_only only')

                    plt.figure()

                    plt.imshow(chainB_whole_pae, label=model_name, cmap="bwr", vmin=0, vmax=30)
                    plt.title('chainB_whole_pae only')

                    plt.figure()

                    plt.imshow(chainB_only, label=model_name, cmap="bwr", vmin=0, vmax=30)
                    plt.title('chainB only')

                    plt.imshow(chainAB_topdiagonal, label=model_name, cmap="bwr", vmin=0, vmax=30)
                    plt.title('chainAB_topdiagonal')

                    plt.figure()
                    plt.imshow(chainAB_bottomdiagonal, label=model_name, cmap="bwr", vmin=0, vmax=30)
                    plt.title('chainAB_bottomdiagonal')

                    plt.show()

        return chainA_domain_avg, chainB_domain_avg, chainABtop_domain_avg

    def selected_res(self, chainA_plddtdom, chainB_plddtdom, chainAdom_pae, chainBdom_pae, comp_dict_clean,
                     chainABtop_domain_avg, pae_cutoff):

        chain_Ares = comp_dict_clean[0]
        chain_Bres = comp_dict_clean[1]

        chainA_selindex = []
        chainA_selres = []
        for i in range(len(chainAdom_pae)):
            if chainAdom_pae[i] <= pae_cutoff:
                chainA_selindex.append(i)

        for i in range(len(chainA_plddtdom)):
            for j in chainA_selindex:
                if i == j:
                    chainA_selres.append(chainA_plddtdom[i])

        # chainABtop_selindex = []
        # chainABtop_selres = []
        # for i in range(len(chainABtop_domain_avg)):
        #     if chainABtop_domain_avg[i] <= 21:
        #         chainABtop_selindex.append(i)
        # for i in range(len(chainA_plddtdom)):
        #     for j in chainABtop_selindex:
        #         if i == j:
        #             chainABtop_selres.append(chainA_plddtdom[i])
        # # print("chainABtop_selindex")
        # print(chainABtop_selindex)

        chainB_selindex = []
        chainB_selres = []

        for i in range(len(chainBdom_pae)):
            # print(i, chainBdom_pae[i])
            # print(chainBdom_pae[i])
            if chainBdom_pae[i] <= pae_cutoff:
                chainB_selindex.append(i)

        for i in range(len(chainB_plddtdom)):
            for j in chainB_selindex:
                if i == j:
                    chainB_selres.append(chainB_plddtdom[i])

        checkA = check_nested_level(chainA_selres)
        checkB = check_nested_level(chainB_selres)

        if not chainA_selindex or not chainB_selindex:
            print(f"pae scores too high for selection, all pae >= {pae_cutoff}")

        else:
            chainA_filres = []
            if checkA:
                for i in range(1, len(chainA_selres)):
                    current_sublist = chainA_selres[i]
                    previous_sublist = chainA_selres[i - 1]
                    res_difference = current_sublist[0] - previous_sublist[-1]
                    if res_difference <= 5 or current_sublist[
                        0] <= chain_Ares:  # or current res is in middle of sequence (<=last res# of chain)
                        # if res_difference <= 15:  # or current res is in middle of sequence (<=last res# of chain)
                        missing_values = list(range(previous_sublist[-1] + 1, current_sublist[0]))
                        combined_sublist = previous_sublist + missing_values + current_sublist
                        chainA_filres.append(combined_sublist)
                    else:
                        # print("else chainA_selres", chainA_selres)
                        chainA_filres = chainA_selres
            else:
                chainA_filres = chainA_selres  # if list only contains 1 sublist

            chainB_filres = []
            if checkB:
                for i in range(1, len(chainB_selres)):
                    current_sublist = chainB_selres[i]
                    previous_sublist = chainB_selres[i - 1]
                    res_difference = current_sublist[0] - previous_sublist[-1]
                    if res_difference <= 5 or current_sublist[
                        0] <= chain_Bres:  # or current res is in middle of sequence (<=last res# of chain)
                        # if res_difference <= 15:  # or current res is in middle of sequence (<=last res# of chain)
                        missing_values = list(range(previous_sublist[-1] + 1, current_sublist[0]))
                        combined_sublist = previous_sublist + missing_values + current_sublist
                        chainB_filres.append(combined_sublist)
                    else:
                        # print("else chainB_selres", chainB_selres)
                        chainB_filres = chainB_selres

            else:
                chainB_filres = chainB_selres

            chainA_reslist_temp = [item for tpl in chainA_filres for item in tpl]
            chainB_reslist_temp = [item for tpl in chainB_filres for item in tpl]

            newpdb_chainA_reslist = chainA_reslist_temp
            newpdb_chainB_reslist = chainB_reslist_temp

            if len(newpdb_chainB_reslist) <= 10:
                print(
                    f'{complex_name}: filtered sequence for chain B too short length:{len(newpdb_chainB_reslist)}, not meaningful, removing example')
            else:
                return newpdb_chainA_reslist, newpdb_chainB_reslist


class Plotting(ParseAFoutputs):

    def __init__(self, ):
        super().__init__(json_file, pkl_file, complex_name, complex_seq_len_dict)

    def plot_pae_per_model(self,
                           pae_per_model,
                           comp_dict_clean,
                           ranked_list,
                           out_path,
                           num_models=5,
                           save_plots=False, plot_top=False):
        """
        Generate PAE plot for all models
        Args:
            pae_per_model: The dictionary of collected PAE/PLDDT scores
            complex_name: The name of current complex
            complex_seq_len_dict: The dictionary of complex and its sequence lengths
            ranked_list: list of models in ranked order
            out_path: path to save plots
            num_models: number of models
            save_plots: If True will save PAE plot in 'out_path' directory
            plot_top: if True, will only plot PAE for top ranked model

        Returns: None

        """
        top_ranked_model_name = ranked_list

        if len(pae_per_model) < 25:  # less than 25 means incomplete runs for default multimer 5x5models
            return

        num_runs_per_model = math.ceil(len(pae_per_model) / num_models)
        fig = plt.figure(figsize=(3 * num_models, 2 * num_runs_per_model), dpi=100)
        fig.suptitle(f"Predicted Alignment Error (PAE) - {self.complex_name}")

        for n, (model_name, value) in enumerate(pae_per_model.items()):
            if plot_top and model_name == top_ranked_model_name + ".pkl":
                print("VALUE PAE!!!!", type(value["pae"]))
                plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=30)
                plt.title(model_name[:-4] + 'pae')
                plt.colorbar()
            elif not plot_top:
                # print(value["pae"])

                plt.subplot(num_runs_per_model, num_models, n + 1)
                subplot_label = model_name[:-4]
                plt.title(subplot_label)
                plt.imshow(value["pae"], label=subplot_label, cmap="bwr", vmin=0, vmax=30)
                plt.colorbar()
            plt.yticks([0, comp_dict_clean[0], sum(comp_dict_clean)], ['0', 'A', 'B'], color='black', fontsize=14)

            # print(comp_dict_clean[0])

            plt.axhline(y=comp_dict_clean[0], color='black', linestyle='--', linewidth=2)
            plt.axvline(x=comp_dict_clean[0], color='black', linestyle='--', linewidth=2)

        fig.tight_layout()

        if save_plots:
            plot_type = "top_" if plot_top else ""
            pae_save_path = f"{out_path}/{self.complex_name}_{plot_type}pae.png"
            # print(f"Saving {'TOP ' if plot_top else ''} pae plot", pae_save_path)
            plt.savefig(pae_save_path)
            plt.close()
        else:
            plt.show()

    def plot_plddt_per_position(self, pae_per_model,
                                ranked_list,
                                comp_dict_clean,
                                out_path,
                                save_plots=False,
                                plot_top=False):
        """
        Plots plddt score
        Args:
            pae_per_model: The dictionary of collected PAE/PLDDT scores
            complex_name: name of current complex
            ranked_list: ranked list of all models, in order 1-5
            comp_dict_clean: list of chain A length and chain B end point
            out_path: save directory path
            save_plots: if true, saves plot to out_path
            plot_top: if true, plots only the top ranked model

        Returns:None

        """
        fig = plt.figure(figsize=(30, 15))
        ax1 = fig.add_subplot(111)

        plt.title(f"Predicted pLDDT per position - {self.complex_name}")
        yticks = np.arange(0, 101, 10)
        plt.yticks(yticks, color='black', fontsize=14)
        plt.ylabel("Predicted LDDT", color='black', fontsize=14)
        plt.axhline(y=90, color='black', linestyle='--', label='Threshold at 90')
        # plt.axhline(y=75, color='black', linestyle='--', label='domains threshold at 75')

        top_ranked_model_name = next(iter(ranked_list))

        for model_name, value in pae_per_model.items():
            if plot_top and model_name == top_ranked_model_name + ".pkl":
                name = top_ranked_model_name
                legend_name = model_name[0:7] + "_" + model_name[20:26]  # short legend name
                ax1.plot(value["plddt"], label=f"{legend_name}, plddts: {round(ranked_list[name], 6)}")
            elif not plot_top:
                name = model_name[:-4]
                legend_name = model_name[0:7] + "_" + model_name[20:26]  # short legend name
                ax1.plot(value["plddt"], label=f"{legend_name}, plddts: {round(ranked_list[name], 6)}")

            clean_list = [0, comp_dict_clean[0], len(value["plddt"])]
        ax2 = ax1.twiny()
        ax1.set_xlabel(r"Residues")
        fig.subplots_adjust(bottom=0.2)
        ax2.set_xticks(clean_list, color='black', fontsize=14)
        ax2.set_xticklabels(['Label A', 'Label A/B', 'Label B'], color='black', fontsize=14)

        # Move twinned axis ticks and label from top to bottom
        ax2.xaxis.set_ticks_position("bottom")
        ax2.xaxis.set_label_position("bottom")

        # Offset the twin axis below the host
        ax2.spines["bottom"].set_position(("axes", -0.15))
        ax2.set_frame_on(True)
        ax2.patch.set_visible(False)
        # ax2.spines['top'].set_position(('outward', 40))
        for sp in ax2.spines.values():
            sp.set_visible(False)
        ax2.spines["bottom"].set_visible(True)
        ax2.set_xlim(ax1.get_xlim())

        ax2.set_xticks(clean_list, color='black', fontsize=14)
        ax2.set_xticklabels(['Label A', 'Label A/B', 'Label B'], color='black', fontsize=14)
        # ax2.set_yticklabels(['Label A', 'Label B'])

        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True)
        # ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True)

        if save_plots:
            plot_type = "top_" if plot_top else ""
            plddt_save_path = f"{out_path}/{self.complex_name}_{plot_type}plddt.png"
            # print(f"Saving {'TOP ' if plot_top else ''}predicted_lddt plot", plddt_save_path)
            plt.savefig(plddt_save_path)
            plt.close()
        else:
            plt.show()

    def plot_iptms(self, ranked_list):
        """
        Plot a bar plot of iptms score for all models
        Args:
            ranked_list: list of ranked models + associated iptms
        Returns: none

        """
        # print(ranked_list)
        model_names = list(ranked_list.keys())
        scores = list(ranked_list.values())
        plt.title('Ranked models iPTM score')
        plt.bar(model_names, scores, color='blue')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()


def check_nested_level(my_list):
    """
    Checks if a nested list contains more than 1 sublist
    Args:
        my_list: list
    Returns:True if more than 1 nested sublist, False if only 1 sublist
    """
    nested_list_count = sum(isinstance(item, list) for item in my_list)
    return nested_list_count > 1


def path_split(cluster):
    if cluster:
        return str(sub).split("/")[-1]
    else:
        return str(sub).split("\\")[-1]


def split_sublists(lst, list_length):
    result = []
    for sublist in lst:
        if len(sublist) > list_length:
            for i in range(0, len(sublist), list_length):
                result.append(sublist[i:i + list_length])
        else:
            result.append(sublist)
    return result


if __name__ == "__main__":
    U = UtilityFunctions()

    if len(sys.argv) > 1:
        cluster = True
    else:
        cluster = False

    if cluster:
        aim1 = True
        aim2 = False

        cluster_path = sys.argv[2]
        clean_completed_files = sys.argv[3]

        print("clean_completed_files", clean_completed_files)
        print(f"Current dataset::{aim1}")
        print("cluster_path", cluster_path)

        p = Path(cluster_path) if aim1 else Path(
            '/projects/ccib/lamoureux/as3190/alphafold_multimer/aim2_output_data/')
        complex_seq_len_path = '/projects/ccib/lamoureux/as3190/alphafold_multimer/complex_seq_lengths_dict.pkl'
        complex_seq_len_dict = U.read_pkl(complex_seq_len_path)
        completed_list_path = (p / clean_completed_files) if aim1 else (p / 'aim2_clean_completed_directory.txt')
        aim1_workindataset = p / "trim_complex_aim1_561examples_1743total.txt"  # total run exp 1723 (NOT 17) (5/22/2024)

    else:
        p = Path('/Users/Anushriya/Documents/Lamoureux_lab/alphafold_multimer/output_data/noninteraction_pairs_AF/')
        complex_seq_len_path = "/Users/Anushriya/Documents/Lamoureux_lab/alphafold_multimer/output_data/PTM_tests/complex_seq_lengths_dict.pkl"
        complex_seq_len_dict = U.read_pkl(complex_seq_len_path)
        completed_list_path = p / "clean_completed_directory.txt"  # all completed files
        aim1_workindataset_1027exp = p / "workingset_complex_aim1_469exp.txt"  # only specific files as needed

    ###################################################################################################################
    "set the PLDDT cutoff threshold"
    plttd_threshold = 50  # initial plddt confidence threshold

    negative_interaction = True  # if True= complexes are negative non-interaction pairs, if False: Positive interaction pairs
    if negative_interaction:
        print("Processing NEGATIVE NON INTERACTION COMPLEX")
    else:
        print("Processing POSITIVE INTERACITON COMPLEX")

    plotting = True
    save_plots = True
    # plot_top if True, only plot for top model
    plot_top = True

    "True = calculates averages for all metrics"
    calculate_metrics = True

    "Flags for saving new trimmed pdb files, True=Write new pdb file"
    save_trimmed_pdbs = False

    "Workingset = Structures that meet all metrics (based on metrics_analysis.py), True=text file with sub-directory names of complexes"
    workingset_filter = False

    aim1 = True
    aim2 = False

    if workingset_filter:
        print("CURRENTLY PARSING WORKING DATASET (EXAMPLES MEETING ALL THRESHOLDS)")
        "Parse through all subdirectories within working set dataset"
        completed_list = []
        with open(aim1_workindataset, "r") as file:
            for line in file:
                line = line.strip()
                completed_list.append(line)
        print("Number of examples in working set", len(completed_list))

    else:
        "Parse through all subdirectories within output_data"
        completed_list = []
        with open(completed_list_path, "r") as file:
            for line in file:
                line = line.strip()
                completed_list.append(line)
        print("Number of newly completed complex", len(completed_list))

    ###################################################################################################################
    current_flags = f'\nCluster= {cluster}' \
                    f'\ncomplex_seq_len_dict= {len(complex_seq_len_dict)}' \
                    f'\nplttd_threshold= {plttd_threshold}' \
                    f'\nplotting= {plotting}' \
                    f'\nsave_plots= {save_plots}' \
                    f'\nplot_top= {plot_top}' \
                    f'\nsave_trimmed_pdbs= {save_trimmed_pdbs}' \
                    f'\nAIM1= {aim1}' \
                    f'\nAIM2= {aim2}'

    print("###############################")
    print('current_flags:::', current_flags)
    print("###############################")
    ###################################################################################################################

    if os.path.exists(p):
        sub_dirs = [sub_dir for sub_dir in p.iterdir() if sub_dir.is_dir()]  # List all subdirectories
    else:
        print("Path not found")

    for sub in sub_dirs:
        out_path = sub
        # print("out_path", out_path)
        complex_name = path_split(cluster)  # path names/splitting differs between windows/linux systems
        if complex_name in completed_list:
            print("Current complex", complex_name)
            json_file = sub
            pkl_file = sub

            parseaf_outputs = ParseAFoutputs(json_file, pkl_file,
                                             complex_name, complex_seq_len_dict)

            ranked_model = parseaf_outputs.get_top_rankings()

            if ranked_model[0]:
                # Dictionary of plddt/pae/ptm/iptm metrics
                pae_plddt_per_model = parseaf_outputs.all_models_pae()

                comp_dict_clean = parseaf_outputs.clean_comp_lenght_dict(negative_interaction=negative_interaction)

                res_over_threshold = parseaf_outputs.get_residues(pae_plddt_per_model,
                                                                  ranked_model[1],
                                                                  threshold=plttd_threshold,
                                                                  )

                if calculate_metrics:
                    print("Getting PAE and PLDDT averages")
                    get_avgs = GetAvgs(res_over_threshold)
                    # """Get averages of residues that have >= confidence_threshold %plddt score,
                    averages_pae_plddt = get_avgs.get_avg_pae_plddt(pae_plddt_per_model,
                                                                    comp_dict_clean,
                                                                    ranked_model[-1],
                                                                    )
                    get_avgs.all_metrics(averages_pae_plddt[0],
                                         averages_pae_plddt[1],
                                         averages_pae_plddt[2],
                                         averages_pae_plddt[3],
                                         pae_plddt_per_model,
                                         ranked_model[-1],
                                         out_path,
                                         threshold=plttd_threshold,
                                         aim1=aim1,
                                         aim2=aim2
                                         )

                if plotting:
                    plot = Plotting()

                    plot.plot_pae_per_model(pae_plddt_per_model,
                                            comp_dict_clean,
                                            ranked_model[-1],
                                            out_path,
                                            num_models=5,
                                            save_plots=save_plots, plot_top=plot_top
                                            )

                    plot.plot_plddt_per_position(pae_plddt_per_model,
                                                 ranked_model[1],
                                                 comp_dict_clean,
                                                 out_path,
                                                 save_plots=save_plots,
                                                 plot_top=plot_top)

                if save_trimmed_pdbs:
                    """Finding plddt domains:"""

                    # Stricter plddt threshold for trimming PDB files
                    internal_plddt_threshold = 60
                    # pae_cutoff = 17
                    pae_cutoff = 18
                    save_png = True

                    print('trimming internal plddt threshold::', internal_plddt_threshold)
                    find_subplddt_dom = FindPLDDTdomain(internal_plddt_threshold)

                    plddt_domains = find_subplddt_dom.find_plddt_domains(pae_plddt_per_model,
                                                                         ranked_model[1],
                                                                         comp_dict_clean,
                                                                         )

                    domain_paes = find_subplddt_dom.domains_avgPAE(plddt_domains[0],
                                                                   plddt_domains[1],
                                                                   comp_dict_clean,
                                                                   pae_plddt_per_model,
                                                                   ranked_model[-1],
                                                                   debug=False
                                                                   )
                    filter_pdb = find_subplddt_dom.selected_res(plddt_domains[0],
                                                                plddt_domains[1],
                                                                domain_paes[0],
                                                                domain_paes[1],
                                                                comp_dict_clean,
                                                                domain_paes[2],
                                                                pae_cutoff=pae_cutoff)
                    if filter_pdb is None:
                        print(f'complex: {complex_name}, did not satisfy pae score for trimming')

                    else:
                        chainA = filter_pdb[0]
                        chainB = filter_pdb[1]

                        cmd.extend('write_filter_pdb', write_filter_pdb)
                        write_filter_pdb(out_path, complex_name, chainA, chainB, save_png=save_png)
