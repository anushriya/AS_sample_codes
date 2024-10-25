import warnings
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import os

warnings.filterwarnings('ignore')

desired_width = 320
pd.set_option('display.width', desired_width)
pd.set_option('display.max_columns', 12)


def filter_afmetrics(path, iptm_cutoff=0.50, diagonal_AB_avg_pae_cutoff=15, save_workingset=False):
    """
    Generate a pairwise relationship plot of all AF-Multimer metrics for all complex
    Args:
        diagonal_AB_avg_pae_cutoff: Desired cut off of avg pae of diagonal
        iptm_cutoff: Desired cut off of iPTM
        path: Directory path for concatenated .txt file of all complex metrics
        save_workingset: if True save workingset IDS to txt file
    Returns:None

    """
    iptm_cutoff = iptm_cutoff
    diagonal_AB_avg_pae_cutoff = diagonal_AB_avg_pae_cutoff

    df = pd.read_csv(path, sep='\t', low_memory=True, warn_bad_lines=True)

    df_no_duplicates = df.drop_duplicates().drop(1)  # (check the drop, if needed)

    print(len(df_no_duplicates))
    column2convert = [col for col in df_no_duplicates.columns if col != 'complex_name' and
                      (col != 'chainA_plddt_qual' and col != 'chainB_plddt_qual')]

    df_no_duplicates[column2convert] = df_no_duplicates[column2convert].astype('float64')

    '''First constraint: at least 1 residue over 90%plddt, 0=Valid structures, 1=Not Valid structures'''
    df_no_duplicates['valid'] = np.where(
        (df_no_duplicates['chainA_plddt_qual'] == 'pass') &
        (df_no_duplicates['chainB_plddt_qual'] == 'pass'),
        0, 1  # 0 = valid (structure accepted), 1 = invalid (structure not accepted)
    )

    df_no_duplicates = df_no_duplicates.assign(color=1)  # assign all in this dataframe=1, only first constraint applied

    "Drop un-necessary columns"
    df_no_duplicates = df_no_duplicates.drop("chainB_plddt_qual", axis=1)
    df_no_duplicates = df_no_duplicates.drop("chainA_plddt_qual", axis=1)
    # df_no_duplicates = df_no_duplicates.drop("max_predicted_aligned_error", axis=1)

    df_no_duplicates = df_no_duplicates.drop("chainA_numres_plddtover90", axis=1)
    df_no_duplicates = df_no_duplicates.drop("chainB_numres_plddtover90", axis=1)

    """ Make a copy of dataframe with only accepted structures """
    df_filtered_structures = df_no_duplicates[df_no_duplicates['valid'] == 0].copy(deep=True)

    """ Add further restrictions on iptm and diagonal_AB_avg_pae """

    condition_iptm = df_filtered_structures['iptm'] >= iptm_cutoff
    condition_diagonal_AB_avg_pae = df_filtered_structures['diagonal_AB_avg_pae'] <= diagonal_AB_avg_pae_cutoff
    combined_condition = condition_iptm & condition_diagonal_AB_avg_pae

    # Drop rows based on the combined condition
    df_filtered_structures.drop(df_filtered_structures[~combined_condition].index, inplace=True)
    df_filtered_structures = df_filtered_structures.assign(color=0)  # assign color code for all valid structures

    """Concat original and filtered dataframes"""
    df_concat = pd.concat([df_no_duplicates, df_filtered_structures], axis=0)

    print(f"Number of examples accepted:(constraints: at least 1 resiude>90% plddt & "
          f"iptm>={iptm_cutoff}, & "
          f"diagonal_AB_avg_pae <={diagonal_AB_avg_pae_cutoff}):::",
          len(df_filtered_structures))

    aim1_workingset_complex = df_filtered_structures["complex_name"].to_list()
    print("Number of complexes meeting ALL thresholds::", len(set(aim1_workingset_complex)))

    if save_workingset:
        print("saving working dataset to txt file")
        print("save path", str(path).split("\\")[:-1])
        path_split = str(path).split("\\")[:-1]
        base_path = "\\".join(path_split[:-1])
        print("base_path", base_path)
        output_name = os.path.join(base_path, f'trim_complex_aim1_{len(set(aim1_workingset_complex))}examples.txt')
        with open(output_name, 'w') as file:
            for item in aim1_workingset_complex:
                file.write(f'{item}\n')

    return df_concat, df_filtered_structures


def plotting(df_concat, save_path, show=False):
    plot_subsections = True
    if plot_subsections:
        print(df_concat.columns)
        selected_cols = ['complex_name', 'iptm', 'ptm', 'diagonal_AB_avg_pae', 'valid','color']
        df_selected = df_concat[selected_cols]
        print(df_selected)
        df_plot = df_selected
    else:
        df_plot = df_concat

    sns.set(rc={'figure.figsize': (5, 5)})
    # sns.set(rc={'figure.figsize': (20, 20)})
    column_names = [col for col in df_plot.columns if col != 'complex_name' and col != 'valid' and col != 'color'
                    ]

    sns.set(style="ticks", color_codes=True)
    plot = sns.pairplot(df_plot, vars=column_names, hue='color', diag_kind="kde", dropna=True,
                        plot_kws={'alpha': 0.5})

    new_title = 'Accept/Reject'
    plot._legend.set_title(new_title)
    new_labels = ['VALID', 'NOTVALID']  # 0=valid, 1=invalid
    for t, l in zip(plot._legend.texts, new_labels):
        t.set_text(l)

    for ax in plot.axes.flatten():
        # rotate x axis labels
        ax.set_xlabel(ax.get_xlabel(), rotation=2)  # 0
        # rotate y axis labels
        ax.set_ylabel(ax.get_ylabel(), rotation=90)  # 90
        # set y labels alignment
        ax.yaxis.get_label().set_horizontalalignment('center')

    plt.subplots_adjust(bottom=0.07)

    new_title = 'Accept/Reject'
    plot._legend.set_title(new_title)
    new_labels = ['Valid', 'Not Valid']
    for t, l in zip(plot._legend.texts, new_labels):
        t.set_text(l)

    if show:
        plt.show()
    else:
        output_name = save_path + 'AFmetrics_plots.png'
        plt.savefig(output_name)


def quality_checks(df_validstructures, save_path, show=False):
    new_column_names = ['complex_name', 'iptm', 'ptm', 'ranking_confidence', 'chainA_avg_pae', 'chainB_avg_pae',
                        'diagonal_AB_avg_pae', 'chainA_plddt_avg', 'chainB_plddt_avg', 'chainA_respct_plddtover90',
                        'chainB_respct_plddtover90', 'valid', 'color']

    df_validstructures.columns = new_column_names

    iptm_bins = [.50, .60, .70, .80, .99]
    df_validstructures['iptm_rank'] = pd.cut(df_validstructures['iptm'],
                                             bins=iptm_bins,
                                             labels=[f'very_low', f'low',
                                                     f'med', f'high'],
                                             include_lowest=True)

    ranking_confidencebins = [.50, .60, .70, .80, .99]
    df_validstructures['ranking_confidence_rank'] = pd.cut(df_validstructures['ranking_confidence'],
                                                           bins=ranking_confidencebins,
                                                           labels=[
                                                               f'very_low',
                                                               f'low',
                                                               f'med',
                                                               f'high'],
                                                           include_lowest=True)


    # iptm_bin_counts = df_validstructures['iptm_rank'].value_counts()

    # Define the color palette based on the 'rank' column
    rank_palette = {'very_low': 'orange', 'low': 'blue', 'med': 'green', 'high': 'red'}
    column_names = [col for col in df_validstructures.columns if
                    col != 'complex_name' and col != 'iptm_rank' and col != 'ptm_rank' and col != 'color' and col != 'valid' and col != 'ranking_confidence_rank'
                    ]

    # hue = 'iptm_rank'
    hue = 'ranking_confidence_rank'
    # hue = 'ptm_rank'

    sns.pairplot(df_validstructures, vars=column_names, hue=hue, diag_kind="kde", dropna=True, palette=rank_palette,
                 plot_kws={'alpha': 0.5})
    if show:
        plt.show()
    else:
        output_name = save_path + 'AFmetrics_474ws_RCranked_plots.png'
        plt.savefig(output_name)

    return df_validstructures


if __name__ == '__main__':

    aim1 = True
    aim2 = False

    if aim1:
        p = Path(
            '/Users/Anushriya/Documents/Lamoureux_lab/alphafold_multimer/output_data/all_metrics/metrics_correctPAEs/')
        afmetrics_txt = p / 'combined_CORRECTPAE_metrics_1448exp.txt' ## 1448 completed example

    if aim2:
        p = Path(
            '/Users/Anushriya/Documents/Lamoureux_lab/alphafold_multimer/output_data/aim2_output/')
        afmetrics_txt = p / 'combined_aim2_metrics_110examples.txt'

    plot_show = True
    save_workingset = False
    save_path = p

    # metrics returns 1) full concat_df, 2)df_filtered_structures: structures meeting all constraints:(plddt>=90, iptm_cutoff = 0.50,diagonal_AB_avg_pae_cutoff = 15)

    iptm_cutoff = 0.50
    diagonal_AB_avg_pae_cutoff = 15
    metrics = filter_afmetrics(afmetrics_txt, iptm_cutoff=iptm_cutoff,
                               diagonal_AB_avg_pae_cutoff=diagonal_AB_avg_pae_cutoff,
                               save_workingset=save_workingset)

    # Plot all metrics pairwise relationship
    plotting(metrics[0], save_path, show=plot_show)

    # Plot all metrics, colored by iptm qualities
    # quality_checks(metrics[1], save_path, show=False)
