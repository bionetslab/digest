#!/usr/bin/python3

import os
import math
import pandas as pd
from .. import config as c
import seaborn as sns
from pathlib import Path
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size': 17, 'axes.titlelocation': "left", 'axes.titleweight': "bold", 'axes.labelsize': 21})
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

replacements = {"related_genes": "related\ngenes", "related_variants": "related\nvariants"}


def create_plots(results, mode, tar, tar_id, out_dir, prefix, file_type: str = "pdf"):
    """

    :param results: results generated from single_validation method
    :param mode: comparison mode [set, id-set, set-set, cluster]
    :param tar: path to the file with the target input
    :param tar_id: id type of target input
    :param out_dir: output directory for results
    :param prefix: prefix for the file name
    :param file_type: file ending the plots should have [Default=pdf]
    :return:
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)  # make sure output dir exists
    if mode == "cluster":
        cluster_plot(results=results, user_input={"clustering": tar, "type": tar_id},
                     out_dir=out_dir, prefix=prefix, file_type=file_type)
    else:
        set_plot(results=results, user_input={"set": tar, "type": tar_id},
                 out_dir=out_dir, prefix=prefix, file_type=file_type)


def cluster_plot(results, user_input, out_dir, prefix, file_type: str = "pdf"):
    in_type = "diseases" if user_input["type"] in c.SUPPORTED_DISEASE_IDS else "genes"
    # ===== Prepare for scatterplot =====
    p_value_df = pd.DataFrame.from_dict(results["p_values"]["values"]).rename_axis('attribute').reset_index()
    p_value_df = p_value_df.replace(replacements).sort_values(['attribute']).reset_index(drop=True)
    for val in results["p_values"]["values"]:
        p_value_df["log_p-values"] = p_value_df[val].apply(lambda x: -math.log10(x))
        # ===== Plot scatterplot =====
        p_value_plot(title="Empirical P-value (" + val + ")", p_value_df=p_value_df, out_dir=out_dir,
                     prefix=prefix + "_" + val, file_type=file_type)
    # ===== Prepare for mappability plot =====
    mapped_df = user_input["clustering"][['id', 'cluster']]
    cluster_sizes = mapped_df['cluster'].value_counts().to_dict()
    for att in results["input_values"]["mapped_ids"]:
        mapped_df[att] = [1 if x in results["input_values"]["mapped_ids"][att] else 0 for x in mapped_df['id']]
    mapped_df = mapped_df.groupby('cluster', as_index=False).agg(sum).melt('cluster', var_name='attribute',
                                                                           value_name='count')
    mapped_df = mapped_df.replace(replacements).sort_values(['attribute']).reset_index(drop=True)
    mapped_df["fraction"] = mapped_df.apply(lambda x: x['count'] / cluster_sizes[x['cluster']], axis=1)
    # ===== Plot mappability plot =====
    mappability_plot(title="Mappability of input", in_type=in_type, mapped_df=mapped_df, out_dir=out_dir,
                     prefix=prefix, cluster=True, file_type=file_type)


def set_plot(results, user_input, out_dir, prefix, file_type: str = "pdf"):
    in_type = "diseases" if user_input["type"] in c.SUPPORTED_DISEASE_IDS else "genes"
    # ===== Prepare for scatterplot =====
    p_value_df = pd.DataFrame.from_dict(
        {'p_values': results["p_values"]["values"][next(iter(results["p_values"]["values"]))]}).rename_axis(
        'attribute').reset_index()
    p_value_df["log_p-values"] = p_value_df["p_values"].apply(lambda x: -math.log10(x))
    p_value_df = p_value_df.replace(replacements).sort_values(['attribute']).reset_index(drop=True)
    # ===== Plot scatterplot =====
    p_value_plot(title="Empirical P-value", p_value_df=p_value_df, out_dir=out_dir, prefix=prefix, file_type=file_type)
    # ===== Prepare for mappability plot =====
    mapped_df = pd.DataFrame()
    for att in results["input_values"]["mapped_ids"]:
        mapped_df[att] = [1 if x in results["input_values"]["mapped_ids"][att] else 0 for x in user_input["set"]]
    mapped_df = mapped_df.T
    mapped_df["count"] = mapped_df.sum(axis=1)
    mapped_df["fraction"] = mapped_df['count'].apply(lambda x: x / len(user_input["set"]))
    mapped_df = mapped_df.rename_axis('attribute').reset_index()
    mapped_df = mapped_df.replace(replacements).sort_values(['attribute']).reset_index(drop=True)
    # ===== Plot mappability plot =====
    mappability_plot(title="Mappability of input", in_type=in_type, mapped_df=mapped_df, out_dir=out_dir,
                     prefix=prefix, cluster=False, file_type=file_type)


def p_value_plot(title, p_value_df, out_dir, prefix, file_type: str = "pdf"):
    fig = plt.figure(figsize=(6, 6), dpi=80)
    ax = sns.scatterplot(x=p_value_df['attribute'], y=p_value_df['log_p-values'], s=150)
    ax.set(title=title, ylabel="-log10(P)", xlabel="", ylim=(0, 3.1))
    ax.axhline(y=-math.log10(0.05), color="red", linestyle='--')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, prefix + '_p-value.' + file_type), bbox_inches='tight')


def mappability_plot(title, in_type, mapped_df, out_dir, prefix, cluster=False, file_type: str = "pdf"):
    if cluster:
        fig = plt.figure(figsize=(7, 6), dpi=80)
        ax = sns.barplot(x="attribute", y='fraction', data=mapped_df, hue="cluster")
    else:
        fig = plt.figure(figsize=(6, 6), dpi=80)
        ax = sns.barplot(x="attribute", y='fraction', data=mapped_df)
    ax.set(title=title, xlabel="", ylabel="Fraction of " + in_type + " with\nnon-empty annotation sets", ylim=(0, 1.1))
    if cluster:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="Cluster")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, prefix + '_mappability.' + file_type), bbox_inches='tight')
