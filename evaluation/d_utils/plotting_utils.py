#!/usr/bin/python3

import os
import math
import numpy as np
import pandas as pd
from .. import config as c
import seaborn as sns
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from collections import defaultdict
from ..mappers.mapper import Mapper, FileMapper
# only in full biodigest
import graph_tool as gt
import graph_tool.util as gtu
import graph_tool.topology as gtt
from graph_tool import GraphView, draw, Graph

sns.set_palette("colorblind")

plt.rcParams.update({'font.size': 17, 'axes.titlelocation': "left", 'axes.titleweight': "bold", 'axes.labelsize': 21})
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

replacements = {"related_genes": "related\ngenes", "related_variants": "related\nvariants"}
eval_terms = {"DI-based": {"single": "Dunn index", "multi": "Dunn indices"},
              "SS-based": {"single": "Sillhouette score", "multi": "Sillhouette scores"},
              "DBI-based": {"single": "Davis-Bouldin index", "multi": "Davis-Bouldin indices"},
              "JI-based": {"single": "Jaccard index", "multi": "Jaccard indices"},
              "OC-based": {"single": "Overlap coefficient", "multi": "Overlap coefficients"}}
annot_terms = {'GO.BP': 'GO.BP-based', 'GO.CC':'GO.CC-based', 'GO.MF':'GO.MF-based',
               'KEGG': 'KEGG-based', "related_genes": "associated genes", "related_variants": "associated variants"}

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
    if mode == "clustering":
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
        p_value_plot(title="Empirical P-value\nbased on " +eval_terms[val]["single"], p_value_df=p_value_df, out_dir=out_dir,
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
    mappability_plot(title="Mappability of input\ninput to annotations", in_type=in_type, mapped_df=mapped_df, out_dir=out_dir,
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
    for val in results["p_values"]["values"]:
        p_value_plot(title="Empirical P-value\nbased on " +eval_terms[val]["single"], p_value_df=p_value_df,
                     out_dir=out_dir, prefix=prefix + "_" + val, file_type=file_type)
    # ===== Prepare for mappability plot =====
    mapped_df = pd.DataFrame()
    for att in results["input_values"]["mapped_ids"]:
        mapped_df[att] = [1 if len(results["input_values"]["mapped_ids"][att][x]) > 0 else 0 for x in user_input["set"]]
    mapped_df = mapped_df.T
    mapped_df["count"] = mapped_df.sum(axis=1)
    mapped_df["fraction"] = mapped_df['count'].apply(lambda x: x / len(user_input["set"]))
    mapped_df = mapped_df.rename_axis('attribute').reset_index()
    mapped_df = mapped_df.replace(replacements).sort_values(['attribute']).reset_index(drop=True)
    # ===== Plot mappability plot =====
    mappability_plot(title="Mappability of\ninput to annotations", in_type=in_type, mapped_df=mapped_df, out_dir=out_dir,
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


def create_extended_plots(results, mode, tar, out_dir, prefix, file_type: str = "pdf", mapper:Mapper = FileMapper()):
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
    value_distribution_plots(results=results, out_dir=out_dir, prefix=prefix, file_type=file_type)
    term_annotation_plots(results=results, out_dir=out_dir, prefix=prefix, file_type=file_type)
    sankey_plot(results=results, mode=mode, mapper=mapper, out_dir=out_dir, prefix=prefix,
                file_type=file_type, tar_cluster=tar)


def value_distribution_plots(results, out_dir, prefix, file_type: str = "pdf"):
    for eval_term in results["input_values"]['values']:
        df = pd.melt(pd.DataFrame(results['random_values']['values'][eval_term]))
        df['value'] = df['value'].astype(float)
        for term_index, term in enumerate(results["input_values"]['values'][eval_term]):
            fig = plt.figure(figsize=(7, 6), dpi=80)
            plt.axvline(float(results["input_values"]['values'][eval_term][term]), color='darkred', lw=10)
            ax = sns.histplot(df[df["variable"] == term], x="value", kde=True, color=sns.color_palette()[term_index],
                              bins=10)
            if term in replacements:
                plt.title("Distribution of " + eval_terms[eval_term]["multi"] + "\nbased on " + annot_terms[term])
                plt.xlabel(eval_terms[eval_term]["single"] + " based\non " + annot_terms[term])
            else:
                plt.title("Distribution of\n" + annot_terms[term] + " " + eval_terms[eval_term]["multi"])
                plt.xlabel(annot_terms[term] + " " + eval_terms[eval_term]["single"])
            plt.ylabel("Number of runs\non randomized data")
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.05, 0.95, "Empirical\nP-value:\n%.3f" % results["p_values"]['values'][eval_term][term],
                    transform=ax.transAxes, fontsize=14,
                    verticalalignment='top', bbox=props)
            fig.tight_layout()
            fig.savefig(os.path.join(out_dir, prefix + '_' + eval_term + '_' + term + '_distribution.' + file_type),
                        bbox_inches='tight')


def term_annotation_plots(results, out_dir, prefix, file_type: str = "pdf"):
    df = pd.DataFrame(results["input_values"]["mapped_ids"]).fillna("").applymap(len)
    for term_index, term in enumerate(df.columns):
        fig = plt.figure(figsize=(7, 6), dpi=80)
        sns.histplot(df, x=term, kde=True, color=sns.color_palette()[term_index], bins=10)
        if term in replacements:
            plt.title("Distribution of numbers\nof " + annot_terms[term])
            plt.xlabel("Number of " + annot_terms[term])
        else:
            plt.title("Distribution of numbers\nof associated " + term + "-terms")
            plt.xlabel("Number of associated\n" + term + "-terms")
        plt.ylabel("Number of IDs\ncontained in query")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, prefix + '_' + term + '_annotation_distribution.' + file_type),
                    bbox_inches='tight')


def sankey_plot(results, mode, out_dir, prefix, file_type: str = "pdf", tar_cluster=None, include_others=False,
                mapper:Mapper = FileMapper()):
    full_df = pd.DataFrame(results["input_values"]["mapped_ids"])
    for term_index, term in enumerate(full_df.columns):
        if len(full_df[full_df[term].str.len() != 0]) == 0: # save empty plot
            fig = plt.figure(dpi=80)
            plt.show()
            fig.savefig(os.path.join(out_dir, prefix + '_' + term + '_sankey.' + file_type),
                        bbox_inches='tight')
        else:
            def col_sort(df, colname, hierarchy):
                df[colname] = df[colname].astype("category")
                df[colname].cat.set_categories(hierarchy, inplace=True)
                df_sorted = df.sort_values([colname], ascending=True)
                return df_sorted

            d = full_df[[term]].dropna().rename_axis('left').reset_index().explode(term)
            if mode == "clustering":
                d = d.replace({"left": tar_cluster.set_index('id')["cluster"].to_dict()})
                df = pd.DataFrame(columns=["id", term])
                for cluster in d['left'].unique():
                    df = pd.concat([df, d[d["left"] == cluster][term].value_counts(normalize=True).rename_axis(
                        'id').reset_index(name=term)])
                ids = df.groupby(['id']).sum().rename_axis('id').reset_index().nlargest(10, term)["id"]
            else:
                ids = d[term].value_counts().nlargest(10).index
            if include_others:
                d.loc[~d[term].isin(ids), term] = "other"
            else:
                d = d[d[term].isin(ids)]
            d[term] = d[term].astype({term: str}, errors='raise')
            if term == "related_genes":
                rename_ids, _ = mapper.get_loaded_mapping(in_set=d[term], id_type="entrezgene", key="gene_ids")
                rename_ids = rename_ids.explode("symbol").set_index('entrezgene')['symbol'].to_dict()
                d[term] = d[term].map(rename_ids).fillna(d[term])
            # ===== Save hierarchy =====
            hierarchy_left = d["left"].value_counts().index.tolist()
            hierarchy_right = d[term].value_counts().index.tolist()
            # ===== Sort =====
            d = d.reset_index(drop=True).value_counts().reset_index()
            d = col_sort(df=d, colname="left", hierarchy=hierarchy_left)
            d = col_sort(df=d, colname=term, hierarchy=hierarchy_right)
            # ===== Assign colors =====
            z = {**dict(zip(d["left"].unique(), ["#808080"] * len(d.index.unique()))),
                 **dict(zip(d[term].unique(), sns.color_palette() * 10))}

            d.rename(columns={term: 'right', 0: 'weight'}, inplace=True)
            # ===== Plot =====
            sankey(data=d, aspect=20, right_color=True, color_dict=z, term=term,
                   out_dir=out_dir, prefix=prefix, file_type=file_type)


def sankey(data, out_dir, prefix, file_type: str = "pdf", color_dict=None, aspect=4, right_color=False, term=None):
    """
    Make Sankey Diagram showing flow from left-->right

    :param color_dict: dictionary of colors to use for each label {'label':'color'}
    :param aspect: vertical extent of the diagram in units of horizontal extent
    :param right_color: if true, each strip in the diagram will be be colored according to its left label

    original: https://github.com/anazalea/pySankey
    """

    # ===== Identify all labels that appear 'left' or 'right' =====
    all_labels = pd.Series(np.r_[data.left.unique(), data.right.unique()]).unique()

    # ===== Identify left and right labels =====
    left_labels = data.left.unique()
    right_labels = data.right.unique()

    # ===== If no color_dict given, make one =====
    if color_dict is None:
        color_dict = {}
        palette = "hls"
        color_palette = sns.color_palette(palette, len(all_labels))
        for i, label in enumerate(all_labels):
            color_dict[label] = color_palette[i]

    # ===== Determine widths of individual strips =====
    ns_l = defaultdict()
    ns_r = defaultdict()
    for left_label in left_labels:
        left_dict = {}
        right_dict = {}
        for right_label in right_labels:
            left_dict[right_label] = data[(data.left == left_label) & (data.right == right_label)].weight.sum()
            right_dict[right_label] = data[(data.left == left_label) & (data.right == right_label)].weight.sum()
        ns_l[left_label] = left_dict
        ns_r[left_label] = right_dict

    # ===== Determine positions of left label patches and total widths =====
    left_widths = defaultdict()
    for i, left_label in enumerate(left_labels):
        myD = {}
        myD['left'] = data[data.left == left_label].weight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['left']
        else:
            myD['bottom'] = left_widths[left_labels[i - 1]]['top'] + 0.02 * data.weight.sum()
            myD['top'] = myD['bottom'] + myD['left']
            top_edge = myD['top']
        left_widths[left_label] = myD

    # ===== Determine positions of right label patches and total widths =====
    right_widths = defaultdict()
    for i, right_label in enumerate(right_labels):
        myD = {}
        myD['right'] = data[data.right == right_label].weight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['right']
        else:
            myD['bottom'] = right_widths[right_labels[i - 1]]['top'] + 0.02 * data.weight.sum()
            myD['top'] = myD['bottom'] + myD['right']
            top_edge = myD['top']
        right_widths[right_label] = myD

    # ===== Total vertical extent of diagram =====
    xMax = top_edge / aspect

    # ===== Draw vertical bars on left and right of each  label's section & print label =====
    fig = plt.figure(dpi=80)
    for left_label in left_labels:
        plt.fill_between(
            [-0.02 * xMax, 0],
            2 * [left_widths[left_label]['bottom']],
            2 * [left_widths[left_label]['bottom'] + left_widths[left_label]['left']],
            color=color_dict[left_label],
            alpha=0.99
        )
        plt.text(
            -0.05 * xMax,
            left_widths[left_label]['bottom'] + 0.5 * left_widths[left_label]['left'],
            left_label,
            {'ha': 'right', 'va': 'center'},
            fontsize=min(left_widths[left_label]['left'] + 9, 14)
        )
    for right_label in right_labels:
        plt.fill_between(
            [xMax, 1.02 * xMax], 2 * [right_widths[right_label]['bottom']],
                                 2 * [right_widths[right_label]['bottom'] + right_widths[right_label]['right']],
            color=color_dict[right_label],
            alpha=0.99
        )
        plt.text(
            1.05 * xMax,
            right_widths[right_label]['bottom'] + 0.5 * right_widths[right_label]['right'],
            right_label,
            {'ha': 'left', 'va': 'center'},
            fontsize=min(right_widths[right_label]['right'] + 9, 14)
        )

    # ===== Plot strips =====
    for left_label in left_labels:
        for right_label in right_labels:
            label_color = left_label
            if right_color:
                label_color = right_label
            if len(data[(data.left == left_label) & (data.right == right_label)]) > 0:
                # Create array of y values for each strip, half at left value,
                # half at right, convolve
                ys_d = np.array(50 * [left_widths[left_label]['bottom']] + 50 * [right_widths[right_label]['bottom']])
                ys_d = np.convolve(ys_d, 0.05 * np.ones(20), mode='valid')
                ys_d = np.convolve(ys_d, 0.05 * np.ones(20), mode='valid')
                ys_u = np.array(50 * [left_widths[left_label]['bottom'] + ns_l[left_label][right_label]] + 50 * [
                    right_widths[right_label]['bottom'] + ns_r[left_label][right_label]])
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')

                # Update bottom edges at each label so next strip starts at the right place
                left_widths[left_label]['bottom'] += ns_l[left_label][right_label]
                right_widths[right_label]['bottom'] += ns_r[left_label][right_label]
                plt.fill_between(
                    np.linspace(0, xMax, len(ys_d)), ys_d, ys_u, alpha=0.65,
                    color=color_dict[label_color]
                )
    plt.gca().axis('off')
    plt.gcf().set_size_inches(6, 6)
    if term in replacements:
        plt.title("Top 10 most frequent\n" + annot_terms[term] + "\nlinked to each ID")
    else:
        plt.title("Top 10 most frequent\n" + term + "-terms\nlinked to each ID")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, prefix + '_' + term + '_sankey.' + file_type),
                bbox_inches='tight')


def contribution_heatmap(df, axis_limit, contribution_type: str, input_type, num,
                         out_dir, prefix, title_ext="",  file_type: str = "pdf"):
    fig = plt.figure(figsize=(7, 6), dpi=80)
    sns.heatmap(data=df, cmap=sns.color_palette("vlag"), yticklabels=1, vmin=-axis_limit, vmax=axis_limit,
                cbar_kws={'label': 'Significance contribution'})
    plt.title("Top " + str(num) + " " + input_type + " with\nlargest " + contribution_type
              + " contributions\nto empirical P-values" + title_ext)
    plt.yticks(rotation=0)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, prefix  + '_contribution_heatmap.' + file_type),
                bbox_inches='tight')


# def create_contribution_graphs(result_sig, tar, network_data, out_dir, prefix, file_type: str = "pdf"):
#     if network_data is None:
#         if input_type == "diseases":
#             network_data = {"network_file": os.path.join(mapper.files_dir, "ddi_graph.graphml"),
#                             "id_type": new_id_type, "prop_name": "id"}
#         else:
#             input_type = {"network_file": os.path.join(mapper.files_dir, "ggi_graph.graphml"),
#                             "id_type": new_id_type, "prop_name": "id"}
#     g = gt.load_graph(network_data["network_file"])
#     nodes = set(tar)
#     vfilt = g.new_vertex_property('bool')
#     for vertex in g.get_vertices():
#         if g.vertex_properties.name[vertex] in nodes:
#             vfilt[vertex] = True
#     sub_g = Graph(GraphView(g, vfilt=vfilt), prune=True)
#
#     for stat_type in result_sig:
#         sig_df = pd.DataFrame(result_sig[stat_type]).T
#         cmap = sns.color_palette("vlag", as_cmap=True)
#         lim = sig_df.abs().max().max()
#         norm = colors.Normalize(vmin=-lim, vmax=lim)
#
#         for col in sig_df.columns:
#             rgba_values = cmap(norm(sig_df[col]))
#             # save colors for heatmap
#             col_dic = dict()
#             for i,g in enumerate(sig_df[col].index):
#                 col_dic[g] = colors.to_hex(rgba_values[i])
#             v_colors = sub_g.new_vertex_property("string")
#             # assign colors to nodes
#             for v in sub_g.vertices():
#                 v_colors[v] = col_dic[sub_g.vp.name[v]]
#
#             plt.switch_backend("cairo")
#             fig = plt.figure(figsize=(12, 8))
#             gs = GridSpec(nrows=1, ncols=2)
#             ax1 = fig.add_subplot(gs[0, :])
#             ax1.axis('off')
#             draw.graph_draw(sub_g, vertex_text = sub_g.vertex_properties['name'], vertex_font_size = 0.4, vertex_size=1.5,
#                         vertex_text_position = -2, vertex_fill_color=v_colors, mplfig=ax1)
#             if col in replacements:
#                 title_extension = "based on " + col + " annotations"
#             else:
#                 title_extension = "based on " + annot_terms[col]
#             ax1.set(title="Subnetwork with nodes colored by\nsignificance contribution w.r.t. P-values\n" + title_extension)
#             ax2 = fig.add_subplot(gs[0, 1])
#             sns.heatmap(data=sig_df, cmap=cmap, yticklabels=1, vmin=-lim, vmax=lim,
#                         cbar=False, ax=ax2).set_visible(False)
#
#             mappable = ax2.get_children()[0]
#             cbar = plt.colorbar(mappable, ax = [ax1,ax2], orientation = 'horizontal', pad=-0.01)
#             cbar.ax.set_xlabel('Significance contribution')
#             fig.savefig(os.path.join(out_dir, prefix + "_" + stat_type+ "_" + col + '_contribution_graph.' + file_type), bbox_inches='tight')


def create_contribution_plots(result_sig, input_type, out_dir, prefix, file_type: str = "pdf"):

    for stat_type in result_sig:
        sig_df = pd.DataFrame(result_sig[stat_type]).T
        sig_df = sig_df.rename(columns=replacements)
        limit = sig_df.abs().max().max()
        sub_df = sig_df.loc[list(sig_df.abs().max(axis=1).sort_values(ascending=False)[0:min(len(sig_df.index), 15)].index)]
        contribution_heatmap(df=sub_df, axis_limit=limit, contribution_type="absolute", input_type=input_type,  num = 15,
                             out_dir=out_dir, prefix=prefix+"_"+stat_type+"_absolute.", file_type=file_type)
        for col in sig_df.columns:
            sub_df = sig_df.loc[list(sig_df[col].sort_values(ascending=False)[0:min(len(sig_df.index), 15)].index)]
            contribution_heatmap(df=sub_df, axis_limit=limit, contribution_type="positive",
                                 input_type=input_type,  num = 15, title_ext=" for "+col,
                                 out_dir=out_dir, prefix=prefix+"_"+stat_type+"_"+col+"_positive.", file_type=file_type)
            sub_df = sig_df.loc[list(sig_df[col].sort_values(ascending=True)[0:min(len(sig_df.index), 15)].index)]
            contribution_heatmap(df=sub_df, axis_limit=limit, contribution_type="negative",
                                 input_type=input_type, num = 15, title_ext=" for "+col,
                                 out_dir=out_dir, prefix=prefix+"_"+stat_type+"_"+col+"_negative.", file_type=file_type)