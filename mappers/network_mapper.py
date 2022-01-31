#!/bin/python3

import os
import time
import graph_tool as gt

from d_utils import runner_utils as ru

start_time = time.time()
iid_files = {"rat": "rat_annotated_PPIs",
             "human": "human_annotated_PPIs",
             "mouse": "mouse_annotated_PPIs",
             "yeast": "yeast_annotated_PPIs",
             "fly": "fly_annotated_PPIs",
             "worm": "worm_annotated_PPIs",
             "chicken": "chicken_annotated_PPIs"}
id_position = {"uniprot": [0, 1],
               "symbol": [2, 3]}


def create_edge_list(organism: str, id_type: str, out_dir: str):
    ru.print_current_usage("Load file from iid database...")
    os.system("mkdir "+out_dir+"tmp")
    os.system("wget 'http://iid.ophid.utoronto.ca/static/download/" + iid_files[organism] + ".txt.gz' -P "+out_dir+"tmp/")
    ru.print_current_usage("Unpack file from iid database...")
    os.system("gzip -d "+out_dir+"tmp/" + iid_files[organism] + ".txt.gz")
    ru.print_current_usage('Create edge list...')
    edge_list = list()
    with open(out_dir+"tmp/" + iid_files[organism] + ".txt") as fp:
        fp.readline()  # ignore header
        for line in fp:
            line_split = line.split("\t")
            edge_list.append([line_split[id_position[id_type][0]], line_split[id_position[id_type][1]]])
    os.system("rm -R "+out_dir+"tmp")
    return edge_list


def create_network(edge_list: list):
    g = gt.Graph(directed=False)
    node_names = g.add_edge_list(edge_list, hashed=True)
    g.vertex_properties["ID"] = node_names
    return g


