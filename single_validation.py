#!/usr/bin/python3

import pandas as pd
from d_utils import runner_utils as ru, config
from evaluation import comparator as comp


def single_validation(args):
    """
    Single validation of a set, cluster a id versus set and set versus set.

    :param args: user given parameters
    :return:
    """
    if args.mode == "set":
        print(comp.compare_set(id_set=pd.read_csv(args.target, header=None, sep="\t")[0], id_type=args.target_id_type))
    if args.mode == "set-set":
        print(comp.compare_set_to_set(ref=pd.read_csv(args.reference, header=None, sep="\t")[0],
                                      ref_id_type=args.reference_id_type,
                                      tar=pd.read_csv(args.target, header=None, sep="\t")[0],
                                      tar_id_type=args.target_id_type))
    if args.mode == "id-set":
        print(comp.compare_id_to_set(ref_id=args.reference, ref_id_type=args.reference_id_type,
                                     tar=pd.read_csv(args.target, header=None, sep="\t")[0],
                                     tar_id_type=args.target_id_type))
    if args.mode == "cluster":
        print(
            comp.compare_clusters(clusters=pd.read_csv(args.target, header=None, sep="\t"),
                                  id_type=args.target_id_type))
    return


def get_random_runs_values(id_len, id_type):
    pd.read_csv(file=config.FILES_DIR + 'random_runs.csv', sep=",")
    return


if __name__ == "__main__":
    desc = "            Evaluation of disease and gene sets and clusters."
    arguments = ru.save_parameters(script_desc=desc, arguments=('r', 'ri', 't', 'ti', 'm', 'o'))
    single_validation(args=arguments)
