#!/usr/bin/python3

import pandas as pd
import numpy as np
import d_utils.runner_utils as du


def single_validation():
    return





if __name__ == "__main__":
    desc = "Evaluation of gene sets against seed gene sets or disease id; disease id sets and disease clusters."
    args = du.save_parameters(script_desc=desc, arguments=('m', 'o'))
    single_validation(mapping_file=args.disease_mapping, out_dir=args.out_dir)