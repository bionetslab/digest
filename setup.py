#!/usr/bin/python3

import d_utils.runner_utils as du
import digest.mappers.mapping_transformer as mt


if __name__ == "__main__":
    desc = "[Setup] Transform disease mapping file."
    args = du.save_parameters(script_desc=desc, arguments=('m', 'o'))
    mt.transform_mapping(mapping_file=args.disease_mapping, out_dir=args.out_dir)