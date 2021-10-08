#!/bin/python3

import os
import argparse
import time
import psutil


def save_parameters(script_desc, arguments):
    """
    Save command line options into local variables.

    :return: values assigned to input arguments
    """
    global start_time
    start_time = time.time()
    descr = "\n############################################################################\n"
    descr += "###################### DiGeSt - %(prog)s ########################\n"
    descr += script_desc
    descr += "\nusage: python3 %(prog)s [required arguments] [optional arguments]\n"
    epilo = "\n############################################################################\n"
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawTextHelpFormatter, epilog=epilo,
                                     usage=argparse.SUPPRESS)
    required_args = parser.add_argument_group("required arguments")
    if 'n' in arguments:
        required_args.add_argument('-n', '--number_of_nodes', type=int, required=True,
                                   help='Number of nodes.')
    if 'g' in arguments:
        required_args.add_argument('-g', '--ggi_network', type=str, required=True,
                                   help='Gene-Gene Interaction network to get diseases from.')
    if 's' in arguments:
        required_args.add_argument('-s', '--gene_to_snp', type=str, required=True,
                                   help='Gene to SNP mapping file.')
    if 'c' in arguments:
        required_args.add_argument('-c', '--comorbidity', type=str, required=True,
                                   help='File with the comorbidity of diseases.')
    if 'd' in arguments:
        required_args.add_argument('-d', '--diseases', type=int, required=True,
                                   help='Number of diseases that element should have.')
    if 'i' in arguments:
        required_args.add_argument('-i', '--in_dir', type=str, required=True,
                                   help='Input directory with needed files.')
    optional_args = parser.add_argument_group("optional arguments")
    if 'm' in arguments:
        optional_args.add_argument('-m', '--number_deviation', type=int, default=0,
                                   help='Acceptable deviation for randomly chosen number of nodes. [Default=0]')
    if 'v' in arguments:
        optional_args.add_argument('-v', '--vertex_identity', type=str, default='',
                                   help='File with elements for possible vertex ids e.g. gene symbols. [Default='']')
    if 'e' in arguments:
        optional_args.add_argument('-e', '--number_elements', type=int, default=1,
                                   help='Number of elements to create from that type. [Default=1]')
    if 'o' in arguments:
        optional_args.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory. [Default=./]')
    args = parser.parse_args()
    return args


def print_current_usage(text):
    memory_usage = '{0:.2f}'.format(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    time_usage = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
    print('[{}|{}MB]'.format(time_usage, memory_usage) + text)

