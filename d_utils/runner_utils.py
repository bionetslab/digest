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
                                     usage=argparse.SUPPRESS, add_help=False)
    required_args = parser.add_argument_group("required arguments")
    if 'm' in arguments:
        required_args.add_argument('-m', '--disease_mapping', type=str, required=True,
                                   help='Disease mapping file.')
    if 'r' in arguments:
        required_args.add_argument('-r', '--reference', type=str, required=True,
                                   help='Disease mapping file.')
    optional_args = parser.add_argument_group("optional arguments")
    if 'o' in arguments:
        optional_args.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory. [Default=./]')
    optional_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
    args = parser.parse_args()
    return args


def print_current_usage(text):
    memory_usage = '{0:.2f}'.format(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    time_usage = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
    print('[{}|{}MB]'.format(time_usage, memory_usage) + text)

