#!/bin/python3

import os
import sys
import argparse
import time
import psutil
import d_utils.config as c

start_time = time.time()


def save_parameters(script_desc: str, arguments):
    """
    Save command line options into local variables.

    :return: values assigned to input arguments
    """
    descr = "\n############################################################################\n"
    descr += "###################### DiGeSt - %(prog)s ########################\n"
    descr += script_desc
    descr += "\n############################################################################\n"
    descr += "\nusage: python3 %(prog)s [required arguments] [optional arguments]\n"
    epilo = _get_epilog(script_name=os.path.basename(sys.argv[0]))
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawTextHelpFormatter, epilog=epilo,
                                     usage=argparse.SUPPRESS, add_help=False)
    required_args = parser.add_argument_group("required arguments")
    if 'r' in arguments:
        required_args.add_argument('-r', '--reference', type=str,
                                   help='[Only for mode id-set and set-set] Reference file or id. ')
    if 'ri' in arguments:
        required_args.add_argument('-ri', '--reference_id_type', type=str,
                                   choices=c.SUPPORTED_GENE_IDS + c.SUPPORTED_DISEASE_IDS, metavar='REFERENCE_ID_TYPE',
                                   help='[Only for mode id-set and set-set] Reference id type. See possible options below.')
    if 't' in arguments:
        required_args.add_argument('-t', '--target', type=str, required=True,
                                   help='Target file with set or clusters.')
    if 'ti' in arguments:
        required_args.add_argument('-ti', '--target_id_type', type=str, required=True,
                                   choices=c.SUPPORTED_GENE_IDS + c.SUPPORTED_DISEASE_IDS, metavar='TARGET_ID_TYPE',
                                   help='Target id type. See possible options below.')
    if 'm' in arguments:
        required_args.add_argument('-m', '--mode', type=str, required=True,
                                   choices=['set', 'set-set', 'id-set', 'cluster'],
                                   help='Desired mode. See possible options below.')
    # =============================================================================
    # for network creation
    # ============================================================================
    if 'n' in arguments:
        required_args.add_argument('-n', '--node_type', type=str, required=True,
                                   choices=c.SUPPORTED_GENE_IDS + c.SUPPORTED_DISEASE_IDS, metavar='NODE_TYPE',
                                   help='ID type of nodes for network.\n(choose from entrez, ensembl, symbol, uniprot, '
                                        'mondo, omim, snomedct, umls, orpha, mesh, doid, ICD-10)')
    if 'd' in arguments:
        required_args.add_argument('-d', '--distance_type', type=str, required=True,
                                   choices=list(c.DISTANCES.values()), metavar='EDGE_TYPE',
                                   help='Distance type for edges.\n(choose from go_BP, go_CC, go_MF, pathway_kegg, '
                                        'related_genes, related_variants, related_pathways)')
    optional_args = parser.add_argument_group("optional arguments")
    if 'o' in arguments:
        optional_args.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory. [Default=./]')
    if 'e' in arguments:
        optional_args.add_argument("-e", "--enriched", action='store_true', default=False,
                                   help="Set flag, if only enriched attributes of the reference should be used.")
    if 'c' in arguments:
        optional_args.add_argument("-c", "--runs", type=int, default=c.NUMBER_OF_RANDOM_RUNS,
                                   help="Number of runs with random target values for p-value calculation.")
    if 'b' in arguments:
        optional_args.add_argument("-b", "--background_model", type=str, default="complete",
                                   choices=['complete', 'term-pres'],
                                   help="Model defining how random values should be picked. See possible options below.")
    if 'pr' in arguments:
        optional_args.add_argument("-pr", "--replace", type=int, default=100,
                                   help="Percentage of how many of the original ids should be replaced with random ids."
                                        " [Default=100]")
    if 'v' in arguments:
        optional_args.add_argument("-v", "--verbose", action='store_true', default=False,
                                   help="Set flag, if additional info like ids without assigned attributes should "
                                        "be printed.")
    optional_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
    args = parser.parse_args()
    return args


def _get_epilog(script_name):
    epilog = ""
    if script_name == 'single_validation.py':
        epilog += "\n----------------------------------------------------------------------------\n"
        epilog += "\nsupported id types\n"
        epilog += "  for genes\t\t" + ', '.join(c.SUPPORTED_GENE_IDS) + "\n"
        epilog += "  for diseases\t\t" + ', '.join(c.SUPPORTED_DISEASE_IDS) + "\n"
        epilog += "\nsupported modes\n"
        epilog += "  set\t\t\tCompare similarity inside the set. Either genes or diseases.\n"
        epilog += "  set-set\t\tCompare target set to reference set. Both either genes or diseases.\n"
        epilog += "  id-set\t\tCompare target set to reference id. Set either genes or diseases, id of disease.\n"
        epilog += "  cluster\t\tCompare cluster quality inside clustering. Either genes or diseases.\n"
        epilog += "\nsupported background models\n"
        epilog += "  complete\t\tRandom ids will be picked completely randomly.\n"
        epilog += "  deg-pres\t\tDegree preserving pick of random ids based on PPI network from IID.\n"
        epilog += "  term-pres\t\tRandom ids will preserve the number of mapped terms for the replaced ids.\n"
    epilog += "\n############################################################################\n"
    return epilog


def print_current_usage(text):
    memory_usage = '{0:.2f}'.format(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    time_usage = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
    print('[{}|{}MB] '.format(time_usage, memory_usage) + text)
