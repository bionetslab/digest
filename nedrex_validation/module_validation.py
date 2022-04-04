import sys
import random
import pandas as pd
import graph_tool as gt
import graph_tool.util as gtu
import graph_tool.topology as gtt
from graph_tool import GraphView
import input_data_generator

# =============================================================================
def print_usage():
    print(' ')
    print('        usage: python3 module_validation.py network_file module_file true_drugs_file permut_num only_approved_drugs outfile_name')
    print('        -----------------------------------------------------------------')
    print('        network_file       : The PPI or GGI network file must be provided as graphml or gt format.')
    print('                             Nodes attribute "primaryDomainId" is used as node id.')
    print('        module_file        : A table containing the proteins or genes of disease module (if table contains')
    print('                             more than one column they must be tab-separated; the first column will be read as module nodes)')
    print('        true_drugs_file    : A table containing the list of drugs indicated (true drugs) for disease (if table contains')
    print('                             more than one column they must be tab-separated; the first column will be read as drugs)')
    print('        permut_num         : The number of random modules to build background distribution')
    print('        only_approved_drugs: If only approved drugs should be considered in the ranking,')
    print('                             it should be set to "Y". If including all approved and unapproved drugs is desired should be set to "N"')
    print('        outfile_name       : P-values will be saved in this file')
    print(' ')


# =============================================================================
def check_input_style(input_list):
    try:
        network_file = input_list[1]
        module_file = input_list[2]
        true_drugs_file = input_list[3]
        permut_num = input_list[4]
        only_approved_drugs = (input_list[5].lower() == 'y')
        outfile_name = input_list[6]
    # if input is not correctly given, print out a usage message and exit
    except:
        print_usage()
        sys.exit(0)
        return
    return network_file, module_file, true_drugs_file, permut_num, only_approved_drugs, outfile_name

# =============================================================================
# def read_input(network_file, seed_file, only_direct_drugs, only_approved_drugs):
def read_input(network_file, module_file, true_drugs_file, only_approved_drugs):
    G = gt.load_graph(network_file)
    module_nodes = list(pd.read_csv(module_file, sep='\t', header=None, index_col=False)[0])
    true_drugs = set(pd.read_csv(true_drugs_file, sep='\t', header=None, index_col=False)[0])
    print('Number of true drugs to use as reference for the validation: %d' %len(true_drugs))
    if 'drugbank.' not in pd.read_csv(true_drugs_file, sep='\t', header=None, index_col=False)[0][1]:
        true_drugs = set('drugbank.'+d for d in true_drugs)
    only_approved = only_approved_drugs
    if 'entrez.' in module_nodes[1] or module_nodes[1].isnumeric():
        n_type = 'gene'
        print('Type of module nodes is: gene')
    elif 'uniprot.' in module_nodes[1] or (6 <= len(module_nodes[1]) <= 10 and module_nodes[1][0].isalpha()):
        n_type = 'protein'
        print('Type of module nodes is: protein')
    else:
        n_type = 'invalid'
        print('Invalid type of module nodes (neither protein nor gene)')
    return G, module_nodes, true_drugs, only_approved, n_type

# =============================================================================
def find_num_ccs(G, module_nodes):
    G.set_directed(False)
    module_nodes_ids = [int(gtu.find_vertex(G, prop=G.vertex_properties['primaryDomainId'], match=node)[0]) for node in module_nodes]
    v_modul_filt = G.new_vertex_property('bool')
    for i in module_nodes_ids:
        v_modul_filt[i] = True

    induced = GraphView(G, v_modul_filt)
    ees = [(G.vertex_index[e.source()], G.vertex_index[e.target()]) for e in induced.edges()]
    induced_module = gt.Graph(directed=False)
    node_id_to_node = {node_id: induced_module.add_vertex() for node_id in module_nodes}
    node_id_property = induced_module.new_vp('string', vals=module_nodes)
    induced_module.vertex_properties['primaryDomainId'] = node_id_property
    edges = []
    for e in ees:
        source = node_id_to_node[G.vertex_properties['primaryDomainId'][e[0]]]
        target = node_id_to_node[G.vertex_properties['primaryDomainId'][e[1]]]
        edges.append((source, target))
    induced_module.add_edge_list(edges)

    comp, hist = gtt.label_components(induced_module, directed=None, attractors=False)
    comp_labels = [comp.a[node] for node in range(induced_module.num_vertices())]
    ccs_num = max(comp_labels) + 1
    return ccs_num

# ===========================================================================
def generate_rand_modules(G, N, ccs_num, MS):
    random_cc_modules_list = []
    rand_module_nodes_list = []
    for r in range(N):
        random_cc_module, rand_module_nodes, iterr = _rand_module(G, N, ccs_num, MS)
        if iterr <= 20:
            random_cc_modules_list.append(random_cc_module)
            rand_module_nodes_list.append(rand_module_nodes)
        else:
            random_cc_module, rand_module_nodes, iterr = _rand_module(G, N, ccs_num, MS)
            random_cc_modules_list.append(random_cc_module)
            rand_module_nodes_list.append(rand_module_nodes)
    return rand_module_nodes_list

# ===========================================================================
def _rand_module(G, N, ccs_num, MS):
    # step 1: randomly select ccs_num number of seeds (which are not neighbors)
    G_nodes = list(range(G.num_vertices()))
    random_cc_module = dict()
    initial_seeds = []
    neighb_initial_seeds = []
    for c in range(ccs_num):
        # random_cc_modules[c] = []
        added = False
        while not added:
            s = random.choice(G_nodes)
            # s = 16810
            if s not in initial_seeds and s not in neighb_initial_seeds:
                added = True
                initial_seeds.append(s)
                neighb_initial_seeds.extend(list(G.get_all_neighbors(s)))
                random_cc_module[s] = []
                # random_cc_modules[s].append(s)
    rand_module_nodes = set()
    rand_module_nodes.update(initial_seeds)
    # step 2: expand the selected seeds by their neighbors enforcing the number of CCs remain the same
    for cs in initial_seeds:
        rn = random.choice(G.get_all_neighbors(cs))
        intersect = False
        for k in random_cc_module.keys():
            if k != cs and (len(set(G.get_all_neighbors(rn)).intersection(set(random_cc_module[k]))) > 0 or rn in G.get_all_neighbors(k)):
                intersect = True
                break
        if not intersect:
            random_cc_module[cs].append(rn)
            rand_module_nodes.add(rn)
        if len(rand_module_nodes) >= MS:
            break
    iterr = 0
    while len(rand_module_nodes) < MS and iterr < 20:
        iterr += 1
        for k in initial_seeds:
            vv = random_cc_module[k]
            # for k, vv in random_cc_modules.items():
            for v in vv:
                rn = random.choice(G.get_all_neighbors(v))
                if rn not in vv:
                    intersect = False
                    for i in random_cc_module.keys():
                        if i != k and (len(set(G.get_all_neighbors(rn)).intersection(set(random_cc_module[i]))) > 0 or rn in G.get_all_neighbors(i)):
                            intersect = True
                            break
                    if not intersect:
                        random_cc_module[k].append(rn)
                        rand_module_nodes.add(rn)
                if len(rand_module_nodes) >= MS:
                    break
            if len(rand_module_nodes) >= MS:
                break
    return random_cc_module, rand_module_nodes, iterr

# ===========================================================================
if __name__ == '__main__':
    # -----------------------------------------------------
    # Checking for input from the command line:
    # -----------------------------------------------------
    #
    # [1] File providing the PPI or GGI network in graphml or gt format, as nodes uniport AC or entrez ID are accepted respectively.
    #     (valid attribute for node id: "primaryDomainId")
    #
    # [2] File containing the disease module proteins/genes, uniport AC or entrez ID are accepted respectively (if table contains more than one
    #     column they must be tab-separated; the first column will be used only)
    #
    # [3] File containing the list of drugs indicated to treat the disease (i.e. the list of drugs that should be considered true positive).
    #     Acceptable ID for drugs is drugbank ID. (if table contains more than one column they must be tab-separated; the first column will be used only)
    #
    # [4] The number of random modules to build background distribution for empirical p-value computation (at least 1000)
    #
    # [5] Letter "Y" to include only approved drugs in p-value computation or letter "N" to include all the drugs
    #
    # [6] name of the output file
    #
    # Returns one empirical p-value for the result (to-be-validated) disease module

    input_list = sys.argv
    # check if input style is correct
    network_file, module_file, true_drugs_file, permut_num, only_approved_drugs, outfile_name = check_input_style(input_list)
    print('Network file: ' + network_file)
    print('Input module file: ' + module_file)
    print('Input true drugs file: ' + true_drugs_file)
    print('Number of permutations: ' + permut_num)
    print('Only approved drugs: ' + str(only_approved_drugs))
    print('Output file: ' + outfile_name)

    # read the network, input module's nodes and input true drugs:
    G, module_nodes, true_drugs, only_approved, n_type = read_input(network_file, module_file, true_drugs_file, only_approved_drugs)
    N = int(permut_num)
    # Get list of all NeDRex' drugs
    dbids = input_data_generator.get_drugs_nedrex(only_approved)
    # Get protein/Gene-drug data
    drug_map = input_data_generator.protein_drug_api(only_approved)
    if n_type == 'gene':
        drug_map = input_data_generator.gene_drug_api(only_approved, drug_map)
    # rename result_drugs to drugs_targeting_module and define a function using API to get this set of drugs, consider only_approved boolean as well
    drugs_targeting_module = input_data_generator.drugs_targeting_module_api(set(module_nodes), n_type, only_approved, dbids)
    print('Number of overlapping drugs between set of drugs targeting the input disease module and reference list of drugs: %d' %len(true_drugs.intersection(drugs_targeting_module)))
    # Find number of CCs for the input module
    ccs_num = find_num_ccs(G, module_nodes)
    # Generate random modules of matched number of connected components
    rand_module_nodes_list = generate_rand_modules(G, N, ccs_num, len(module_nodes))
    node_ids = [G.vertex_properties["primaryDomainId"][node] for node in range(G.num_vertices())]
    ### get drugs targeting random modules
    drugs_rand_modules_list = [set() for i in range(N)]
    for i in range(N):
        for n in rand_module_nodes_list[i]:
            if node_ids[n] in drug_map.keys():
                drugs_rand_modules_list[i].update(drug_map[node_ids[n]])

    ### computing empirical p-values:
    result_true_freq = len(drugs_targeting_module.intersection(true_drugs))
    rand_true_freq = []
    for i in range(N):
        rand_true_freq.append(len(drugs_rand_modules_list[i].intersection(true_drugs)))
    greater = 0
    for i in range(N):
        if rand_true_freq[i] >= result_true_freq:
            greater += 1
    empiric_pval = float(greater+1)/float(N+1)

    ### computing empirical p-values for precision:
    result_true_prec = float(len(drugs_targeting_module.intersection(true_drugs)))/float(len(drugs_targeting_module))
    rand_true_prec = []
    for i in range(N):
        if len(drugs_rand_modules_list[i]) != 0:
            rand_true_prec.append(float(len(drugs_rand_modules_list[i].intersection(true_drugs)))/float(len(drugs_rand_modules_list[i])))
        elif len(drugs_rand_modules_list[i]) == 0:
            rand_true_prec.append(0.0)
    greater_prec = 0
    for i in range(N):
        if rand_true_prec[i] >= result_true_prec:
            greater_prec += 1
    empiric_pval_prec = float(greater_prec+1)/float(N+1)

    with open(outfile_name, 'w') as f:
        f.write(f"The computed empirical p-value for the input disease module based on {N} permutations: %.6f" % empiric_pval)
        f.write('\n')
        f.write(f"The computed empirical p-value (precision-based) for the input disease module based on {N} permutationswith: %.6f" % empiric_pval_prec)

    print("\n results have been saved to '%s' \n" % outfile_name)
    print("empirical p-value: '%.6f' \n" % empiric_pval)
    print("empirical precision-based p-value: '%.6f' \n" % empiric_pval_prec)
