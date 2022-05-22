import pandas as pd
from .mappers import mapping_utils as mu
from .mappers.mapper import Mapper
from . import config as c
from abc import abstractmethod
import random
# only in full biodigest
import graph_tool as gt
import graph_tool.util as gtu
import graph_tool.topology as gtt
from graph_tool import GraphView


class BackgroundModel():

    @abstractmethod
    def get_module(self, **kwargs):
        pass


class CompleteModel(BackgroundModel):

    def __init__(self, prev_id_type, full_id_map):
        self.full_id_set = full_id_map[full_id_map[c.ID_TYPE_KEY[prev_id_type]] != ""][
            c.ID_TYPE_KEY[prev_id_type]].tolist()

    def get_module(self, to_replace):
        random_sample = set(random.sample(self.full_id_set, len(to_replace)))
        return random_sample


class TermPresModel(BackgroundModel):

    def __init__(self, mapper: Mapper, prev_id_type, new_id_type, map_id_type, map_att_type, term):
        # prepare candidates
        att_map = mu.map_to_prev_id(main_id_type=c.ID_TYPE_KEY[new_id_type],
                                    id_type=c.ID_TYPE_KEY[prev_id_type],
                                    id_mapping=mapper.loaded_mappings[map_id_type],
                                    att_mapping=mapper.loaded_mappings[map_att_type])
        self.atts_to_size(pd_map=att_map)
        self.size_mapping_to_dict(pd_size_map=self.att_len, id_col=c.ID_TYPE_KEY[prev_id_type], term_col=term,
                                  threshold=100)

    def get_module(self, to_replace, term, prev_id_type):
        random_sample = set()
        for replace_id in to_replace:
            if replace_id in self.size_mapping:  # only if id is mappable to other ids
                random_sample.add(
                    self.att_len[self.att_len[term].isin(self.size_mapping[replace_id])][
                        c.ID_TYPE_KEY[prev_id_type]].sample(
                        n=1).values[0])
        return random_sample

    def atts_to_size(self, pd_map: pd.DataFrame):
        att_len = pd_map.copy()
        att_len[att_len.columns[1:]] = att_len[att_len.columns[1:]].applymap(mu.set_to_len)
        att_len['sum'] = att_len[att_len.columns[1:]].sum(axis=1)
        self.att_len = att_len

    def size_mapping_to_dict(self, pd_size_map: pd.DataFrame, id_col: str, term_col: str, threshold: int = 100):
        size_to_occ = pd.DataFrame(pd_size_map[term_col].value_counts()).sort_index().to_dict()[term_col]
        pd_size_map = pd_size_map.sort_values(by=[id_col]).reset_index(drop=True)
        new_dict = dict()
        term_sizes = pd_size_map[term_col].unique().tolist()
        for index, key in enumerate(term_sizes):
            curr_keys = [key]
            if size_to_occ[key] < threshold:
                sum_tmp, add_top, add_bottom = size_to_occ[key], index, index
                while sum_tmp < threshold:
                    if add_top - 1 >= 0:
                        add_top = add_top - 1
                        sum_tmp = sum_tmp + size_to_occ[term_sizes[add_top]]
                        curr_keys.append(term_sizes[add_top])
                    if add_bottom + 1 < len(term_sizes):
                        add_bottom = add_bottom + 1
                        sum_tmp = sum_tmp + size_to_occ[term_sizes[add_bottom]]
                        curr_keys.append(term_sizes[add_bottom])
            for cur_id in pd_size_map[pd_size_map[term_col] == key][id_col]:
                new_dict[cur_id] = curr_keys
        self.size_mapping = new_dict


class NetworkModel(BackgroundModel):

    def __init__(self, network_data: dict, to_replace, N):
        self.network_type = network_data['id_type']
        if network_data["network_file"].endswith(".sif"):
            df = pd.read_csv(network_data["network_file"], sep="\t", header=None)
            G = gt.Graph(directed=False)
            v_ids = G.add_edge_list(df[[0, 2]].values, hashed=True)
            G.vertex_properties['id'] = v_ids
        else:
            G = gt.load_graph(network_data["network_file"])
            G.vertex_properties['id'] = G.vertex_properties[network_data['prop_name']]
        # Reduce to ids that are present in the network
        to_replace_filtered = to_replace.intersection(set(G.vertex_properties['id']))
        # Find number of CCs for the input module
        ccs_num = self.find_num_ccs(G, to_replace_filtered)
        # Generate random modules of matched number of connected components
        self.generate_rand_modules(G, N, ccs_num, len(to_replace_filtered))
        self.node_ids = {node: G.vertex_properties["id"][node] for node in range(G.num_vertices())}

    def get_module(self, index):
        return set([self.node_ids[node] for node in self.rand_module_nodes_list[index]])

    def find_num_ccs(self, G, module_nodes):
        G.set_directed(False)
        module_nodes_ids = [int(gtu.find_vertex(G, prop=G.vertex_properties['id'], match=node)[0]) for node
                            in module_nodes]
        v_modul_filt = G.new_vertex_property('bool')
        for i in module_nodes_ids:
            v_modul_filt[i] = True

        induced = GraphView(G, v_modul_filt)
        ees = [(G.vertex_index[e.source()], G.vertex_index[e.target()]) for e in induced.edges()]
        induced_module = gt.Graph(directed=False)
        node_id_to_node = {node_id: induced_module.add_vertex() for node_id in module_nodes}
        node_id_property = induced_module.new_vp('string', vals=module_nodes)
        induced_module.vertex_properties['id'] = node_id_property
        edges = []
        for e in ees:
            source = node_id_to_node[G.vertex_properties['id'][e[0]]]
            target = node_id_to_node[G.vertex_properties['id'][e[1]]]
            edges.append((source, target))
        induced_module.add_edge_list(edges)

        comp, hist = gtt.label_components(induced_module, directed=None, attractors=False)
        comp_labels = [comp.a[node] for node in range(induced_module.num_vertices())]
        ccs_num = max(comp_labels) + 1
        return ccs_num

    def generate_rand_modules(self, G, N, ccs_num, MS):
        random_cc_modules_list = []
        rand_module_nodes_list = []
        for r in range(N):
            random_cc_module, rand_module_nodes, iterr = self._rand_module(G, ccs_num, MS)
            if iterr <= 20:
                random_cc_modules_list.append(random_cc_module)
                rand_module_nodes_list.append(rand_module_nodes)
            else:
                random_cc_module, rand_module_nodes, iterr = self._rand_module(G, ccs_num, MS)
                random_cc_modules_list.append(random_cc_module)
                rand_module_nodes_list.append(rand_module_nodes)
            if len(rand_module_nodes) == 0:
                print("tja...")
        self.rand_module_nodes_list = rand_module_nodes_list


    def _rand_module(self, G, ccs_num, MS):
        # step 1: randomly select ccs_num number of seeds (which are not neighbors)
        G_nodes = list(range(G.num_vertices()))
        random_cc_module = dict()
        initial_seeds = []
        neighb_initial_seeds = []
        for c in range(ccs_num):
            added = False
            while not added:
                s = random.choice(G_nodes)

                if s not in initial_seeds and s not in neighb_initial_seeds and len(G.get_all_neighbors(s)>0):
                    added = True
                    initial_seeds.append(s)
                    neighb_initial_seeds.extend(list(G.get_all_neighbors(s)))
                    random_cc_module[s] = []

        rand_module_nodes = set()
        rand_module_nodes.update(initial_seeds)
        # step 2: expand the selected seeds by their neighbors enforcing the number of CCs remain the same
        for cs in initial_seeds:
            #print(G.get_all_neighbors(cs))
            #print(cs)
            rn = random.choice(G.get_all_neighbors(cs))
            intersect = False
            for k in random_cc_module.keys():
                if k != cs and (len(set(G.get_all_neighbors(rn)).intersection(
                        set(random_cc_module[k]))) > 0 or rn in G.get_all_neighbors(k)):
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
                            if i != k and (len(set(G.get_all_neighbors(rn)).intersection(
                                    set(random_cc_module[i]))) > 0 or rn in G.get_all_neighbors(i)):
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