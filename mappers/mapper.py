#!/usr/bin/python3

import pandas as pd
import numpy as np
from abc import abstractmethod
from d_utils import config, mapping_utils as mu
import scipy.sparse as sp
import pickle


class Mapper:
    loaded_mappings = {'gene_ids': pd.DataFrame(), 'gene_atts': pd.DataFrame(), 'disorder_ids': pd.DataFrame(),
                       'disorder_atts': pd.DataFrame()}

    loaded_distance_ids = {'gene_mat_ids': dict(), 'disease_mat_ids': dict()}
    loaded_distances = {'go_BP': sp.csr_matrix(None), 'go_CC': sp.csr_matrix(None), 'go_MF': sp.csr_matrix(None),
                        'pathway_kegg': sp.csr_matrix(None), 'related_genes': sp.csr_matrix(None),
                        'related_variants': sp.csr_matrix(None), 'related_pathways': sp.csr_matrix(None)}

    @abstractmethod
    def load_mappings(self, set_type: str):
        pass

    @abstractmethod
    def load_distances(self, set_type: str):
        pass

    @abstractmethod
    def load_file(self, key: str, in_type: str):
        pass

    def update_mappings(self, in_df: pd.DataFrame(), key: str):
        if not self.loaded_mappings[key].empty:
            self.loaded_mappings[key] = pd.concat([self.loaded_mappings[key], in_df], ignore_index=True)
        else:
            self.loaded_mappings[key] = in_df

    def get_loaded_mapping(self, in_set, id_type: str, key: str):
        if not self.loaded_mappings[key].empty:
            hit_mapping = self.loaded_mappings[key].loc[self.loaded_mappings[key][id_type].isin(in_set)]
            if not hit_mapping.empty:
                return hit_mapping, set(in_set) - set(hit_mapping[id_type])
        return pd.DataFrame(), in_set

    def update_distance_ids(self, in_series: pd.Series, key: str) -> pd.Series:
        if self.loaded_distance_ids[key]:  # is not empty
            for index, value in enumerate(iterable=in_series[len(self.loaded_distance_ids[key].keys()):],
                                          start=len(self.loaded_distance_ids[key].keys())):
                self.loaded_distance_ids[key][value] = index
            return in_series[len(self.loaded_distance_ids[key].keys()):]
        else:
            for index, value in enumerate(iterable=in_series, start=0):
                self.loaded_distance_ids[key][value] = index
            return in_series

    def update_distances(self, in_mat: sp.coo_matrix, key: str):
        if self.loaded_distances[key].nnz > 0:
            old_mat = self.loaded_distances[key].tocoo()
            row = np.concatenate((old_mat.row, in_mat.row), axis=None)
            col = np.concatenate((old_mat.col, in_mat.col), axis=None)
            data = np.concatenate((old_mat.data, in_mat.data), axis=None)
            self.loaded_distances[key] = sp.csr_matrix((data, (row, col)), shape=(
                len(self.loaded_distance_ids[key].keys()), len(self.loaded_distance_ids[key].keys())))
        else:
            self.loaded_distances[key] = in_mat.tocsr()

    def get_loaded_distance(self, id1, id2, id_type: str, key: str):
        if id1 in self.loaded_distance_ids[id_type] and id2 in self.loaded_distance_ids[id_type]:
            return self.loaded_distances[key][self.loaded_distance_ids[id_type][id1],
                                              self.loaded_distance_ids[id_type][id2]]
        return None

    def get_loaded_distances(self, in_set, id_type: str, key: str):
        if self.loaded_distance_ids[id_type]:  # is not empty
            index_to_id, hit_values = dict(), list()
            for value in in_set:
                if value in self.loaded_distance_ids[id_type]:
                    index_to_id[self.loaded_distance_ids[id_type][value]] = value
            cur_mat = self.loaded_distances[key].tocoo()
            for index, value in enumerate(cur_mat.row):
                if value in index_to_id:
                    if cur_mat.col[index] in index_to_id:
                        hit_values.append(cur_mat.data[index])
            return hit_values, set(in_set) - set(self.loaded_distance_ids[id_type].keys())
        return list(), in_set

    def get_full_set(self, id_type: str, mapping_name: str) -> set:
        return set(self.loaded_mappings[mapping_name][id_type])

    @abstractmethod
    def save_mappings(self):
        pass

    @abstractmethod
    def save_distances(self):
        pass

    @abstractmethod
    def save_file(self, in_object, key: str, in_type: str):
        pass

    def drop_mappings(self):
        self.loaded_mappings = {'gene_ids': pd.DataFrame(), 'gene_atts': pd.DataFrame(), 'disorder_ids': pd.DataFrame(),
                                'disorder_atts': pd.DataFrame()}


class FileMapper(Mapper):
    file_names = {'gene_ids': config.FILES_DIR + 'gene_id_mapping.csv',
                  'gene_atts': config.FILES_DIR + 'gene_att_mapping.csv',
                  'disorder_ids': config.FILES_DIR + 'disease_id_mapping.csv',
                  'disorder_atts': config.FILES_DIR + 'disease_att_mapping.csv',
                  'gene_mat_ids': config.FILES_DIR + 'gene_mat_ids.pkl',
                  'disease_mat_ids': config.FILES_DIR + 'disease_mat_ids.pkl',
                  'go_BP': config.FILES_DIR + 'gene_dist_go_BP.npz',
                  'go_CC': config.FILES_DIR + 'gene_dist_go_CC.npz',
                  'go_MF': config.FILES_DIR + 'gene_dist_go_MF.npz',
                  'pathway_kegg': config.FILES_DIR + 'gene_dist_pathway_kegg.npz',
                  'related_genes': config.FILES_DIR + 'disease_dist_rel_genes.npz',
                  'related_variants': config.FILES_DIR + 'disease_dist_rel_variants.npz',
                  'related_pathways': config.FILES_DIR + 'disease_dist_rel_pathways.npz'}

    def load_mappings(self, set_type):
        if set_type in config.SUPPORTED_GENE_IDS:
            for mapping_key in ['gene_ids', 'gene_atts']:
                self.load_file(key=mapping_key, in_type='mapping')
        else:  # if set_type in config.SUPPORTED_DISEASE_IDS
            for mapping_key in ['disorder_ids', 'disorder_atts']:
                self.load_file(key=mapping_key, in_type='mapping')

    def _load_file_mapping(self, file, sep, mapping_name):
        """
        Get previous mappings from file.

        :param file: file with previous distances
        :param sep: separator of values in file
        :param mapping_name: key name of local dictionary for saving data
        :return: dataframe with previous mapping
        """
        # ===== Get mapping from local mapping file =====
        mapping = pd.read_csv(file, sep=sep, header=0, dtype=str)
        if mapping_name == "disorder_ids":
            icd_unstack = mu.split_and_expand_column(data=mapping, split_string=",", column_name="ICD-10")
            mapping = pd.concat([icd_unstack, mapping[mapping['ICD-10'] != '']])
        # ===== Save mapping to local dictionary =====
        self.loaded_mappings[mapping_name] = mapping

    def load_distances(self, set_type):
        if set_type in config.SUPPORTED_GENE_IDS:
            self.load_file(key='gene_mat_ids', in_type='distance_id')
            for distance_key in ['go_BP', 'go_CC', 'go_MF', 'pathway_kegg']:
                self.load_file(key=distance_key, in_type='distance')
        else:  # if set_type in config.SUPPORTED_DISEASE_IDS
            self.load_file(key='disease_mat_ids', in_type='distance_id')
            for distance_key in ['related_genes', 'related_variants', 'related_pathways']:
                self.load_file(key=distance_key, in_type='distance')

    def load_file(self, key: str, in_type: str):
        if in_type == "mapping":
            self._load_file_mapping(file=self.file_names[key], sep=",", mapping_name=key)
        elif in_type == "distance":
            self.loaded_distances[key] = sp.load_npz(self.file_names[key]).tocrs()
        else:  # in_type == "distance_id"
            with open(self.file_names[key], 'rb') as f:
                self.loaded_distance_ids[key] = pickle.load(f)

    def save_mappings(self):
        for mapping_key in ['gene_ids', 'gene_atts', 'disorder_ids', 'disorder_atts']:
            self.save_file(in_object=self.loaded_mappings[mapping_key], key=mapping_key, in_type='mapping')

    def save_distances(self):
        for distance_id_key in ['gene_mat_ids', 'disease_mat_ids']:
            self.save_file(in_object=self.loaded_distance_ids[distance_id_key], key=distance_id_key,
                           in_type='distance_id')
        for distance_key in ['go_BP', 'go_CC', 'go_MF', 'pathway_kegg', 'related_genes', 'related_variants',
                             'related_pathways']:
            self.save_file(in_object=self.loaded_distances[distance_key], key=distance_key, in_type='distance')

    def save_file(self, in_object, key: str, in_type: str):
        if in_type == "mapping":
            in_object.to_csv(self.file_names[key], index=True)
        elif in_type == "distance":
            sp.save_npz(self.file_names[key].tocoo(), in_object)
        else:  # in_type == "distance_id"
            with open(self.file_names[key], 'wb+') as f:
                pickle.dump(in_object, f)
