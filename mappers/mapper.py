#!/usr/bin/python3

import pandas as pd
from abc import abstractmethod
from d_utils import config, mapping_utils as mu
import scipy.sparse as sp
import pickle


class Mapper:
    loaded_mappings = {'gene_ids': pd.DataFrame(), 'gene_atts': pd.DataFrame(), 'disorder_ids': pd.DataFrame(),
                       'disorder_atts': pd.DataFrame()}

    loaded_distance_ids = {'gene_mat_ids': dict(), 'disease_mat_ids': dict()}
    loaded_distances = {'go_BP': sp.coo_matrix(None), 'go_CC': sp.coo_matrix(None), 'go_MF': sp.coo_matrix(None),
                        'pathway_kegg': sp.coo_matrix(None), 'related_genes': sp.coo_matrix(None),
                        'related_variants': sp.coo_matrix(None), 'related_pathways': sp.coo_matrix(None)}

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

    def update_distance_ids(self, in_list: list, key: str):
        if self.loaded_distance_ids[key]:  # is not empty
            print("todo")
        else:
            for index, value in enumerate(in_list):
                self.loaded_distance_ids[key][index] = value

    def update_distances(self, in_mat: sp.coo_matrix, key: str):
        if self.loaded_distances[key].nnz > 0:
            print("todo")
        else:
            self.loaded_distances[key] = in_mat

    def get_loaded_distances(self, in_set, key: str):
        if not self.loaded_mappings[key].empty:
            hit_mapping = self.loaded_mappings[key]['mat'][self.loaded_mappings[key].index.isin(in_set),
                                                           self.loaded_mappings[key].columns.isin(in_set)]
            if not hit_mapping.empty:
                return hit_mapping, set(in_set) - set(hit_mapping.columns)
        return pd.DataFrame(), in_set

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
            self.loaded_distances[key] = sp.load_npz(self.file_names[key])
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
            sp.save_npz(self.file_names[key], in_object)
        else:  # in_type == "distance_id"
            with open(self.file_names[key], 'wb+') as f:
                pickle.dump(in_object, f)
