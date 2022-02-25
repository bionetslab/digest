#!/usr/bin/python3

import pandas as pd
import numpy as np
from pathlib import Path
from abc import abstractmethod
from .. import config
from . import mapping_utils as mu
import scipy.sparse as sp
import pickle
import os


class Mapper:
    loaded_mappings = {'gene_ids': pd.DataFrame(), 'gene_atts': pd.DataFrame(), 'disorder_ids': pd.DataFrame(),
                       'disorder_atts': pd.DataFrame()}

    loaded_distance_ids = {'jaccard': {'gene_mat_ids': dict(), 'disease_mat_ids': dict()},
                           'overlap': {'gene_mat_ids': dict(), 'disease_mat_ids': dict()}}
    loaded_distances = {
        'jaccard': {'go_BP': sp.csr_matrix(None), 'go_CC': sp.csr_matrix(None), 'go_MF': sp.csr_matrix(None),
                    'pathway_kegg': sp.csr_matrix(None), 'related_genes': sp.csr_matrix(None),
                    'related_variants': sp.csr_matrix(None), 'related_pathways': sp.csr_matrix(None)},
        'overlap': {'go_BP': sp.csr_matrix(None), 'go_CC': sp.csr_matrix(None), 'go_MF': sp.csr_matrix(None),
                    'pathway_kegg': sp.csr_matrix(None), 'related_genes': sp.csr_matrix(None),
                    'related_variants': sp.csr_matrix(None), 'related_pathways': sp.csr_matrix(None)}}

    changed_mappings = set()

    load: bool = True

    def __init__(self, preload: bool = False):
        if preload:
            self.load_mappings()
            self.load_distances(set_type="entrez", distance_measure="jaccard")
            self.load_distances(set_type="mondo", distance_measure="jaccard")
            self.load_distances(set_type="entrez", distance_measure="overlap")
            self.load_distances(set_type="mondo", distance_measure="overlap")
            self.load = False

    @abstractmethod
    def load_mappings(self):
        pass

    @abstractmethod
    def load_distances(self, set_type: str, distance_measure: str):
        pass

    @abstractmethod
    def load_file(self, key: str, in_type: str):
        pass

    def get_loaded_mapping(self, in_set, id_type: str, key: str):
        if not self.loaded_mappings[key].empty:
            current_mapping = self.loaded_mappings[key].copy()
            if id_type not in ['entrezgene', 'mondo']:
                current_mapping = current_mapping.explode(id_type)
            hit_mapping = current_mapping.loc[current_mapping[id_type].isin(in_set)]
            if not hit_mapping.empty:
                return hit_mapping, set(in_set) - set(hit_mapping[id_type])
        return pd.DataFrame(), in_set

    def update_mappings(self, in_df: pd.DataFrame(), key: str):
        self.changed_mappings.add(key)
        if not self.loaded_mappings[key].empty:
            self.loaded_mappings[key] = pd.concat([self.loaded_mappings[key], in_df], ignore_index=True)
        else:
            self.loaded_mappings[key] = in_df

    def get_loaded_mapping_ids(self, in_ids, id_type: str) -> pd.DataFrame:
        if id_type in config.SUPPORTED_GENE_IDS:
            mapping, _ = self.get_loaded_mapping(in_set=in_ids, id_type=config.ID_TYPE_KEY[id_type], key="gene_ids")
            return mapping
        else:  # if set_type in config.SUPPORTED_DISEASE_IDS
            mapping, _ = self.get_loaded_mapping(in_set=in_ids, id_type=config.ID_TYPE_KEY[id_type], key="disorder_ids")
            return mapping

    def update_distance_ids(self, in_series: pd.Series, key: str, distance_measure: str) -> pd.Series:
        self.changed_mappings.add(distance_measure + "_" + key)
        if self.loaded_distance_ids[distance_measure][key]:  # is not empty
            if len(self.loaded_distance_ids[distance_measure][key].keys()) < len(in_series):
                new_ids = in_series[len(self.loaded_distance_ids[distance_measure][key].keys()):]
                for index, value in enumerate(iterable=new_ids,
                                              start=len(self.loaded_distance_ids[distance_measure][key].keys())):
                    self.loaded_distance_ids[distance_measure][key][value] = index
                return new_ids
            else:
                return pd.Series([],dtype=pd.StringDtype())
        else:
            for index, value in enumerate(iterable=in_series, start=0):
                self.loaded_distance_ids[distance_measure][key][value] = index
            return in_series

    def get_loaded_distances(self, in_series: pd.Series, id_type: str, key: str, distance_measure: str,
                             to_series: pd.Series = None) -> sp.csr_matrix:
        if self.loaded_distance_ids[distance_measure][id_type]:  # is not empty
            indices = list()
            for element in in_series:
                indices.append(self.loaded_distance_ids[distance_measure][id_type][element])
            if to_series is not None:
                to_indices = list()
                for element in to_series:
                    to_indices.append(self.loaded_distance_ids[distance_measure][id_type][element])
                return self.loaded_distances[distance_measure][key][indices, :][:, to_indices]
            else:
                return self.loaded_distances[distance_measure][key][indices, :][:, indices]
        else:
            return sp.csr_matrix(None)

    def update_distances(self, in_mat: sp.coo_matrix, id_type: str, key: str, distance_measure: str):
        self.changed_mappings.add(distance_measure + "_" + key)
        if self.loaded_distances[distance_measure][key].nnz > 0:
            old_mat = self.loaded_distances[distance_measure][key].tocoo()
            row = np.concatenate((old_mat.row, in_mat.row), axis=None)
            col = np.concatenate((old_mat.col, in_mat.col), axis=None)
            data = np.concatenate((old_mat.data, in_mat.data), axis=None)
            self.loaded_distances[distance_measure][key] = sp.csr_matrix((data, (row, col)), shape=(
                len(self.loaded_distance_ids[distance_measure][id_type].keys()),
                len(self.loaded_distance_ids[distance_measure][id_type].keys())))
        else:
            self.loaded_distances[distance_measure][key] = in_mat.tocsr()

    def get_full_set(self, id_type: str, mapping_name: str) -> pd.DataFrame:
        current_mapping = self.loaded_mappings[mapping_name]
        if id_type not in ['entrezgene', 'mondo']:
            current_mapping = current_mapping.explode(id_type).fillna("")
        return current_mapping

    @abstractmethod
    def save_mappings(self):
        pass

    @abstractmethod
    def save_distances(self):
        pass

    @abstractmethod
    def save_file(self, in_object, key: str, in_type: str):
        pass

    @abstractmethod
    def check_for_setup_sources(self):
        pass

    def drop_mappings(self):
        self.loaded_mappings = {'gene_ids': pd.DataFrame(), 'gene_atts': pd.DataFrame(), 'disorder_ids': pd.DataFrame(),
                                'disorder_atts': pd.DataFrame()}


class FileMapper(Mapper):
    file_names = {'gene_ids': 'gene_id_mapping.csv',
                  'gene_atts': 'gene_att_mapping.csv',
                  'disorder_ids': 'disease_id_mapping.csv',
                  'disorder_atts': 'disease_att_mapping.csv',
                  'gene_mat_ids': 'gene_mat_ids.pkl',
                  'disease_mat_ids': 'disease_mat_ids.pkl',
                  'go_BP': 'gene_dist_go_BP.npz',
                  'go_CC': 'gene_dist_go_CC.npz',
                  'go_MF': 'gene_dist_go_MF.npz',
                  'pathway_kegg': 'gene_dist_pathway_kegg.npz',
                  'related_genes': 'disease_dist_rel_genes.npz',
                  'related_variants': 'disease_dist_rel_variants.npz',
                  'related_pathways': 'disease_dist_rel_pathways.npz'}

    def __init__(self, preload: bool = False, files_dir=config.FILES_DIR):
        self.files_dir = files_dir
        super().__init__(preload=preload)

    def load_mappings(self):
        if self.load:
            for mapping_key in ['gene_atts', 'disorder_atts', 'gene_ids', 'disorder_ids']:
                self.load_file(key=mapping_key, in_type='mapping')
                self.loaded_mappings[mapping_key][self.loaded_mappings[mapping_key].columns[1:]] = \
                    self.loaded_mappings[mapping_key][self.loaded_mappings[mapping_key].columns[1:]].fillna(
                        '').applymap(mu.string_to_set)

    def _load_file_mapping(self, file, sep, mapping_name):
        """
        Get previous mappings from file.

        :param file: file with previous distances
        :param sep: separator of values in file
        :param mapping_name: key name of local dictionary for saving data
        :return: dataframe with previous mapping
        """
        # ===== Get mapping from local mapping file =====
        mapping = pd.read_csv(file, sep=sep, header=0, dtype=str).fillna('')
        if mapping_name == "disorder_ids":
            icd_unstack = mu.split_and_expand_column(data=mapping, split_string=",", column_name="ICD-10")
            mapping = pd.concat([icd_unstack, mapping[mapping['ICD-10'] != '']])
        # ===== Save mapping to local dictionary =====
        self.loaded_mappings[mapping_name] = mapping

    def load_distances(self, set_type: str, distance_measure: str):
        if self.load:
            if set_type in config.SUPPORTED_GENE_IDS:
                self.load_file(key='gene_mat_ids', in_type='gene_id', distance_measure=distance_measure)
                for distance_key in ['go_BP', 'go_CC', 'go_MF', 'pathway_kegg']:
                    self.load_file(key=distance_key, in_type='distance', distance_measure=distance_measure)
            else:  # if set_type in config.SUPPORTED_DISEASE_IDS
                self.load_file(key='disease_mat_ids', in_type='distance_id', distance_measure=distance_measure)
                for distance_key in ['related_genes', 'related_variants', 'related_pathways']:
                    self.load_file(key=distance_key, in_type='distance', distance_measure=distance_measure)

    def load_file(self, key: str, in_type: str, distance_measure: str = None):
        if in_type == "mapping":
            self._load_file_mapping(file=os.path.join(self.files_dir, self.file_names[key]), sep=",", mapping_name=key)
        elif in_type == "distance":
            self.loaded_distances[distance_measure][key] = sp.load_npz(
                os.path.join(self.files_dir, distance_measure, self.file_names[key])).tocsr()
        else:  # in_type == "distance_id"
            with open(os.path.join(self.files_dir, distance_measure, self.file_names[key]), 'rb') as f:
                self.loaded_distance_ids[distance_measure][key] = pickle.load(f)

    def save_mappings(self):
        for mapping_key in ['gene_ids', 'disorder_ids']:
            if not self.loaded_mappings[mapping_key].empty and mapping_key in self.changed_mappings:
                self.save_file(in_object=self.loaded_mappings[mapping_key], key=mapping_key, in_type='mapping')
        for mapping_key in ['gene_atts', 'disorder_atts']:
            if not self.loaded_mappings[mapping_key].empty and mapping_key in self.changed_mappings:
                df = self.loaded_mappings[mapping_key]
                df[df.columns[1:]] = df[df.columns[1:]].fillna('').applymap(mu.set_to_string)
                self.save_file(in_object=df, key=mapping_key, in_type='mapping')

    def save_distances(self):
        for distance_measure in ['jaccard', 'overlap']:
            os.system("mkdir -p " + os.path.join(self.files_dir, distance_measure))
            for distance_id_key in ['gene_mat_ids', 'disease_mat_ids']:
                if distance_measure + "_" + distance_id_key in self.changed_mappings:
                    self.save_file(in_object=self.loaded_distance_ids[distance_measure][distance_id_key],
                                   key=distance_id_key, distance_measure=distance_measure,
                                   in_type='distance_id')
            for distance_key in ['go_BP', 'go_CC', 'go_MF', 'pathway_kegg', 'related_genes', 'related_variants',
                                 'related_pathways']:
                if distance_measure + "_" + distance_key in self.changed_mappings:
                    self.save_file(in_object=self.loaded_distances[distance_measure][distance_key], key=distance_key,
                                   in_type='distance', distance_measure=distance_measure)

    def save_file(self, in_object, key: str, in_type: str, distance_measure: str = None):
        if in_type == "mapping":
            if not self.loaded_mappings[key].empty:
                in_object.to_csv(os.path.join(self.files_dir, self.file_names[key]), index=False)
        elif in_type == "distance":
            if self.loaded_distances[distance_measure][key].nnz > 0:
                sp.save_npz(os.path.join(self.files_dir, distance_measure, self.file_names[key]), in_object.tocoo())
        else:  # in_type == "distance_id"
            if self.loaded_distance_ids[distance_measure][key]:
                with open(os.path.join(self.files_dir, distance_measure, self.file_names[key]), 'wb+') as f:
                    pickle.dump(in_object, f)

    def check_for_setup_sources(self):
        for key in ['gene_atts', 'disorder_atts', 'gene_ids', 'disorder_ids']:
            if not Path(os.path.join(self.files_dir, self.file_names[key])).is_file():
                raise Exception(self.file_names[key] + "does not exist. Please run setup.")
        for distance_measure in ['jaccard', 'overlap']:
            for key in ['gene_mat_ids', 'disease_mat_ids', 'go_BP', 'go_CC', 'go_MF', 'pathway_kegg',
                        'related_genes', 'related_variants', 'related_pathways']:
                if not Path(os.path.join(self.files_dir, distance_measure, self.file_names[key])).is_file():
                    raise Exception(self.file_names[key] + "does not exist. Please run setup.")
