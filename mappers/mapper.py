#!/usr/bin/python3

import pandas as pd
from abc import abstractmethod
from d_utils import config, mapping_utils as mu
from pathlib import Path


class Mapper:
    loaded_mappings = {'gene_ids': pd.DataFrame(), 'gene_atts': pd.DataFrame(), 'disorder_ids': pd.DataFrame(),
                       'disorder_atts': pd.DataFrame()}
    loaded_distances = {'go_BP': pd.DataFrame(), 'go_CC': pd.DataFrame(), 'go_MF': pd.DataFrame(),
                        'pathway_kegg': pd.DataFrame(), 'related_genes': pd.DataFrame(),
                        'related_variants': pd.DataFrame(), 'related_pathways': pd.DataFrame()}

    @abstractmethod
    def load_mappings(self, set_type: str):
        pass

    @abstractmethod
    def load_distance(self, key: str):
        pass

    @abstractmethod
    def load_distances(self, set_type: str):
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

    def update_distances(self, in_df: pd.DataFrame(), key: str):
        if not self.loaded_distances[key].empty:
            new_columns = list(in_df.columns.difference(self.loaded_distances[key].columns))
            self.loaded_distances[key] = self.loaded_distances[key].join(in_df)
            # todo
        else:
            self.loaded_distances[key] = in_df

    def get_loaded_distances(self, in_set, key: str):
        if not self.loaded_mappings[key].empty:
            hit_mapping = self.loaded_mappings[key].loc[self.loaded_mappings[key].index.isin(in_set),
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
    def save_directly(self, in_df: pd.DataFrame(), key: str):
        pass

    def drop_mappings(self):
        self.loaded_mappings = {'gene_ids': pd.DataFrame(), 'gene_atts': pd.DataFrame(), 'disorder_ids': pd.DataFrame(),
                                'disorder_atts': pd.DataFrame()}


class FileMapper(Mapper):
    file_names = {'gene_ids': config.FILES_DIR + 'gene_id_mapping.csv',
                  'gene_atts': config.FILES_DIR + 'gene_att_mapping.csv',
                  'disorder_ids': config.FILES_DIR + 'disease_id_mapping.csv',
                  'disorder_atts': config.FILES_DIR + 'disease_att_mapping.csv',
                  'go_BP': config.FILES_DIR + 'gene_dist_go_BP.csv',
                  'go_CC': config.FILES_DIR + 'gene_dist_go_CC.csv',
                  'go_MF': config.FILES_DIR + 'gene_dist_go_MF.csv',
                  'pathway_kegg': config.FILES_DIR + 'gene_dist_pathway_kegg.csv',
                  'related_genes': config.FILES_DIR + 'disease_dist_rel_genes.csv',
                  'related_variants': config.FILES_DIR + 'disease_dist_rel_variants.csv',
                  'related_pathways': config.FILES_DIR + 'disease_dist_rel_pathways.csv'}

    def load_mappings(self, set_type):
        if set_type in config.SUPPORTED_GENE_IDS:
            self._load_file_mapping(file=self.file_names['gene_ids'], sep=",", mapping_name='gene_ids')
            self._load_file_mapping(file=self.file_names['gene_atts'], sep=",", mapping_name='gene_atts')
        else:  # if set_type in config.SUPPORTED_DISEASE_IDS
            self._load_file_mapping(file=self.file_names['disorder_ids'], sep=",", mapping_name='disorder_ids')
            self._load_file_mapping(file=self.file_names['disorder_atts'], sep=",", mapping_name='disorder_atts')

    def load_distance(self, key: str):
        self._load_file_distances(file=self.file_names[key], sep=",", mapping_name=key)

    def load_distances(self, set_type):
        if set_type in config.SUPPORTED_GENE_IDS:
            self._load_file_distances(file=self.file_names['go_BP'], sep=",", mapping_name='go_BP')
            self._load_file_distances(file=self.file_names['go_CC'], sep=",", mapping_name='go-CC')
            self._load_file_distances(file=self.file_names['go_MF'], sep=",", mapping_name='go_MF')
            self._load_file_distances(file=self.file_names['pathway_kegg'], sep=",", mapping_name='pathway_kegg')
        else:  # if set_type in config.SUPPORTED_DISEASE_IDS
            self._load_file_distances(file=self.file_names['related_genes'], sep=",", mapping_name='related_genes')
            self._load_file_distances(file=self.file_names['related_variants'], sep=",",
                                      mapping_name='related_variants')
            self._load_file_distances(file=self.file_names['related_pathways'], sep=",",
                                      mapping_name='related_pathways')

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
        if not self.loaded_mappings[mapping_name].empty:
            self.loaded_mappings[mapping_name] = pd.concat(
                [self.loaded_mappings[mapping_name], mapping]).drop_duplicates()
        else:
            self.loaded_mappings[mapping_name] = mapping

    def _load_file_distances(self, file, sep, mapping_name):
        """
        Get previous distances from file.

        :param file: file with previous distances
        :param sep: separator of values in file
        :param mapping_name: key name of local dictionary for saving data
        :return: dataframe with previous mapping
        """
        # ===== Get distances from local distance file =====
        cols = pd.read_csv(file, sep=sep, nrows=1).columns.tolist()
        distances = pd.read_csv(file, sep=sep, header=0, dtype={cols[0]: str})
        distances.set_index(cols[0], inplace=True)
        # ===== Save mapping to local dictionary =====
        if not self.loaded_distances[mapping_name].empty:
            print("todo")
        else:
            self.loaded_distances[mapping_name] = distances

    def save_mappings(self):
        self._save_file_mapping(mapping_name='gene_ids')
        self._save_file_mapping(mapping_name='gene_atts')
        self._save_file_mapping(mapping_name='disorder_ids')
        self._save_file_mapping(mapping_name='disorder_atts')

    def _save_file_mapping(self, mapping_name):
        if not self.loaded_mappings[mapping_name].empty:
            if Path(self.file_names[mapping_name]).is_file():
                mapping = pd.read_csv(self.file_names[mapping_name], sep=",", header=0, dtype=str)
                mapping = pd.concat([self.loaded_mappings[mapping_name], mapping]).drop_duplicates()
                mapping.to_csv(self.file_names[mapping_name], index=False)
            else:
                self.loaded_mappings[mapping_name].to_csv(self.file_names[mapping_name], index=False)

    def save_distances(self):
        self._save_file_distances(mapping_name='go_BP')
        self._save_file_distances(mapping_name='go_CC')
        self._save_file_distances(mapping_name='go_MF')
        self._save_file_distances(mapping_name='pathway_kegg')
        self._save_file_distances(mapping_name='related_genes')
        self._save_file_distances(mapping_name='related_variants')
        self._save_file_distances(mapping_name='related_pathways')

    def _save_file_distances(self, mapping_name):
        self.loaded_distances[mapping_name].to_csv(self.file_names[mapping_name], index=True)

    def save_directly(self, in_df: pd.DataFrame(), key: str):
        in_df.to_csv(self.file_names[key], index=True)
