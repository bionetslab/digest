#!/usr/bin/python3

import pandas as pd
from abc import abstractmethod
from d_utils import config, mapping_utils as mu
from pathlib import Path


class Mapper:
    loaded_mappings = {'gene_ids': pd.DataFrame(), 'gene_atts': pd.DataFrame(), 'disorder_ids': pd.DataFrame(),
                       'disorder_atts': pd.DataFrame()}

    @abstractmethod
    def load_mappings(self, in_set: set, set_type: str):
        pass

    def update_mappings(self, in_df: pd.DataFrame(), key: str):
        if not self.loaded_mappings[key].empty:
            self.loaded_mappings[key] = pd.concat(
                [self.loaded_mappings[key], in_df]).drop_duplicates()
        else:
            self.loaded_mappings[key] = in_df

    def get_loaded_mapping(self, in_set, id_type, key: str):
        if not self.loaded_mappings[key].empty:
            hit_mapping = self.loaded_mappings[key][self.loaded_mappings[key][id_type].isin(in_set)]
            if not hit_mapping.empty:
                return hit_mapping, set(in_set) - set(hit_mapping[id_type])
        return pd.DataFrame(), in_set

    @abstractmethod
    def save_mappings(self):
        pass

    def drop_mappings(self):
        self.loaded_mappings = {'gene_ids': pd.DataFrame(), 'gene_atts': pd.DataFrame(), 'disorder_ids': pd.DataFrame(),
                                'disorder_atts': pd.DataFrame()}


class FileMapper(Mapper):
    file_names = {'gene_ids': config.FILES_DIR + 'gene_id_mapping.csv',
                  'gene_atts': config.FILES_DIR + 'gene_att_mapping.csv',
                  'disorder_ids': config.FILES_DIR + 'disorders.map',
                  'disorder_atts': config.FILES_DIR + 'disease_disgenet_mapping.csv'}

    def load_mappings(self, in_set, set_type):
        if set_type in config.SUPPORTED_GENE_IDS:
            hit_mapping = self._load_file_mapping(in_set=in_set, id_type=set_type, file=self.file_names['gene_ids'],
                                                  sep=",", mapping_name='gene_ids')
            if hit_mapping is not None:
                self._load_file_mapping(in_set=set(hit_mapping['entrezgene']), id_type='entrez',
                                        file=self.file_names['gene_atts'], sep=",", mapping_name='gene_atts')

        else:  # if set_type in config.SUPPORTED_DISEASE_IDS
            hit_mapping = self._load_file_mapping(in_set=in_set, id_type=set_type, file=self.file_names['disorder_ids'],
                                                  sep="\t", mapping_name='disorder_ids')
            self._load_file_mapping(in_set=list(set('MONDO:' + hit_mapping['mondo'])), id_type='mondo',
                                    file=self.file_names['disorder_atts'], sep=",", mapping_name='disorder_atts')

    def _load_file_mapping(self, in_set, id_type, file, sep, mapping_name):
        """
        Get previous mapping from file sub setting it to the
        element in in_set. Also returning the list of elements
        from the set that are not present in the previous mapping.

        :param in_set: set of genes or diseases
        :param id_type: id type of set
        :param file: file with previous mapping
        :param sep: separator of values in file
        :param mapping_name: key name of local dictionary for saving data
        :return: dataframe with previous mapping of input set or subset;
        """
        # ===== Get mapping from local mapping file =====
        mapping = pd.read_csv(file, sep=sep, header=0, dtype=str)
        if id_type == "ICD-10":
            mapping = mu.split_and_expand_column(data=mapping, split_string=",", column_name="ICD-10")
        id_type = config.ID_TYPE_KEY[id_type] if id_type in config.ID_TYPE_KEY else id_type
        # ===== Get mapping subset matching input set =====
        print(mapping)
        hit_mapping = mapping[mapping[id_type].isin(in_set)].copy()
        # ===== Save mapping to local dictionary =====
        if len(hit_mapping.index) > 0:
            if self.loaded_mappings[mapping_name] is not None:
                self.loaded_mappings[mapping_name] = pd.concat(
                    [self.loaded_mappings[mapping_name], hit_mapping]).drop_duplicates()
            else:
                self.loaded_mappings[mapping_name] = hit_mapping
            return hit_mapping
        else:
            return None

    def save_mappings(self):
        self._save_file_mapping(mapping_name='gene_ids')
        self._save_file_mapping(mapping_name='gene_atts')
        self._save_file_mapping(mapping_name='disorder_atts')

    def _save_file_mapping(self, mapping_name):
        if not self.loaded_mappings[mapping_name].empty:
            if Path(self.file_names[mapping_name]).is_file():
                mapping = pd.read_csv(self.file_names[mapping_name], sep=",", header=0, dtype=str)
                mapping = pd.concat([self.loaded_mappings[mapping_name], mapping]).drop_duplicates()
                mapping.to_csv(self.file_names[mapping_name], index=False)
            else:
                self.loaded_mappings[mapping_name].to_csv(self.file_names[mapping_name], index=False)
