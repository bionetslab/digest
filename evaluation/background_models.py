import pandas as pd
from .mappers import mapping_utils as mu
from .mappers.mapper import Mapper
from . import config as c
from abc import abstractmethod
import random


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