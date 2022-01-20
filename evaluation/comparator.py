#!/usr/bin/python3

from d_utils import eval_utils as eu, config as c
from mappers import gene_mapper as gm, disease_mapper as dm
from mappers.mapper import Mapper
from evaluation import score_calculator as sc
from abc import abstractmethod


class Comparator:
    def __init__(self, mapper: Mapper, verbose: bool = False):
        self.mapper = mapper
        self.verbose = verbose
        self.mapping = None
        self.id_set = None
        self.id_type = None
        self.att_id, self.att_key, self.sparse_key = None, None, None

    def load_target(self, id_set, id_type):
        self.id_set = id_set
        self.id_type = id_type
        if id_type in c.SUPPORTED_DISEASE_IDS:
            self.mapping = dm.get_disease_to_attributes(disease_set=id_set, id_type=id_type, mapper=self.mapper)
            self.sparse_key, self.att_key, self.att_id = 'disease_mat_ids', 'disorder_atts', 'mondo'
        else:  # if id_type in c.SUPPORTED_GENE_IDS:
            self.mapping = gm.get_gene_to_attributes(gene_set=id_set, id_type=id_type, mapper=self.mapper)
            self.sparse_key, self.att_key, self.att_id = 'gene_mat_ids', 'gene_atts', 'entrezgene'

    @abstractmethod
    def compare(self, threshold: float = 0.0):
        pass


class SetComparator(Comparator):
    """
    Compare the set on itself based on connected attributes. See config for more info.
    """

    def compare(self, threshold: float = 0.0):
        result = dict()
        new_ids = self.mapper.update_distance_ids(in_series=self.mapper.loaded_mappings[self.att_key][self.att_id],
                                                  key=self.sparse_key)
        for attribute in self.mapping.columns[1:]:
            subset_df = self.mapping[self.mapping[attribute].str.len() > 0]
            missing_values = len(self.mapping) - len(subset_df)
            if missing_values > 0:
                print("Missing values for " + attribute + " :" + str(missing_values) + "/" + str(
                    len(self.id_set))) if self.verbose else None
            if len(new_ids) > 0:
                comp_mat = eu.get_distance_matrix(full_att_series=self.mapper.loaded_mappings[self.att_key][attribute],
                                                  from_ids=self.mapper.loaded_mappings[self.att_key][self.att_id],
                                                  id_to_index=self.mapper.loaded_distance_ids[self.sparse_key],
                                                  to_ids=new_ids)
                self.mapper.update_distances(in_mat=comp_mat, key=c.DISTANCES[attribute], id_type=self.sparse_key)
            ids = self.mapper.get_loaded_mapping_ids(in_ids=set(subset_df[subset_df.columns[0]]),
                                                     id_type=self.id_type)
            sub_mat = self.mapper.get_loaded_distances(in_series=ids[self.att_id], id_type=self.sparse_key,
                                                       key=c.DISTANCES[attribute])
            result[attribute] = sub_mat.sum() / ((len(self.mapping) * (len(self.mapping) - 1)) / 2)
        return result


class SetSetComparator(Comparator):
    """
    Compare two sets of the same type (either genes or diseases) with each other.
    The tar set is evaluated how good it matches to the ref set.
    """

    def __init__(self, mapper: Mapper, enriched: bool = False, verbose: bool = False):
        super().__init__(mapper=mapper, verbose=verbose)
        self.enriched = enriched
        self.ref_dict = None

    def load_reference(self, ref, ref_id_type):
        if ref_id_type in c.SUPPORTED_DISEASE_IDS:
            reference_mapping = dm.get_disease_to_attributes(disease_set=ref, id_type=ref_id_type, mapper=self.mapper)
        else:  # if ref_id_type in c.SUPPORTED_GENE_IDS:
            if self.enriched:
                reference_mapping = gm.get_enriched_attributes(gene_set=ref, id_type=ref_id_type, mapper=self.mapper)
            else:
                reference_mapping = gm.get_gene_to_attributes(gene_set=ref, id_type=ref_id_type, mapper=self.mapper)
        if self.enriched:
            self.ref_dict = eu.create_ref_dict(mapping=reference_mapping, keys=c.ENRICH_KEY.keys(), enriched=True)
        else:
            self.ref_dict = eu.create_ref_dict(mapping=reference_mapping, keys=reference_mapping.columns[1:],
                                               enriched=False)

    def compare(self, threshold: float = 0.0):
        return eu.evaluate_values(mapping=self.mapping, ref_dict=self.ref_dict, threshold=threshold,
                                  keys=self.mapping.columns[1:])


class IDSetComparator(Comparator):
    """
    Compare disease id to either gene or disease set.
    The tar set is evaluated how good it matches to the ref id.
    """

    def __init__(self, mapper: Mapper, verbose: bool = False):
        super().__init__(mapper=mapper, verbose=verbose)
        self.ref_dict = None

    def load_reference(self, ref, ref_id_type):
        id_mapping = dm.get_disease_to_attributes(disease_set={ref}, id_type=ref_id_type, mapper=self.mapper)
        if self.id_type in c.SUPPORTED_DISEASE_IDS:
            self.ref_dict = eu.create_ref_dict(mapping=id_mapping, keys=id_mapping.columns[1:])
        else:  # if targets_id_type in c.SUPPORTED_GENE_IDS:
            id_mapping = id_mapping.rename(columns={'ctd.pathway_related_to_disease': 'pathway.kegg'})
            self.ref_dict = eu.create_ref_dict(mapping=id_mapping, keys={'pathway.kegg'})

    def compare(self, threshold: float = 0.0):
        return eu.evaluate_values(mapping=self.mapping, ref_dict=self.ref_dict, threshold=threshold,
                                  keys={'pathway.kegg'})


class ClusterComparator(Comparator):
    """
    Evaluate the quality of clustering of given set with diseases or genes and
    assigned clusters. Additionally calculate statistical values with
    silhouette score and dunn index.
    """

    def __init__(self, mapper: Mapper, verbose: bool = False):
        super().__init__(mapper, verbose)
        self.clustering = None

    def load_target(self, id_set, id_type):
        id_set['cluster_index'] = id_set.groupby(1).ngroup()
        super().load_target(id_set=id_set[0], id_type=id_type)
        self.clustering = id_set[[0, 'cluster_index']]

    def compare(self, threshold: float = 0.0):
        result_di, result_ss = dict(), dict()
        new_ids = self.mapper.update_distance_ids(in_series=self.mapper.loaded_mappings[self.att_key][self.att_id],
                                                  key=self.sparse_key)
        for attribute in self.mapping.columns[1:]:
            subset_df = self.mapping[self.mapping[attribute].str.len() > 0]
            subset_clusters = self.clustering[self.clustering[0].isin(subset_df[self.id_type])]
            missing_values = len(self.mapping) - len(subset_df)
            if missing_values > 0:
                print("Missing values for " + attribute + " :" + str(missing_values) + "/" + str(
                    len(self.mapping))) if self.verbose else None

            if len(new_ids) > 0:
                comp_mat = eu.get_distance_matrix(full_att_series=self.mapper.loaded_mappings[self.att_key][attribute],
                                                  from_ids=self.mapper.loaded_mappings[self.att_key][self.att_id],
                                                  id_to_index=self.mapper.loaded_distance_ids[self.sparse_key],
                                                  to_ids=new_ids)
                self.mapper.update_distances(in_mat=comp_mat, key=c.DISTANCES[attribute], id_type=self.sparse_key)

            ids = self.mapper.get_loaded_mapping_ids(in_ids=set(subset_df[subset_df.columns[0]]),
                                                     id_type=self.id_type)
            distances = self.mapper.get_loaded_distances(in_series=ids[self.att_id], id_type=self.sparse_key,
                                                         key=c.DISTANCES[attribute])
            distances = dict(distances.todok().items())

            ss_score = sc.silhouette_score(ids_cluster=subset_clusters, ids_mapping=ids, distances=distances,
                                           mapper=self.mapper,
                                           ids={'id_type': c.ID_TYPE_KEY[self.id_type], 'sparse_key': self.sparse_key,
                                                'att_id': self.att_id, 'attribute': c.DISTANCES[attribute]})
            di_score = sc.dunn_index(ids_cluster=subset_clusters, ids_mapping=ids, distances=distances,
                                     mapper=self.mapper,
                                     ids={'id_type': c.ID_TYPE_KEY[self.id_type], 'sparse_key': self.sparse_key,
                                          'att_id': self.att_id, 'attribute': c.DISTANCES[attribute]})
            result_di[attribute] = di_score
            result_ss[attribute] = ss_score[0]  # ss_score[1] all intermediate results
        return result_di, result_ss
