#!/usr/bin/python3

from d_utils import runner_utils as ru, config
from mappers.mapper import Mapper, FileMapper
import graph_tool as gt

def create_network(node_type, edge_type, out_dir, mapper:Mapper=FileMapper()):
    ru.print_current_usage('Starting generation of network ...')
    ru.print_current_usage('Load distance file ...')
    mapper.load_distance(key=edge_type)
    if node_type in config.SUPPORTED_DISEASE_IDS:
        if node_type != 'mondo':
            ru.print_current_usage('Map mondo ID to '+node_type+' ...')
    else: # node_type in config.SUPPORTED_GENE_IDS:
        if node_type != 'entrez':
            ru.print_current_usage('Map entrez ID to ' + node_type + ' ...')
    ru.print_current_usage('Create edge list ...')
    edge_list = list()
    for i in range(0, len(mapper.loaded_distances[edge_type].columns)-1):
        for j in range(i+1, len(mapper.loaded_distances[edge_type].columns)):
            weight = mapper.loaded_distances[edge_type].loc[mapper.loaded_distances[edge_type].columns[i],
                                                            mapper.loaded_distances[edge_type].columns[j]]
            if weight > 0.0:
                edge_list.append((mapper.loaded_distances[edge_type].columns[i],
                                  mapper.loaded_distances[edge_type].columns[j],
                                  weight))
    ru.print_current_usage('Create graph ...')
    g = gt.Graph(directed=False)
    edge_weights = g.new_edge_property('double')
    node_names = g.add_edge_list(edge_list, hashed=True, eprops=[edge_weights])
    g.edge_properties["weight"] = edge_weights
    g.vertex_properties["ID"] = node_names
    ru.print_current_usage('Save data with '+g.num_vertices()+' vertices ...')
    g.save(out_dir+node_type+"_"+edge_type+"_network.graphml", fmt="graphml")
    ru.print_current_usage('Finished ...')

if __name__ == "__main__":
    desc = "            Transform precalculated distances to networks."
    args = ru.save_parameters(script_desc=desc, arguments=('n', 'd', 'o'))
    create_network(node_type=args.node_type, edge_type=args.distance_type, out_dir=args.out_dir)