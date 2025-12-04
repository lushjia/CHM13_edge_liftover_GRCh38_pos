"""
edge - graph var relationship

# usage: python edge_graphvar.py -edge edge.txt -var graphvar.vcf -o output_prefix

# output:
1) output_prefix.edge_var.txt: edge - graph var index relationship
    format: node1  node2  var_index1,var_index2,...
2) output_prefix.var_index.txt: graph var index - var id relationship
    format: var_index  var_id

"""

import argparse
import re 
from pysam import VariantFile
from collections import defaultdict

def parse_edge_file(edge_file):
    # record edges in a set
    edge_set = set()
    original_edge_format = {} # (node1, node2) : edge_str
    with open(edge_file, 'r') as ef:
        for line in ef:
            edge = line.strip() # >130839136>130839168
            nodes = re.split(">|<", edge)[1:] # ['130839136', '130839168']
            # smaller node first 
            if nodes[0] > nodes[1]:
                nodes[0], nodes[1] = nodes[1], nodes[0]
            # >130839136>130839168
            edge_set.add((nodes[0], nodes[1]))
            original_edge_format[(nodes[0], nodes[1])] = edge

    return edge_set, original_edge_format

def process_vcf(vcf_filename, edge_set):
    # record 
    pangenie_decompose_vcf = VariantFile(vcf_filename)
    var_index = 0
    var_dict = {} # index: var_id
    edge_var_dict = defaultdict(list) # edge: index ... 
    # output 
    # index: var_id 
    for variant in pangenie_decompose_vcf.fetch():
        var_edges = set()
        var_id = variant.id
        alt_paths = variant.info['AT'] # ref and alt path # >141461110>141461111>141461114,>141461110>141461113>141461114
        # get edges from alt paths
        # since biallelic file, only two path in AT field
        record = False
        for num_path in range(2):
            path = alt_paths[num_path]
            nodes = re.split(">|<", path)[1:]
            for ref_node_index in range(len(nodes)-1):
                node_pair = [nodes[ref_node_index], nodes[ref_node_index+1]]
                # smaller node first
                if node_pair[0] > node_pair[1]:
                    node_pair[0], node_pair[1] = node_pair[1], node_pair[0]
                var_edges.add((node_pair[0], node_pair[1]))
        # record edges that are in edge_set
        for node_pair in var_edges:
            if node_pair in edge_set:
                edge_var_dict[node_pair].append(var_index)
                record = True
        if record:
            var_dict[var_index] = var_id 
            var_index += 1
    return var_dict, edge_var_dict



if __name__== "__main__":
    parser = argparse.ArgumentParser(prog='edge_graphvar.py', description=__doc__)
    parser.add_argument('-edge', metavar='txt', required=True, help='edge list file, one edge per line, e.g. >130839136>130839168')
    parser.add_argument('-var', metavar='VCF', required=True, help='graph variant VCF file')
    parser.add_argument('-o', metavar='OUTPREFIX', required=True, help='Prefix of the output files.')
    args = parser.parse_args()
    
    edge_set, original_edge_format = parse_edge_file(args.edge) # set of (node1, node2) ... 
    var_dict, edge_var_dict = process_vcf(args.var, edge_set) # 

    # output files
    with open(f'{args.o}.edge_var.txt', 'w') as out_edge_var:
        for edge in edge_var_dict:
            var_indices = edge_var_dict[edge]
            var_indices_str = ','.join([str(idx) for idx in var_indices])
            edge_str = original_edge_format[edge]
            out_edge_var.write(f'{edge_str}\t{var_indices_str}\n')
    with open(f'{args.o}.var_index.txt','w') as out_var_index:
        for var_index in var_dict:
            var_id = var_dict[var_index]
            out_var_index.write(f'{var_index}\t{var_id}\n')
