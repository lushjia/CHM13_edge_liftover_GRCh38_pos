"""
edge (chm13 graph) - GRCh38 positions

Usage 
python edge_b38pos.py -edge_var edge_var.txt -var_index var_index.txt -liftover lifted_over.vcf -o output_prefix

# output 
1) output_prefix_edge_b38pos.txt: edge - GRCh38 positions of associated variants
    format: edge    b38:pos1_varid,b38:pos2_varid,...
2) output_prefix_edge_b38pos_leftmost.txt: edge - left most GRCh38
    format: edge    var_id    chr    pos

"""

import argparse
from pysam import VariantFile

def load_edge_var(edge_var_file):
    edge_var_dict = {}
    with open(edge_var_file, 'r') as f:
        for line in f:
            line = line.strip()
            edge, var_indices = line.split('\t', 1)
            var_index_list = var_indices.split(',')
            edge_var_dict[edge] = var_index_list
    return edge_var_dict

def load_var_index(var_index_file):
    var_index_dict = {}
    with open(var_index_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                index, var_id = line.split('\t', 1)
                var_index_dict[index] = var_id
    return var_index_dict

def find_var_b38_pos(liftover_vcf, var_index):
    var_set = set(var_index.values())
    vcf_in = VariantFile(liftover_vcf)
    var_b38_pos = {}
    for variant in vcf_in.fetch():
        var_id = variant.id
        if var_id in var_set:
            chrom = variant.chrom
            pos = variant.pos
            var_b38_pos[var_id] = f"{chrom}:{pos}"
    return var_b38_pos

if __name__== "__main__":
    parser = argparse.ArgumentParser(prog='edge_b38pos.py', description=__doc__)
    parser.add_argument('-edge_var', metavar='txt', required=True, help='edge variant correlation file, e.g. >130839136>130839168 1,3,10')
    parser.add_argument('-var_index', metavar='txt', required=True, help='variant index file, e.g. 1 >162724237>162724240-1')
    parser.add_argument('-liftover', metavar='vcf', required=True, help='lifted over VCF file with GRCh38 positions.')
    parser.add_argument('-o', metavar='OUTPREFIX', required=True, help='Prefix of the output files.')
    args = parser.parse_args()

    edge_var_dict = load_edge_var(args.edge_var)
    var_index = load_var_index(args.var_index)
    var_b38_pos = find_var_b38_pos(args.liftover, var_index) # Get GRCh38 positions for variants; var: b38_pos
    # output edge with GRCh38 positions
    # 1. # edge    b38_pos1,b38_pos2,...
    # 2. # edge    b38_pos (left most/smallest)
    with open(f'{args.o}_edge_b38pos.txt', 'w') as out_f, open(f'{args.o}_edge_b38pos_leftmost.txt', 'w') as out_f_left: 
        for edge, var_indices in edge_var_dict.items():
            b38_pos_list = []
            for var_i in var_indices:
                var_id = var_index.get(var_i)
                if var_id and var_id in var_b38_pos:
                    b38_pos_list.append(var_b38_pos[var_id] + "_" + var_id) # chr:pos_varid chr1:297653_>125818493>125818520-1
            # if successfully found in GRCh38 by liftover
            if b38_pos_list:
                # write all position 
                out_f.write(f"{edge}\t{','.join(b38_pos_list)}\n") # 125818493 125818494 chr1:297653_>125818493>125818520-1
                # write left most position
                ## if only one pos 
                if len(b38_pos_list) == 1:
                    chr = b38_pos_list[0].rsplit("_",1)[0].split(":")[0]
                    pos = b38_pos_list[0].rsplit("_",1)[0].split(":")[1]
                    id = b38_pos_list[0].rsplit("_",1)[1]
                    out_f_left.write(f"{edge}\t{id}\t{chr}\t{pos}\n")
                else:
                    # find the most common chr and smallest pos
                    chr_count = {}
                    chr_pos_dict = {}
                    for b38_pos in b38_pos_list:
                        chr = b38_pos.rsplit("_",1)[0].split(":")[0]
                        pos = int(b38_pos.rsplit("_",1)[0].split(":")[1])
                        id = b38_pos.rsplit("_",1)[1]
                        if chr not in chr_count:
                            chr_count[chr] = 0
                            chr_pos_dict[chr] = []
                        chr_count[chr] += 1
                        chr_pos_dict[chr].append((pos, id))
                    # find chr with max count
                    max_chr = max(chr_count, key=chr_count.get) # chr1
                    # find smallest pos in that chr
                    pos_id_list = chr_pos_dict[max_chr]
                    pos_id_list.sort() # sort by pos
                    smallest_pos, id = pos_id_list[0]
                    out_f_left.write(f"{edge}\t{id}\t{max_chr}\t{smallest_pos}\n")
            else:
                # no position found, skip this edge
                continue








