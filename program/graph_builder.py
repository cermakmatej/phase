from tqdm import tqdm
from collections import defaultdict
import math
import re
import sys
import input_parser

# Calculate weight of allele pair on single read
def calculate_weight(map_qual, base_qual_x, base_qual_y, epsilon = 0.01):
    p_wrong_mapping = math.pow(10, - map_qual/10)
    p_wrong_x = math.pow(10, -base_qual_x/10)
    p_wrong_y = math.pow(10, -base_qual_y/10)
    return (1-p_wrong_mapping)*(1-p_wrong_x)*(1-p_wrong_y) + epsilon

# Add edge of given weight to graph
def add_edge(graph, u, v, w):
        graph[u][v] = w
        graph[v][u] = w

# Build variant graph for every sample in input
def build_variant_graph(input_file, variants):
    if input_file == "stdin":
        sample_entries, sample_allele_ids = input_parser.parse_input(sys.stdin, variants)
    else:
        with open(input_file) as f:
            sample_entries, sample_allele_ids = input_parser.parse_input(f, variants) 

    graphs = []


    for sample in range(len(sample_entries)):
        print(f"Building variant graph for sample {sample}", file=sys.stderr)

        entries = sample_entries[sample]
        allele_ids = sample_allele_ids[sample]
        
        G = defaultdict(dict)  # dict[vertex] -> dict[neighbor] = weight

        for read in tqdm(entries.keys()):

            alleles_on_read, conflict_alleles, base_quals, map_qual  = input_parser.get_alleles_on_read(read, entries, allele_ids)              

            # Calculate weight adjustments for each pair of alleles on the read
            for a1 in alleles_on_read:
                for a2 in alleles_on_read:
                    if a1 < a2 and a1 in allele_ids[parse_allele_id(a1)] and a2 in allele_ids[parse_allele_id(a2)]:
                        w = calculate_weight(map_qual, base_quals[a1], base_quals[a2])
                        old = G[a1].get(a2, 0.0)
                        add_edge(G, a1, a2, old + w)
                        
                        
        # Add edges between alternative alleles
        for var, alleles in allele_ids.items():
            if len(alleles) == 2:
                a1,a2 = list(alleles)
                add_edge(G, a1, a2, 0)      
            # Remove alleles without alternative alleles in the graph
            if len(alleles) == 1:
                loner = alleles.pop()
                for neighbor in G[loner].keys():
                    G[neighbor].pop(loner)
                G.pop(loner)                  

        graphs.append(G)
    return graphs

# Parse id of the variant of given allele
def parse_allele_id(allele_id):
    match = re.match(r"(.+):(\d+)_([ACGTN]+)", allele_id)
    if match:
        chrom, pos, base = match.groups()
        return chrom, pos
    return None
