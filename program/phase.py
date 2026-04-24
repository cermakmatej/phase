import re
import sys
from tqdm import tqdm
from collections import defaultdict
import dynamic_programming as phaser
from graph_builder import build_variant_graph
import output_parser 
import argparse

# Final version of the main phase program
parser = argparse.ArgumentParser()
parser.add_argument("--variant_file", default=None, type=str, help="Variant file for VCF output.")

# Parse allele_id represented as a string to (chrom, pos, base) representation
def parse_allele_id(allele_id):
    match = re.match(r"(.+):(\d+)_([ACGTN]+)", allele_id)
    if match:
        chrom, pos, base = match.groups()
        return chrom, int(pos), base
    return None

# Find conected components in a graph
def connected_components(graph):
    visited = set()
    components = []

    for start in graph:
        if start not in visited:
            stack = [start]
            component = []

            while stack:
                u = stack.pop()
                if u not in visited:
                    visited.add(u)
                    component.append(u)
                    for v in graph[u]:
                        if v not in visited:
                            stack.append(v)

            components.append(component)

    return components

# Find inducted subgraph on a subset of vertices
def subgraph(graph, vertices):
    vertices = set(vertices)
    subgraph = defaultdict(dict)

    for u in vertices:
        if u in graph: 
            for v, w in graph[u].items():
                if v in vertices:  
                    subgraph[u][v] = w

    return subgraph

def main(args:argparse.Namespace):
    #input = "data/22.txt"

    input = "stdin"
    #input = "sim_input_art_30.txt"
    #variants = vcf_output_parser.load_vcf("sim_sample_22_art_30.vcf")

    # Load variants and build the graphs for each sample
    variants = output_parser.load_vcf(args.variant_file)
    graphs = build_variant_graph(input, variants)

    # For each sample
    for i in range(len(graphs)):
        print(f"Phasing sample {i}", file=sys.stderr)

        G = graphs[i]
        print(f"##SAMPLE: {i}")

        # Find conected components and sort them based on first variant
        components = connected_components(G)
        components = sorted(components, key = lambda x: int(parse_allele_id(sorted(list(subgraph(G,x).keys()), key = lambda y: int(parse_allele_id(y)[1]))[0])[1]))
        
        # For each component = block
        for j, component_nodes in tqdm(enumerate(components, 1)):
            
            # Get the inducted subgraph on the contained variants
            H = subgraph(G, component_nodes)
            
            # If block has < 4 alleles, pass
            if len(H.keys()) < 4:
                continue

            divergence, partition, masked = phaser.phase(H)

            print(f"#BLOCK: variants: {len(H.keys())//2} KL divergence: {divergence:.4f}")

            vcf = output_parser.dict_to_out(partition, variants, masked)            
            print(vcf)
        
        print(f"Completed phasing sample {i}", file=sys.stderr)

if __name__ == "__main__":
    main_args = parser.parse_args([] if "__file__" not in globals() else None)
    main(main_args)
