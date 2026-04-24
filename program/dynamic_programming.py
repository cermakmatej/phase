import re
from collections import defaultdict
import numpy as np


def parse_allele_id(allele_id):
    match = re.match(r"(.+):(\d+)_([ACGTN]+)", allele_id)
    if match:
        chrom, pos, base= match.groups()
        return chrom, int(pos), base
    return None

# Group alleles represented as vertices of a graph into variants
def group_alleles_by_variant(G: defaultdict):
    variant_groups = defaultdict(list)
    for node in G.keys():
        parsed = parse_allele_id(node)
        if parsed:
            chrom, pos, _ = parsed
            key = (chrom, pos)
            variant_groups[key].append(node)
    return variant_groups

# Compute total weight of a given group of alleles in a graph
def compute_group_weight(G: defaultdict, group: list):
    sum = 0
    for u in group:
        for v in group:
            if v in G[u].keys() and u<v:
                if G[u][v]!=0:
                    #sum += (G[u][v]) ** 0.5
                    sum += G[u][v]
    return sum

# Compute KL divergence of two distributions
def KL(a, b):
    if np.any((a > 0) & (b == 0)):
        return np.inf
    mask = (a > 0) & (b > 0)
    return np.sum(a[mask] * np.log(a[mask] / b[mask]))


# Find mask as the difference of the two strings with minimum Hamming distance
def binary_diff(a: int, b: int, c: int, l) -> list[int]:
    def bit_diff(x: int, y: int) -> list[int]:
        
        return [i for i in range(l) if ((x >> i) & 1) != ((y >> i) & 1)]

    pairs = [
        bit_diff(a, b),
        bit_diff(a, c),
        bit_diff(b, c),
    ]

    return min(pairs, key=len)


# Mask phased variants based on KL divergence
def mask(weights, S, l, epsilon = 0.05): 

    # inicialize ideal distribution
    expected_distr = np.zeros(S)
    expected_distr[0] = 0.5 - epsilon
    expected_distr[1] = 0.5 - epsilon
    for i in range(2, S):
        expected_distr[i]=(2*epsilon)/(S-2)

    masked = []

    for start_pos, row in enumerate(weights):
        sorted_row = -np.sort(-np.asarray(row))
        
        sorted_distr = sorted_row/np.sum(row)
        real_distr = row/np.sum(row)
        
        KLdiv = KL(sorted_distr, expected_distr)

        while KLdiv > 2.4 and len(masked)<l-1:
            first_indices = np.where(row==sorted_row[0])[0]
            second_indices = np.where(row==sorted_row[1])[0]
            third_indices = np.where(row==sorted_row[2])[0]
            indices = np.concatenate((first_indices, second_indices, third_indices))
            first = (indices[0])
            second = (indices[1])
            third = (indices[2])
            if third == first:
                third = indices[4]
            mask = binary_diff(first, second, third, l)


            for m in mask:
                bit = 1 << m
                for i in range(len(real_distr)):
                    if i & bit:               
                        j = i & ~bit         
                        real_distr[j] += real_distr[i]
                        real_distr[i] = 0

            row = real_distr
            sorted_row = -np.sort(-np.asarray(row)) 
            KLdiv = KL(sorted_row, expected_distr)

            for m in mask:
                masked.append([start_pos, m])
    return masked, KLdiv


def dynamic_programming(G: defaultdict, variants):
    K = 8

    vertices = list(G.keys())
    edges = []
    for u in G.keys():
        for v in G[u].keys():
            if u<v:
                edges.append((u,v,G[u][v]))
    edges.sort(key=lambda e: e[2], reverse=True)
    variant_keys = list(variants.keys())
    variant_keys.sort(key= lambda pair: pair[1])
    
    window_size = min(K, len(variant_keys))


    sub_haplos = [format(i, f'0{window_size}b') for i in range(2**window_size)]
    
    weights = []
    for i in range(len(variant_keys)-window_size+1):
        row = []
        for h in sub_haplos:
            group = []
            zero = False
            for k in range(window_size):
                allele_options = variants[variant_keys[i + k]]
                a = h[k]
                try:
                    group.append(allele_options[int(a)])
                except:
                    zero = True
                    break
            if not zero:
                row.append(compute_group_weight(G, group))
            else:

                row.append(0)

        weights.append(row)

    S = len(sub_haplos)
    n_windows = len(weights)


    dp = [[0] * S for _ in range(n_windows)]
    parent = [[None] * S for _ in range(n_windows)]

    for s in range(S):
        dp[0][s] = weights[0][s]

    # transitions
    for i in range(1, n_windows):
        for s_new in range(S):
            new_h = sub_haplos[s_new]
            for s_old in range(S):
                old_h = sub_haplos[s_old]

                # window compatibility:
                # old window:    [a b c d]
                # new window:      [b c d e]
                if old_h[1:] != new_h[:-1]:
                    continue
                
                complement = ""
                for char in new_h:
                    if char == "0":
                        complement +="1"
                    else:
                        complement +="0"
                complement=sub_haplos.index(complement)

                cand = dp[i-1][s_old] + weights[i][s_new] + weights[i][complement]

                if cand > dp[i][s_new]:
                    dp[i][s_new] = cand
                    parent[i][s_new] = s_old

    # backtrack optimal haplotype from dynamic programming matrix
    def backtrack_haplotype(state):
        states = [state]
        for i in range(n_windows-1, 0, -1):
            states.append(parent[i][states[-1]])
        
        states = states[::-1]

        haplotype_bits = sub_haplos[states[0]]
        for st in states[1:]:
            haplotype_bits += sub_haplos[st][-1]

        result = []
        for i, bit in enumerate(haplotype_bits):
            result.append(variants[variant_keys[i]][int(bit)])
        return result


    masks, KLdiv = mask(weights, S, window_size)

    masked_vars = set()
    for start_pos, pos_in_window in masks:
        masked_vars.add(variant_keys[start_pos + window_size-pos_in_window-1])


    last_state = max(range(S), key=lambda s: dp[-1][s])

    result = backtrack_haplotype(last_state)

    A = result


    B = []
    for allele in vertices:
        if allele not in A:
            B.append(allele)

    partition = (A, B)
    return KLdiv, partition, masked_vars



def phase(G):
    variants = group_alleles_by_variant(G)
    return dynamic_programming(G, variants)

