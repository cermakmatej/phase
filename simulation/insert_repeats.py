#!/usr/bin/env python3
import argparse
import random
import numpy.random
import gzip
import bisect
from collections import defaultdict

random.seed(42)
numpy.random.seed(42)

def read_fasta_one_seq(path: str) -> str:
    seq_lines = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"): continue
            seq_lines.append(line.strip().upper())
    return "".join(seq_lines)

def write_fasta(seq: str, name: str, path: str, width: int = 60):
    with open(path, "w") as out:
        out.write(f">{name}\n")
        for i in range(0, len(seq), width):
            out.write(seq[i:i+width] + "\n")


def parse_vcf_to_dict(file_path):
    variants_by_pos = defaultdict(list)
    meta_lines = []
    header_columns = []
    
    open_func = gzip.open if file_path.endswith('.gz') else open
    read_mode = 'rt' if file_path.endswith('.gz') else 'r'
    
    print(f"Parsing VCF file: {file_path}...")
    try:
        with open_func(file_path, read_mode, encoding='utf-8') as f:
            for line in f:
                if line.startswith('##'):
                    meta_lines.append(line.strip())
                    continue
                    
                if line.startswith('#CHROM'):
                    meta_lines.append(line.strip())
                    header_columns = line.strip().lstrip('#').split('\t')
                    continue
                    
                values = line.strip().split('\t')
                
                if not header_columns:
                    header_columns = [f"Column_{i}" for i in range(len(values))]
                
                variant_data = dict(zip(header_columns, values))
                chrom = variant_data.get('CHROM', 'unknown')
                pos = int(variant_data.get('POS', '0'))
                
                variants_by_pos[(chrom, pos)].append(variant_data)
                
        return dict(variants_by_pos), meta_lines, header_columns
    except Exception as e:
        print(f"[ERROR] An error occurred while loading VCF: {e}")
        return {}, [], []

def write_vcf(vcf_records, headers, header_columns, out_path, final_seq_len):
    with open(out_path, 'w') as f:
        for h in headers:
            if h.startswith("##contig=") or h.startswith("#CHROM"):
                continue
            f.write(h + "\n")
        
        f.write(f"##contig=<ID=ref_with_repeats,length={final_seq_len}>\n")
        
        f.write("#" + "\t".join(header_columns).lstrip("#") + "\n")
        
        vcf_records.sort(key=lambda x: int(x['POS']))
        
        for rec in vcf_records:
            row = [str(rec.get(col, '.')) for col in header_columns]
            f.write("\t".join(row) + "\n")

def get_variants_in_region(vcf_dict, target_chrom, start_pos, end_pos):
    region_variants = []
    target_chrom = str(target_chrom)
    
    for pos in range(start_pos, end_pos + 1):
        if (target_chrom, pos) in vcf_dict:
            distance_from_start = pos - start_pos
            for variant in vcf_dict[(target_chrom, pos)]:
                variant_with_dist = variant.copy()
                variant_with_dist['DISTANCE'] = distance_from_start
                region_variants.append(variant_with_dist)
                
    return sorted(region_variants, key=lambda x: int(x.get('DISTANCE', 0)))

def apply_variants_to_motif(motif_seq, vars_in_motif, prob_alt=0.5):
    mutated_motif = list(motif_seq)
    applied_vars = []
    
    for var in vars_in_motif:
        if random.random() < prob_alt:
            dist = var['DISTANCE']
            ref = var['REF']
            
            alt = var['ALT'].split(',')[0] 
            
            if len(ref) == 1 and len(alt) == 1:
                
                if 0 <= dist < len(mutated_motif) and mutated_motif[dist].upper() == ref.upper():

                    
                    applied_vars.append({
                        'orig_var': var,
                        'motif_offset': dist,
                        'alt_len': 1
                    })
                    
    return "".join(mutated_motif), applied_vars


def get_rep_seq_variants(mutated_motif, applied_vars, block_len):
    rep_seq_vars = []
    motif_len = max(len(mutated_motif), 1)
    num_repeats = (block_len // motif_len) + 1
    
    for i in range(num_repeats):
        block_offset = i * motif_len
        for var_info in applied_vars:
            final_offset = block_offset + var_info['motif_offset']
            
            if final_offset + var_info['alt_len'] <= block_len:
                rep_seq_vars.append({
                    'orig_var': var_info['orig_var'],
                    'rep_offset': final_offset
                })
                
    return rep_seq_vars


def get_motif_from_seq(seq: str, motif_len: int):
    if len(seq) < motif_len:
        return seq, 0, len(seq)
    start = random.randint(0, len(seq) - motif_len)
    return seq[start : start + motif_len], start, start + motif_len - 1

def make_repeat_block(motif: str, block_len: int) -> str:
    if not motif: return "A" * block_len
    r = (motif * ((block_len // len(motif)) + 1))[:block_len]
    return r

def insert_repeats(seq: str, target_percent, vcf_data):
    variants_dict, meta_lines, header_columns = vcf_data
    
    chromosomes_in_vcf = list(set(k[0] for k in variants_dict.keys()))
    target_chrom = chromosomes_in_vcf[0] if chromosomes_in_vcf else "22"
    chrom_name = target_chrom
    
    genome_len = len(seq)
    target = target_percent / 100
    target_total = int(genome_len * target / (1 - target))
    print(f"Target repetitive length: {target_total:,} bp ({target_percent:.2f}%)")

    seq_list = list(seq)
    inserted = 0
    repeats = []
    variants_found_total = 0

    while inserted < target_total:
        motif_len = max(numpy.random.poisson(10), 1)
        motif, start, end = get_motif_from_seq(seq, motif_len)
        
        vcf_start = start + 1
        vcf_end = end + 1
        
        vars_in_motif = get_variants_in_region(variants_dict, "22", vcf_start, vcf_end)
        if vars_in_motif:
            variants_found_total += len(vars_in_motif)


        if vars_in_motif:
            block_len = len(motif)
            insertions = random.randint(1, 100)
        else:
            block_len = random.choice([10, 30, 100, 300, 500]) * len(motif)
            insertions = random.randint(1, 50)
            passed_insertions = 0
        
        for i in range(insertions):
            if inserted >= target_total:
                break
                
            if vars_in_motif:
                mutated_motif, applied_vars = apply_variants_to_motif(motif, vars_in_motif, prob_alt=0.5)
            else:
                mutated_motif = motif
                applied_vars = []

            rep_seq = make_repeat_block(mutated_motif, block_len)
            rep_seq_vars = get_rep_seq_variants(mutated_motif, applied_vars, block_len)
                
            pos = random.randint(0, len(seq_list))
            
            idx = bisect.bisect_right(repeats, (pos, "", []))
            
            prekryv = False
            if idx > 0:
                prev_pos, prev_seq, _ = repeats[idx - 1]
                if prev_pos + len(prev_seq) > pos:
                    prekryv = True
                    
            if not prekryv and idx < len(repeats):
                next_pos, next_seq, _ = repeats[idx]
                if pos + len(rep_seq) > next_pos:
                    prekryv = True
                    
            if not prekryv:
                passed_insertions +=1

                repeats.insert(idx, (pos, rep_seq, rep_seq_vars))
                inserted += len(rep_seq)+motif_len
        print(f"Motif lenth: {motif_len}, Tandem repeats: {block_len/motif_len}, Insertions: {passed_insertions}")


    print("Assembling final sequence and VCF records...")
    final_seq_blocks = []
    final_vcf_records = []
    
    current_orig_pos = 0
    cumulative_inserted_len = 0
    new_vars = 0
    
    for pos, rep_seq, rep_seq_vars in repeats:
        final_seq_blocks.append("".join(seq_list[current_orig_pos:pos]))
        current_new_seq_pos = pos + cumulative_inserted_len
        
        for var_info in rep_seq_vars:
            final_vcf_pos = current_new_seq_pos + var_info['rep_offset'] + 1 
            new_record = var_info['orig_var'].copy()
            new_record['POS'] = str(final_vcf_pos)
            new_record['CHROM'] = "ref_with_repeats"
            final_vcf_records.append(new_record)
            new_vars +=1
            
        final_seq_blocks.append(rep_seq)
        cumulative_inserted_len += len(rep_seq)
        current_orig_pos = pos
        
    final_seq_blocks.append("".join(seq_list[current_orig_pos:]))
    final_sequence = "".join(final_seq_blocks)
    
    print("Shifting and appending original variants...")
    target_chrom_str = str(chrom_name)
    target_chrom_chr = f"chr{chrom_name}"
    
    ins_positions = [r[0] for r in repeats]
    cum_shifts = []
    current_shift = 0
    for r in repeats:
        current_shift += len(r[1])
        cum_shifts.append(current_shift)
        
    for (chrom, pos), var_list in variants_dict.items():
        if chrom in (target_chrom_str, target_chrom_chr):
            orig_pos_0based = pos - 1
            
            idx = bisect.bisect_right(ins_positions, orig_pos_0based)
            shift = cum_shifts[idx - 1] if idx > 0 else 0
            
            new_pos_1based = pos + shift
            
            for var in var_list:
                shifted_var = var.copy()
                shifted_var['POS'] = str(new_pos_1based)
                shifted_var['CHROM'] = "ref_with_repeats"
                final_vcf_records.append(shifted_var)
                
    print("Sorting final VCF records by POS...")
    final_vcf_records.sort(key=lambda x: int(x['POS']))
    
    print(f"Final VCF contains {len(final_vcf_records):,} records (original + newly generated {new_vars}).")
    return final_sequence, final_vcf_records


def main():
    p = argparse.ArgumentParser(description="Insert repeats into reference sequence.")
    p.add_argument("--ref", required=True, help="Input FASTA")
    p.add_argument("--vcf", required=True, help="Input VCF")
    p.add_argument("--target_percent", type=float, default=0.0, help="Target percent of genome length to be repetitive")
    p.add_argument("--out_fa", required=True, help="Output FASTA filename")
    p.add_argument("--out_vcf", required=True, help="Output VCF filename")

    args = p.parse_args()

    seq = read_fasta_one_seq(args.ref)
    print(f"Loaded reference length: {len(seq):,} bp")
    
    vcf_data = parse_vcf_to_dict(args.vcf)

    new_seq, new_vcf_records = insert_repeats(
        seq,
        target_percent=args.target_percent,
        vcf_data=vcf_data
    )

    write_fasta(new_seq, "ref_with_repeats", args.out_fa)
    print(f"Wrote new FASTA: {args.out_fa}")
    
    write_vcf(new_vcf_records, vcf_data[1], vcf_data[2], args.out_vcf, len(new_seq))
    print("Done")

if __name__ == "__main__":
    main()
