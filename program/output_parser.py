from collections import defaultdict

def parse_key(key: str):
    try:
        chrom, pos_tag = key.split(":")
        pos_str, base = pos_tag.split("_", 1)
        pos = int(pos_str)
    except Exception as e:
        raise ValueError(f"Invalid key '{key}': {e}")
    return chrom, pos, base

def load_vcf(guide_vcf: str):
    variants = defaultdict(dict)
    with open(guide_vcf) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            chrom, pos, vid, ref, alt, qual, flt, info, fmt, *rest = line.strip().split("\t")
            pos = int(pos)
            variants[(chrom, pos)]["ref"]=ref

            variants[(chrom, pos)]["alt"] = alt.split(",") if alt != "." else []
    return variants


def dict_to_out(haplotypes: tuple, variants: dict, masked) -> str:
    lines = []
    sites = defaultdict(dict)

    for haplotype in range(len(haplotypes)):
        group = haplotypes[haplotype]
        for allele in group:
            chrom, pos, base = parse_key(allele)
            sites[(chrom, pos)][haplotype] = base
    
    for (chrom, pos) in sorted(sites.keys(), key=lambda x: (x[0], x[1])):
        
        site = sites[(chrom, pos)]

        ref = variants[(chrom, pos)]["ref"]

        alts = variants[(chrom, pos)]["alt"]

        try:
            alleles = [ref]+alts

            gt = str(alleles.index(site[0]))+"|"+str(alleles.index(site[1]))

            if len(alts)>1 and gt.find('0')!=-1:
                print(f"#Alleles in predicted haplotype inconsistent with supplied vcf on variant {chrom},{pos}")
                continue
        except:
            print(f"#Alleles in predicted haplotype inconsistent with supplied vcf on variant {chrom},{pos}")
            continue
            
        alts = ",".join(alts)
        
                
        if (chrom, pos) in masked:
                lines.append(f"{chrom}\t{pos}\t{ref}\t{alts}\tUNPHASED\t{gt}")
        else:
            lines.append(f"{chrom}\t{pos}\t{ref}\t{alts}\t.\t{gt}")
        
        

    return "\n".join(lines)