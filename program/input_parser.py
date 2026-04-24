from collections import defaultdict

# Parse input in format  (chr    pos     ismpl   base    BQ      MQ      read)
def parse_input(input, valid_variants):
    entries_on_reads = defaultdict(lambda: defaultdict(set))
    variants = defaultdict(lambda: defaultdict(set))
    
    for line in input:
        if not line.startswith("#"):
            chr, pos, ismpl, base, bq, mq, readname = line.strip().split()
            
            # convert to 1-based
            pos_int = int(pos) + 1
            pos_str = str(pos_int)
            ismpl = int(ismpl)
            
            if (chr, pos_int) in valid_variants:
                vcf_site = valid_variants[(chr, pos_int)]

                
                if len(vcf_site["alt"])>1:
                    valid_bases = vcf_site["alt"]
                else:
                    valid_bases = [vcf_site["ref"]] + vcf_site["alt"]
                
                if base in valid_bases:
                    allele_id = f"{chr}:{pos_str}_{base}"
                    
                    variants[ismpl][(chr, pos_str)].add(allele_id)
                    entries_on_reads[ismpl][readname].add((chr, pos_str, base, bq, mq))
    
    # remove variants without at least two alleles
    for sample in list(variants.keys()):
        to_remove = []
        for variant_id, alleles in variants[sample].items():
            if len(alleles) < 2:
                to_remove.append(variant_id)

        for variant_id in to_remove:
            variants[sample].pop(variant_id)
            
    return entries_on_reads, variants

# Search for alleles covered by given read
def get_alleles_on_read(read, entries_on_reads, variants):
    entries = entries_on_reads[read]
    alleles_on_read = set()
    conflict_alleles = set()
    base_qualities = dict()
    
    for chr, pos, allele, bq, mq in entries:
        allele_id = chr + ":" + pos + "_" + allele
        alleles_on_read.add(allele_id)
        map_qual = int(mq)
        base_qualities[allele_id] = int(bq)
        
        # conflict alleles are not used in current implementation
        for al in variants[(chr, pos)]:
            if al!=allele_id:
                conflict_alleles.add(al)
                base_qualities[al] = int(bq)

    return alleles_on_read, conflict_alleles, base_qualities, map_qual