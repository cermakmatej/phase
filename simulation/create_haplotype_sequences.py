import pysam
import random
import argparse

# Create haplotype sequences given reference sequence and variant file
def make_haplotype_fastas(ref_fasta, phased_vcf, out_prefix, bias = 0):
    ref = pysam.FastaFile(ref_fasta)
    vcf = pysam.VariantFile(phased_vcf)

    hap1_seqs = {}
    hap2_seqs = {}


    for reference in ref.references:
        chrom = reference.split(':')[0]
        seq = list(ref.fetch(reference))
        hap1_seq = seq.copy()
        hap2_seq = seq.copy()

        first = True

        for rec in vcf:
            


            if rec.chrom != chrom:
                continue

            if "GT" not in rec.samples[0]:
                continue
                
            gt = rec.samples[0]["GT"]

            if gt is None or not rec.samples[0].phased:
                continue
        
            
            alleles = [rec.ref] + list(rec.alts)

            
            pos = rec.pos - 1  


            # if heterozygous with phase (0|1 or 1|0)
            if len(gt) == 2:
                a1 = alleles[gt[0]]
                a2 = alleles[gt[1]]

                # insert alleles into sequence
                hap1_seq[pos:pos+len(rec.ref)] = list(a1)
                hap2_seq[pos:pos+len(rec.ref)] = list(a2)

        hap1_seqs[reference] = "".join(hap1_seq)
        hap2_seqs[reference] = "".join(hap2_seq)

    with open(f"{out_prefix}_hap1.fa", "w") as f1, open(f"{out_prefix}_hap2.fa", "w") as f2:
        for chrom in ref.references:
            f1.write(f">{chrom}\n{hap1_seqs[chrom]}\n")
            f2.write(f">{chrom}\n{hap2_seqs[chrom]}\n")

    ref.close()
    vcf.close()


def main():
    p = argparse.ArgumentParser(description="Create haplotype sequences from reference and variants.")
    p.add_argument("--ref", required=True, help="Input FASTA ")
    p.add_argument("--vcf", required=True, help="Input VCF ")
    p.add_argument("--out", required=True, help="Output prefix")

    args = p.parse_args()

    make_haplotype_fastas(args.ref, args.vcf, args.out)

if __name__ == "__main__":
    main()
