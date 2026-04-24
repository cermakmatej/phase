import pysam
import random

input_vcf = "gnomad.genomes.r2.1.1.sites.22.vcf.bgz"
output_vcf = "simulated_indiv_chr22.vcf"

vcf_in = pysam.VariantFile(input_vcf)

vcf_in.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

vcf_in.header.add_sample("Simulated_Indiv")

vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)

for rec in vcf_in:
    try:
        af = rec.info['AF'][0]
    except KeyError:
        continue

    p_ref_ref = (1 - af) ** 2
    p_ref_alt = 2 * af * (1 - af)

    rand_val = random.random()

    if rand_val < p_ref_ref:
        gt = (0, 0)
    elif rand_val < p_ref_ref + p_ref_alt:
        gt = (0, 1)
    else:
        gt = (1, 1)

    if gt != (0, 0):
        new_rec = vcf_out.new_record(
            contig=rec.contig,
            start=rec.start,
            stop=rec.stop,
            alleles=rec.alleles,
            id=rec.id,
            qual=rec.qual
        )

        new_rec.samples["Simulated_Indiv"]["GT"] = gt
        vcf_out.write(new_rec)

vcf_in.close()
vcf_out.close()

print("Simulation succesful")
