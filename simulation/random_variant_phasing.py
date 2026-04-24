import pysam
import random
import argparse

# From unphased vcf create randomly phased vcf
def random_phase_vcf(input_vcf, output_vcf):
    vcf_in = pysam.VariantFile(input_vcf, "r")
    vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

    for rec in vcf_in.fetch():
        for sample in rec.samples:
            gt = rec.samples[sample]["GT"]
            if gt is None:
                continue

            if len(gt) == 2 and gt[0] != gt[1]:
                if random.random() < 0.5:
                    phased_gt = (gt[0], gt[1])
                else:
                    phased_gt = (gt[1], gt[0])
                rec.samples[sample]["GT"] = phased_gt
                rec.samples[sample].phased = True
            else:
                rec.samples[sample]["GT"] = gt
        vcf_out.write(rec)

    vcf_in.close()
    vcf_out.close()


def main():
    p = argparse.ArgumentParser(description="Randomly phase variants.")
    p.add_argument("--input", required=True, help="Input VCF (unphased)")
    p.add_argument("--out", required=True, help="Outpu VCF (randomly phased)")

    args = p.parse_args()

    random_phase_vcf(args.input, args.out)

if __name__ == "__main__":
    main()



