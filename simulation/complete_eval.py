import argparse
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="Evaluate phased variants against a ground truth VCF.")
    parser.add_argument("--variant_file", required=True, type=str, help="Predicted variant file to evaluate (VCF or custom format).")
    parser.add_argument("--true_vcf", default="phased_hets_22.vcf", type=str, help="Ground truth VCF file.")
    return parser.parse_args()

def load_vcf(vcf_path: str, sample_col: int = -1):
    data = []
    block = {}
    block_id = 0

    with open(vcf_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#SET"):
                if block_id != 0:
                    data.append(block)
                block_id += 1
                block = {}
                continue
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            filt = fields[6]
            fmt = fields[8]
            sample = fields[sample_col]

            if filt == "UNPHASED":
                block[(chrom, pos)] = "UNPHASED"
            else:
                fmt_fields = fmt.split(":")
                try:
                    gt_idx = fmt_fields.index("GT")
                    gt_val = sample.split(":")[gt_idx]
                    hap0, hap1 = gt_val.replace("/", "|").split("|")
                    block[(chrom, pos)] = (hap0, hap1)
                except (ValueError, IndexError):
                    if ":" not in sample:
                        try:
                            hap0, hap1 = sample.replace("/", "|").split("|")
                            block[(chrom, pos)] = (hap0, hap1)
                        except ValueError:
                            pass

        if block or block_id == 0:
            data.append(block)

    return data

def load_output_file(file_path: str):
    data = []
    block = {}

    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("##SAMPLE:"):
                continue

            if line.startswith("#BLOCK:"):
                if block:
                    data.append(block)
                    block = {}
                continue

            fields = line.split()
            if len(fields) >= 6:
                chrom = fields[0]
                pos = int(fields[1])
                filt = fields[4]
                gt = fields[5]

                if filt == "UNPHASED" or gt == "UNPHASED":
                    block[(chrom, pos)] = "UNPHASED"
                else:
                    hap0, hap1 = gt.replace("/", "|").split("|")
                    block[(chrom, pos)] = (hap0, hap1)

    if block:
        data.append(block)

    return data

def load_pred_file(file_path: str):
    with open(file_path) as f:
        first_line = f.readline().strip()

    if first_line.startswith("##SAMPLE:") or first_line.startswith("#BLOCK:"):
        return load_output_file(file_path)
    else:
        return load_vcf(file_path)

def compare_sets(v1, v2):
    matches, mismatches, unphased, only2 = [], [], [], []

    for key in v2.keys():
        if key in v1:
            hap1 = v1[key]
            hap2 = v2[key]

            if hap2 == "UNPHASED" or hap1 == "UNPHASED":
                unphased.append((key, hap1, hap2))
            elif hap1 == hap2:
                matches.append((key, hap1, hap2))
            else:
                mismatches.append((key, hap1, hap2))
        else:
            only2.append(key)

    if len(matches) < len(mismatches):
        matches, mismatches = mismatches, matches

    return matches, mismatches, unphased, only2

def calculate_switch_errors(true_dict, pred_block):
    switches = 0
    opportunities = 0

    valid_keys = [
        k for k in sorted(pred_block.keys())
        if k in true_dict and pred_block[k] != "UNPHASED" and true_dict[k] != "UNPHASED"
    ]

    if len(valid_keys) < 2:
        return 0, 0

    prev_match = None
    for key in valid_keys:
        is_match = (pred_block[key] == true_dict[key])

        if prev_match is not None:
            opportunities += 1
            if is_match != prev_match:
                switches += 1

        prev_match = is_match

    return switches, opportunities

def main(args):
    try:
        true_data = load_vcf(args.true_vcf)
        if not true_data:
            print("Error: Reference VCF contains no data.", file=sys.stderr)
            sys.exit(1)
        true = true_data[0]

        pred = load_pred_file(args.variant_file)
        if not pred:
            print("Error: predicted file contains no data.", file=sys.stderr)
            sys.exit(1)

    except FileNotFoundError as e:
        print(f"Error opening file: {e}", file=sys.stderr)
        sys.exit(1)

    total_match = 0
    total_mismatch = 0
    total_unphased = 0
    total_only2 = 0
    sum_mismatch_ratio = 0
    fully_phased = 0

    total_switches = 0
    total_switch_opportunities = 0

    all_predicted_keys = set()

    for entry in pred:
        all_predicted_keys.update(entry.keys())

        matches, mismatches, unphased, only2 = compare_sets(true, entry)
        switches, opportunities = calculate_switch_errors(true, entry)


        total_match += len(matches)
        total_mismatch += len(mismatches)
        total_unphased += len(unphased)
        total_only2 += len(only2)

        total_switches += switches
        total_switch_opportunities += opportunities

        total_variants = len(matches) + len(mismatches)

        if total_variants > 0:
            sum_mismatch_ratio += len(mismatches) / total_variants

        if len(mismatches) == 0 and len(unphased) == 0 and len(only2) == 0 and total_variants > 0:
            fully_phased += 1

    num_blocks = len(pred)
    total_only1 = len(set(true.keys()) - all_predicted_keys)
    global_ser = (total_switches / total_switch_opportunities) if total_switch_opportunities > 0 else 0

    print("================ EVALUATION SUMMARY ================")
    print(f"Total blocks evaluated:  {num_blocks}")
    print(f"Total matches:           {total_match}")
    print(f"Total mismatches:        {total_mismatch}")
    print(f"Total unphased:          {total_unphased}")
    print(f"Overall accuracy:        {total_match/(total_match+total_mismatch+total_unphased)}")
    print(f"Only in true (Missed):   {total_only1}")
    print(f"Only in pred (False):    {total_only2}")
    print("----------------------------------------------------")
    print(f"Switch errors:           {total_switches} (out of {total_switch_opportunities} opportunities)")
    print(f"Global Switch Error Rate:{global_ser:.4f}")
    print("----------------------------------------------------")
    print(f"Average mismatch ratio:  {sum_mismatch_ratio/num_blocks:.4f}" if num_blocks > 0 else "No valid blocks to calc ratio")
    print(f"Fully phased sets:       {fully_phased/num_blocks:.4f}" if num_blocks > 0 else "No valid blocks to calc fully phased")
    print("====================================================")

if __name__ == "__main__":
    main_args = parse_arguments()
    main(main_args)

