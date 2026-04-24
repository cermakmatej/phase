#!/usr/bin/env bash
set -euo pipefail

bcftools view -v snps -g het simulated_indiv_chr22.vcf -Ov -o simulated_hets_chr22.vcf

for i in {0..40..5}
do

python ./insert_repeats.py --ref chr22.fa --vcf simulated_hets_chr22.vcf --out_fa references/reference_"$i".fa --out_vcf tmp_variants_"$i".vcf --target_percent "$i" > logs/log"$i".txt

python ./random_variant_phasing.py --in tmp_variants_"$i".vcf --out vars/phased_variants_"$i".vcf

python ./create_haplotype_sequences.py --ref references/reference_"$i".fa --vcf vars/phased_variants_"$i".vcf --out tmp_sim_22_repeats_"$i"

../utils/art_bin_MountRainier/art_illumina -ss HS25 -i tmp_sim_22_repeats_"$i"_hap1.fa -p -l 150 -m 550 -s 150 -f 15 -o tmp_hap1_R
../utils/art_bin_MountRainier/art_illumina -ss HS25 -i tmp_sim_22_repeats_"$i"_hap2.fa -p -l 150 -m 550 -s 150 -f 15 -o tmp_hap2_R

# renaming reads to prevent collisions of names caused by simulation
awk 'NR%4==1 {print "@" "hap1_" substr($0,2); next} {print}' tmp_hap1_R1.fq > tmp_hap1_R1.renamed.fq
awk 'NR%4==1 {print "@" "hap1_" substr($0,2); next} {print}' tmp_hap1_R2.fq > tmp_hap1_R2.renamed.fq

awk 'NR%4==1 {print "@" "hap2_" substr($0,2); next} {print}' tmp_hap2_R1.fq > tmp_hap2_R1.renamed.fq
awk 'NR%4==1 {print "@" "hap2_" substr($0,2); next} {print}' tmp_hap2_R2.fq > tmp_hap2_R2.renamed.fq


cat tmp_hap1_R1.renamed.fq tmp_hap2_R1.renamed.fq > tmp_sim_R1.fq
cat tmp_hap1_R2.renamed.fq tmp_hap2_R2.renamed.fq > tmp_sim_R2.fq

bwa index references/reference_"$i".fa
bwa mem -t 12 references/reference_"$i".fa tmp_sim_R1.fq tmp_sim_R2.fq | samtools sort -o sim.bam
samtools index sim.bam

bcftools mpileup -f references/reference_"$i".fa -Ou sim.bam | bcftools call   -mv   -Ov   -o sim_reps"$i".vcf

bcftools +phase -a sim.bam -f references/reference_"$i".fa -s sim_reps"$i".vcf | python ../phase.py --variant_file sim_reps"$i".vcf > reps"$i"_output.txt

python ./complete_eval.py --variant_file reps"$i"_output.txt --true_vcf vars/phased_variants_"$i".vcf> evaluated_reps"$i".txt
rm tmp_*
done

