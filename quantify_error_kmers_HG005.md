# Analysis to account for the error kmers in QV that we aren't polishing

### Quantify number of error kmers falling in regions of hifi coverage dropout in HG005 assembly

1. run mosdepth quantize on HiFi bam
```
#!/bin/bash
#SBATCH --job-name=HG005_mosdepth
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/polishing_dropouts/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth --threads 4 \
    -f /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --quantize 0:1:10:150: \
    /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized \
    /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_PHv6_DPmm2model1/toil_hprc_deepPolisher_out/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.PHARAOHv6.bam

zcat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.tsv

awk 'BEGIN {FS=OFS="\t"} {
        if ($4 == "0:1") {
            $4 = "NO_COVERAGE"
        } else if ($4 == "1:10") {
            $4 = "LOW_COVERAGE"
        } else if ($4 == "10:150") {
            $4 = "CALLABLE"
        } else if ($4 == "150:inf") {
            $4 = "HIGH_COVERAGE"
        }
        print }' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.tsv > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.bed

```
2. extract only "no coverage" regions
```
grep NO_COVERAGE HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.bed | bedtools merge -i - > HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed

grep LOW_COVERAGE HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.bed | bedtools merge -i - > HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.LOW_COVERAGE.bed
```
3. Overlap with raw FP kmer bed files

```
# total error kmers in no coverage regions
cat HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.asm_only.bed HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.altHap_only.bed > HG005_dip_fp_kmers.bed

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/raw_results/HG005_dip_fp_kmers.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed -wa | sort | uniq | wc -l

# 147419 / 1617645

# error kmers unchanged by polishing in no coverage regions

# FP kmers unchanged by polishing
bedtools intersect \
    -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.dip_only.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/annotate_fp_kmers/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ.polished.merqury.dip_only.bed \
    | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed -wa | sort | uniq | awk '{sum+=$4;} END{print sum;}'

156351 / 1575969

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_and_LOW_COVERAGE.bed -wa | sort | uniq | awk '{sum+=$4;} END{print sum;}'

893947 / 1575969
```

```
# total error kmers in low (1-10) and no coverage regions
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/raw_results/HG005_dip_fp_kmers.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_and_LOW_COVERAGE.bed -wa | sort | uniq | wc -l

# 147419 / 1617645

855079 / 1617645

# error kmers unchanged by polishing in no coverage regions

# FP kmers unchanged by polishing
bedtools intersect \
    -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.dip_only.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/annotate_fp_kmers/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ.polished.merqury.dip_only.bed \
    | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed

```

Repeat for 5x coverage
```
#!/bin/bash
#SBATCH --job-name=HG005_mosdepth
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/polishing_dropouts/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth --threads 4 \
    -f /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --quantize 0:5:10:150: \
    /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized \
    /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_PHv6_DPmm2model1/toil_hprc_deepPolisher_out/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.PHARAOHv6.bam

zcat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.tsv

awk 'BEGIN {FS=OFS="\t"} {
        if ($4 == "0:5") {
            $4 = "NO_COVERAGE"
        } else if ($4 == "5:10") {
            $4 = "LOW_COVERAGE"
        } else if ($4 == "10:150") {
            $4 = "CALLABLE"
        } else if ($4 == "150:inf") {
            $4 = "HIGH_COVERAGE"
        }
        print }' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.tsv > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.bed

grep NO_COVERAGE /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed
```
intersect
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/raw_results/HG005_dip_fp_kmers.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed -wa | sort | uniq | wc -l

# 409805

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed -wa | sort | uniq | awk '{sum+=$4;} END{print sum;}'
# 443860
```


Try with 3

```
#!/bin/bash
#SBATCH --job-name=HG005_mosdepth
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/polishing_dropouts/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth --threads 4 \
    -f /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --quantize 0:3:10:150: \
    /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized \
    /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_PHv6_DPmm2model1/toil_hprc_deepPolisher_out/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.PHARAOHv6.bam

zcat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.tsv

awk 'BEGIN {FS=OFS="\t"} {
        if ($4 == "0:3") {
            $4 = "NO_COVERAGE"
        } else if ($4 == "3:10") {
            $4 = "LOW_COVERAGE"
        } else if ($4 == "10:150") {
            $4 = "CALLABLE"
        } else if ($4 == "150:inf") {
            $4 = "HIGH_COVERAGE"
        }
        print }' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.tsv > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.bed

grep NO_COVERAGE /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed
```
intersect
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/raw_results/HG005_dip_fp_kmers.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed -wa | sort | uniq | wc -l

# 289570

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed -wa | sort | uniq | awk '{sum+=$4;} END{print sum;}'
# 307487
```

how many polishing edits are in 5,10x coverage?

```
# subset to windows > 100bp
awk -v OFS="\t" '{if ($3 -$2 >=100) print $1,$2,$3 }' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.gt.100bp.bed

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GQ20_INS1_GQ12_DEL1_GQ5_else/HG005.mm2_model1.polisher_output.GQ20_INS1_GQ12_DEL1_GQ5_else.vcf.gz -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.gt.100bp.bed | sort | uniq | wc -l

6441

awk -v OFS="\t" '{if ($3 -$2 >=100) print $1,$2,$3 }' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.100bp.bed

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GQ20_INS1_GQ12_DEL1_GQ5_else/HG005.mm2_model1.polisher_output.GQ20_INS1_GQ12_DEL1_GQ5_else.vcf.gz -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.100bp.bed | sort | uniq | wc -l


```


Checking coverage dropouts greater than 50bp so we dont pick up errors in insertions
```
# 3 x
awk -v OFS="\t" '{if ($3 -$2 >100) print $1,$2,$3 }' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.3x.gt.100bp.bed

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_3/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.3x.gt.100bp.bed -wa | sort | uniq | awk '{sum+=$4;} END{print sum;}'
# 259121

# 5 x
awk -v OFS="\t" '{if ($3 -$2 >100) print $1,$2,$3 }' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed  > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed -wa | sort | uniq | wc -l

359392/1,574,857

# 10 x
awk -v OFS="\t" '{if ($3 -$2 >100) print $1,$2,$3 }' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_and_LOW_COVERAGE.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_and_LOW_COVERAGE.10x.gt100bp.bed

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_and_LOW_COVERAGE.10x.gt100bp.bed -wa | sort | uniq | awk '{sum+=$4;} END{print sum;}'

# 0x

awk -v OFS="\t" '{if ($3 -$2 >100) print $1,$2,$3 }' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.gt100bp.bed

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.gt100bp.bed -wa | sort | uniq | awk '{sum+=$4;} END{print sum;}'
```

## Project CHM13 simple repeat annotations to HG005 raw assembly and intersect with unchanged error kmers

Downloaded from:
https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/

```
cd /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation
grep Low_complexity chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed > chm13_simple_low_complexity.bed
grep Simple_repeat chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed >>chm13_simple_low_complexity.bed
```

dipcall HG005 to chm13
Dipcall raw assembly against GRCh38
```

{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/chm13v2.0.fa.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/chm13v2.0.fa",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/dipcall_chm13

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24,18]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem single_machine \
    --maxCores 32 \
    --batchLogsDir ./toil_logs \
    /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/dipcall.wdl \
    dipcall_inputs.json \
    --outputDirectory ./dipcall_outfiles \
    --outputFile dipcall_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```

Run projection with Mobin's script
```
# 2. removed unmapped records in paf files
cd /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/dipcall_chm13/dipcall_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall/

grep "tp:A:P" HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.paf > ../HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf
grep "tp:A:P" HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.paf > ../HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf

bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/chm13_simple_low_complexity.srt.bed | bedtools merge -i - > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/chm13_simple_low_complexity.srt.mrg.bed

# 3. Project repeat annotation bedfile to y2 assembly, each haplotype separately
# Hap1
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/dipcall_chm13/dipcall_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf \
  --blocks /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/chm13_simple_low_complexity.srt.mrg.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap1.chm13_simple_low_complexity.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap1.chm13_simple_low_complexity.projection.bed

# Hap 2
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/dipcall_chm13/dipcall_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf \
  --blocks /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/chm13_simple_low_complexity.srt.mrg.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap2.chm13_simple_low_complexity.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap2.chm13_simple_low_complexity.projection.bed

cat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap2.chm13_simple_low_complexity.projection.bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap1.chm13_simple_low_complexity.projection.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed

# Number of bases in repeat annotation
bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed | bedtools merge -i - | awk '{sum += $3-$2}END{print sum}'

# 143921510 / 5972405005 = 2.4% of genome
```

Intersect projected annotations with error kmers
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed -wa | sort | uniq | wc -l

# 569002 / 1,574,857 = 36 %

# How many bases of repeat annotation overlap with <5x coverage

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed -wa | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_unchanged_error_kmers_5x_cov.bed

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed -wa | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_unchanged_error_kmers_repeats.bed

bedtools intersect -f 1 -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_unchanged_error_kmers_repeats.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_unchanged_error_kmers_5x_cov.bed | sort | uniq | wc -l

116200 / 1,574,857 = 7.3%

(569002 - 116200) / 1,574,857 = 28%
(359392 - 116200) /1,574,857 = 15%
```


Look in IGV at examples of 47% remaining

```
cat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_unchanged_error_kmers_5x_cov.bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_unchanged_error_kmers_repeats.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_unchanged_error_kmers_repeats_5x_coverage.bed

bedtools subtract -wa -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_unchanged_error_kmers_repeats_5x_coverage.bed | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/remaining_47_percent.bed
```
Projecting all chm13 repeats over

```
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/dipcall_chm13/dipcall_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf \
  --blocks /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/chm13_collapsed_RM.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap1.chm13_RM.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap1.chm13_RM.projection.bed

# Hap 2
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/dipcall_chm13/dipcall_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf \
  --blocks /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/chm13_collapsed_RM.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap2.chm13_RM.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap2.chm13_RM.projection.bed

cat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap2.chm13_RM.projection.bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap1.chm13_RM.projection.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_RM.projection.bed
```

```
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'asm2ref' \
  --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.hap1PolToRaw.paf \
  --blocks /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/stop_codon.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/stop_codon_projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/stop_codon_projection.bed
```

### Check correspondence with element data

Create new element HG005 meryl dbs for k=21 and k=31


Check coverage for HG005 element data
```
#!/bin/bash

#SBATCH --job-name=element_mosdepth
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=200gb
#SBATCH --threads-per-core=1
#SBATCH --output=slurm_logs/downsample_%x_%j_%A_%a.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -n --fast-mode -t 4 /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/mosdepth/element_HG005 \
    /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/HG005.element_cloudbreak_350bp_ins.grch38.bam

source /private/home/mmastora/progs/miniconda3/etc/profile.d/conda.sh
conda activate analysis

python3 ~/progs/mosdepth/scripts/plot-dist.py /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/mosdepth/element_HG005.mosdepth.global.dist.txt &> /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/mosdepth/coverages.txt
```

Downsampling to 30x for meryl to match illumina
```
#!/bin/bash
#SBATCH --job-name=meryl
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --output=%x.%j.log
#SBATCH --time=6:00:00

export PATH=/private/home/mmastora/progs/meryl-1.4.1/bin:$PATH

samtools view -s 0.65 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/HG005.element_cloudbreak_350bp_ins.grch38.bam > /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/HG005.element_cloudbreak_350bp_ins.grch38.30x.bam

samtools fastq -@32 /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/HG005.element_cloudbreak_350bp_ins.grch38.30x.bam > /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/HG005.element_cloudbreak_350bp_ins.grch38.30x.fastq

meryl count threads=32 k=31 /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/HG005.element_cloudbreak_350bp_ins.grch38.30x.fastq output /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.k31.element.30x.meryl

meryl count threads=32 k=21 /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/HG005.element_cloudbreak_350bp_ins.grch38.30x.fastq output /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.k21.element.30x.meryl
```
Create yak DB
```
#!/bin/bash
#SBATCH --job-name=yak_HG5
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

time docker run --rm -u `id -u`:`id -g` -v /private/groups:/private/groups \
    juklucas/hpp_yak:latest yak count \
    -k21 -b37 -t32 \
    -o /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/yak_files/HG005.k21.30x.element.yak /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/HG005.element_cloudbreak_350bp_ins.grch38.30x.fastq

time docker run --rm -u `id -u`:`id -g` -v /private/groups:/private/groups \
    juklucas/hpp_yak:latest yak count \
    -k31 -b37 -t32 \
    -o /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/yak_files/HG005.k31.30x.element.yak /private/groups/patenlab/mira/hprc_polishing/data/element_HG005/HG005.element_cloudbreak_350bp_ins.grch38.30x.fastq
```

Get illumina kmers unchanged by polishing
```
# -f 1 requires a be completely contained in by to be reported
bedtools intersect -f 1 -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.dip_only.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/annotate_fp_kmers/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ.polished.merqury.dip_only.bed \
    | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed

# 1574857 /
```

Get element FP kmers unchanged by polishing
```
# k31
cat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_DP_manuscript_element_k31/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.hap1PolToRaw_asm_only.projection.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_DP_manuscript_element_k31/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.hap2PolToRaw_altHap_only.projection.bed > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_DP_manuscript_element_k31/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/HG005_polished_dip_FP_kmers_elementk31_projection.bed
```

```
bedtools intersect -f 1 -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_DP_manuscript_element_k31/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/raw_results/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.dip_only.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_DP_manuscript_element_k31/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/HG005_polished_dip_FP_kmers_elementk31_projection.bed \
    | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed

3097256 / 3142727


# How many error kmers unchanged by polishing in element overlap with illumina?

# in element, in illumina
bedtools intersect \
    -wa -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed \
    | sort | uniq | wc -l

1514000 / 3098598 = 48%

# in illumina, in element
bedtools intersect -wa \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed \
    | sort | uniq | wc -l

1132843 / 1575969 = 71%

# require 100% reciprocal overlap

bedtools intersect -f 1 -r \
    -wa -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed \
    | sort | uniq | wc -l

948307 / 3097256 = 30%

# in illumina, in element
bedtools intersect -f 1 -r \
    -wa -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed \
    | sort | uniq | wc -l

948307 / 1574857 = 60%
```

How many kmers validated by element are in cov dropout or repeats?
```
# illumina kmers validated by element

bedtools intersect -f 1 -r \
    -wa -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed


# 5x coverage or less, for gt 100 bp

bedtools intersect -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed | sort | uniq | wc -l

294555 / 948307 = 31 %

# overlap simple repeat or low complexity

bedtools intersect -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed | sort | uniq | wc -l

315053 / 948307 = 33%

# both simple repeat and low complexity

bedtools intersect -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed | sort | uniq >  /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.5x.gt100bp.bed

bedtools intersect -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.low_complexity.bed

s
bedtools intersect -f 1 -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.low_complexity.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.5x.gt100bp.bed | sort | uniq | wc -l

# 82955
```


Counts so far:
```
948307 / 1574857 error kmers unchanged by polishing validated by element

294555-82955 = 211600/948307 = 22% in low coverage regions
315053-82955 = 232098/948307 = 24% in low complexity repeats
82955 / 948307 = 8.7% in low coverage and low complexity repeat
```

Checking overlap of non-validated error kmers with low complexity and coverage regions

```
bedtools intersect -v -f 1 -r \
    -wa -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_NOTin_elementk31.bed

# 5x coverage or less, for gt 100 bp

bedtools intersect -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_NOTin_elementk31.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed | sort | uniq | wc -l

(64837-33245) / 626550 = 10 %

# overlap simple repeat or low complexity

bedtools intersect -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_NOTin_elementk31.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed | sort | uniq | wc -l

253949 / 626550 = 33%

# both simple repeat and low complexity

bedtools intersect -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_NOTin_elementk31.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed | sort | uniq >  /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_NOTin_elementk31.5x.gt100bp.bed

bedtools intersect -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_NOTin_elementk31.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_NOTin_elementk31.low_complexity.bed

s
bedtools intersect -f 1 -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_NOTin_elementk31.low_complexity.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_NOTin_elementk31.5x.gt100bp.bed | sort | uniq | wc -l

33245 / 626550 =
253949- 33245

```
### Check for higher GC content


Merge error kmers of the three categories:
```
# illumina unchanged by polishing
bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed | bedtools merge -i - > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.bed

# illumina and element
bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed | bedtools merge -i - > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.bed

# just element unchanged
bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed | bedtools merge -i - > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.bed
```

Calculate GC and AT content for each
```
bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.nuc.bed

bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.nuc.bed

bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.nuc.bed  
```

Get txt file of counts
```
ls | grep .nuc.bed | while read line
    do cut -f4 $line | grep -v "_pct" > ${line}.AT.txt
    cut -f5 $line | grep -v "_pct" > ${line}.GC.txt
    done  
```

randomly permute same size regions for each file across genome
```
bedtools shuffle -seed 2400 -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.bed -g /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.fai.genome > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.RANDOM.bed

bedtools shuffle -seed 2400 -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.bed -g /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.fai.genome > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.RANDOM.bed

bedtools shuffle -seed 2400 -i  /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.bed -g /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.fai.genome > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.RANDOM.bed
```

Calculate GC and AC content for the random shuffled regions
```
bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.RANDOM.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.RANDOM.nuc.bed

bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.RANDOM.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.RANDOM.nuc.bed

bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.RANDOM.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.RANDOM.nuc.bed
```
Plot histograms: https://colab.research.google.com/drive/18kY97xhAcgFaaGU6tH2L-IYEnmbbvC0E?authuser=1#scrollTo=dtQs42ojHGsY

#### Also just check number of error kmers in each category
```
# calculate GC content for all error kmers unchanged by polishing
bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed

# count error kmers with > n % GC content

awk '$5 > 0.9' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 23215 / 1574857 = 1.4%

awk '$5 > 0.8' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 156986 / 1574857 = 10%

awk '$5 > 0.7' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 651829 / 1574857 = 41%

awk '$5 > 0.6' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 1132050 / 1574857 = 71%

# count error kmers with > n % AT content

awk '$4 > 0.9' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 63799 / 1574857 = 4%

awk '$4 > 0.8' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 109984 / 1574857 = 6.9%

awk '$4 > 0.7' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 154277 / 1574857 = 9.7%

awk '$4 > 0.6' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 204499 / 1574857 = 12.9%
```

Checking the ones validated by element
```
bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed

awk '$5 > 0.9' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | wc -l
# 12821 / 948307 = 1.3 %

awk '$5 > 0.8' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | wc -l
# 92141 / 948307 = 9.7%

awk '$5 > 0.7' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | wc -l
# 419247 / 948307 = 44.2%

awk '$5 > 0.6' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | wc -l
# 726149 / 948307 = 76%

# count error kmers with > n % AT content

awk '$4 > 0.9' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | wc -l
# 34815 / 948307 = 2%

awk '$4 > 0.8' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | wc -l
# 45809 / 948307 = 2.9%

awk '$4 > 0.7' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | wc -l
# 59799 / 948307 = 3.7%

awk '$4 > 0.6' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | wc -l
# 89316 / 948307 = 9.4%
```

GC / AT content Permutation histogram plots:

Get GC content of 100,000 random k=31 windows in genome
```
awk -v OFS='\t' {'print $1,$2'} /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.fai > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.fai.genome

# illumina unchanged n=1574857
bedtools random -l 31 -n 100000 -seed 71346 \
    -g /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.fai.genome \
    > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/100k_random_k31.bed

bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/100k_random_k31.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/100k_random_k31.nuc.bed

cut -f 7 /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/100k_random_k31.nuc.bed | grep -v "pct" > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/100k_random_k31.nuc.perc_AT.txt

cut -f 8 /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/100k_random_k31.nuc.bed | grep -v "pct" > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/100k_random_k31.nuc.perc_GC.txt

# illumina x element n = 948307
```

k=31 error kmers unchanged by polishing, illumina
```
cut -f 4 /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | grep -v "pct" > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/HG005_ilm_k31_error_kmers_unchanged_by_polishing.perc_AT.txt

cut -f 5 /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | grep -v "pct" > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/HG005_ilm_k31_error_kmers_unchanged_by_polishing.perc_GC.txt
```

k=31 error kmers unchanged by polishing, element
```
bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.nuc.bed

cut -f 4 /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.nuc.bed | grep -v "pct" >/private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/HG005_error_kmers_unchanged_by_polishing.elementk31.nuc.perc_AT.txt

cut -f 5 /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.nuc.bed | grep -v "pct" >/private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/HG005_error_kmers_unchanged_by_polishing.elementk31.nuc.perc_GC.txt
```

k=31 error kmers unchanged by polishing in illumina and element
```
cut -f4 /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | grep -v "pct" > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.AT.txt

cut -f5 /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | grep -v "pct" > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/histograms/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.GC.txt
```


Number of error kmers with > 70% GC content and how they overlap with the other categories

```
awk '$5 > 0.7' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | grep -v "pct" > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.70_perc_gc.bed

cat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_unchanged_error_kmers_repeats.bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_unchanged_error_kmers_5x_cov.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_unchanged_error_kmers_5x_cov.repeats.bed

bedtools subtract -wa -f 1 -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.70_perc_gc.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_unchanged_error_kmers_5x_cov.repeats.bed | sort | uniq | wc -l

# 213710 / 1,574,857

/private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.70_perc_gc.bed | sort | uniq | wc -l
```

#### Annotate homopolymers

```
from Bio import SeqIO

def find_homopolymers(sequence, min_length=5):
    """Find homopolymer regions in a sequence."""
    homopolymers = []
    n = len(sequence)
    i = 0
    while i < n:
        start = i
        while i < n and sequence[i] == sequence[start]:
            i += 1
        length = i - start
        if length >= min_length:
            homopolymers.append((start, i, sequence[start]))
    return homopolymers

def fasta_to_bed(fasta_file, bed_file, min_length=5):
    """Convert FASTA file to BED file with homopolymers."""
    with open(bed_file, 'w') as out_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequence = str(record.seq)
            homopolymers = find_homopolymers(sequence, min_length)
            for start, end, base in homopolymers:
                # Output BED format: chromosome, start, end, name, score, strand
                out_file.write(f"{record.id}\t{start}\t{end}\t{base}_homopolymer\t0\t+\n")

# Example usage
fasta_file = '/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa'  # Replace with your FASTA file path
bed_file = '/private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_homopolymers.gt5.bed'     # Replace with desired output BED file path
fasta_to_bed(fasta_file, bed_file)
```

```
/private/home/mmastora/progs/scripts/find_homopolymers.py
```


```
# overlaps with homopolymer > 5bp
bedtools intersect -wa -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_homopolymers.gt5.bed | sort | uniq | wc -l

1050681 / 1574857 = 66%

# homopolymer > 10 bp
bedtools intersect -wa -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_homopolymers.gt10.bed | sort | uniq | wc -l

166080 / 1574857 = 10 %
```


Combine with other categories
```
cut -f 1-3 /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.70_perc_gc.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.70_perc_gc.cut.bed

bedtools intersect -wa -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_homopolymers.gt5.bed | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_error_kmers_unchanged_by_polishing.homopolymers.bed

bedtools intersect -wa -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_homopolymers.gt10.bed | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_error_kmers_unchanged_by_polishing.homopolymers.10.bed

cat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_unchanged_error_kmers_5x_cov.repeats.bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.70_perc_gc.cut.bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_error_kmers_unchanged_by_polishing.homopolymers.10.bed | sort | uniq | wc -l

1412654 / 1574857 = 89% homopolymers > 5

1118379 / 1574857 = 71 % 
```
