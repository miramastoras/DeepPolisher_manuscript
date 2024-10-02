### Defining QV stratifications for DeepPolisher verkko model 2

#### HG005

Run mosdepth
```
#!/bin/bash
#SBATCH --job-name=mosdepth_HG005
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
    -v /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -Q 1 --threads 4 \
    -f /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --quantize 0:5:10:150: \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth/HG005_mosdepth \
    /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_PHv6_DPmm2model1/toil_hprc_deepPolisher_out/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.PHARAOHv6.bam

zcat /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth/HG005_mosdepth.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth/HG005_mosdepth_quantized.tsv

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
        print }' /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth/HG005_mosdepth_quantized.tsv > /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth/HG005_mosdepth_quantized.bed

grep NO_COVERAGE /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth/HG005_mosdepth_quantized.bed > /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth/HG005_mosdepth_quantized.lt5x_cov.bed

awk -v OFS="\t" '{if ($3 -$2 >100) print $1,$2,$3 }' /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth/HG005_mosdepth_quantized.lt5x_cov.bed > /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth/HG005_mosdepth_quantized.lt5x_cov.gt100bp.MAPQ1.bed
```
Subtract "hifi dropout" regions from the GIAB confidence regions projections.

Workflow for aligning assemblies to GRCh38 and projecting over the confidence regions: https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/align_asm_project_blocks/DP_manuscript_merqury_stratifications

CSV file of all file locations for this analysis: https://docs.google.com/spreadsheets/d/17HVGUJ7cvf8SvxkLNVg4cGcIvDTaDCJvnGf56CT27d4/edit?gid=900221925#gid=900221925

Location of bed files:
```
/private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/Merqury_stratifications/subtract_dropout_from_GIAB
```

Subtract:
```
cd /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/Merqury_stratifications/subtract_dropout_from_GIAB

mkdir diploid_beds

# combine bed files
cut -f1 -d"," Merqury_stratifications.csv | grep -v "sample_id" | while read line ; do
    Hap1Bed=`grep $line Merqury_stratifications.csv | cut -f10 -d","`
    Hap2Bed=`grep $line Merqury_stratifications.csv | cut -f11 -d","`
    cat $Hap1Bed $Hap2Bed > diploid_beds/${line}_dip.HG002_intersect_HG005_GIAB_v4.2.1.projection.bed
  done

mkdir QV_beds
# subtract from GIAB
cut -f1 -d"," Merqury_stratifications.csv | grep -v "sample_id" | while read line ; do
    mosdepth=`grep $line Merqury_stratifications.csv | cut -f9 -d","`
    bedtools subtract -a diploid_beds/${line}_dip.HG002_intersect_HG005_GIAB_v4.2.1.projection.bed -b ${mosdepth} \
    > QV_beds/${line}_dip.GIAB_conf_projection.gtMAPQ1_5x.bed
  done

# list for entry in CSV
cut -f1 -d"," Merqury_stratifications.csv | grep -v "sample_id" | while read line ; do
    realpath QV_beds/${line}_dip.GIAB_conf_projection.gtMAPQ1_5x.bed
  done
```

Now, run merqury stratifications wdl:
