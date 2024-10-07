## Defining QV stratifications for DeepPolisher verkko model 2

### HG005

#### Run QV in the GIAB regions, without the ONT gaps

```
# project ont only regions to polished asm coords
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest \
    python3 /home/programs/src/project_blocks_multi_thread.py --threads 32 --mode 'ref2asm' \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_hifi_verkko_model_2_hprc_filters/analysis/hprc_polishing_QC_no_meryl_outputs/HG005_hifi_verkko_model_2_hprc_filters.dipPolToRaw.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/final_ont.gaps.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/final_ont.gaps.polished_projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/final_ont.gaps.polished_projection.bed
```
Now subtract ONT filling regions from GIAB conf refions
```
# combine hap1 and hap2 for polished asm
cat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/Merqury_stratifications/align_asm_project_blocks/HG005_verkko_polished/analysis/align_asm_project_blocks_outputs/HG005_hifi_verkko_model_2_hprc_filters_hap1.asmToRefHap1_HG002_intersect_HG005_GIAB_v4.2.1.projection.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/Merqury_stratifications/align_asm_project_blocks/HG005_verkko_polished/analysis/align_asm_project_blocks_outputs/HG005_hifi_verkko_model_2_hprc_filters_hap2.asmToRefHap2_HG002_intersect_HG005_GIAB_v4.2.1.projection.bed > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/Merqury_stratifications/align_asm_project_blocks/HG005_verkko_polished/analysis/align_asm_project_blocks_outputs/HG005_hifi_verkko_model_2_hprc_filters_dip.asmToRef_HG002_intersect_HG005_GIAB_v4.2.1.projection.bed

bedtools subtract \
    -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/Merqury_stratifications/align_asm_project_blocks/HG005_verkko_polished/analysis/align_asm_project_blocks_outputs/HG005_hifi_verkko_model_2_hprc_filters_dip.asmToRef_HG002_intersect_HG005_GIAB_v4.2.1.projection.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/final_ont.gaps.polished_projection.bed \
    > /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/HG005.verkko.polished.dip.GIAB_conf.projection.no_ont_gaps.bed

cat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/Merqury_stratifications/align_asm_project_blocks/HG005_verkko_raw/analysis/align_asm_project_blocks_outputs/assembly.asmToRefHap1_HG002_intersect_HG005_GIAB_v4.2.1.projection.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/Merqury_stratifications/align_asm_project_blocks/HG005_verkko_raw/analysis/align_asm_project_blocks_outputs/assembly.asmToRefHap2_HG002_intersect_HG005_GIAB_v4.2.1.projection.bed > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/Merqury_stratifications/align_asm_project_blocks/HG005_verkko_raw/analysis/align_asm_project_blocks_outputs/assembly.asmToRefDip_HG002_intersect_HG005_GIAB_v4.2.1.projection.bed

bedtools subtract \
    -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/Merqury_stratifications/align_asm_project_blocks/HG005_verkko_raw/analysis/align_asm_project_blocks_outputs/assembly.asmToRefDip_HG002_intersect_HG005_GIAB_v4.2.1.projection.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/final_ont.gaps.polished_projection.bed \
    > /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/HG005.verkko.raw.dip.GIAB_conf.projection.no_ont_gaps.bed
```

subset FA
```
bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_verkko/assembly.fasta -bed /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/HG005.verkko.raw.dip.GIAB_conf.projection.no_ont_gaps.bed -fo /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/HG005.verkko.raw.dip.GIAB_conf.projection.no_ont_gaps.fa

bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_hifi_verkko_model_2_hprc_filters/applyPolish_dipcall_outputs/HG005_hifi_verkko_model_2_hprc_filters.dip.polished.fasta -bed /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/HG005.verkko.polished.dip.GIAB_conf.projection.no_ont_gaps.bed -fo /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/HG005.verkko.polished.dip.GIAB_conf.projection.no_ont_gaps.fa
```
Run merqury
```
#!/bin/bash
#SBATCH --job-name=merqury
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

time docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    -v /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/merqury_raw_GIAB/:/data \
    juklucas/hpp_merqury:latest merqury.sh \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k31.30x.meryl \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/HG005.verkko.raw.dip.GIAB_conf.projection.no_ont_gaps.fa \
    merqury_raw_verkko

time docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    -v /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/merqury_pol_GIAB/:/data \
    juklucas/hpp_merqury:latest merqury.sh \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k31.30x.meryl \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/HG005.verkko.polished.dip.GIAB_conf.projection.no_ont_gaps.fa \
    merqury_polished_verkko
```
Percent of genome
```
awk '{sum += $3-$2}END{print sum}'  /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/exclude_ONT_gaps/HG005.verkko.raw.dip.GIAB_conf.projection.no_ont_gaps.bed
# 4935864693

awk '{sum += $3-$2}END{print sum}' /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/Merqury_stratifications/align_asm_project_blocks/HG005_verkko_raw/analysis/align_asm_project_blocks_outputs/assembly.asmToRefDip_HG002_intersect_HG005_GIAB_v4.2.1.projection.bed

# 4936281804 - 4935864693 = 417111  bases of ONT gaps in confidence regions

awk '{sum += $3-$2}END{print sum}' /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/final_ont.gaps.bed
# 2194000 bases of ONT gaps whole genome
```

#### Coverage dropouts stratifications for HG005 and primates

Run mosdepth for hifiasm
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

Run mosdepth for verkko
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
    -v /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth_verkko/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -Q 1 --threads 4 \
    -f /private/groups/migalab/juklucas/polishing/HG005/verkko/trio_asm/assembly.fasta \
    --quantize 0:5:10:150: \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth_verkko/HG005_mosdepth \
    /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/phoenix_batch_submissions_manuscript/HG005_verkko_model1/analysis/hprc_DeepPolisher_outputs/HG005_verkko_model1.hifi.to.diploid.asm.PHARAOH.bam

zcat /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth_verkko/HG005_mosdepth.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth_verkko/HG005_mosdepth_quantized.tsv

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
        print }' /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth_verkko/HG005_mosdepth_quantized.tsv > /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth_verkko/HG005_mosdepth_quantized.bed

grep NO_COVERAGE /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth_verkko/HG005_mosdepth_quantized.bed > /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth_verkko/HG005_mosdepth_quantized.lt5x_cov.bed

awk -v OFS="\t" '{if ($3 -$2 >100) print $1,$2,$3 }' /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth_verkko/HG005_mosdepth_quantized.lt5x_cov.bed > /private/groups/patenlab/mira/hprc_polishing/verkko_model_debug/HG005_QV/mosdepth_verkko/HG005_mosdepth_quantized.lt5x_cov.gt100bp.MAPQ1.bed
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
cut -f1 -d"," /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/merqury_stratifications/DeepPolisher_manuscript/Merqury_stratifications.csv | grep -v "sample_id" | while read line ; do
    echo $line
    Hap1Bed=`grep $line Merqury_stratifications.csv | cut -f10 -d","`
    Hap2Bed=`grep $line Merqury_stratifications.csv | cut -f11 -d","`
    cat $Hap1Bed $Hap2Bed > diploid_beds/${line}_dip.HG002_intersect_HG005_GIAB_v4.2.1.projection.bed
  done

mkdir QV_beds
# subtract from GIAB
cut -f1 -d"," /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/merqury_stratifications/DeepPolisher_manuscript/Merqury_stratifications.csv | grep -v "sample_id" | while read line ; do
    echo $line
    mosdepth=`grep $line /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/merqury_stratifications/DeepPolisher_manuscript/Merqury_stratifications.csv | cut -f9 -d","`
    bedtools subtract -a diploid_beds/${line}_dip.HG002_intersect_HG005_GIAB_v4.2.1.projection.bed -b ${mosdepth} > QV_beds/${line}_dip.GIAB_conf_projection.gtMAPQ1_5x.bed
  done

# list for entry in CSV
cut -f1 -d"," /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/merqury_stratifications/DeepPolisher_manuscript/Merqury_stratifications.csv | grep -v "sample_id" | while read line ; do
    realpath QV_beds/${line}_dip.GIAB_conf_projection.gtMAPQ1_5x.bed
  done
```

Check number of bases subtracted for each
```
cut -f1 -d"," /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/merqury_stratifications/DeepPolisher_manuscript/Merqury_stratifications.csv | grep -v "sample_id" | while read line ; do
    echo $line " number of bases subtracted: "
    mosdepth=`grep $line /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/merqury_stratifications/DeepPolisher_manuscript/Merqury_stratifications.csv | cut -f9 -d","`
    bedtools intersect -a diploid_beds/${line}_dip.HG002_intersect_HG005_GIAB_v4.2.1.projection.bed -b ${mosdepth} | awk '{sum += $3-$2}END{print sum}'
  done
```
