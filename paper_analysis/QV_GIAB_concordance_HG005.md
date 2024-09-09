## Manual investigation of HG005 GIAB FP / FN variants alongside QV

In this analysis, we are gathering files for viewing in IGV to perform a manual investigation into the cause of the vastly different error rates suggested by the GIAB variants and the QV score 

Directory for analysis:
```
/private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/
```

### 1. Polishing edits annotated by inducing fp kmers, neutral, unchanging, fixing fp kmers,
as well as bed files for the projected FP kmers from merqury
```
cp /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/HG005_y2_DCv1.2_PHv6_DPmm2_model1_docker_v0.0.8_12122023/annotate_edit_with_fp_kmers/HG005_mm2_model1_k31/annotate_edit_with_fp_kmers_outputs/* /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/
```
Location of run on phoenix:
https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/annotate_edit_with_fp_kmers/HG005_mm2_model1_dockerv0.8

### 2. Get GIAB variants projections to raw assembly

Intersect raw hap.py vcf and polished hap.py vcf
```
mkdir -p /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy

bcftools isec /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/HG005_y2_DCv1.2_PHv6_DPmm2_model1_docker_v0.0.8_12122023/applyPolish_dipcall_happy/happy_out.vcf.gz /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallhapDotPyVCF/HG005_y2_raw_happy.vcf.gz -p /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/HG005_pol_raw_happy_isec
```

Get GIAB FP FN in polished asm not in raw asm
```
grep "^#" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/HG005_pol_raw_happy_isec/0000.vcf > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unique_to_pol_asm.vcf

grep -v "^#" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/HG005_pol_raw_happy_isec/0000.vcf | grep :F > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unique_to_pol_asm.vcf

# 2149 variants
```
Get GIAB TP induced by polishing, projections from hg38 to raw assembly
```
grep "^#" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/HG005_pol_raw_happy_isec/0000.vcf > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/TP.happy.unique_to_pol_asm.vcf

grep -v "^#" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/HG005_pol_raw_happy_isec/0000.vcf | grep :TP > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/TP.happy.unique_to_pol_asm.vcf

# 941
```
Get GIAB FP FN unchanged in the assembly
```
grep "^#" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/HG005_pol_raw_happy_isec/0002.vcf > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unchanged_by_polishing.vcf

grep -v "^#" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/HG005_pol_raw_happy_isec/0002.vcf | grep :F >> /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unchanged_by_polishing.vcf

# 8132
```
Project these over
```
# convert to bed
export PATH=$PATH:/private/home/mmastora/progs/bin/

/private/home/mmastora/progs/bin/vcf2bed < /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unique_to_pol_asm.vcf | awk '{print $1"\t"$2-10"\t"$3+10"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12}'> /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unique_to_pol_asm.vcf.10bp.bed

/private/home/mmastora/progs/bin/vcf2bed < /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unchanged_by_polishing.vcf | awk '{print $1"\t"$2-10"\t"$3+10"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12}'> /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unchanged_by_polishing.vcf.10bp.bed

/private/home/mmastora/progs/bin/vcf2bed < /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/TP.happy.unique_to_pol_asm.vcf | awk '{print $1"\t"$2-10"\t"$3+10"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12}'> /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/TP.happy.unique_to_pol_asm.vcf.10bp.bed

# hap 1 fp fn polished only
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unique_to_pol_asm.vcf.10bp.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unique_to_pol_asm.hap1.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unique_to_pol_asm.hap1.projection.bed \
    --threads 16

# hap 2 fp fn polished only
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unique_to_pol_asm.vcf.10bp.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unique_to_pol_asm.hap2.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unique_to_pol_asm.hap2.projection.bed \
    --threads 16

# hap 1 fp fn unchanged
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unchanged_by_polishing.vcf.10bp.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unchanged_by_polishing.hap1.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unchanged_by_polishing.hap1.projection.bed \
    --threads 16

# hap2 FP FN unchanged
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/FPFN.happy.unchanged_by_polishing.vcf.10bp.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unchanged_by_polishing.hap2.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unchanged_by_polishing.hap2.projection.bed \
    --threads 16

#
# hap 1 fp fn unchanged
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/TP.happy.unique_to_pol_asm.vcf.10bp.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/TP.happy.unique_to_pol_asm.hap1.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/TP.happy.unique_to_pol_asm.hap1.projection.bed \
    --threads 16

# hap2 FP FN unchanged
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/intersect_happy/TP.happy.unique_to_pol_asm.vcf.10bp.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/TP.happy.unique_to_pol_asm.hap2.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/TP.happy.unique_to_pol_asm.hap2.projection.bed \
    --threads 16
```

Combined projectable and projected
```
# TP unique to polished
paste -d "\t" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/TP.happy.unique_to_pol_asm.hap2.projection.bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/TP.happy.unique_to_pol_asm.hap2.projectable.bed > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/TP.happy.unique_to_pol_asm.hap2.projected.bed

paste -d "\t" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/TP.happy.unique_to_pol_asm.hap1.projection.bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/TP.happy.unique_to_pol_asm.hap1.projectable.bed > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/TP.happy.unique_to_pol_asm.hap1.projected.bed

# FP FN unchanged
paste -d "\t" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unchanged_by_polishing.hap2.projection.bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unchanged_by_polishing.hap2.projectable.bed > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/FPFN.happy.unchanged_by_polishing.hap2.projected.bed

paste -d "\t" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unchanged_by_polishing.hap1.projection.bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unchanged_by_polishing.hap1.projectable.bed > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/FPFN.happy.unchanged_by_polishing.hap1.projected.bed

#
paste -d "\t" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unique_to_pol_asm.hap2.projection.bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unique_to_pol_asm.hap2.projectable.bed > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/FPFN.happy.unique_to_pol_asm.hap2.projected.bed

paste -d "\t" /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unique_to_pol_asm.hap1.projection.bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/project_blocks_GIAB_to_raw/FPFN.happy.unique_to_pol_asm.hap1.projectable.bed > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/FPFN.happy.unique_to_pol_asm.hap1.projected.bed

# combined haplotypes for IGV track:
cat FPFN.happy.unchanged_by_polishing.hap1.projected.bed FPFN.happy.unchanged_by_polishing.hap2.projected.bed > FPFN.happy.unchanged_by_polishing.dip.projected.bed

cat FPFN.happy.unique_to_pol_asm.hap1.projected.bed FPFN.happy.unique_to_pol_asm.hap2.projected.bed > FPFN.happy.unique_to_pol_asm.dip.projected.bed

cat TP.happy.unique_to_pol_asm.hap1.projected.bed TP.happy.unique_to_pol_asm.hap2.projected.bed > TP.happy.unique_to_pol_asm.dip.projected.bed
```

project HG005 confidence regions as a track
```
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf \
    --blocks /private/groups/patenlab/mira/data/HG005_GRCh38_1_22_v4.2.1_benchmark.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.hap2.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.hap2.projection.bed \
    --threads 16

#
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf \
    --blocks /private/groups/patenlab/mira/data/HG005_GRCh38_1_22_v4.2.1_benchmark.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.hap1.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.hap1.projection.bed \
    --threads 16

cat /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.hap1.projection.bed /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.hap2.projection.bed > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.dip.projection.bed


# project dipcall x conf regions
bedtools intersect -a /private/groups/patenlab/mira/data/HG005_GRCh38_1_22_v4.2.1_benchmark.bed -b /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.bed  > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.raw_dipcall.isec.conf_regions.bed

docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.raw_dipcall.isec.conf_regions.srt.mrg.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.dipcall.bed.hap2.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.dipcall.bed.hap2.projection.bed \
    --threads 16

#
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.raw_dipcall.isec.conf_regions.srt.mrg.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.dipcall.bed.hap1.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.dipcall.bed.hap1.projection.bed \
    --threads 16

cat /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.dipcall.bed.hap1.projection.bed /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.dipcall.bed.hap2.projection.bed > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.dipcall.bed.dip.projection.bed
```

Project all FP/FN in GIAB polished to raw asm
```
# pull FP FN variants
zcat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/HG005_y2_DCv1.2_PHv6_DPmm2_model1_docker_v0.0.8_12122023/applyPolish_dipcall_happy/happy_out.vcf.gz | grep "^#" > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/polished.happy.FPFN.vcf

zcat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/HG005_y2_DCv1.2_PHv6_DPmm2_model1_docker_v0.0.8_12122023/applyPolish_dipcall_happy/happy_out.vcf.gz | grep -v "^#" | grep :F > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/polished.happy.FPFN.vcf

# convert to bed
export PATH=$PATH:/private/home/mmastora/progs/bin/
/private/home/mmastora/progs/bin/vcf2bed < /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/polished.happy.FPFN.vcf | awk '{print $1"\t"$2-10"\t"$3+10"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12}'> /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/polished.happy.FPFN.vcf.10bp.bed

# project over
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/polished.happy.FPFN.vcf.10bp.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/PolFPFN.projectable.hap1.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/PolFPFN.projection.hap1.bed \
    --threads 16

# hap2 FP FN unchanged
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/polished.happy.FPFN.vcf.10bp.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/PolFPFN.projectable.hap2.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/PolFPFN.projection.hap2.bed \
    --threads 16

cat /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/PolFPFN.projection.hap2.bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/PolFPFN.projection.hap1.bed > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect/PolFPFN.projection.dip.bed
```

```
cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/all_variants_intersect

# 11151 / 11816 are projectable from grch38

bedtools intersect -a PolFPFN.projection.dip.bed -b ../HG005_mm2_model1_k31.fixed_fp_kmer_blocks.bed | sort | uniq | wc -l

# 126

bedtools intersect -a PolFPFN.projection.dip.bed -b ../HG005_mm2_model1_k31.induced_fp_kmer_blocks.bed | sort | uniq | wc -l
# 2228

bedtools intersect -a PolFPFN.projection.dip.bed -b ../HG005_mm2_model1_k31.unchanged_fp_kmer_blocks.bed | sort | uniq | wc -l
# 736

# Assuming that we double each count because only one haplotype is wrong, we still have 4970 edits that are either neutral or don't project right to intersect
```

### 3. Manual investigation in IGV:

Select 10 FP/FN induced giab variants randomly
```
cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/randoms

shuf -n 10 ../FPFN.happy.unique_to_pol_asm.hap1.projected.bed > FPFN.happy.unique_to_pol_asm.hap1.rand.bed

cut -f 10 FPFN.happy.unique_to_pol_asm.hap1.rand.bed | while read line ; do grep $line ../FPFN.happy.unique_to_pol_asm.hap2.projected.bed ; done >> FPFN.happy.unique_to_pol_asm.hap2.rand.bed

# reorganize so hap1 and hap2 coords are next to each other in the file , in order
cat FPFN.happy.unique_to_pol_asm.hap1.rand.bed FPFN.happy.unique_to_pol_asm.hap2.rand.bed | sort -V -k9,9 -k10,10 | cut -f1-5,7-13,15- > FPFN.happy.unique_to_pol_asm.dip.rand.bed
```
Select 10 FP kmers induced inside GIAB confidence regions
```
bedtools intersect -a ../HG005_mm2_model1_k31.deeppolisher.annotated_induced_fp_kmers.vcf -b  /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.dip.projection.bed | shuf -n 10 > induced_fp_kmers.rand.vcf

grep
```
Select 10 fixed GIAB variants
```
cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/randoms

shuf -n 10 ../TP.happy.unique_to_pol_asm.hap1.projected.bed > TP.happy.unique_to_pol_asm.hap1.rand.bed

cut -f 10 TP.happy.unique_to_pol_asm.hap1.rand.bed | while read line ; do grep $line ../TP.happy.unique_to_pol_asm.hap2.projected.bed ; done >> TP.happy.unique_to_pol_asm.hap2.rand.bed

# reorganize so hap1 and hap2 coords are next to each other in the file , in order
cat TP.happy.unique_to_pol_asm.hap1.rand.bed TP.happy.unique_to_pol_asm.hap2.rand.bed | sort -V -k9,9 -k10,10 | cut -f1-5,7-13,15- >TP.happy.unique_to_pol_asm.dip.rand.bed
```

Select 10 fixed FP kmers inside conf regions
```
cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_DCv1.2_PHv6_DPmm2_model1/QV_GIAB_concordance/randoms


bedtools intersect -a ../HG005_mm2_model1_k31.deeppolisher.annotated_fixed_fp_kmers.vcf -b  /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG005_GRCh38_1_22_v4.2.1_benchmark.dipcall.bed.dip.projection.bed | shuf -n 10 > fixed_fp_kmers.rand.vcf

# case 1 where GIAB disagrees with ilm and dp
grep 1089789 FPFN.happy.unchanged_by_polishing.hap2.projected.bed | grep h2tg000012l

chr5	107516163	107516184
h2tg000012l     108978994
h1tg000010l	108581954	108581986

# case 2 where GIAB disagrees with ilm and dp

h1tg000010l     124524385

grep 107516163 FPFN.happy.unchanged_by_polishing.hap1.projected.bed
```


#### Preparing IGV files

Look at illumina all to one alignments and both haplotypes and hifi alignments
```
cd /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one
```

Download files here:
```
/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/reads

cat *.fastq.gz ../HG005.ilm.fastq.gz
```
BWA index
```
docker run --rm -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira quay.io/masri2019/hpp_bwa:latest bwa index -p /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa

docker run --rm -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira quay.io/masri2019/hpp_bwa:latest bwa index -p /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa
```

BWA mem alignment
```
#!/bin/bash
#SBATCH --job-name=HG005_ilm_all2mat
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira quay.io/masri2019/hpp_bwa:latest bwa mem -t32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/reads/HG005.ilm.fastq.gz | samtools view -b -h > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.50x.all_to_mat.trio_hifiasm_0.19.5.DC_1.2_40x.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.50x.all_to_mat.trio_hifiasm_0.19.5.DC_1.2_40x.bam
```
```
#!/bin/bash
#SBATCH --job-name=HG005_ilm_all2pat
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira quay.io/masri2019/hpp_bwa:latest bwa mem -t32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/reads/HG005.ilm.fastq.gz | samtools view -b -h > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/paternal/HG005.ilm.50x.all_to_pat.trio_hifiasm_0.19.5.DC_1.2_40x.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/paternal/HG005.ilm.50x.all_to_pat.trio_hifiasm_0.19.5.DC_1.2_40x.bam
```
