## Annotating FP kmers, fig2

Combine HG002 and HG005 hap FP kmer bedfiles from merqury to use
```
cat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG002_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.altHap_only.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG002_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.asm_only.bed | bedtools sort -i - | bedtools merge -i - -c 1 -o count > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG002_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.dip_only.bed

cat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG002_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/Raw.insideConf.subBed.merqury.altHap_only.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG002_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/Raw.insideConf.subBed.merqury.asm_only.bed > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG002_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/Raw.insideConf.subBed.merqury.dip.bed

cat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG002_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/Raw.outsideConf.subBed.merqury.altHap_only.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG002_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/Raw.outsideConf.subBed.merqury.asm_only.bed > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG002_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/Raw.outsideConf.subBed.merqury.dip.bed

cat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv5_winnowmap_model5_dockerv0.8/analysis/hprc_polishing_QC_no_meryl_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.asm_only.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv5_winnowmap_model5_dockerv0.8/analysis/hprc_polishing_QC_no_meryl_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.altHap_only.bed | bedtools sort -i - | bedtools merge -i - -c 1 -o count > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv5_winnowmap_model5_dockerv0.8/analysis/hprc_polishing_QC_no_meryl_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.dip_only.bed

cat
```

project inside and outside conf regions to raw assemblies
```
# HG002 hap1 inside

docker run -it --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/flagger@sha256:5d738412b56bac5a64227569c1d6e57e7920e3d3e5724c17ab233f92279bcff6 \
    python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode "ref2asm" \
    --paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap1.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap1.projection.bed \
    --threads 16

# HG002 hap2 inside
docker run -it --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/flagger@sha256:5d738412b56bac5a64227569c1d6e57e7920e3d3e5724c17ab233f92279bcff6 \
    python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode "ref2asm" \
    --paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap2.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap2.projection.bed \
    --threads 16

cat /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap2.projection.bed /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap1.projection.bed > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002_intersect_HG005_GIAB_v4.2.1.dip.projection.bed

# HG002 hap1 outside
docker run -it --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/flagger@sha256:5d738412b56bac5a64227569c1d6e57e7920e3d3e5724c17ab233f92279bcff6 \
    python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode "ref2asm" \
    --paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/outside_HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap1.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap1.projection.bed \
    --threads 16

# HG002 hap2 outside
docker run -it --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/flagger@sha256:5d738412b56bac5a64227569c1d6e57e7920e3d3e5724c17ab233f92279bcff6 \
    python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode "ref2asm" \
    --paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/outside_HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap2.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap2.projection.bed \
    --threads 16

cat /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap2.projection.bed /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap1.projection.bed > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.dip.projection.bed

# HG005 inside hap1

docker run -it --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/flagger@sha256:5d738412b56bac5a64227569c1d6e57e7920e3d3e5724c17ab233f92279bcff6 \
    python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode "ref2asm" \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap1.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap1.projection.bed \
    --threads 16

# HG005 inside hap2
docker run -it --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/flagger@sha256:5d738412b56bac5a64227569c1d6e57e7920e3d3e5724c17ab233f92279bcff6 \
    python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode "ref2asm" \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap2.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap2.projection.bed \
    --threads 16

cat /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap2.projection.bed /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG002_intersect_HG005_GIAB_v4.2.1.hap1.projection.bed > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG002_intersect_HG005_GIAB_v4.2.1.dip.projection.bed


# HG005 outside hap1

docker run -it --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/flagger@sha256:5d738412b56bac5a64227569c1d6e57e7920e3d3e5724c17ab233f92279bcff6 \
    python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode "ref2asm" \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf \
    --blocks /private/groups/hprc/ref_files/giab/outside_HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap1.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap1.projection.bed \
    --threads 16

# HG005 outside hap2

docker run -it --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/flagger@sha256:5d738412b56bac5a64227569c1d6e57e7920e3d3e5724c17ab233f92279bcff6 \
    python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode "ref2asm" \
    --paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf \
    --blocks /private/groups/hprc/ref_files/giab/outside_HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap2.projectable.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap2.projection.bed \
    --threads 16

cat /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap1.projection.bed /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.hap2.projection.bed > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.dip.projection.bed
```

```
mkdir -p annotate_fp_kmers

cd /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl

inside_conf=/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002_intersect_HG005_GIAB_v4.2.1.dip.projection.bed
outside_conf=/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.dip.projection.bed
```

```sh
# /bin/bash

sample=$1
raw_asm_fp_kmers=$2
inside_conf=$3
outside_conf=$4

mkdir -p ./annotate_fp_kmers/${sample}
###  whole genome

# get total FP kmers
tar -zxvf ${sample}/analysis/hprc_polishing_QC_no_meryl_outputs/*.polished.merqury.tar.gz -C ./annotate_fp_kmers/${sample}/

total_fp_kmers=`cat ./annotate_fp_kmers/${sample}/*.polished.merqury.altHap_only.bed ./annotate_fp_kmers/${sample}/*.polished.merqury.asm_only.bed | wc -l`

# get # FP kmers projectable to raw
projected_fp_kmers=`cat ${sample}/analysis/hprc_polishing_QC_no_meryl_outputs/*hap1PolToRaw_asm_only.projection.bed ${sample}/analysis/hprc_polishing_QC_no_meryl_outputs/*hap2PolToRaw_altHap_only.projection.bed | awk '{sum+=$4;} END{print sum;}'`

# merge haplotype projected files
cat ${sample}/analysis/hprc_polishing_QC_no_meryl_outputs/*hap1PolToRaw_asm_only.projection.bed ${sample}/analysis/hprc_polishing_QC_no_meryl_outputs/*hap2PolToRaw_altHap_only.projection.bed > ./annotate_fp_kmers/${sample}/${sample}.polished.merqury.dip_only.bed

# get FP kmer blocks unchanged by polishing
fp_kmers_unchanged_wg=`bedtools intersect -f 1 \
    -a ./annotate_fp_kmers/${sample}/${sample}.polished.merqury.dip_only.bed \
    -b ${raw_asm_fp_kmers} \
    | sort | uniq | awk '{sum+=$4;} END{print sum;}'`

###  inside conf

# subset FP kmers polished projected to the confidence regions  
bedtools intersect -f 1 -a ./annotate_fp_kmers/${sample}/${sample}.polished.merqury.dip_only.bed -b ${inside_conf} | sort | uniq > ./annotate_fp_kmers/${sample}/${sample}.polished.merqury.dip_only.insideConf.bed

# subset total FP kmers
tar -zxvf ${sample}/analysis/hprc_polishing_QC_no_meryl_outputs/Polished.insideConf.subBed.merqury.tar.gz -C ./annotate_fp_kmers/${sample}/

total_fp_kmers_conf=`cat ./annotate_fp_kmers/${sample}/Polished.insideConf.subBed.merqury.altHap_only.bed ./annotate_fp_kmers/${sample}/Polished.insideConf.subBed.merqury.asm_only.bed | wc -l`

# get projected FP kmers

projected_fp_kmers_conf=`awk '{sum+=$4;} END{print sum;}' ./annotate_fp_kmers/${sample}/${sample}.polished.merqury.dip_only.insideConf.bed`

# FP kmers unchanged by polishing, inside conf
fp_kmers_unchanged_conf=`bedtools intersect -f 1 \
    -a ./annotate_fp_kmers/${sample}/${sample}.polished.merqury.dip_only.insideConf.bed \
    -b ${raw_asm_fp_kmers} \
    | sort | uniq | awk '{sum+=$4;} END{print sum;}'`

echo ${sample},${total_fp_kmers},${projected_fp_kmers},${fp_kmers_unchanged_wg},${total_fp_kmers_conf},${projected_fp_kmers_conf},${fp_kmers_unchanged_conf} >> ./annotate_fp_kmers/all_results.csv
```

```
cd /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl

echo "sample,total_fp_kmers_wg,projected_fp_kmers_wg,unchanged_fp_kmers_wg,total_fp_kmers_conf,projected_fp_kmers_conf,unchanged_fp_kmers_conf" > ./annotate_fp_kmers/all_results.csv

# HG002 samples
for sample in HG002_nextPolish2 HG002_deepvariant HG002_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ HG002_t2t_polish
    do bash annotate_fp_kmers/annotate_fp_kmers.sh ${sample} /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG002_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.dip_only.bed /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002_intersect_HG005_GIAB_v4.2.1.dip.projection.bed /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.dip.projection.bed
  done

# HG005 samples
for sample in HG005_nextPolish2 HG005_deepvariant HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ HG005_t2t_polish
    do bash annotate_fp_kmers/annotate_fp_kmers.sh ${sample} /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv5_winnowmap_model5_dockerv0.8/analysis/hprc_polishing_QC_no_meryl_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.dip_only.bed /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/HG002_intersect_HG005_GIAB_v4.2.1.dip.projection.bed /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/grch38/outside_HG002_intersect_HG005_GIAB_v4.2.1.dip.projection.bed
    done
```

https://genomeark.s3.amazonaws.com/index.html?prefix=species/Gorilla_gorilla/mGorGor1/genomic_data/pacbio_hifi/previous-versions/
https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/mPanPan1.MT.cur.20231122.fasta.gz

s3://genomeark/species/Gorilla_gorilla/mGorGor1/genomic_data/pacbio_hifi/m64076_210215_140546.hifi_reads.bam
https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/
https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/
