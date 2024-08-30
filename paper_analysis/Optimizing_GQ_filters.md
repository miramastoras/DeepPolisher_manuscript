## Optimizing GQ filters on DeepPolisher vcf for HPRC intermediate assembly

2/16/2024

#### Step 1: Selecting samples with varying QVs

- Select 4 samples with varying QV improvements after polishing with unfiltered DeepPolisher vcf

Samples chosen: HG01975,HG01993,HG04115,HG02129

#### Step 2: Get Merqury FP kmer bed files for raw assemblies, polished assemblies, using k=21 and k=31

Run hprc polishing QC workflow: [Located here](https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/hprc_polishing_QC/optimize_GQ_filters_HPRC)

#### Step 3: Run Mobin's wdl for annotating the polishing vcf with FP kmers

[Located here](https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/hprc_polishing_QC/optimize_GQ_filters_HPRC/annotate_edit_with_fp_kmers)

#### Step 4: Find optimal GQ filters across 4 samples

For k=21 and k=31 samples

 Add all the FP kmers from each sample together, then calculate the optimal GQ filter for all variants, INS-1 and the rest, and INS-1, DEL-1 and the rest. Calculate the total QV points improvement


Get total FP kmers in raw assembly for Mobin's calculation:
```
cd /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/optimize_GQ_filters

for line in HG01975 HG01993 HG04115 HG02129 ; do tar -zxvf ${line}_k21/hprc_polishing_QC_outputs/${line}.merqury.tar.gz -C ./tmp/ ; totalK=`cat ./tmp/${line}.merqury.altHap_only.bed ./tmp/${line}.merqury.asm_only.bed | wc -l` ; echo counts,${line}_k21,${totalK}; done  | grep counts | cut -f2-3 -d"," >> raw_fp_kmer_counts.csv

for line in HG01975 HG01993 HG04115 HG02129 ; do tar -zxvf ${line}_k31/hprc_polishing_QC_outputs/${line}.merqury.tar.gz -C ./tmp/ ; totalK=`cat ./tmp/${line}.merqury.altHap_only.bed ./tmp/${line}.merqury.asm_only.bed | wc -l` ; echo counts,${line}_k31,${totalK}; done  | grep counts | cut -f2-3 -d"," >> raw_fp_kmer_counts.csv
```

Python notebook used to calculate optimal QV: https://github.com/miramastoras/DeepPolisher_manuscript/blob/main/scripts/Optimizing_GQ_filters_HPRC_int_asm.ipynb
