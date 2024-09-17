## Establishing a high confidence subset of the GIABv4.2.1 benchmark for benchmarking DeepPolisher

### 1. Download HG002_T2T_v1.0.1 assemblies

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/HG002/assemblies/

Location on phoenix cluster:
```
ls /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/data

hg002v1.0.1.mat.fasta.gz
hg002v1.0.1.mat.fasta.gz.fa.gz
hg002v1.0.1.pat.fasta.gz
hg002v1.0.1.pat.fasta.gz.fai
```
### 2. Run dipcall on HG002_T2T_v1.0.1 vs GRCh38


Use custom minimap2 parameters suggested by Justin Zook (-z200000,10000) to align across larger SVs and more divergent regions like the MHC

`dipcall_inputs.json`:

```json
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/data/hg002v1.0.1.mat.fasta.gz",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/data/hg002v1.0.1.pat.fasta.gz",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_Q100_GRCh38_z2k

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem single_machine \
    --maxCores 32 \
    --batchLogsDir ./toil_logs \
    /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/dipcall.z2k.wdl \
    dipcall_inputs.json \
    --outputDirectory ./dipcall_outfiles \
    --outputFile dipcall_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    --restart \
    2>&1 | tee log.txt
```




### 3. Intersect high confidence regions for Q100 and GIAB v4.2.1

```
bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/data/GRCh38_HG002-T2TQ100-V1.0_smvar.benchmark.bed > /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.bed

# intersect dipcall bed, to only use places alignable between the two assemblies

bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_Q100_GRCh38_z2k/dipcall_outfiles/hg002v1.0.1.dipcall.bed > /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_Q100_GRCh38_z2k/HG002_GIABv4.2.1_Q100_confidence_intersected.dipcall_z2k.bed

```

### 4. Run happy on the output of step 2 against GIAB v4.2.1 VCF, in high confidence regions for Q100 and GIAB

Run hap.py
```
#!/bin/bash
#SBATCH --job-name=happy
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00

docker run --rm -u `id -u`:`id -g` \
-v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
/private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_Q100_GRCh38_z2k/dipcall_outfiles/hg002v1.0.1.dipcall.vcf.gz \
-r /private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
-f /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_Q100_GRCh38_z2k/HG002_GIABv4.2.1_Q100_confidence_intersected.dipcall_z2k.bed \
-o /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_Q100_GRCh38_z2k/happy_out \
--pass-only --no-roc --no-json --engine=vcfeval --threads=32
```

Results:
https://docs.google.com/spreadsheets/d/11l6R63-JQf10AbeJrZQk52U2HaIA5BdHZb3HYNxRkU8/edit#gid=0


### 4. From hap.py output, subtract all regions with FP/FN from GIAB high confidence bed, +/- 50bp

```
cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_Q100_GRCh38_z2k

# extract FP/FN from happy vcf
zcat happy_out.vcf.gz | grep "^#" > happy_out.FPFN.vcf
zcat happy_out.vcf.gz | grep -v "^#" | grep :F >> happy_out.FPFN.vcf

# count variants
zcat happy_out.vcf.gz | grep -v "^#" | wc -l
# 4828845

# total hap.py records that are FP or FN
grep -v "^#" happy_out.FPFN.vcf | wc -l
# 2317 hap.py records

# convert FP/FN hap.py variants to bed, add 50bp on each side of the variant to avoid representation issues
export PATH=$PATH:/private/home/mmastora/progs/bin/
/private/home/mmastora/progs/bin/vcf2bed --do-not-split < happy_out.FPFN.vcf | awk '{print $1"\t"$2-50"\t"$3+50"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12}' > happy_out.FPFN.vcf.50bp.bed

# subtract FP/FN hap.py variants from GIAB confidence bed
bedtools subtract -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_Q100_GRCh38_z2k/happy_out.FPFN.vcf.50bp.bed > /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed

# count number of GIAB 4.2.1 variants in the vcf that were subtracted
bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -b /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed | wc -l
# 3890596

bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed | wc -l
# 3889759

# 3890596 - 3889766 = 830 variants are subtracted

# how much of subtracted sequence is in segdups?

bedtools subtract -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed | awk '{sum += $3-$2}END{print sum}'
# 140098

bedtools subtract -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed | bedtools intersect -a - -b /private/groups/patenlab/mira/data/GRCh38_stratifications_v3.3/SegmentalDuplications/GRCh38_segdups.bed | awk '{sum += $3-$2}END{print sum}'
# 89327 / 140098 = 63.7%
```

New bed file for hap.py calculations in polishing manuscript:
```
/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed
```
