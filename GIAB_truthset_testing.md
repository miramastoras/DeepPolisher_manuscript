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
cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=4:00:00 --partition=high_priority"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
    /private/groups/hprc/polishing/hpp_production_workflows/QC/wdl/tasks/dipcall.wdl \
    dipcall_inputs.json \
    --outputDirectory ./dipcall_outfiles \
    --outputFile dipcall_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```
### 3. Intersect high confidence regions for Q100 and GIAB v4.2.1

```
bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/data/GRCh38_HG002-T2TQ100-V1.0_smvar.benchmark.bed > /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.bed

# intersect dipcall bed, to only use places alignable between the two assemblies

bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/dipcall_outfiles/hg002v1.0.1.dipcall.bed > /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.dipcall.bed
```

### 4. Run happy on the output of step 2 against GIAB v4.2.1 VCF, in high confidence regions for Q100 and GIAB

Run hap.py
```
cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/happy/happy_out
```

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
/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/dipcall_outfiles/hg002v1.0.1.dipcall.vcf.gz \
-r /private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
-f /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.dipcall.bed \
-o /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_T2T_GRCh38/happy_dipcall/happy_out \
--pass-only --no-roc --no-json --engine=vcfeval --threads=32
```
Results:
https://docs.google.com/spreadsheets/d/11l6R63-JQf10AbeJrZQk52U2HaIA5BdHZb3HYNxRkU8/edit#gid=0


### 4. From hap.py output, subtract all regions with FP/FN from GIAB high confidence bed, +/- 50bp

```
cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_T2T_GRCh38/happy_dipcall/

# extract FP/FN from happy vcf
zcat happy_out.vcf.gz | grep "^#" > happy_out.FPFN.vcf
zcat happy_out.vcf.gz | grep -v "^#" | grep :F >> happy_out.FPFN.vcf

# count variants
zcat happy_out.vcf.gz | grep -v "^#" | wc -l
# 4828845

# total hap.py records that are FP or FN
grep -v "^#" happy_out.FPFN.vcf | wc -l
# 2151 records will be subtracted

# convert FP/FN hap.py variants to bed, add 50bp on each side of the variant
export PATH=$PATH:/private/home/mmastora/progs/bin/
/private/home/mmastora/progs/bin/vcf2bed --do-not-split < happy_out.FPFN.vcf | awk '{print $1"\t"$2-50"\t"$3+50"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12}' > happy_out.FPFN.vcf.50bp.bed

# subtract FP/FN hap.py variants from GIAB confidence bed
bedtools subtract -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_T2T_GRCh38/happy_dipcall/happy_out.FPFN.vcf.50bp.bed > /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.bed

# count number of GIAB 4.2.1 variants in the vcf that were subtracted
bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -b /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed | wc -l
# 3890596

bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.bed | wc -l
# 3889759

# 3890596 - 3889759 = 837 variants are subtracted

# how much of subtracted sequence is in segdups?

bedtools subtract -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.bed | awk '{sum += $3-$2}END{print sum}'
# 137412

bedtools subtract -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.bed | bedtools intersect -a - -b /private/groups/patenlab/mira/data/GRCh38_stratifications_v3.3/SegmentalDuplications/GRCh38_segdups.bed | awk '{sum += $3-$2}END{print sum}'
# 85766 / 137412 = 62%
```
### 5. Checking hap.py results with Justin Zook's

dipcall vcf: https://drive.google.com/drive/u/1/folders/1etIgAyu-WB4Wn0LmvPsj7erUjUsx4SFm

Rerunning hap.py with Justin's dipcall vcf

Dipcall bed:
https://drive.google.com/drive/u/1/folders/1LULEWMdIS9XYKU8kFl8R9hQxLpB-P33L
```
# Intersect beds with his dipcall bed:
bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GRCh38_HG2-T2TQ100-V1.0_dipcall-z2k.dip.bed > /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.dipcall_JustinZ.bed

awk '{sum += $3-$2}END{print sum}'  /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.dipcall_JustinZ.bed
# 2517725219

awk '{sum += $3-$2}END{print sum}' /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.dipcall.
# 2512576364

# 5,148,855 bases missing from my dipcall bed

bedtools subtract -a /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GRCh38_HG2-T2TQ100-V1.0_dipcall-z2k.dip.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/dipcall_outfiles/hg002v1.0.1.dipcall.bed | awk '{sum += $3-$2}END{print sum}'

# 23888872

bedtools subtract -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GRCh38_HG2-T2TQ100-V1.0_dipcall-z2k.dip.bed -a /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/dipcall_outfiles/hg002v1.0.1.dipcall.bed | awk '{sum += $3-$2}END{print sum}'

# 11111036
```


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
/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GRCh38_v4.2.1_HG2-T2TQ100-V1.0_smvar_dipcall-z2k.vcf.gz \
-r /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz \
-f /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.dipcall_JustinZ.bed \
-o /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/Justin_dipcall/happy_out \
--pass-only --no-roc --no-json --engine=vcfeval --threads=32
```

```
cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_T2T_GRCh38/Justin_dipcall/

# extract FP/FN from happy vcf
zcat happy_out.vcf.gz | grep "^#" > happy_out.FPFN.vcf
zcat happy_out.vcf.gz | grep -v "^#" | grep :F >> happy_out.FPFN.vcf

# count variants
zcat happy_out.vcf.gz | grep -v "^#" | wc -l
# 4828845

# total hap.py records that are FP or FN
grep -v "^#" happy_out.FPFN.vcf | wc -l
# 2151 records will be subtracted

# convert FP/FN hap.py variants to bed, add 50bp on each side of the variant
export PATH=$PATH:/private/home/mmastora/progs/bin/
/private/home/mmastora/progs/bin/vcf2bed --do-not-split < happy_out.FPFN.vcf | awk '{print $1"\t"$2-50"\t"$3+50"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12}' > happy_out.FPFN.vcf.50bp.bed

# subtract FP/FN hap.py variants from GIAB confidence bed
bedtools subtract -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_T2T_GRCh38/Justin_dipcall/happy_out.FPFN.vcf.50bp.bed > /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.Justin_dipcall.bed

# count number of GIAB 4.2.1 variants in the vcf that were subtracted
bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -b /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed | wc -l
# 3890596

bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.Justin_dipcall.bed | wc -l
# 3889759

# 3890596 - 3889759 = 837 variants are subtracted

# how much of subtracted sequence is in segdups?

bedtools subtract -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.bed | awk '{sum += $3-$2}END{print sum}'
# 137412

bedtools subtract -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.bed | bedtools intersect -a - -b /private/groups/patenlab/mira/data/GRCh38_stratifications_v3.3/SegmentalDuplications/GRCh38_segdups.bed | awk '{sum += $3-$2}END{print sum}'
# 85766 / 137412 = 62%
```


Rerun dipcall with -z200000,10000

`dipcall_inputs.json`:

```json
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/data/hg002v1.0.1.mat.fasta.gz",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/data/hg002v1.0.1.pat.fasta.gz",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_Q100_GRCh38_decoys_z2k

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
    2>&1 | tee log.txt
```
Rerun dipcall with -z200000,10000 and the other fasta file

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
