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

intersect dipcall bed

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
-f /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.bed \
-o /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_T2T_GRCh38/happy/happy_out \
--pass-only --no-roc --no-json --engine=vcfeval --threads=16

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
cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_T2T_GRCh38/happy/

# extract FP/FN from happy vcf
zcat happy_out.vcf.gz | grep "^#" > happy_out.FPFN.vcf
zcat happy_out.vcf.gz | grep -v "^#" | grep :F >> happy_out.FPFN.vcf

# count variants subtracted
zcat happy_out.vcf.gz | grep -v "^#" | wc -l
# 4828845

grep -v "^#" happy_out.FPFN.vcf | wc -l
# 18088 to be subtracted

bedtools intersect -a happy_out.FPFN.vcf -b /private/groups/patenlab/mira/data/HG002_GIABv4.2.1_Q100_confidence_intersected.bed | wc -l

# convert to bed, add 50bp on each side of the variant
export PATH=$PATH:/private/home/mmastora/progs/bin/

# all GIAB
/private/home/mmastora/progs/bin/vcf2bed --do-not-split < happy_out.FPFN.vcf | awk '{print $1"\t"$2-50"\t"$3+50"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12}' > happy_out.FPFN.vcf.50bp.bed

# subtract FP/FN from GIAB confidence set
bedtools subtract -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/dipcall_T2T_GRCh38/happy/happy_out.FPFN.vcf.50bp.bed > /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.bed

# count number of GIAB 4.2.1 variants that were subtracted
bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -b /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed | wc -l
# 3890596

bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_truthset_testing/GIAB_T2T_Q100_conf_beds_concordant_50bp.bed | wc -l
# 3873818

16778 variants are subtracted
```
