## Creating a new truth set for training DeepPolisher, using only places where GIAB v4.2.1 and T2T v1.0.1 agree

02/22/2024

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
### 3. Run happy on the output of step 2 against GIAB v4.2.1 VCF, in high confidence regions

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
#SBATCH --cpus-per-task=16
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00

docker run --rm -u `id -u`:`id -g` \
-v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
/private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/dipcall_outfiles/hg002v1.0.1.dipcall.vcf.gz \
-r /private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
-f /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
-o /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/happy/happy_out \
--pass-only --no-roc --no-json --engine=vcfeval --threads=16
```

Results: [Rows 2-3](https://docs.google.com/spreadsheets/d/1V4aKo-NwwafFELxX0-LcCRxG9hmJ03t3jmxiLfdmEbU/edit#gid=0)
### 4. From hap.py output, subtract all regions with FP/FN from GIAB high confidence bed

```
cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/happy

# extract FP/FN from happy vcf
zcat happy_out.vcf.gz | grep "^#" > happy_out.FPFN.vcf
zcat happy_out.vcf.gz | grep -v "^#" | grep :F >> happy_out.FPFN.vcf

# count variants subtracted
zcat happy_out.vcf.gz | grep -v "^#" | wc -l
# 4828845

grep -v "^#" happy_out.FPFN.vcf | wc -l
# 32503 to be subtracted

# subtract FP/FN from GIAB confidence set
bedtools subtract -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/happy/happy_out.FPFN.vcf > /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/HG002_GRCh38_T2T_platinum.bed
```

### 5. Run happy again on the output of step 2 against GIAB v4.2.1 but with  HG002_GRCh38_T2T_platinum.bed. F1-score should be close to 1

```
#!/bin/bash
#SBATCH --job-name=happy
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=16
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00

cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/happy_platinum_bed/

docker run --rm -u `id -u`:`id -g` \
-v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
/private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/dipcall_outfiles/hg002v1.0.1.dipcall.vcf.gz \
-r /private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
-f /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/HG002_GRCh38_T2T_platinum.bed \
-o /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_GRCh38/happy_platinum_bed/happy_out \
--pass-only --no-roc --no-json --engine=vcfeval --threads=16
```
Results: [Rows 4-5](https://docs.google.com/spreadsheets/d/1V4aKo-NwwafFELxX0-LcCRxG9hmJ03t3jmxiLfdmEbU/edit#gid=0)

### 6. Project HG002_GRCh38_T2T_platinum.bed to hg002_hprc_y2 assembly

Using alignments of GRCh38 to hprc y2 assemblies generated by a previous dipcall run:

```
/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.grch38.paf

/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.grch38.paf
```
Project platinum bed to hprc y2 assemblies
```
# project to pat
docker run --rm -u `id -u`:`id -g` \
-v /private/groups:/private/groups \
mobinasri/flagger:latest python3 \
/home/programs/src/project_blocks_multi_thread.py \
--threads 16 --mode 'ref2asm' \
--paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.grch38.paf \
--blocks /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/HG002_GRCh38_T2T_platinum.bed \
--outputProjectable /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/project_platinum_to_hprc_y2/HG002_GRCh38_T2T_platinum.projectable.y2.pat.bed \
--outputProjection /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/project_platinum_to_hprc_y2/HG002_GRCh38_T2T_platinum.projection.y2.pat.bed

# project to mat
docker run --rm -u `id -u`:`id -g` \
-v /private/groups:/private/groups \
mobinasri/flagger:latest python3 \
/home/programs/src/project_blocks_multi_thread.py \
--threads 16 --mode 'ref2asm' \
--paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/dipcall_grch38/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.grch38.paf \
--blocks /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/HG002_GRCh38_T2T_platinum.bed \
--outputProjectable /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/project_platinum_to_hprc_y2/HG002_GRCh38_T2T_platinum.projectable.y2.mat.bed \
--outputProjection /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/project_platinum_to_hprc_y2/HG002_GRCh38_T2T_platinum.projection.y2.mat.bed

```
Count projectable bases
```
# total bases
awk '{sum += $3-$2}END{print sum}' /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/HG002_GRCh38_T2T_platinum.bed

2542201256

# total projectable bases to pat
awk '{sum += $3-$2}END{print sum}' /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/project_platinum_to_hprc_y2/HG002_GRCh38_T2T_platinum.projectable.y2.pat.bed

2539442469

# total projectable bases to mat
awk '{sum += $3-$2}END{print sum}' /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/project_platinum_to_hprc_y2/HG002_GRCh38_T2T_platinum.projectable.y2.mat.bed

2541655670
```

### 7. Run dipcall on HG002_T2T_v1.0.1 vs HG002 hprc y2 assembly



`dipcall_inputs.json`:

```json
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/data/hg002v1.0.1.mat.fasta.gz",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/data/hg002v1.0.1.pat.fasta.gz",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}
```

Dipcall T2T against HPRC y2 Pat

`dipcall_inputs.json`:

```json
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/data/hg002v1.0.1.mat.fasta.gz",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/data/hg002v1.0.1.pat.fasta.gz",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}
```
toil locations:
```
cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_HPRC_y2_mat

cd /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_HPRC_y2_pat
```
toil command:
```
export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=4:00:00 --partition=medium"
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

### 8. Intersect dipcall outputs (VCFs) with HG002_GRCh38_T2T_platinum_projected.bed

```
# pat to pat
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/dipcall_T2T_HPRC_y2_pat/dipcall_outfiles/hg002v1.0.1.dipcall.vcf.gz -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/project_platinum_to_hprc_y2/HG002_GRCh38_T2T_platinum.projection.y2.pat.bed
```

### 9. Apply