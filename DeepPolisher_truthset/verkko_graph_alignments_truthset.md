## Creating a truthset for training a DeepPolisher verkko model from read alignments produced during the assembly process

### 1. Dipcall verkko assembly against truthset (T2T Q100 v1.1)

Sergey provided the assembly and alignments here: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=backup/Sergey/HG002/

Separate hap1 and hap2
Sergey confirmed hap2 = paternal, hap1 = maternal
```
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1

cat verkko_assembly.fasta | seqkit grep -r -p '.*haplotype1' > HG002_verkko2.2.haplotype1.maternal.fasta
cat verkko_assembly.fasta | seqkit grep -r -p '.*haplotype2' > HG002_verkko2.2.haplotype2.paternal.fasta
```
Set up input json files

```
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/HG002_verkko2.2.haplotype1.maternal.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/HG002_verkko2.2.haplotype1.maternal.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.mat_Y_EBV_MT.copy.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.mat_Y_EBV_MT.fasta",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}
```
```
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/HG002_verkko2.2.haplotype2.paternal.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/HG002_verkko2.2.haplotype2.paternal.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.pat_X_EBV_MT.copy.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.pat_X_EBV_MT.fasta",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_pat
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_mat

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24,18]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
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
### 2. Project GIAB confidence regions to verkko assembly

Dipcall raw assembly against GRCh38
```
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/HG002_verkko2.2.haplotype2.paternal.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/HG002_verkko2.2.haplotype1.maternal.fasta",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_grch38

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24,18]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
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
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_grch38/dipcall_outfiles

grep "tp:A:P" HG002_verkko2.2.haplotype1.hap1.paf > ../HG002_verkko2.2.hap1.AP.paf
grep "tp:A:P" HG002_verkko2.2.haplotype1.hap2.paf > ../HG002_verkko2.2.hap2.AP.paf

# 3. Project GIAB confidence bedfile to y2 assembly, each haplotype separately
# Hap1
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_grch38/dipcall_outfiles/HG002_verkko2.2.hap1.AP.paf \
  --blocks /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_grch38/HG002_verkko_2.2_hap1_GIAB_conf.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_grch38/HG002_verkko_2.2_hap1_GIAB_conf.projection.bed
# Hap 2
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_grch38/dipcall_outfiles/HG002_verkko2.2.hap2.AP.paf \
  --blocks /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_grch38/HG002_verkko_2.2_hap2_GIAB_conf.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_grch38/HG002_verkko_2.2_hap2_GIAB_conf.projection.bed
```

Pre-processing truthset files:
```
# merge alignment to truth bams:
samtools merge -o /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_Q100_v1.1_dipcall.merged.bam /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_mat/dipcall_outfiles/hg002v1.1_Y_EBV_MT.dipcall/hg002v1.1_Y_EBV_MT.hap1.bam /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_pat/dipcall_outfiles/hg002v1.1_X_EBV_MT.dipcall/hg002v1.1_X_EBV_MT.hap1.bam
```

```
# merge confidence projection bed files
cat /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_hap1_GIAB_confident_regions.projection.bed /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_hap2_GIAB_confident_regions.projection.bed > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_GIAB_confident_regions.diploid.projection.bed
```
```
# filter for primary alignments to grch38

cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_grch38/dipcall_outfiles/HG002_verkko2.2.haplotype1.dipcall/

samtools merge -o HG002_verkko2.2.dip.GRCH38.merged.bam HG002_verkko2.2.haplotype1.hap1.bam HG002_verkko2.2.haplotype1.hap2.bam

samtools view -F0x900 HG002_verkko2.2.dip.GRCH38.merged.bam -o HG002_verkko2.2.dip.GRCH38.merged.primary.bam

# convert to bed
bedtools bamtobed -i HG002_verkko2.2.dip.GRCH38.merged.primary.bam > HG002_verkko2.2.dip.GRCH38.merged.primary.bed

awk '{print $4 , $1}' HG002_verkko2.2.dip.GRCH38.merged.primary.bed | sort | uniq | wc -l
awk '{print $4 , $1}' HG002_verkko2.2.dip.GRCH38.merged.primary.bed > chrom_maps

# separate fasta into chr 1-19 and chr 20 based on maps
grep -v "chr21" chrom_maps | grep -v "chr22" | grep -v "chr20" | grep -v "chrUn" > chrom_maps_1_19
grep "chr20" chrom_maps > chrom_maps_20

cut -f1 -d " " chrom_maps_1_19 | while read line ; do samtools faidx /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_assembly.fasta $line >> /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_chr1_19.fasta ; done

cut -f1 -d " " chrom_maps_20 | while read line ; do samtools faidx /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_assembly.fasta $line >> /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_chr20.fasta ; done
```
Uploaded for kishwar:
```
gsutil -o GSUtil:parallel_composite_upload_threshold=50GB -m cp -r verkko_graph_alignments_Q100v1.1_truthset gs://pepper-deepvariant/mira/
```

#### Test on chr 20

```
samtools view -bh verkko_assembly.bam haplotype2-0000081 haplotype1-0000017 > verkko_assembly.chr20.bam
samtools index verkko_assembly.chr20.bam
```
https://docs.google.com/spreadsheets/d/1zQ1ICWvyYBqy_7L7vS6QsvlwGJyMIbwQrcZxH_Cbu1g/edit?gid=1300297667#gid=1300297667

Polish assembly
```
bcftools consensus -H2 -f /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_chr20.fasta /private/groups/patenlab/mira/phoenix_batch_executions/workflows/DeepPolisher/HG002_verkko_graph_alignment_chr20/analysis/DeepPolisher_outputs/polisher_output.no_filters.vcf.gz > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/HG002_verkko_2.2_chr20.checkpoint34.fasta
```
Dipcall:
```
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/HG002_verkko_2.2_chr20.checkpoint34.hap1.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/HG002_verkko_2.2_chr20.checkpoint34.hap2.fasta",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_chr20.hap1.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_chr20.hap2.fasta",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
```
#!/bin/bash
#SBATCH --job-name=dipcall_truthset
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --exclude=phoenix-[09,10,22,23,24,18]

cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_raw

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24,18]"
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
hap.py
```
bash /private/home/mmastora/progs/scripts/GIAB_happy_chr20.sh \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/dipcall_outfiles/HG002_verkko_2.2_chr20.checkpoint34.hap2.dipcall.vcf.gz \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG002_y2_raw/dipCallTar/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall.GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/happy/happy \
    HG002

s
bash /private/home/mmastora/progs/scripts/GIAB_happy_chr20.sh \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_raw/dipcall_outfiles/HG002_verkko_2.2_chr20.hap2.dipcall.vcf.gz \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG002_y2_raw/dipCallTar/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall.GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_raw/happy/happy \
    HG002
```
Test on straight v4.2.1 without intersection with Q100
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/whole_genome_dipcall/dipcall_raw/dipcall_outfiles/HG002_verkko2.2.haplotype2.dipcall.bed -b /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/whole_genome_dipcall/dipcall_raw/dipcall_outfiles/HG002_verkko2.2.haplotype2.dipcall.GIABv4.2.1_benchmark_noinconsistent.bed

bash /private/home/mmastora/progs/scripts/GIAB_happy_chr20.sh \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/dipcall_outfiles/HG002_verkko_2.2_chr20.checkpoint34.hap2.dipcall.vcf.gz \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/whole_genome_dipcall/dipcall_raw/dipcall_outfiles/HG002_verkko2.2.haplotype2.dipcall.GIABv4.2.1_benchmark_noinconsistent.bed \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/happy_GIABv4.2.1/happy \
    HG002

bash /private/home/mmastora/progs/scripts/GIAB_happy_chr20.sh \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_raw/dipcall_outfiles/HG002_verkko_2.2_chr20.hap2.dipcall.vcf.gz \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/whole_genome_dipcall/dipcall_raw/dipcall_outfiles/HG002_verkko2.2.haplotype2.dipcall.GIABv4.2.1_benchmark_noinconsistent.bed \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_raw/happy_GIABv4.2.1/happy \
    HG002
```

Run dipcall on whole genome raw
```
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/HG002_verkko2.2.haplotype1.maternal.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/HG002_verkko2.2.haplotype2.paternal.fasta",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
```

```

```
#!/bin/bash
#SBATCH --job-name=dipcall_truthset
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --exclude=phoenix-[09,10,22,23,24,18]

cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/whole_genome_dipcall/dipcall_raw

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24,18]"
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
```
bash /private/home/mmastora/progs/scripts/GIAB_happy.sh \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/whole_genome_dipcall/dipcall_raw/dipcall_outfiles/HG002_verkko2.2.haplotype2.dipcall.vcf.gz \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG002_y2_raw/dipCallTar/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall.GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/whole_genome_dipcall/dipcall_raw/happy/happy \
    HG002
```

Look at some FP FN calls in IGV

Convert FP / FN variants to bed
```
grep "^#" /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/happy/happy.vcf > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/happy/happy.FPFN.vcf

grep -v "^#" /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/happy/happy.vcf | grep ":F" >> /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/happy/happy.FPFN.vcf

export PATH=$PATH:/private/home/mmastora/progs/bin/
/private/home/mmastora/progs/bin/vcf2bed --do-not-split < /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/happy/happy.FPFN.vcf | awk '{print $1"\t"$2-10"\t"$3+10"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12}' > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/happy/happy.FPFN.vcf.10bp.bed
```

Project polished asm TP/FP/FN to raw assembly
```
# 3. Project GIAB confidence bedfile to y2 assembly, each haplotype separately
# Hap1
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_raw/dipcall_outfiles/HG002_verkko_2.2_chr20.hap1.paf \
  --blocks /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/happy/happy.FPFN.vcf.10bp.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/polished_FP_FN.projectable.hap1.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/polished_FP_FN.projection.hap1.bed

s
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_raw/dipcall_outfiles/HG002_verkko_2.2_chr20.hap2.paf \
  --blocks /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/dipcall_pol/happy/happy.FPFN.vcf.10bp.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/polished_FP_FN.projectable.hap2.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/polished_FP_FN.projection.hap2.bed

cat /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/polished_FP_FN.projection.hap2.bed /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/polished_FP_FN.projection.hap1.bed > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/polished_FP_FN.projection.dip.bed
```

```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint34/polished_FP_FN.projection.dip.bed -b /private/groups/patenlab/mira/phoenix_batch_executions/workflows/DeepPolisher/HG002_verkko_graph_alignment_chr20/analysis/DeepPolisher_outputs/polisher_output.no_filters.vcf.gz
```
#### Testing checkpoint 505

```
bcftools consensus -H2 -f /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_chr20.hap1.fasta /private/groups/patenlab/mira/phoenix_batch_executions/workflows/DeepPolisher/HG002_verkko_graph_alignment_chr20_checkpoint505/analysis/DeepPolisher_outputs/polisher_output.no_filters.vcf.gz > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint505/HG002_verkko_2.2_chr20.checkpoint505.hap1.fasta

bcftools consensus -H2 -f /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_chr20.hap2.fasta /private/groups/patenlab/mira/phoenix_batch_executions/workflows/DeepPolisher/HG002_verkko_graph_alignment_chr20_checkpoint505/analysis/DeepPolisher_outputs/polisher_output.no_filters.vcf.gz > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint505/HG002_verkko_2.2_chr20.checkpoint505.hap2.fasta
```
dipcall
```
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint505/HG002_verkko_2.2_chr20.checkpoint505.hap1.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint505/HG002_verkko_2.2_chr20.checkpoint505.hap2.fasta",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
```
#!/bin/bash
#SBATCH --job-name=dipcall_truthset
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --exclude=phoenix-[09,10,22,23,24,18]

cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint505/dipcall_pol

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24,18]"
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
run happy
```
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint505/dipcall_pol

bash /private/home/mmastora/progs/scripts/GIAB_happy_chr20.sh \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint505/dipcall_pol/dipcall_outfiles/HG002_verkko_2.2_chr20.checkpoint505.hap2.dipcall.vcf.gz \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG002_y2_raw/dipCallTar/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall.GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint505/dipcall_pol/happy \
    HG002

s
bash /private/home/mmastora/progs/scripts/GIAB_happy_chr20.sh \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint505/dipcall_pol/dipcall_outfiles/HG002_verkko_2.2_chr20.checkpoint505.hap2.dipcall.vcf.gz \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/whole_genome_dipcall/dipcall_raw/dipcall_outfiles/HG002_verkko2.2.haplotype2.dipcall.GIABv4.2.1_benchmark_noinconsistent.bed \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/test_chr20_checkpoint505/dipcall_pol_GIAB_4.2.1/happy \
    HG002

export PATH=$PATH:/private/home/mmastora/progs/bin/
/private/home/mmastora/progs/bin/vcf2bed --do-not-split < /private/groups/patenlab/mira/phoenix_batch_executions/workflows/DeepPolisher/HG002_verkko_graph_alignment_chr20_checkpoint505/analysis/DeepPolisher_outputs/polisher_output.no_filters.vcf > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/DeepPolisher/HG002_verkko_graph_alignment_chr20_checkpoint505/analysis/DeepPolisher_outputs/polisher_output.no_filters.vcf.bed


bedtools intersect -b /private/groups/patenlab/mira/phoenix_batch_executions/workflows/DeepPolisher/HG002_verkko_graph_alignment_chr20_checkpoint505/analysis/DeepPolisher_outputs/polisher_output.no_filters.vcf -a /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_graph_alignments_Q100v1.1_truthset/HG002_verkko_2.2_GIAB_confident_regions.diploid.projection.fix.bed | sort | uniq | wc -l
```
#### Downsample verkko bam to 20x

First check avg coverage with mosdepth
```
docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -n --fast-mode -t 4 /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/mosdepth/verkko_graph_alignment_HG002 \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_assembly.bam

source /private/home/mmastora/progs/miniconda3/etc/profile.d/conda.sh
conda activate analysis

python3 ~/progs/mosdepth/scripts/plot-dist.py verkko_graph_alignment_HG002.mosdepth.global.dist.txt &> coverages.txt

# 31 x total
```
Downsample
```
#!/bin/bash
#SBATCH --job-name=downsample_verkko
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=200gb
#SBATCH --threads-per-core=1
#SBATCH --time=12:00:00

samtools view -s .67 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_assembly.bam > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_assembly.20x.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/verkko_assembly.20x.bam
```

Get bedfile excluding the suspicious regions from the Q100 genome to test training DP with more examples.

Downloaded this file:
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.1/benchmark/resources/hg002v1.1.resources.tar.gz

Subtract suspicious bed from genome bed
```sh
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups \
    pegi3s/bedtools bedtools subtract \
    -a /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/resources/v1.1.genome.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/resources/v1.1.excluded_regions.bed \
    > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/Q100v1.1_minus_excluded_regions.genome.bed
```
Project bed file to verkko 2.2 assembly coordinates
```sh
# Project Q100 confidence bedfile to verkko assembly, each haplotype separately
# Hap1
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'asm2ref' \
  --paf /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_pat/dipcall_outfiles/hg002v1.1_X_EBV_MT.hap1.AP.paf \
  --blocks /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/Q100v1.1_minus_excluded_regions.genome.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/Q100_minus_excluded_regions/Q100_wo_excluded_regions_verkko2.2.pat.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/Q100_minus_excluded_regions/Q100_wo_excluded_regions_verkko2.2.pat.projection.bed

# Hap 2
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'asm2ref' \
  --paf /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/dipcall_mat/dipcall_outfiles/hg002v1.1_Y_EBV_MT.hap1.AP.paf \
  --blocks /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/Q100v1.1_minus_excluded_regions.genome.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/Q100_minus_excluded_regions/Q100_wo_excluded_regions_verkko2.2.mat.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/Q100_minus_excluded_regions/Q100_wo_excluded_regions_verkko2.2.mat.projection.bed
```
Merge projection beds
```
cat /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/Q100_minus_excluded_regions/Q100_wo_excluded_regions_verkko2.2.mat.projection.bed /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/Q100_minus_excluded_regions/Q100_wo_excluded_regions_verkko2.2.pat.projection.bed > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1/Q100_minus_excluded_regions/Q100_wo_excluded_regions_verkko2.2.dip.projection.bed
```

### Verkko truthset for updated version on 8/5/2025

### 1. Dipcall verkko assembly against truthset (T2T Q100 v1.1)

Sergey provided the assembly and alignments here: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=backup/Sergey/HG002/

Assembly is HiC, so need to rearrange mat and pat contigs

https://s3-us-west-2.amazonaws.com/human-pangenomics/backup/Sergey/HG002/verkko_assemby.trio.bed

Removing the contigs with the large switch from both bed files.
```
haplotype2-0000075	16	79111402	mat,pat
haplotype1-0000007	1362	80581086	mat,pat
```
```
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1_08052025

grep mat verkko_assemby.trio.bed | grep -v "haplotype2-0000075" | grep -v "haplotype1-0000007" | cut -f1 > verkko_assemby.mat_contigs.txt
grep pat verkko_assemby.trio.bed | grep -v "haplotype2-0000075" | grep -v "haplotype1-0000007" | cut -f1 > verkko_assemby.pat_contigs.txt

samtools faidx verkko_assembly.fasta -r verkko_assemby.mat_contigs.txt > HG002_verkko_08052025.maternal.fasta

samtools faidx verkko_assembly.fasta -r verkko_assemby.pat_contigs.txt > HG002_verkko_08052025.paternal.fasta
```
Set up input json files

```
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1_08052025/HG002_verkko_08052025.maternal.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1_08052025/HG002_verkko_08052025.maternal.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.mat_Y_EBV_MT.copy.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.mat_Y_EBV_MT.fasta",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}
```
```
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1_08052025/HG002_verkko_08052025.paternal.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1_08052025/HG002_verkko_08052025.paternal.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.pat_X_EBV_MT.copy.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.pat_X_EBV_MT.fasta",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1_08052025/dipcall_pat
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/verkko_graph_alignments_Q100v1.1_08052025/dipcall_mat

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24,18]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
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
Get bed file of regions with low element mappability:
```

```
