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
