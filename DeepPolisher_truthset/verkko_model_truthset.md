## Creating a truthset for training a DeepPolisher hifi verkko model

### 1. Dipcall verkko assembly against truthset (T2T Q100 v0.9)

```
#!/bin/bash
#SBATCH --job-name=dipcall_truthset
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
    humanpangenomics/hpp_dipcall_v0.3:latest /opt/dipcall/dipcall.kit/run-dipcall \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_pat \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_pat_verkko_2.0_2024_06_05.fa \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.pat_X_EBV_MT.fasta \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.pat_X_EBV_MT.fasta > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_pat/HG002_t2t_v0.9_2_hprc_pat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
        humanpangenomics/hpp_dipcall_v0.3:latest make -j2 -f /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_pat/HG002_t2t_v0.9_2_hprc_pat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
    humanpangenomics/hpp_dipcall_v0.3:latest /opt/dipcall/dipcall.kit/run-dipcall \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_mat \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_mat_verkko_2.0_2024_06_05.fa \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.mat_Y_EBV_MT.fasta \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.mat_Y_EBV_MT.fasta > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_mat/HG002_t2t_v0.9_2_hprc_mat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
        humanpangenomics/hpp_dipcall_v0.3:latest make -j2 -f /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_mat/HG002_t2t_v0.9_2_hprc_mat.mak
```

### 2. Project GIAB confidence regions to verkko assembly

Dipcall raw assembly against GRCh38
```

{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_mat_verkko_2.0_2024_06_05.fa",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_pat_verkko_2.0_2024_06_05.fa",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem single_machine \
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
cd /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall/dipcall_outfiles/HG002_verkko_2.0_2024_06_05.dipcall

grep "tp:A:P" HG002_verkko_2.0_2024_06_05.hap1.paf > HG002_verkko_2.0_2024_06_05.hap1.AP.paf
grep "tp:A:P" HG002_verkko_2.0_2024_06_05.hap2.paf > HG002_verkko_2.0_2024_06_05.hap2.AP.paf

# 3. Project GIAB confidence bedfile to y2 assembly, each haplotype separately
# Hap1
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall/dipcall_outfiles/HG002_verkko_2.0_2024_06_05.dipcall/HG002_verkko_2.0_2024_06_05.hap1.AP.paf \
  --blocks /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap1_GIAB_conf.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap1_GIAB_conf.projection.bed

# Hap 2
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall/dipcall_outfiles/HG002_verkko_2.0_2024_06_05.dipcall/HG002_verkko_2.0_2024_06_05.hap2.AP.paf \
  --blocks /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap2_GIAB_conf.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap2_GIAB_conf.projection.bed
```

Uploaded for Kishwar here:
```
gsutil -o GSUtil:parallel_composite_upload_threshold=50GB -m cp -r verkko_model_truthset_HG002_v0.9 gs://pepper-deepvariant/mira/

gs://pepper-deepvariant/mira/verkko_model_truthset_HG002_v0.9
```

Run GIAB numbers for unpolished assembly
```
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall.bed -b /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG002_y2_raw/dipCallTar/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall.GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/happy/dipcall_bed_isec_conf.bed

bash /private/home/mmastora/progs/scripts/GIAB_happy.sh \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall.vcf.gz \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/happy/dipcall_bed_isec_conf.bed \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/happy/verkko_ont_porec_happy_out \
HG002

bash /private/home/mmastora/progs/scripts/GIAB_happy_chr20.sh \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall.vcf.gz \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/happy/dipcall_bed_isec_conf.bed \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/happy/verkko_ont_porec_happy_out \
HG002
```


## Creating a truthset for training a DeepPolisher ONT verkko model

### 0. Re-arrange haplotypes so that all PAT contigs are in one fasta and all MAT contigs are in the other.

Align ONT assembly to Q100 v0.9
```
#!/bin/bash
#SBATCH --job-name=minimap2_verkko_ONT_asm_Q100v0.9
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

docker run -u `id -u`:`id -g` -v /private:/private mobinasri/long_read_aligner:v0.4.0 minimap2 -ax asm5 -I8g -t32 /private/groups/patenlab/mira/data/hg002v0.9.fasta /private/groups/migalab/kkyriaki/experiments/verkko_assemblies/VERKKO_LC2024_POREC/assembly.fasta  | samtools view -bh - > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/VERKKO_LC2024_POREC.Q100v0.9.dip.bam
```

Subset to primary alignments only
```
samtools view -F0x900 /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/VERKKO_LC2024_POREC.Q100v0.9.dip.bam -o /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/VERKKO_LC2024_POREC.Q100v0.9.dip.primary.bam
```

Convert to bed file
```
bedtools bamtobed -i /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/VERKKO_LC2024_POREC.Q100v0.9.dip.primary.bam > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/VERKKO_LC2024_POREC.Q100v0.9.dip.primary.bam.bed

#
awk '{print $4 , $1}' /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/VERKKO_LC2024_POREC.Q100v0.9.dip.primary.bam.bed | sort | uniq > contig_maps

 | wc -l
grep MATERNAL contig_maps | wc -l
```

Re-partition fasta files according to contig mapping:
```
grep PATERNAL contig_maps | cut -f1 -d" " | while read line ; do samtools faidx /private/groups/patenlab/mira/ONT_DeepPolisher/assemblies/verkko_ont_porec/assembly.fasta $line >> /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/assembly.pat.fa ; done

grep MATERNAL contig_maps | cut -f1 -d" "  | cut -f1 -d"," | while read line ; do samtools faidx /private/groups/patenlab/mira/ONT_DeepPolisher/assemblies/verkko_ont_porec/assembly.fasta $line >> /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/assembly.mat.fa ; done
```


### 1. Dipcall verkko ONT-only + POREC assembly against truthset (T2T Q100 v0.9)

```
#!/bin/bash
#SBATCH --job-name=dipcall_truthset
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
    humanpangenomics/hpp_dipcall_v0.3:latest /opt/dipcall/dipcall.kit/run-dipcall \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall/dipcall_verkko_porec_pat \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/assembly.pat.fa \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.pat_X_EBV_MT.fasta \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.pat_X_EBV_MT.fasta > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall/dipcall_verkko_porec_pat/HG002_t2t_v0.9_2_verkko_porec_pat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
        humanpangenomics/hpp_dipcall_v0.3:latest make -j2 -f /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall/dipcall_verkko_porec_pat/HG002_t2t_v0.9_2_verkko_porec_pat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
    humanpangenomics/hpp_dipcall_v0.3:latest /opt/dipcall/dipcall.kit/run-dipcall \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall/dipcall_verkko_porec_mat \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/assembly.mat.fa \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.mat_Y_EBV_MT.fasta \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.mat_Y_EBV_MT.fasta > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall/dipcall_verkko_porec_mat/HG002_t2t_v0.9_2_verkko_porec_mat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
        humanpangenomics/hpp_dipcall_v0.3:latest make -j2 -f /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall/dipcall_verkko_porec_mat/HG002_t2t_v0.9_2_verkko_porec_mat.mak
```
### 2. Project GIAB confidence regions to verkko assembly

Dipcall raw assembly against GRCh38
```

{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/assembly.mat.fa",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/assembly.pat.fa",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem single_machine \
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
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall

grep "tp:A:P" assembly.hap1.paf > ../assembly.hap1.AP.paf
grep "tp:A:P" assembly.hap2.paf > ../assembly.hap2.AP.paf

# 3. Project GIAB confidence bedfile to y2 assembly, each haplotype separately
# Hap1
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.hap1.AP.paf \
  --blocks /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/HG002_verkko_porec_GIAB_conf.hap1.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/HG002_verkko_porec_GIAB_conf.hap1.projection.bed

# Hap 2
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.hap2.AP.paf \
  --blocks /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/HG002_verkko_porec_GIAB_conf.hap2.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/HG002_verkko_porec_GIAB_conf.hap2.projection.bed
```

Pre-processing truthset files:
```
# merge alignment to truth bams:
samtools merge -o /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/ont_model_truthset_verkko_porec_Q100v0.9/dipcall_verkko_porec.dipcall.Q100v0.9.merged.bam /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall/dipcall_verkko_porec_pat.hap1.bam /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall/dipcall_verkko_porec_mat.hap2.bam

# merge confidence projection bed files
cat /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/HG002_verkko_porec_GIAB_conf.hap1.projection.bed /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/HG002_verkko_porec_GIAB_conf.hap2.projection.bed > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall/HG002_verkko_porec_GIAB_conf.dip.merged.projection.bed

# filter for primary alignments to grch38

cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall
samtools merge -o assembly.dip.merged.bam assembly.hap1.bam assembly.hap2.bam

samtools view -F0x900 /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall/assembly.dip.merged.bam -o /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall/assembly.dip.merged.primary.bam

# convert to bed
bedtools bamtobed -i assembly.dip.merged.primary.bam > assembly.dip.merged.primary.bam.bed

awk '{print $4 , $1}' assembly.dip.merged.primary.bam.bed | sort | uniq | wc -l
awk '{print $4 , $1}' assembly.dip.merged.primary.bam.bed > chrom_maps

# separate fasta into chr 1-19 and chr 20 based on maps
grep -v "chr21" chrom_maps | grep -v "chr22" | grep -v "chr20" | grep -v "chrUn" > chrom_maps_1_19
grep "chr20" chrom_maps > chrom_maps_20

cut -f1 -d " " /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall/chrom_maps_1_19 | while read line ; do samtools faidx /private/groups/patenlab/mira/ONT_DeepPolisher/assemblies/verkko_ont_porec/assembly.fasta $line >> /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/ont_model_truthset_verkko_porec_Q100v0.9/assembly_verkko_ont_porec_chr1_19.fasta ; done

cut -f1 -d " " /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall/chrom_maps_20 | while read line ; do samtools faidx /private/groups/patenlab/mira/ONT_DeepPolisher/assemblies/verkko_ont_porec/assembly.fasta $line >> /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/ont_model_truthset_verkko_porec_Q100v0.9/assembly_verkko_ont_porec_chr20.fasta ; done
```

Uploaded for Kishwar here:
```
gsutil -o GSUtil:parallel_composite_upload_threshold=50GB -m cp -r ont_model_truthset_verkko_porec_Q100v0.9 gs://pepper-deepvariant/mira/
```

Run GIAB numbers for unpolished assembly
```
cd /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall.bed -b /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG002_y2_raw/dipCallTar/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall.GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/happy/dipcall_bed_isec_conf.bed

bash /private/home/mmastora/progs/scripts/GIAB_happy.sh \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall.vcf.gz \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/happy/dipcall_bed_isec_conf.bed \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/happy/verkko_ont_porec_happy_out \
HG002

bash /private/home/mmastora/progs/scripts/GIAB_happy_chr20.sh \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/dipcall_outfiles/assembly.dipcall.vcf.gz \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/happy/dipcall_bed_isec_conf.bed \
/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/dipcall_grch38/happy/verkko_ont_porec_happy_out \
HG002
```

## Creating a truthset for training a DeepPolisher hifi verkko model based off of latest Q100 version 1.1

### 1. Dipcall verkko assembly against truthset (T2T Q100 v0.9)

```
#!/bin/bash
#SBATCH --job-name=dipcall_truthset
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
    humanpangenomics/hpp_dipcall_v0.3:latest /opt/dipcall/dipcall.kit/run-dipcall \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model_v1.1/dipcall_verkkov2.0_pat \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_pat_verkko_2.0_2024_06_05.fa \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.pat_X_EBV_MT.fasta \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.pat_X_EBV_MT.fasta > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model_v1.1/dipcall_verkkov2.0_pat/HG002_t2t_v1.1_hprc_pat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
        humanpangenomics/hpp_dipcall_v0.3:latest make -j2 -f /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model_v1.1/dipcall_verkkov2.0_pat/HG002_t2t_v1.1_hprc_pat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
    humanpangenomics/hpp_dipcall_v0.3:latest /opt/dipcall/dipcall.kit/run-dipcall \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model_v1.1/dipcall_verkkov2.0_mat \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_mat_verkko_2.0_2024_06_05.fa \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.mat_Y_EBV_MT.fasta \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/truth/hg002v1.1.mat_Y_EBV_MT.fasta > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_mat/HG002_t2t_v1.1_hprc_mat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
        humanpangenomics/hpp_dipcall_v0.3:latest make -j2 -f /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_mat/HG002_t2t_v1.1_hprc_mat.mak
```

### 2. Project GIAB confidence regions to verkko assembly

Dipcall raw assembly against GRCh38
```

{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_mat_verkko_2.0_2024_06_05.fa",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_pat_verkko_2.0_2024_06_05.fa",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem single_machine \
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
cd /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall/dipcall_outfiles/HG002_verkko_2.0_2024_06_05.dipcall

grep "tp:A:P" HG002_verkko_2.0_2024_06_05.hap1.paf > HG002_verkko_2.0_2024_06_05.hap1.AP.paf
grep "tp:A:P" HG002_verkko_2.0_2024_06_05.hap2.paf > HG002_verkko_2.0_2024_06_05.hap2.AP.paf

# 3. Project GIAB confidence bedfile to y2 assembly, each haplotype separately
# Hap1
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall/dipcall_outfiles/HG002_verkko_2.0_2024_06_05.dipcall/HG002_verkko_2.0_2024_06_05.hap1.AP.paf \
  --blocks /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap1_GIAB_conf.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap1_GIAB_conf.projection.bed

# Hap 2
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall/dipcall_outfiles/HG002_verkko_2.0_2024_06_05.dipcall/HG002_verkko_2.0_2024_06_05.hap2.AP.paf \
  --blocks /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap2_GIAB_conf.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap2_GIAB_conf.projection.bed
```

Uploaded for Kishwar here:
```
gsutil -o GSUtil:parallel_composite_upload_threshold=50GB -m cp -r verkko_model_truthset_HG002_Q100v1.1 gs://pepper-deepvariant/mira/
```
