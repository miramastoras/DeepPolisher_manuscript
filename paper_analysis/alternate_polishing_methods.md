## Benchmarking DeepPolisher against alternate polishing methods

Alternate polishing approaches:

- [T2T polishing pipeline](https://github.com/arangrhie/T2T-Polish/blob/master/doc/T2T_polishing_case_study.md)
    35x hifi, 100x ilm, 120x ont
- [HiFi-DeepVariant](https://github.com/google/deepvariant)
- [NextPolish2](https://github.com/Nextomics/NextPolish2)

### 1. T2T polishing pipeline

#### 1.1 HG002

**Generate illumina alignments**

BWA index
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa index -p \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa
```

BWA mem alignment
```
#!/bin/bash
#SBATCH --job-name=HG002_ilm
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa mem -t32 \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R1.fastq.gz \
    /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R2.fastq.gz \
    | samtools view -b -h \
    > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/illumina_all_to_dip/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.HiSeq30x.bwa.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/illumina_all_to_dip/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.HiSeq30x.bwa.bam
```


**Generate R9 winnowmap all to dip alignments**

```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/flagger/wdls/workflows/long_read_aligner_scattered.wdl
```

```
{
  "longReadAlignmentScattered.preset": "map-ont",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "longReadAlignmentScattered.enableRunningSecphase": false,
  "longReadAlignmentScattered.aligner": "winnowmap",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -I8g -Y -L -p0.5",
  "longReadAlignmentScattered.kmerSize": 15,
  "longReadAlignmentScattered.sampleName": "HG002",
  "longReadAlignmentScattered.suffix": "trio_hifiasm_0.19.5.DC_1.2_40x.R941_Guppy5.winnowmap",
  "longReadAlignmentScattered.readFiles": ["/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/UL_R941_Guppy6/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.mm2.R941_Guppy6_UL.srt.bam"],
  "longReadAlignmentScattered.enableAddingMDTag": false
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/UL_R941_Guppy6/all_to_diploid

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=long"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
    /private/home/mmastora/progs/flagger/wdls/workflows/long_read_aligner_scattered.wdl \
    long_read_aligner_scattered_inputs.json \
    --outputDirectory ./long_read_aligner_scattered_outfiles \
    --outputFile long_read_aligner_scattered_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    --restart \
    2>&1 | tee log.txt
```

**HiFi winnowmap alignments**

Check coverage
```
#!/bin/bash
#SBATCH --job-name=HG002_mosdepth
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -n --fast-mode -t 4 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/mosdepth/HG002 \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam

python3 ~/progs/mosdepth/scripts/plot-dist.py /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/mosdepth/HG002.mosdepth.global.dist.txt
```

Confirmed 40x coverage

Merge HiFi and illumina bam files
```
#!/bin/bash
#SBATCH --job-name=merge_ilm_hifi_HG002
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=6:00:00


samtools merge -@32 \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hybrid_hifi_ilm/HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.Bwa.HiSeq_30x_hybrid.bam \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/illumina_all_to_dip/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.HiSeq30x.bwa.srt.bam
```

**Run SNV indel assembly wdl for variant calling**

```
{
  "snv_indel_assembly.sample": "HG002",
  "snv_indel_assembly.deepVariant_t.inputReadsIdx": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hybrid_hifi_ilm/HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.Bwa.HiSeq_30x_hybrid.bam.bai",
  "snv_indel_assembly.deepVariant_t.inputReads": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hybrid_hifi_ilm/HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.Bwa.HiSeq_30x_hybrid.bam",
  "snv_indel_assembly.assemblyIndex": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.fai",
  "snv_indel_assembly.pmdv_t.inputReadsIdx": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/UL_R941_Guppy6/all_to_diploid/long_read_aligner_scattered_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.R941_Guppy5.winnowmap.bam.bai",
  "snv_indel_assembly.assembly": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "snv_indel_assembly.deepVariant_t.dockerImage": "google/deepvariant:1.5.0",
  "snv_indel_assembly.pmdv_t.dockerImage": "kishwars/pepper_deepvariant:r0.8",
  "snv_indel_assembly.pmdv_t.inputReads": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/UL_R941_Guppy6/all_to_diploid/long_read_aligner_scattered_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.R941_Guppy5.winnowmap.bam"
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG002_y2_rerun_without_pharaoh_04052024

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=48:00:00 --partition=high_priority"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
    /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/workflows/snv_indel_assembly.wdl \
    snv_indel_assembly.inputs.json \
    --outputDirectory ./snv_indel_assembly.outfiles \
    --outputFile snv_indel_assembly.outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```

Filter with merfin
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/merfin.wdl
```

```
{
  "runMerfin.readmerDBTarball": "/private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y1_terra_tables/meryl_dbs/ilm.k21.meryl.tar",
  "runMerfin.Merfin.mode": "-polish",
  "runMerfin.Merfin.vcfFile": "/private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG002_y2_rerun_without_pharaoh_04052024/snv_indel_assembly.outfiles/HG002_MERGED_SMALL_VARIANTS.vcf.gz",
  "runMerfin.Merfin.refFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa"
}

```

```
#!/bin/bash
#SBATCH --job-name=merfin_HG002
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --threads-per-core=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=8
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

cd /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG002_y2_rerun_without_pharaoh_04052024/merfin

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
    /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/merfin.wdl \
    merfin.inputs.json \
    --outputDirectory ./merfin.outfiles \
    --outputFile merfin.outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```

#### 1.2 HG005

**Generate illumina alignments**

BWA index
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa index -p \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa
```

Bwa mem align
```
#!/bin/bash
#SBATCH --job-name=HG005_ilm
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa mem -t32 \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/reads/HG005.ilm.fastq.gz \
    | samtools view -b -h \
    > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_dip/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.bwa.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_dip/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.bwa.bam
```

Check coverage
```
#!/bin/bash
#SBATCH --job-name=HG005_mosdepth
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

samtools sort -@32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_dip/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.bwa.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_dip/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.bwa.srt.bam


docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_dip/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -n --fast-mode -t 4 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_dip/mosdepth/HG005 \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_dip/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.bwa.bam

python3 ~/progs/mosdepth/scripts/plot-dist.py /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/mosdepth/HG002.mosdepth.global.dist.txt
```
**Generate R9 ONT alignments**

```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/flagger/wdls/workflows/long_read_aligner_scattered.wdl
```

```
{
  "longReadAlignmentScattered.preset": "map-ont",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "longReadAlignmentScattered.enableRunningSecphase": false,
  "longReadAlignmentScattered.aligner": "winnowmap",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -I8g -Y -L -p0.5",
  "longReadAlignmentScattered.kmerSize": 15,
  "longReadAlignmentScattered.sampleName": "HG005",
  "longReadAlignmentScattered.suffix": "trio_hifiasm_0.19.5.DC_1.2_40x.R941_Guppy5.winnowmap",
  "longReadAlignmentScattered.readFiles": ["/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/ONT/R941_Guppy5/05_25_21_R941_GM24631_3X_Guppy_5.0.7_prom_sup.fastq.gz","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/ONT/R941_Guppy5/05_25_21_R941_GM24631_9X_2_Guppy_5.0.7_prom_sup.fastq.gz","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/ONT/R941_Guppy5/05_25_21_R941_GM24631_9X_3_Guppy_5.0.7_prom_sup.fastq.gz","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/ONT/R941_Guppy5/05_25_21_R941_GM24631_9X_Guppy_5.0.7_prom_sup.fastq.gz"],
  "longReadAlignmentScattered.enableAddingMDTag": false
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/UL_R941_Guppy6/all_to_diploid

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=24:00:00 --partition=long"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
    /private/home/mmastora/progs/flagger/wdls/workflows/long_read_aligner_scattered2.wdl \
    long_read_aligner_scattered_inputs.json \
    --outputDirectory ./long_read_aligner_scattered_outfiles \
    --outputFile long_read_aligner_scattered_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```


**HiFi winnowmap alignments**


Merge HiFi and illumina bam files
```
#!/bin/bash
#SBATCH --job-name=merge_ilm_hifi_HG005
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=6:00:00


samtools merge -@32 \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/hybrid_hifi_ilm/HG005.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.Bwa.HiSeq_20x_hybrid.bam /private/groups/patenlab/masri/hprc/polishing/HG005/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/DC_1.2_alignments/diploid/slurm_run/HG005.diploid.trio_hifiasm_0.19.5.DC_1.2.winnowmap_2.03.DC_1.2_40x.bam \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_dip/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.bwa.srt.bam
```
**Run SNV indel assembly wdl for variant calling**

```
{
  "snv_indel_assembly.sample": "HG005",
  "snv_indel_assembly.deepVariant_t.inputReadsIdx": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/hybrid_hifi_ilm/HG005.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.Bwa.HiSeq_20x_hybrid.bam.bai",
  "snv_indel_assembly.deepVariant_t.inputReads": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/hybrid_hifi_ilm/HG005.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.Bwa.HiSeq_20x_hybrid.bam",
  "snv_indel_assembly.assemblyIndex": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.fai",
  "snv_indel_assembly.pmdv_t.inputReadsIdx": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/UL_R941_Guppy6/all_to_diploid/long_read_aligner_scattered_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.R941_Guppy5.winnowmap.bam.bai",
  "snv_indel_assembly.assembly": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "snv_indel_assembly.deepVariant_t.dockerImage": "google/deepvariant:1.5.0",
  "snv_indel_assembly.pmdv_t.dockerImage": "kishwars/pepper_deepvariant:r0.8",
  "snv_indel_assembly.pmdv_t.inputReads": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/UL_R941_Guppy6/all_to_diploid/long_read_aligner_scattered_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.R941_Guppy5.winnowmap.bam"
}
```



```
#!/bin/bash            
#SBATCH --job-name=HG005_t2t_polish
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --time=72:00:00
#SBATCH --partition=long
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=slurm_logs/submission_%x_%j_%A_%a.log


cd /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=72:00:00 --partition=long"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
    /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/workflows/snv_indel_assembly.wdl \
    snv_indel_assembly.inputs.json \
    --outputDirectory ./snv_indel_assembly.outfiles \
    --outputFile snv_indel_assembly.outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    --restart \
    2>&1 | tee log.txt
```

Was using wrong flow cells for this initial run, so re-ran just hybrid DCv1.2 variant calling   
```
#!/bin/bash
#SBATCH --job-name=HG005_DV_hybrid
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00

docker run --rm -u `id -u`:`id -g` \
-v /private/groups:/private/groups \
google/deepvariant:1.5.0 \
/opt/deepvariant/bin/run_deepvariant \
--model_type HYBRID_PACBIO_ILLUMINA \
--ref /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
--reads /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/hybrid_hifi_ilm/HG005.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.Bwa.HiSeq_20x_hybrid.bam \
--output_vcf /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/HG005_y2_hybrid.variants.vcf.gz \
--num_shards 32 \
--intermediate_results_dir /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/intermediate_results_dir
```

Merge variant callsets
```
# filter new hybrid callset
bcftools view -f "PASS" -e 'FORMAT/VAF<=0.5 | FORMAT/GQ<=30' -Oz /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/HG005_y2_hybrid.variants.vcf.gz > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/HG005_y2_hybrid.variants.filtered.vcf.gz

# merge with PMDV callset

docker run --rm -u `id -u`:`id -g` \
-v /private/groups/patenlab:/private/groups/patenlab \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
/private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/HG005_y2_hybrid.variants.filtered.vcf.gz \
/private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/snv_indel_assembly.outfiles/HG005_PEPPER_DeepVariant.filtered.vcf.gz \
-r /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
-o /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/hybrid_ont_happy \
--pass-only \
--engine=vcfeval \
--threads=16

docker run --rm -u `id -u`:`id -g` \
-v /private/groups/patenlab:/private/groups/patenlab \
kishwars/t2t_polishing:0.1 \
python3 vcf_merge_t2t.py \
-v1 /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/HG005_y2_hybrid.variants.filtered.vcf.gz \
-v2 /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/snv_indel_assembly.outfiles/HG005_PEPPER_DeepVariant.filtered.vcf.gz \
-hv /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/hybrid_ont_happy.vcf.gz \
-o /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/HG005_MERGED_SMALL_VARIANTS.vcf.gz
```

Downsample HG005 ilm  bamfile for meryl and yak dbs for merfin
```
#!/bin/bash
#SBATCH --job-name=downsampling_HG005_ilm
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

# 50x to 30x
samtools view -s 0.6 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.50x.all_to_mat.trio_hifiasm_0.19.5.DC_1.2_40x.srt.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.50x.all_to_mat.trio_hifiasm_0.19.5.DC_1.2_40x.srt.downsampled_30x.bam
```
```
#!/bin/bash
#SBATCH --job-name=samtools_fastq_HG5_ilm
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00

samtools fastq -@32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.50x.all_to_mat.trio_hifiasm_0.19.5.DC_1.2_40x.srt.downsampled_30x.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.fastq
```
Rebuild meryl for merfin
```
#!/bin/bash
#SBATCH --job-name=meryl_HG5
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

export PATH=/private/home/mmastora/progs/meryl-1.4.1/bin:$PATH

meryl count threads=32 k=21 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.fastq output /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k21.30x.meryl

tar -cvf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k21.30x.meryl.tar /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k21.30x.meryl

meryl count threads=32 k=31 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.fastq output /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k31.30x.meryl
```

Rebuild yak

```
#!/bin/bash
#SBATCH --job-name=yak_HG5
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

time docker run --rm -u `id -u`:`id -g` -v /private/groups:/private/groups \
    juklucas/hpp_yak:latest yak count \
    -k21 -b37 -t32 \
    -o /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/yak_files/HG005.k21.30x.yak /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.fastq

time docker run --rm -u `id -u`:`id -g` -v /private/groups:/private/groups \
    juklucas/hpp_yak:latest yak count \
    -k31 -b37 -t32 \
    -o /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/yak_files/HG005.k31.30x.yak /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.fastq
```

Run genomescope for merfin
```
cd /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/merfin_manual

docker run --rm -u `id -u`:`id -g` \
  -v "/private/groups/patenlab/mira":"/private/groups/patenlab/mira" \
  juklucas/hpp_merqury:latest \
  meryl histogram /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k21.30x.meryl > meryl.hist

docker run -it --rm -v /private/groups:/private/groups dmolik/genomescope2:latest xvfb-run /home/genomics/genomescope2.0/genomescope.R -i /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/merfin_manual/meryl.hist -k 21 -o /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/merfin_manual/genomescope_outfiles -p 1 --fitted_hist &> genomescope.stdout

    GenomeScope analyzing /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/merfin_manual/meryl.hist p=1 k=21 outdir=genomescope_outfiles


    a:100%
    Model converged het:0 kcov:44.9 err:0.00548 model fit:1.59 len:2878576465
```


Run Merfin
```
docker run -it --rm -v /private/groups:/private/groups \
    miramastoras/merfin:latest \
    merfin -polish  \
    -vcf /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/HG005_MERGED_SMALL_VARIANTS.vcf.gz \
    -threads 32 \
    -sequence /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    -readmers /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k21.30x.meryl \
    -prob /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/merfin_manual/genomescope_outfiles/lookup_table.txt \
    -peak 44.9 \
    -output /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/HG005_MERGED_SMALL_VARIANTS.filtered.k21.merfin
```

GIAB analysis here: https://github.com/miramastoras/phoenix_batch_submissions/blob/main/polishing/applyPolish_dipcall_happy/GIAB_samples_manuscript/applyPolish_dipcall_input_jsons/HG005_t2t_polish_applyPolish_dipcall.json

QV pipeline here: https://github.com/miramastoras/phoenix_batch_submissions/blob/main/polishing/hprc_polishing_QC_no_meryl/GIAB_samples_DP_manuscript/hprc_polishing_QC_no_meryl_input_jsons/HG005_t2t_polish_hprc_polishing_QC_no_meryl.json

### 2. DeepVariant polishing on winnowmap HiFi alignments

HG002

```
#!/bin/bash
#SBATCH --job-name=dv-hifi_HG002
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

BIN_VERSION="1.6.1"
docker run -u `id -u`:`id -g` \
  -v "/private/groups/patenlab/mira":"/private/groups/patenlab/mira" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
  --reads=/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam \
  --output_vcf=/private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/HG002_winnowmap_DV/HG002.dip.y2_hifiasm_0.19.5.DC_1.2.HiFi_winnowmap.deepvariant.vcf.gz \
  --output_gvcf=/private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/HG002_winnowmap_DV/HG002.dip.y2_hifiasm_0.19.5.DC_1.2.HiFi_winnowmap.deepvariant.gvcf.gz \
  --num_shards=32
```

HG005

```
#!/bin/bash
#SBATCH --job-name=dv-hifi_HG005
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

BIN_VERSION="1.6.1"
docker run -u `id -u`:`id -g` \
  -v "/private/groups/":"/private/groups/" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
  --reads=/private/groups/patenlab/masri/hprc/polishing/HG005/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/DC_1.2_alignments/diploid/slurm_run/HG005.diploid.trio_hifiasm_0.19.5.DC_1.2.winnowmap_2.03.DC_1.2_40x.bam \
  --output_vcf=/private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/HG005_winnowmap_DV/HG005.dip.y2_hifiasm_0.19.5.DC_1.2.40x.HiFi_winnowmap.deepvariant.vcf.gz \
  --output_gvcf=/private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/HG005_winnowmap_DV/HG005.dip.y2_hifiasm_0.19.5.DC_1.2.40x.HiFi_winnowmap.deepvariant.gvcf.gz \
  --num_shards=32
```

Filter to PASS only
```
bcftools view -Oz -f "PASS" /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/HG005_winnowmap_DV/HG005.dip.y2_hifiasm_0.19.5.DC_1.2.40x.HiFi_winnowmap.deepvariant.vcf.gz > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/HG005_winnowmap_DV/HG005.dip.y2_hifiasm_0.19.5.DC_1.2.40x.HiFi_winnowmap.deepvariant.PASS.vcf.gz

bcftools view -Oz -f "PASS" /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/HG002_winnowmap_DV/HG002.dip.y2_hifiasm_0.19.5.DC_1.2.HiFi_winnowmap.deepvariant.vcf.gz > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/HG002_winnowmap_DV/HG002.dip.y2_hifiasm_0.19.5.DC_1.2.HiFi_winnowmap.deepvariant.PASS.vcf.gz
```
Get count of homalt polishing variants
```
zcat /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/HG005_winnowmap_DV/HG005.dip.y2_hifiasm_0.19.5.DC_1.2.40x.HiFi_winnowmap.deepvariant.PASS.vcf.gz | grep -v "^#" | grep "1/1" | wc -l
# 17210

zcat /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/HG002_winnowmap_DV/HG002.dip.y2_hifiasm_0.19.5.DC_1.2.HiFi_winnowmap.deepvariant.PASS.vcf.gz | grep -v "^#" | grep "1/1" | wc -l

```

polishing assembly and GIAB variants analysis: https://github.com/miramastoras/phoenix_batch_submissions/blob/main/polishing/applyPolish_dipcall_happy/GIAB_samples_manuscript/applyPolish_dipcall_input_jsons/HG005_deepvariant_applyPolish_dipcall.json

QV pipeline: https://github.com/miramastoras/phoenix_batch_submissions/blob/main/polishing/hprc_polishing_QC_no_meryl/GIAB_samples_DP_manuscript/hprc_polishing_QC_no_meryl_input_jsons/HG005_deepvariant_hprc_polishing_QC_no_meryl.json

### 3. NextPolish2 polishing on winnowmap HiFi alignments

Following instructions from https://github.com/Nextomics/NextPolish2/blob/main/doc/benchmark5.md

#### HG005
```
srun \
  --job-name "yak" \
  --cpus-per-task 32 \
  --partition medium \
  --mem 250G \
  --time 12:00:00 \
  --pty bash \
  -i

# Quality control and filtering for illumina
fastp -t 32 -f 5 -t 5 --cut_front --cut_tail -i /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.fastq -o /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.clean.fastq

Read1 before filtering:
total reads: 377058477
total bases: 93790588803
Q20 bases: 85580716657(91.2466%)
Q30 bases: 80629417261(85.9675%)

Read1 after filtering:
total reads: 376644600
total bases: 86827705843
Q20 bases: 80120119365(92.2748%)
Q30 bases: 75896593201(87.4106%)

Filtering result:
reads passed filter: 376644600
reads failed due to low quality: 230924
reads failed due to too many N: 110939
reads failed due to too short: 72014
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
```

Make yak files
```
#!/bin/bash
#SBATCH --job-name=yak_HG5
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

time docker run --rm -u `id -u`:`id -g` -v /private/groups:/private/groups \
    juklucas/hpp_yak:latest yak count \
    -k21 -b37 -t32 \
    -o /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.clean.k21.yak /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.clean.fastq

time docker run --rm -u `id -u`:`id -g` -v /private/groups:/private/groups \
    juklucas/hpp_yak:latest yak count \
    -k31 -b37 -t32 \
    -o /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.clean.k31.yak /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.clean.fastq
```

Run nextpolish2 on HG005
```
#!/bin/bash
#SBATCH --job-name=NP_HG5
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

source /private/home/mmastora/progs/miniconda3/etc/profile.d/conda.sh
conda activate nextpolish

time nextPolish2 -t 32 /private/groups/patenlab/masri/hprc/polishing/HG005/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/DC_1.2_alignments/diploid/slurm_run/HG005.diploid.trio_hifiasm_0.19.5.DC_1.2.winnowmap_2.03.DC_1.2_40x.bam /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.gz /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.clean.k21.yak /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/maternal/HG005.ilm.30x.clean.k31.yak > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG005_y2_winnowmap_HiFi/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.NP2.clean.polished.fa
```
Separate by Haplotype
```
cat /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG005_y2_winnowmap_HiFi/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.NP2.clean.polished.fa | seqkit grep -r -p '.*h1tg' > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG005_y2_winnowmap_HiFi/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.NP2.clean.polished.fa

cat /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG005_y2_winnowmap_HiFi/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.NP2.clean.polished.fa | seqkit grep -r -p '.*h2tg' > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG005_y2_winnowmap_HiFi/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.NP2.clean.polished.fa
```

Manually run dipcall and happy evaluation

`dipcall_inputs.json`:

```json
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG005_y2_winnowmap_HiFi/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.NP2.clean.polished.fa",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG005_y2_winnowmap_HiFi/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.NP2.clean.polished.fa",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_clean_nextPolish2

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24]"
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

Run happy
```
bash /private/home/mmastora/progs/scripts/GIAB_happy.sh /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_clean_nextPolish2/dipcall_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.NP2.clean.polished.dipcall.vcf.gz /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall.GIABv4.2.1.confidence.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_clean_nextPolish2/happy_out HG005
```

#### HG002


```
# Download NGS data
srun \
  --job-name "clean" \
  --cpus-per-task 32 \
  --partition medium \
  --mem 250G \
  --time 2:00:00 \
  --pty bash \
  -i

# Quality control and filtering for illumina
fastp -t 32 -f 5 --cut_front --cut_tail -i /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R1.fastq.gz -I /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R2.fastq.gz -o /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R1.clean.fastq.gz -O /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R2.clean.fastq.gz
```

Make yak files
```
#!/bin/bash
#SBATCH --job-name=yak_HG2
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

time docker run --rm -u `id -u`:`id -g` -v /private/groups:/private/groups \
    juklucas/hpp_yak:latest yak count \
    -k21 -b37 -t32 \
    -o /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/yak_files/HG002.ilm.30x.clean.k21.yak /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R1.clean.fastq.gz /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R2.clean.fastq.gz

time docker run --rm -u `id -u`:`id -g` -v /private/groups:/private/groups \
    juklucas/hpp_yak:latest yak count \
    -k31 -b37 -t32 \
    -o /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/yak_files/HG002.ilm.30x.clean.k31.yak /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R1.clean.fastq.gz /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R2.clean.fastq.gz
```
Run nextpolish2 on HG005
```
#!/bin/bash
#SBATCH --job-name=NP_HG2
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=300gb
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

conda activate nextpolish

time nextPolish2 -t 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.gz /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/yak_files/HG002.ilm.30x.clean.k21.yak /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/yak_files/HG002.ilm.30x.clean.k31.yak > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG002_y2_winnowmap_HiFi/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.NP2.clean.polished.fa
```
Separate by Haplotype
```
cat /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG002_y2_winnowmap_HiFi/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.NP2.clean.polished.fa | seqkit grep -r -p '.*h1tg' > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG002_y2_winnowmap_HiFi/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.NP2.clean.polished.fa

cat /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG002_y2_winnowmap_HiFi/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.NP2.clean.polished.fa | seqkit grep -r -p '.*h2tg' > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG002_y2_winnowmap_HiFi/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.NP2.clean.polished.fa
```

```json
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG002_y2_winnowmap_HiFi/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.NP2.clean.polished.fa",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/nextpolish2/HG002_y2_winnowmap_HiFi/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.NP2.clean.polished.fa",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_clean_nextPolish2

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24]"
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

Run happy chr 20
```
bash /private/home/mmastora/progs/scripts/GIAB_happy_chr20.sh /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_clean_nextPolish2/dipcall_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.NP2.clean.polished.dipcall.vcf.gz /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG002_y2_raw/dipCallTar/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall.GIAB_T2T_Q100_conf_beds_concordant_50bp.dipcall_z2k.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_clean_nextPolish2/happy_out HG002
```

### Running deepvariant on pharaoh alignments for hprc samples

DeepVariant:
https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/deepvariant/hprc_manuscript

Polishing:
https://github.com/miramastoras/phoenix_batch_submissions/blob/main/polishing/deepvariant/hprc_manuscript/run_deepvariant.sh

QV pipeline: https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/hprc_polishing_QC/hprc_dv_polishing

### Polishing with DV-hybrid model without filters for HG005 and HG002

```
# HG002
bcftools consensus -H1 -f /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG002_y2_rerun_without_pharaoh_04052024/snv_indel_assembly.outfiles/HG002_deepvariant.vcf.gz > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/DV_hybrid_only/HG002/HG002_y2_DV_hybrid_polished.pat.fa
# 11574

bcftools consensus -H1 -f /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG002_y2_rerun_without_pharaoh_04052024/snv_indel_assembly.outfiles/HG002_deepvariant.vcf.gz > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/DV_hybrid_only/HG002/HG002_y2_DV_hybrid_polished.mat.fa
# 12006

# HG005
bcftools consensus -H1 -f /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/HG005_y2_hybrid.variants.vcf.gz > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/DV_hybrid_only/HG005/HG005_Y2_DV_hybrid_polished.pat.fa
# 9593

bcftools consensus -H1 -f /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/t2t_polish/HG005_y2/variant_calling_manually/HG005_y2_hybrid.variants.vcf.gz > /private/groups/patenlab/mira/hprc_polishing/y2_alt_polishers/DV_hybrid_only/HG005/HG005_Y2_DV_hybrid_polished.mat.fa
# 10288
```

```
bash ~/progs/scripts/GIAB_happy_chr20_21_22.sh /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_deepvariant/applyPolish_dipcall_outputs/HG005_deepvariant_hap1.polished.dipcall.vcf.gz /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall.GIABv4.2.1.confidence.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_deepvariant/happy_outputs_chr20_22/happy_outs HG005
```

```
bash /private/home/mmastora/progs/scripts/GIAB_happy_stratifications.sh /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.dipcall.vcf.gz /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall.GIABv4.2.1.confidence.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/happy_stratifications HG005


bash /private/home/mmastora/progs/scripts/GIAB_happy_stratifications.sh /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.vcf.gz /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall.GIABv4.2.1.confidence.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/happy_stratifications HG005
```
