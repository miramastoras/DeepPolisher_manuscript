### Benchmarking performance of PHARAOH pipeline

Running DeepPolisher on minimap2 alignments

HG002
```
{
  "runDeepPolisher.sampleName": "HG002",
  "runDeepPolisher.ModelFilesTarGZ": "/private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-665.tar.gz",
  "runDeepPolisher.Bai": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/toil_hifi_dv1.2_mm2_all2dip_out/HG002.DCv1.2.minimap2v2.26.bam.bai",
  "runDeepPolisher.useOptimalGQFilter": true,
  "runDeepPolisher.Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "runDeepPolisher.dockerImage": "google/deepconsensus:polisher_v0.0.8_12122023",
  "runDeepPolisher.Bam": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/toil_hifi_dv1.2_mm2_all2dip_out/HG002.DCv1.2.minimap2v2.26.bam"
}
```
Remove reads with de > 0.02
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl
```

```
{
  "runCorrectBam.correctBam.dockerImage": "mobinasri/secphase:dev-v0.2.0-hom",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/toil_hifi_dv1.2_mm2_all2dip_out/HG002.DCv1.2.minimap2v2.26.bam",
  "runCorrectBam.correctBam.options": "--maxDiv 0.02 --minReadLen 0 --minAlignmentLen 0",
  "runCorrectBam.correctBam.suffix": "maxDiv.02"
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid


export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-00:00 --partition=long"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
    /private/groups/hprc/polishing/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl \
    correct_bam_inputs.json \
    --outputDirectory ./correct_bam_outfiles \
    --outputFile correct_bam_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```

Run DeepPolisher
```
cd /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG002_y2_DCv1.2_mm2

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-00:00 --partition=long"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
    /private/groups/hprc/polishing/hpp_production_workflows/QC/wdl/tasks/DeepPolisher.wdl \
    DeepPolisher.inputs.json \
    --outputDirectory ./DeepPolisher.outfiles \
    --outputFile DeepPolisher.outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```

HG005

first need to align DCv1.2 reads to assembly with minimap2
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/flagger/wdls/workflows/long_read_aligner_scattered.wdl
```

```
{
  "longReadAlignmentScattered.preset": "map-hifi",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "longReadAlignmentScattered.enableRunningSecphase": false,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -L -Y -I8g",
  "longReadAlignmentScattered.kmerSize": 19,
  "longReadAlignmentScattered.sampleName": "HG005",
  "longReadAlignmentScattered.suffix": "trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26",
  "longReadAlignmentScattered.readFiles": ["/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200723_190224.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200801_011415.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200802_073944.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64109_200304_195708.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64109_200309_192110.dc.q20.fastq.gz"],
  "longReadAlignmentScattered.enableAddingMDTag": false
}
```

```
#!/bin/bash
#SBATCH --job-name=HG005_DCv1.2_minimap2
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --threads-per-core=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=8
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

cd /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2

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
    2>&1 | tee log.txt
```

Remove reads with de > 0.02
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl
```

```
{
  "runCorrectBam.correctBam.dockerImage": "mobinasri/secphase:dev-v0.2.0-hom",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/long_read_aligner_scattered_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.bam",
  "runCorrectBam.correctBam.options": "--maxDiv 0.02 --minReadLen 0 --minAlignmentLen 0",
  "runCorrectBam.correctBam.suffix": "maxDiv.02"
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/correct_bam

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-00:00 --partition=long"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
    /private/groups/hprc/polishing/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl \
    correct_bam_inputs.json \
    --outputDirectory ./correct_bam_outfiles \
    --outputFile correct_bam_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```

Run DeepPolisher

```
{
  "runDeepPolisher.sampleName": "HG005",
  "runDeepPolisher.ModelFilesTarGZ": "/private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-665.tar.gz",
  "runDeepPolisher.Bai": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/long_read_aligner_scattered_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.bam.bai",
  "runDeepPolisher.useOptimalGQFilter": true,
  "runDeepPolisher.Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "runDeepPolisher.dockerImage": "google/deepconsensus:polisher_v0.0.8_12122023",
  "runDeepPolisher.Bam": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/long_read_aligner_scattered_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.bam"
}
```

Run DeepPolisher
```
cd /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_mm2

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-00:00 --partition=long"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
    /private/groups/hprc/polishing/hpp_production_workflows/QC/wdl/tasks/DeepPolisher.wdl \
    DeepPolisher.inputs.json \
    --outputDirectory ./DeepPolisher.outfiles \
    --outputFile DeepPolisher.outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```
