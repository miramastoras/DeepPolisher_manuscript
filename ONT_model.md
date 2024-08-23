## Training a DeepPolisher ONT model

### ONT Q25 data, shasta v0.12.1

Phased assemblies by GFAse provided by Konstantinos:
```
/private/groups/migalab/kkyriaki/experiments/shasta_assemblies/Q25_400speed_ML05/gfase/homology_new_chainer/phase_0.fasta
/private/groups/migalab/kkyriaki/experiments/shasta_assemblies/Q25_400speed_ML05/gfase/homology_new_chainer/phase_1.fasta
/private/groups/migalab/kkyriaki/experiments/shasta_assemblies/Q25_400speed_ML05/gfase/homology_new_chainer/unphased.fasta
```

Combine assemblies to same fasta file:
```
cd /private/groups/patenlab/mira/ONT_DeepPolisher/assemblies

cat phase_0.fasta phase_1.fasta unphased.fasta > HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.fasta
```
### Prepare alignments for training

Mobin's workflow:
https://github.com/mobinasri/flagger/blob/main/wdls/workflows/long_read_aligner_scattered.wdl

Get long read aligner wdl inputs
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs ~/progs/flagger/wdls/workflows/long_read_aligner_scattered.wdl
```

```
{
  "longReadAlignmentScattered.correctBamOptions": "--primaryOnly -m0 -a0 --maxDiv 0.09",
  "longReadAlignmentScattered.secphaseDockerImage": "mobinasri/secphase:v0.4.3",
  "longReadAlignmentScattered.preset": "map-ont",
  "longReadAlignmentScattered.secphaseVersion": "v0.4.3",
  "longReadAlignmentScattered.alignment.dockerImage": "mobinasri/long_read_aligner:v0.4.0",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/ONT_DeepPolisher/assemblies/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.fasta",
  "longReadAlignmentScattered.secphaseOptions": "--ont",
  "longReadAlignmentScattered.enableRunningSecphase": true,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -L -Y",
  "longReadAlignmentScattered.kmerSize": 15,
  "longReadAlignmentScattered.sampleName": "HG002",
  "longReadAlignmentScattered.suffix": "mm2v2.26.secphasev0.4.3",
  "longReadAlignmentScattered.readFiles": ["/private/nanopore/basecalled/q27/q25_400speed/PAW42666andPAW42495_gt10q10k_doradoTrimmed.fastq.gz"]
}
```

Submit the alignment
```
#!/bin/bash
#SBATCH --job-name=Q25_shasta_alignment_ONT_model
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=4
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

cd /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24]"
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

Run manually due to weird samtools error
```
#!/bin/bash
#SBATCH --job-name=Q25_shasta_alignment_ONT_model
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

docker run -u `id -u`:`id -g` -v /private:/private mobinasri/long_read_aligner:v0.4.0 minimap2 -k 15 --cs --eqx -L -Y -ax map-ont -t32 /private/groups/patenlab/mira/ONT_DeepPolisher/assemblies/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.fasta /private/nanopore/basecalled/q27/q25_400speed/PAW42666andPAW42495_gt10q10k_doradoTrimmed.fastq.gz  > /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.mm2.sam

samtools view -hb /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.mm2.sam > /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.mm2.bam
```

Run secphase
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs ~/progs/hpp_production_workflows/QC/wdl/tasks/secphase.wdl
```
```
{
  "runSecPhase.secphaseOptions": "--ont",
  "runSecPhase.version": "v0.4.3",
  "runSecPhase.secphaseDockerImage": "mobinasri/secphase:v0.4.3",
  "runSecPhase.diploidAssemblyFastaGz": "/private/groups/patenlab/mira/ONT_DeepPolisher/assemblies/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.fasta.gz",
  "runSecPhase.inputBam": "/private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.mm2.bam"
}

```

```
#!/bin/bash
#SBATCH --job-name=Q25_shasta_secphase_ONT_model
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=4
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

cd /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/secphase

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
     ~/progs/hpp_production_workflows/QC/wdl/tasks/secphase.wdl \
    secphase_inputs.json \
    --outputDirectory ./secphase_outfiles \
    --outputFile secphase_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```

Sort and index bam
```
#!/bin/bash
#SBATCH --job-name=ONT_model_sort
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

ulimit -Sn 5000

mkdir -p /data/tmp/mira/

samtools sort -@32 -T /data/tmp/mira/ /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.mm2.bam > /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.mm2.srt.bam

samtools index /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.mm2.srt.bam
```
Run correct bam
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs ~/progs/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl
```

```
{
  "runCorrectBam.correctBam.options": "--primaryOnly -m0 -a0 --maxDiv 0.09",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.mm2.srt.bam",
  "runCorrectBam.correctBam.phasingLogText": "/private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/secphase/secphase_outfiles/secphase_v0.4.3.out.log",
  "runCorrectBam.correctBam.suffix": "secphase.maxDiv.09"
}
```

```
#!/bin/bash
#SBATCH --job-name=Q25_shasta_correctbam_ONT_model
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=4
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

cd /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/correct_bam

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
     ~/progs/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl \
    correct_bam_inputs.json \
    --outputDirectory ./correct_bam_outfiles \
    --outputFile correct_bam_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```
### Spot check alignments for training model

Dipcall against grch38

```json
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/ONT_DeepPolisher/assemblies/phase_1.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/ONT_DeepPolisher/assemblies/phase_0.fasta",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/ONT_DeepPolisher/dipcall/unpolished_HG002_R10_45x_Q25.shasta_v0.12.1.gfase

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

### ONT Herro corrected data, verkko ont+poreC assembly

Run alignment manually
```
#!/bin/bash
#SBATCH --job-name=herro_ont_model_verkko_alignment
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

docker run -u `id -u`:`id -g` -v /private:/private mobinasri/long_read_aligner:v0.4.0 minimap2 -k 15 --cs --eqx -L -Y -ax map-ont -t32 /private/groups/migalab/kkyriaki/experiments/verkko_assemblies/VERKKO_LC2024_POREC/assembly.fasta /private/groups/migalab/kkyriaki/experiments/data/ONT_LC2024/assm/all_correct.fasta  > /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/HG002_verkko_porec_herro_ONT.mm2.sam

samtools view -hb /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/HG002_verkko_porec_herro_ONT.mm2.sam > /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/HG002_verkko_porec_herro_ONT.mm2.bam
```

Run secphase
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs ~/progs/hpp_production_workflows/QC/wdl/tasks/secphase.wdl
```
```
{
  "runSecPhase.secphaseOptions": "--ont",
  "runSecPhase.version": "v0.4.3",
  "runSecPhase.secphaseDockerImage": "mobinasri/secphase:v0.4.3",
  "runSecPhase.diploidAssemblyFastaGz": "/private/groups/patenlab/mira/ONT_DeepPolisher/assemblies/verkko_ont_porec/assembly.fasta.gz",
  "runSecPhase.inputBam": "/private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/HG002_verkko_porec_herro_ONT.mm2.bam"
}

```

```
#!/bin/bash
#SBATCH --job-name=secphase_ONT_verkko_herro
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=4
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

cd /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/secphase

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
     ~/progs/hpp_production_workflows/QC/wdl/tasks/secphase.wdl \
    secphase_inputs.json \
    --outputDirectory ./secphase_outfiles \
    --outputFile secphase_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```

Sort and index bam
```
#!/bin/bash
#SBATCH --job-name=ONT_model_sort
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

ulimit -Sn 5000

mkdir -p /data/tmp/mira/

samtools sort -@32 -T /data/tmp/mira/ /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/HG002_verkko_porec_herro_ONT.mm2.bam > /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/HG002_verkko_porec_herro_ONT.mm2.srt.bam

samtools index /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/HG002_verkko_porec_herro_ONT.mm2.srt.bam
```
Run correct bam
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs ~/progs/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl
```

```
{
  "runCorrectBam.correctBam.options": "--primaryOnly -m0 -a0 --maxDiv 0.09",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/HG002_verkko_porec_herro_ONT.mm2.srt.bam",
  "runCorrectBam.correctBam.phasingLogText": "/private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/secphase/secphase_outfiles/secphase_v0.4.3.out.log",
  "runCorrectBam.correctBam.suffix": "secphase.maxDiv.09"
}
```

```
#!/bin/bash
#SBATCH --job-name=correctbam_ONT_model
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=4
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

cd /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/correct_bam

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
     ~/progs/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl \
    correct_bam_inputs.json \
    --outputDirectory ./correct_bam_outfiles \
    --outputFile correct_bam_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```


Run mosdepth to check coverage
```
#!/bin/bash
#SBATCH --job-name=mosdepth_verkko_ont
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
    -v /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/correct_bam/correct_bam_outfiles/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -n --fast-mode -t 4 /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/correct_bam/correct_bam_outfiles/mosdepth/HG002_verkko_porec_herro_ONT \
    /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/correct_bam/correct_bam_outfiles/ONT_model_verkko_porec_HG002_alignments/HG002_verkko_porec_herro_ONT.mm2.srt.secphase.maxDiv.09.bam

python3 ~/progs/mosdepth/scripts/plot-dist.py /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_herro_ONT/correct_bam/correct_bam_outfiles/mosdepth/HG002_verkko_porec_herro_ONT.mosdepth.global.dist.txt
```

### ONT 6b4, verkko ont+poreC assembly

Get long read aligner wdl inputs
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs ~/progs/flagger/wdls/workflows/long_read_aligner_scattered.wdl
```

```
{
  "longReadAlignmentScattered.correctBamOptions": "--primaryOnly -m0 -a0 --maxDiv 0.09",
  "longReadAlignmentScattered.secphaseDockerImage": "mobinasri/secphase:v0.4.3",
  "longReadAlignmentScattered.preset": "map-ont",
  "longReadAlignmentScattered.secphaseVersion": "v0.4.3",
  "longReadAlignmentScattered.alignment.dockerImage": "mobinasri/long_read_aligner:v0.4.0",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/migalab/kkyriaki/experiments/verkko_assemblies/VERKKO_LC2024_POREC/assembly.fasta",
  "longReadAlignmentScattered.secphaseOptions": "--ont",
  "longReadAlignmentScattered.enableRunningSecphase": true,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -L -Y",
  "longReadAlignmentScattered.kmerSize": 15,
  "longReadAlignmentScattered.sampleName": "HG002",
  "longReadAlignmentScattered.suffix": "mm2v2.26.secphasev0.4.3",
  "longReadAlignmentScattered.readFiles": ["/private/groups/migalab/kkyriaki/experiments/data/ONT_LC2024/basecalling/apk/PAW41746.fastq"]
}
```

Submit the alignment
```
#!/bin/bash
#SBATCH --job-name=verkko_poreC_6b4
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=4
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

cd /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_6b4_ONT

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24]"
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

Run mosdepth to check coverage
```
#!/bin/bash
#SBATCH --job-name=mosdepth_6b4_verkko_ont
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
    -v /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_6b4_ONT/long_read_aligner_scattered_outfiles/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -n --fast-mode -t 4 /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_6b4_ONT/long_read_aligner_scattered_outfiles/mosdepth/HG002_verkko_porec_6b4_ONT \
    /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_6b4_ONT/long_read_aligner_scattered_outfiles/HG002_ONT_verkko_porec_6b4_alignments/HG002.mm2v2.26.secphasev0.4.3.secphase_v0.4.3.bam

python3 ~/progs/mosdepth/scripts/plot-dist.py /private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_verkko_porec_6b4_ONT/long_read_aligner_scattered_outfiles/mosdepth/HG002_verkko_porec_6b4_ONT.mosdepth.global.dist.txt
```

21x coverage per haplotype, ~42x total.

Uploaded here for training:
```
gsutil -o GSUtil:parallel_composite_upload_threshold=50GB -m cp -r ONT_model_verkko_porec_HG002_alignments gs://pepper-deepvariant/mira/

```


### Test polishing ONT verkko assembly with hifi reads

```
{
  "hprc_DeepPolisher.ONTReads": ["/private/groups/patenlab/masri/hprc/polishing/HG002/reads/R941_Guppy6/gt_100k/downsampled_40x/03_08_22_R941_HG002_1_Guppy_6.0.6_prom_sup.gt_100kb.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/R941_Guppy6/gt_100k/downsampled_40x/03_08_22_R941_HG002_4_Guppy_6.0.6_prom_sup.gt_100kb.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/R941_Guppy6/gt_100k/downsampled_40x/03_08_22_R941_HG002_3_Guppy_6.0.6_prom_sup.gt_100kb.fastq.gz"],
  "hprc_DeepPolisher.HifiReads": ["/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190714_120746.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190830_220126.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190901_095311.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64012_190920_173625.dc.q20.fastq.gz"],
  "hprc_DeepPolisher.Hap1RawFasta": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/assembly.pat.fa",
  "hprc_DeepPolisher.Hap2RawFasta": "/private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/ont_model/assign_contigs/assembly.mat.fa",
  "hprc_DeepPolisher.DeepPolisherModelFilesTarGZ": "/private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-665.tar.gz",
  "hprc_DeepPolisher.sampleName": "HG002"
}
```


```
mkdir -p slurm_logs

sbatch \
     --job-name=verkko_asm_ont_hifi_reads_hprc-DeepPolisher \
     --array=[1]%1 \
     --partition=long \
     --cpus-per-task=32 \
     --mem=400gb \
     --mail-type=FAIL,END \
     --exclude=phoenix-[09,10,22,23,24] \
     --mail-user=mmastora@ucsc.edu \
     /private/groups/hprc/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
     --wdl /private/groups/hprc/polishing/hpp_production_workflows/QC/wdl/workflows/hprc_DeepPolisher.wdl \
     --sample_csv ./samples.csv \
     --input_json_path '/private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG002_verkko_porec_ONT_asm_HiFi_polish_mm2_model1/hprc_DeepPolisher.inputs.json'
```
