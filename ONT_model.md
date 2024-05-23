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

Run correct bam
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs ~/progs/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl
```

```
{
  "runCorrectBam.correctBam.options": "--primaryOnly -m0 -a0 --maxDiv 0.09",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/ONT_DeepPolisher/alignments/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.minimap2/HG002_R10_45x_Q25.shasta_v0.12.1.gfase.dip.mm2.bam",
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
