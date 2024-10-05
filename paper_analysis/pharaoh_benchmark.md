### Benchmarking performance of PHARAOH pipeline

HG002

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
Running DeepPolisher on minimap2 alignments

```
{
  "runDeepPolisher.sampleName": "HG002",
  "runDeepPolisher.ModelFilesTarGZ": "/private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-665.tar.gz",
  "runDeepPolisher.Bai": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/correct_bam/correct_bam_outfiles/HG002.DCv1.2.minimap2v2.26.maxDiv.02.bam.bai",
  "runDeepPolisher.useOptimalGQFilter": true,
  "runDeepPolisher.Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "runDeepPolisher.dockerImage": "google/deepconsensus:polisher_v0.0.8_12122023",
  "runDeepPolisher.Bam": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/correct_bam/correct_bam_outfiles/HG002.DCv1.2.minimap2v2.26.maxDiv.02.bam"
}
```

Run DeepPolisher
```
cd /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG002_y2_DCv1.2_mm2

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-00:00 --partition=high_priority"
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
    --restart \
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
    --restart \
    2>&1 | tee log.txt
```

Run DeepPolisher

```
{
  "runDeepPolisher.sampleName": "HG005",
  "runDeepPolisher.ModelFilesTarGZ": "/private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-665.tar.gz",
  "runDeepPolisher.Bai": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/correct_bam/correct_bam_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.maxDiv.02.bam.bai",
  "runDeepPolisher.useOptimalGQFilter": true,
  "runDeepPolisher.Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "runDeepPolisher.dockerImage": "google/deepconsensus:polisher_v0.0.8_12122023",
  "runDeepPolisher.Bam": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/correct_bam/correct_bam_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.maxDiv.02.bam"
}
```

Run DeepPolisher
```
cd /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_mm2

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-00:00 --partition=high_priority"
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
    --restart \
    2>&1 | tee log.txt
```

clush -w phoenix-[00-09,11-20] "rm -rvf /data/tmp/mmastora/"

debug 20x HG005
```
#!/bin/bash
#SBATCH --job-name=HG005_DP_20x
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --time=7-00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    google/deepconsensus:polisher_v0.0.8_12122023 \
    polisher make_images \
    --bam /private/groups/patenlab/mira/hprc_polishing/data/phoenix_batch_submissions/correct_bam/HG005_20x_dp_no_phasing/analysis/correct_bam_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.20x.maxDiv.02.bam \
    --fasta /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --output /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/phoenix_batch_submissions_manuscript/HG005_20x_dp_no_phasing/manual/images/images \
    --cpus 64

# Inference on images to generate VCFs
docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    google/deepconsensus:polisher_v0.0.8_12122023 \
    polisher inference \
    --input_dir /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/phoenix_batch_submissions_manuscript/HG005_20x_dp_no_phasing/manual/images/images \
    --out_dir /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/phoenix_batch_submissions_manuscript/HG005_20x_dp_no_phasing/manual/vcf/ \
    --checkpoint /private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-665/checkpoint-665 \
    --reference_file /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --sample_name HG005 \
    --cpus 64
```

Running the rest through this script and not in toil because I can't figure out the bug

```
#!/bin/bash
#SBATCH --job-name=HG005_DP
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --time=7-00:00

BAM=$1
OUTDIR=$2

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    google/deepconsensus:polisher_v0.0.8_12122023 \
    polisher make_images \
    --bam ${BAM} \
    --fasta /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --output ${OUTDIR}/images/images \
    --cpus 64

# Inference on images to generate VCFs
docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    google/deepconsensus:polisher_v0.0.8_12122023 \
    polisher inference \
    --input_dir ${OUTDIR}/images/images \
    --out_dir ${OUTDIR}/vcf/ \
    --checkpoint /private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-665/checkpoint-665 \
    --reference_file /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --sample_name HG005 \
    --cpus 64
```

```
cd /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/phoenix_batch_submissions_manuscript

cut -f1 -d"," GIAB_samples_deepPolisher_manuscript.csv \
    | grep -v "sample_id" | grep "HG005" | while read line ; do \
    bam=`grep $line GIAB_samples_deepPolisher_manuscript.csv | cut -f 8 -d","`
    mkdir -p /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/phoenix_batch_submissions_manuscript/${line}/manual/
    echo "sbatch launch_HG005_manually.sh $bam /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/phoenix_batch_submissions_manuscript/${line}/manual/"
    done
```

```
sbatch launch_HG005_manually.sh /private/groups/patenlab/mira/hprc_polishing/data/phoenix_batch_submissions/correct_bam/HG005_60x_dp_no_phasing/analysis/correct_bam_outputs/HG005.DCv1.2_60x.minimap2v2.26.maxDiv.02.bam /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/phoenix_batch_submissions_manuscript/HG005_60x_dp_no_phasing/manual/

sbatch launch_HG005_manually.sh /private/groups/patenlab/mira/hprc_polishing/data/phoenix_batch_submissions/correct_bam/HG005_50x_dp_no_phasing/analysis/correct_bam_outputs/HG005.DCv1.2_50x.minimap2v2.26.maxDiv.02.bam /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/phoenix_batch_submissions_manuscript/HG005_50x_dp_no_phasing/manual/

sbatch launch_HG005_manually.sh /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/correct_bam/correct_bam_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.maxDiv.02.bam /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/phoenix_batch_submissions_manuscript/HG005_40x_dp_no_phasing/manual/

sbatch launch_HG005_manually.sh /private/groups/patenlab/mira/hprc_polishing/data/phoenix_batch_submissions/correct_bam/HG005_30x_dp_no_phasing/analysis/correct_bam_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.30x.maxDiv.02.bam /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/phoenix_batch_submissions_manuscript/HG005_30x_dp_no_phasing/manual/

sbatch launch_HG005_manually.sh /private/groups/patenlab/mira/hprc_polishing/data/phoenix_batch_submissions/correct_bam/HG005_10x_dp_no_phasing/analysis/correct_bam_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.10x.maxDiv.02.bam /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/phoenix_batch_submissions_manuscript/HG005_10x_dp_no_phasing/manual/
```

### Running minimap2 for 40x to test DP without correct bam or secphase

haven't been able to figure out the bug so need to redo alignments
```
{
  "longReadAlignmentScattered.correctBamOptions": "--maxDiv 0.02",
  "longReadAlignmentScattered.secphaseDockerImage": "mobinasri/secphase:v0.4.3",
  "longReadAlignmentScattered.preset": "map-hifi",
  "longReadAlignmentScattered.secphaseVersion": "v0.4.3",
  "longReadAlignmentScattered.alignment.dockerImage": "mobinasri/long_read_aligner:v0.3.3",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "longReadAlignmentScattered.secphaseOptions": "--hifi",
  "longReadAlignmentScattered.enableRunningSecphase": false,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -Y -L",
  "longReadAlignmentScattered.kmerSize": 19,
  "longReadAlignmentScattered.sampleName": "HG005",
  "longReadAlignmentScattered.suffix": "mm2v2.26.secphasev0.4.3.maxDiv.02",
  "longReadAlignmentScattered.readFiles": [
    "/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200723_190224.dc.q20.fastq.gz",
    "/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200801_011415.dc.q20.fastq.gz",
    "/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200802_073944.dc.q20.fastq.gz",
    "/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64109_200304_195708.dc.q20.fastq.gz",
    "/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64109_200309_192110.dc.q20.fastq.gz"
  ]
}
```

Submit the alignment
```
#!/bin/bash
#SBATCH --job-name=HG005_DCv1.2_40x_minimap2_only
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=4
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

cd /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_40x_minimap2_only

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


    time toil-wdl-runner \
        --jobStore ./jobstore \
        --stats \
        --clean=never \
        --batchSystem slurm \
        --batchLogsDir ./toil_logs \
         ~/progs/hpp_production_workflows/QC/wdl/tasks/DeepPolisher.wdl \
        DeepPolisher_inputs.json \
        --outputDirectory ./DeepPolisher_outfiles \
        --outputFile DeepPolisher_outputs.json \
        --runLocalJobsOnWorkers \
        --retryCount 1 \
        --disableProgress \
        --logDebug \
        2>&1 | tee log.txt
```

Run DP
```
#!/bin/bash
#SBATCH --job-name=HG005_DP
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --time=7-00:00

BAM=$1
OUTDIR=$2

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    google/deepconsensus:polisher_v0.0.8_12122023 \
    polisher make_images \
    --bam ${BAM} \
    --fasta /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --output ${OUTDIR}/images/images \
    --cpus 64

# Inference on images to generate VCFs
docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    google/deepconsensus:polisher_v0.0.8_12122023 \
    polisher inference \
    --input_dir ${OUTDIR}/images/images \
    --out_dir ${OUTDIR}/vcf/ \
    --checkpoint /private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-665/checkpoint-665 \
    --reference_file /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --sample_name HG005 \
    --cpus 64
```

```
sbatch slurm.sh /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_40x_minimap2_only/long_read_aligner_scattered_outfiles/HG005.mm2v2.26.secphasev0.4.3.maxDiv.02.bam /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/HG005_DCv1.2_40x_minimap2_only/
```

Test DP on old pharaoh bam
```
sbatch /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/HG005_DCv1.2_40x_minimap2_only/slurm.sh /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_PHv6_DPmm2model1/toil_hprc_deepPolisher_out/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.PHARAOHv6.bam /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/HG005_pharaohv6/
```

Run DP in wdl
```
{
  "runDeepPolisher.sampleName": "HG005",
  "runDeepPolisher.ModelFilesTarGZ": "/private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-665.tar.gz",
  "runDeepPolisher.Bai": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_40x_minimap2_only/long_read_aligner_scattered_outfiles/HG005.mm2v2.26.secphasev0.4.3.maxDiv.02.bam.bai",
  "runDeepPolisher.Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "runDeepPolisher.dockerImage": "google/deepconsensus:polisher_v0.0.8_12122023",
  "runDeepPolisher.Bam": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_40x_minimap2_only/long_read_aligner_scattered_outfiles/HG005.mm2v2.26.secphasev0.4.3.maxDiv.02.bam",
  "runDeepPolisher.memSize":512,
  "runDeepPolisher.diskSize":512
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/HG005_DCv1.2_40x_minimap2_only

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24,18]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
     ~/progs/hpp_production_workflows/QC/wdl/tasks/DeepPolisher.wdl \
    DeepPolisher_inputs.json \
    --outputDirectory ./DeepPolisher_outfiles \
    --outputFile DeepPolisher_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```

Run with subset only
```
samtools view -bh /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_40x_minimap2_only/long_read_aligner_scattered_outfiles/HG005.mm2v2.26.secphasev0.4.3.maxDiv.02.bam h1tg000001l > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_40x_minimap2_only/long_read_aligner_scattered_outfiles/HG005.mm2v2.26.secphasev0.4.3.maxDiv.02.h1tg000001l.bam

samtools faidx /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa h1tg000001l > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.h1tg000001l.fa
```

Run DP in wdl
```
{
  "runDeepPolisher.sampleName": "HG005",
  "runDeepPolisher.ModelFilesTarGZ": "/private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-665.tar.gz",
  "runDeepPolisher.Bai": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_40x_minimap2_only/long_read_aligner_scattered_outfiles/HG005.mm2v2.26.secphasev0.4.3.maxDiv.02.h1tg000001l.bam.bai",
  "runDeepPolisher.Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.h1tg000001l.fa",
  "runDeepPolisher.dockerImage": "google/deepconsensus:polisher_v0.0.8_12122023",
  "runDeepPolisher.Bam": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_40x_minimap2_only/long_read_aligner_scattered_outfiles/HG005.mm2v2.26.secphasev0.4.3.maxDiv.02.h1tg000001l.bam"
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/HG005_DCv1.2_40x_minimap2_only/h1tg

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24,18]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
     ~/progs/hpp_production_workflows/QC/wdl/tasks/DeepPolisher.wdl \
    DeepPolisher_inputs.json \
    --outputDirectory ./DeepPolisher_outfiles \
    --outputFile DeepPolisher_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```

Take pharaoh reads, re-align to assembly with minimap2

```
cd /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_40x_minimap2_only/PHARAOH_realign

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24,18]"
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
    2>&1 | tee log.txt

```

```
{
  "longReadAlignmentScattered.correctBamOptions": "--maxDiv 0.02",
  "longReadAlignmentScattered.secphaseDockerImage": "mobinasri/secphase:v0.4.3",
  "longReadAlignmentScattered.preset": "map-hifi",
  "longReadAlignmentScattered.secphaseVersion": "v0.4.3",
  "longReadAlignmentScattered.alignment.dockerImage": "mobinasri/long_read_aligner:v0.3.3",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "longReadAlignmentScattered.secphaseOptions": "--hifi",
  "longReadAlignmentScattered.enableRunningSecphase": false,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -Y -L",
  "longReadAlignmentScattered.kmerSize": 19,
  "longReadAlignmentScattered.sampleName": "HG005",
  "longReadAlignmentScattered.suffix": "mm2v2.26.secphasev0.4.3.maxDiv.02",
  "longReadAlignmentScattered.readFiles": ["/private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_PHv6_DPmm2model1/toil_hprc_deepPolisher_out/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.PHARAOHv6.bam"]
}
```


```
{
  "runDeepPolisher.sampleName": "HG005",
  "runDeepPolisher.ModelFilesTarGZ": "/private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-665.tar.gz",
  "runDeepPolisher.Bai": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_40x_minimap2_only/PHARAOH_realign/long_read_aligner_scattered_outfiles/HG005.mm2v2.26.secphasev0.4.3.maxDiv.02.bam.bai",
  "runDeepPolisher.Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "runDeepPolisher.dockerImage": "google/deepconsensus:polisher_v0.0.8_12122023",
  "runDeepPolisher.Bam": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_40x_minimap2_only/PHARAOH_realign/long_read_aligner_scattered_outfiles/HG005.mm2v2.26.secphasev0.4.3.maxDiv.02.bam"
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/deepPolisher_runs/HG005_DCv1.2_40x_minimap2_only/PHARAOH_realign

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=high_priority --exclude=phoenix-[09,10,22,23,24,18]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
     ~/progs/hpp_production_workflows/QC/wdl/tasks/DeepPolisher.wdl \
    DeepPolisher_inputs.json \
    --outputDirectory ./DeepPolisher_outfiles \
    --outputFile DeepPolisher_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```


### Get number of bases in homozygous regions for HG002 and HG005 assemblies

Align hap1 to hap2
```
#!/bin/bash
#SBATCH --job-name=homozygous_regions
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

source /private/home/mmastora/progs/miniconda3/etc/profile.d/conda.sh
conda activate analysis

HAP1=/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa
HAP2=/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa
DIR=HG002_y2

HAP1=/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa
HAP2=/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa
DIR=HG005_y2

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs \
    ${HAP1} \
    ${HAP2} \
    -o /private/groups/patenlab/mira/hprc_polishing/PHARAOH_benchmark/find_homozygous_regions/${DIR}/${DIR}.raw.mat_to_pat.paf

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/secphase:dev-v0.2.0 \
    python3 /home/programs/src/find_homozygous_regions.py \
    -p /private/groups/patenlab/mira/hprc_polishing/PHARAOH_benchmark/find_homozygous_regions/${DIR}/${DIR}.raw.mat_to_pat.paf \
    -m 20000 -e 50000 -o /private/groups/patenlab/mira/hprc_polishing/PHARAOH_benchmark/find_homozygous_regions/${DIR}/${DIR}.homozygous.20kb

```

bases in HG005 homozygous regions
```
601,850,247

487,955,715
```

Rerun PHARAOH for HG005 and HG002 to get secphase modified blocks bed files


java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/polishing_branch_hpp_master/hpp_production_workflows/polishing/wdl/workflows/PHARAOH.wdl


```
{
  "PHARAOH.allHifiToDiploidBai": "File",
  "PHARAOH.allONTToHap2Bam": "File",
  "PHARAOH.allONTToHap1Bam": "File",
  "PHARAOH.diploidFaGz": "File",
  "PHARAOH.Hap1Fasta": "File",
  "PHARAOH.allONTToHap2Bai": "File",
  "PHARAOH.Hap2Fasta": "File",
  "PHARAOH.Hap1FastaIndex": "File",
  "PHARAOH.allHifiToDiploidBam": "File",
  "PHARAOH.allONTToHap1Bai": "File",
  "PHARAOH.Hap2FastaIndex": "File",
  "PHARAOH.sampleName": "String"
}
```

### Reoptimize maxDiv for CCS reads

```
{
  "longReadAlignmentScattered.preset": "map-hifi",
  "longReadAlignmentScattered.alignment.dockerImage": "mobinasri/long_read_aligner:v0.3.3",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "longReadAlignmentScattered.enableRunningSecphase": false,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -Y -L",
  "longReadAlignmentScattered.kmerSize": 19,
  "longReadAlignmentScattered.sampleName": "HG002",
  "longReadAlignmentScattered.suffix": "mm2v2.26",
  "longReadAlignmentScattered.readFiles": [
    "/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/m64015_190920_185703.Q20.fastq",
    "/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/m64015_190922_010918.Q20.fastq",
    "/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/m64011_190830_220126.Q20.fastq",
    "/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/m64011_190901_095311.Q20.fastq",
    "/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/m64012_190920_173625.Q20.fastq"
  ]
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/CCS_40x

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=24:00:00 --partition=long --exclude=phoenix-[09,10,22,23,24,18,15,13,06,01,20,21,07]"
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
    2>&1 | tee log.txt

```

Pull out one contig:
```
# overlap with error k-ers to get the erroneous reads


```

plot de tag distribution:

```
import pysam
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

bamfile = pysam.AlignmentFile("/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/CCS_40x/long_read_aligner_scattered_outfiles/HG002.mm2v2.26.h1tg000001l.bam", "rb")

de_list=[]
for read in bamfile.fetch():
    tag=str(read.get_tag("de"))

    de_list.append(tag)

data=np.array(de_list)
data = data.astype(float)

# Set the bin size and number of labels on the x-axis
bin_size = 3e-4

# Create the histogram data
hist_values, hist_bins = np.histogram(data, bins=int(1 / bin_size))

# Calculate the center of each bin
bin_centers = (hist_bins[:-1] + hist_bins[1:]) / 2 + bin_size*2

# Plot the histogram as lines
plt.plot(bin_centers, hist_values, drawstyle='steps-post', linewidth=0.5)

# Set the x-axis tick labels
x_ticks = np.linspace(min(data), .02, 10)

plt.xticks(x_ticks, fontsize=5)
plt.xlim([min(data),.02])

# Set labels and title
plt.xlabel('de Value')
plt.ylabel('Frequency')
plt.title('de dist of reads overlapping FP kmers, binsize 3e-3')


# Save the plot to an external file (e.g., PNG format)
plt.savefig('/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/CCS_40x/long_read_aligner_scattered_outfiles/de_dist_3e-4.png', dpi=600)

# Close the plot
plt.close()
```
