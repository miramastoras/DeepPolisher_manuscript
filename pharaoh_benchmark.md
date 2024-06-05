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
