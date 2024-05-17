# Experimenting with different read coverages for DeepPolisher

## Titrating HiFi coverage for the PHARAOH + DeepPolisher pipeline

### HG002 coverage titrations

Downsample HiFi read alignment file using samtools - using reads already aligned to assembly with winnowmap from previous experiment.
```
#!/bin/bash
#SBATCH --job-name=downsampling_HG002_DCv1.2
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

# 60x
samtools view -s 0.75 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/downsampled_bams/HG002.DCv1.2.60x.winnowmapv2.03.bam

# 50x
samtools view -s 0.625 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/downsampled_bams/HG002.DCv1.2.50x.winnowmapv2.03.bam

# 30 x - using minimap2 bam so it can be used for DP without pharaoh
samtools view -s 0.75 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/toil_hifi_dv1.2_mm2_all2dip_out/HG002.DCv1.2.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/downsampled_bams/HG002.DCv1.2.30x.minimap2v2.26.bam

# 20x
samtools view -s 0.5 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/toil_hifi_dv1.2_mm2_all2dip_out/HG002.DCv1.2.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/downsampled_bams/HG002.DCv1.2.20x.minimap2v2.26.bam

# 10x
samtools view -s 0.25 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/toil_hifi_dv1.2_mm2_all2dip_out/HG002.DCv1.2.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/downsampled_bams/HG002.DCv1.2.10x.minimap2v2.26.bam
```

Run HPRC_deepPolisher workflow on different coverage titrations

### HG005 coverage titrations

Downsample HiFi read alignment file using samtools - using reads already aligned to assembly with minimap2 from 40x experiment.
```
#!/bin/bash
#SBATCH --job-name=downsampling_HG005_DCv1.2
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

# 30 x
samtools view -s 0.75 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/long_read_aligner_scattered_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/downsampled_bams/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.30x.bam

# 20x
samtools view -s 0.5 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/long_read_aligner_scattered_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/downsampled_bams/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.20x.bam

# 10x
samtools view -s 0.25 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/long_read_aligner_scattered_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/downsampled_bams/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.10x.bam
```

50x - all flow cells we have
```
["/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200723_190224.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200801_011415.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200802_073944.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64109_200304_195708.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64109_200309_192110.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200730_190124.dc.q20.fastq.gz"]
```

## Titrating HiFi coverage for DeepPolisher without PHARAOH

### HG002

samtools view -s 0.75 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam >

# 50x
samtools view -s 0.625 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam >

New alignment of 60 and 50x files minimap2
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/flagger/wdls/workflows/long_read_aligner_scattered.wdl
```

60x
```
{
  "longReadAlignmentScattered.preset": "map-hifi",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "longReadAlignmentScattered.enableRunningSecphase": false,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -L -Y -I8g",
  "longReadAlignmentScattered.kmerSize": 19,
  "longReadAlignmentScattered.sampleName": "HG002",
  "longReadAlignmentScattered.suffix": "DCv1.2_60x.minimap2v2.26",
  "longReadAlignmentScattered.readFiles": ["/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/downsampled_bams/HG002.DCv1.2.60x.winnowmapv2.03.bam"],
  "longReadAlignmentScattered.enableAddingMDTag": false
}
```
50x

```
{
  "longReadAlignmentScattered.preset": "map-hifi",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "longReadAlignmentScattered.enableRunningSecphase": false,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -L -Y -I8g",
  "longReadAlignmentScattered.kmerSize": 19,
  "longReadAlignmentScattered.sampleName": "HG002",
  "longReadAlignmentScattered.suffix": "DCv1.2_50x.minimap2v2.26",
  "longReadAlignmentScattered.readFiles": ["/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/downsampled_bams/HG002.DCv1.2.50x.winnowmapv2.03.bam"],
  "longReadAlignmentScattered.enableAddingMDTag": false
}
```


```
#!/bin/bash
#SBATCH --job-name=HG002_mm2_50x
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --time=7-00:00

cd /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/60x_bam
cd /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/50x_bam

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


Already have 30,20,10x bams. Submit to deepPolisher

### HG005

Only have 50x coverage
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/flagger/wdls/workflows/long_read_aligner_scattered.wdl
```

50x (48x)
```
{
  "longReadAlignmentScattered.preset": "map-hifi",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "longReadAlignmentScattered.enableRunningSecphase": false,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -L -Y -I8g",
  "longReadAlignmentScattered.kmerSize": 19,
  "longReadAlignmentScattered.sampleName": "HG005",
  "longReadAlignmentScattered.suffix": "DCv1.2_60x.minimap2v2.26",
  "longReadAlignmentScattered.readFiles": ["/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200723_190224.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200801_011415.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64017_200802_073944.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64109_200304_195708.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG005/reads/HiFi_DC_v1.2/downsampled_40x/m64109_200309_192110.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200730_190124.dc.q20.fastq.gz"],
  "longReadAlignmentScattered.enableAddingMDTag": false
}
```

```
#!/bin/bash
#SBATCH --job-name=align_HG005_50x
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=4
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-0:00

cd /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/50x_bam

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-0:00 --partition=long --exclude=phoenix-[09,10,22,23,24]"
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

Already have 30,20,10x bams. Submit to deepPolisher


## Titrating HiFi coverage for CCS and Sequel 2 DCv1.1 data

HG002 CCS:

53x max coverage:
```
m64011_190830_220126.Q20.fastq	8.80676
m64011_190901_095311.Q20.fastq	8.3522
m64012_190920_173625.Q20.fastq	9.81073
m64012_190921_234837.Q20.fastq	9.98385
m64015_190920_185703.Q20.fastq	9.50113
m64015_190922_010918.Q20.fastq	9.10393
sample.fastq	53.0788
```
30x (29.2x) coverage:
```
m64012_190920_173625.Q20.fastq	9.81073
m64012_190921_234837.Q20.fastq	9.98385
m64015_190920_185703.Q20.fastq	9.50113
```
20x (19.7)
```
m64012_190920_173625.Q20.fastq	9.81073
m64012_190921_234837.Q20.fastq	9.98385
```
10x
```
m64012_190921_234837.Q20.fastq	9.98385
```
HG005 CCS :

50x (48x) coverage: all files
```
m64017_200723_190224.fastq	6.59322
m64017_200730_190124.fastq	6.69309
m64017_200801_011415.fastq	5.57823
m64017_200802_073944.fastq	4.31361
m64109_200304_195708.fastq	8.82954
m64109_200309_192110.fastq	8.58724
m64109_200311_013444.fastq	8.11014
sample.fastq	48.3507
```

HG005 CCS :

50x (48x) coverage: all files
```
m64017_200723_190224.fastq	6.59322
m64017_200730_190124.fastq	6.69309
m64017_200801_011415.fastq	5.57823
m64017_200802_073944.fastq	4.31361
m64109_200304_195708.fastq	8.82954
m64109_200309_192110.fastq	8.58724
m64109_200311_013444.fastq	8.11014
sample.fastq	48.3507
```

30x (31.89)
```
m64017_200723_190224.fastq	6.59322
m64017_200730_190124.fastq	6.69309
m64017_200801_011415.fastq	5.57823
m64017_200802_073944.fastq	4.31361
m64109_200304_195708.fastq	8.82954
```

20x
```
m64109_200309_192110.fastq	8.58724
m64109_200311_013444.fastq	8.11014
m64017_200802_073944.fastq	4.31361
```

10x ()
```
m64017_200730_190124.fastq	6.69309
m64017_200802_073944.fastq	4.31361
```

HG002 DCv1.1
```
66x bam:
/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.bam
```
Downsample other coverages:
```
#!/bin/bash
#SBATCH --job-name=downsampling_HG002_DCv1.1
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

samtools view -s 0.75 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.50x.bam

samtools view -s 0.5 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.30x.bam

samtools view -s 0.3 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.20x.bam

samtools view -s 0.15 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.10x.bam
```
