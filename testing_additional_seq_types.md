### Testing the mm2 model 1 on different types of HiFi data other than 1.2

#### CCS HiFi data

HG002: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/
HG005: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/HG005/raw_data/PacBio_HiFi/

From the HG002 readme, select data to match what we have for HG005:
```
Sequel II System with Chemistry 2.0
    Size selection  15 kb or 20 kb selected on SageELF
    Run time        30 hrs per SMRT Cell 8M
    CCS             "Circular Consensus Sequencing" analysis in SMRT Link v8.0
    SRA             https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA586863
    15 kb Library, ~36-fold coverage
        m64012_190920_173625
        m64012_190921_234837
        m64015_190920_185703
        m64015_190922_010918
    20 kb Library, ~16-fold coverage
        m64011_190830_220126
        m64011_190901_095311
```

Download data:

```
#!/bin/bash
#SBATCH --job-name=download_CCS
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

cd /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/15kb/m64012_190920_173625.Q20.fastq
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/15kb/m64012_190921_234837.Q20.fastq
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/15kb/m64015_190920_185703.Q20.fastq
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/15kb/m64015_190922_010918.Q20.fastq
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/20kb/m64011_190830_220126.Q20.fastq
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/20kb/m64011_190901_095311.Q20.fastq
```
Check read coverage for HG002 data
```
{
  "runReadStats.dockerImage": "mobinasri/fai_read_stats:latest",
  "runReadStats.reads": ["/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/m64011_190830_220126.Q20.fastq","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/m64011_190901_095311.Q20.fastq","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/m64012_190920_173625.Q20.fastq","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/m64012_190921_234837.Q20.fastq","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/m64015_190920_185703.Q20.fastq","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/m64015_190922_010918.Q20.fastq"]
}
```

```
#!/bin/bash
#SBATCH --job-name=read_stats
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --exclude=phoenix-06
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

cd /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium"
export TOIL_COORDINATION_DIR=/data/tmp


LOCAL_FOLDER=/data/tmp/$(whoami)/read_stats_HG002
mkdir -p ${LOCAL_FOLDER}

toil-wdl-runner \
    --jobStore ${LOCAL_FOLDER}/jobstore \
    --batchSystem single_machine \
    --maxCores 32 \
    --batchLogsDir /private/groups/patenlab/mira/WashU_pedigree/data/read_stats/HiFi_DCv1.2/toil_logs \
    ~/progs/master_hpp_production_workflows/QC/wdl/tasks/read_stats.wdl \
    /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/read_stats.inputs.json  \
    --outputDirectory /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/outputs \
    --outputFile /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel2_CCS/outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    2>&1 | tee log.txt

toil clean ${LOCAL_FOLDER}/jobstore

```

```
grep "total_reads" sample_all.report.tsv > coverage_per_file.tsv
grep "N50" sample_all.report.tsv | paste -d "\t" coverage_per_file.tsv - > tmp; mv tmp coverage_per_file.tsv
awk '{ print $1"\t"$3*$6 / 3000000000 }' coverage_per_file.tsv


m64011_190830_220126.Q20.fastq	8.80676
m64011_190901_095311.Q20.fastq	8.3522
m64012_190920_173625.Q20.fastq	9.81073
m64012_190921_234837.Q20.fastq	9.98385
m64015_190920_185703.Q20.fastq	9.50113
m64015_190922_010918.Q20.fastq	9.10393
sample.fastq	53.0788


Selecting these files for 45x coverage
m64015_190920_185703.Q20.fastq	9.50113
m64015_190922_010918.Q20.fastq	9.10393
m64011_190830_220126.Q20.fastq	8.80676
m64011_190901_095311.Q20.fastq	8.3522
m64012_190920_173625.Q20.fastq	9.81073
```

Download HG005 data
```
#!/bin/bash
#SBATCH --job-name=download_CCS_HG5
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG005/raw_data/PacBio_HiFi/PBmixSequel788_1_A01_PBXX_30hours_19kbV2PD_70pM_HumanHG005_CCS/m64109_200304_195708.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG005/raw_data/PacBio_HiFi/PBmixSequel789_1_A01_PBXW_30hours_15kbV2PD_70pM_HumanHG005_CCS/m64109_200309_192110.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG005/raw_data/PacBio_HiFi/PBmixSequel789_2_B01_PBXX_30hours_19kbV2PD_70pM_HumanHG005_CCS/m64109_200311_013444.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG005/raw_data/PacBio_HiFi/PBmixSequel840_1_A01_PCCD_30hours_15kbV2PD_70pM_HumanHG005_CCS/m64017_200723_190224.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG005/raw_data/PacBio_HiFi/PBmixSequel842_1_A01_PCCD_30hours_15kbV2PD_70pM_HumanHG005_CCS/m64017_200730_190124.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG005/raw_data/PacBio_HiFi/PBmixSequel842_2_B01_PCCD_30hours_15kbV2PD_70pM_HumanHG005_CCS/m64017_200801_011415.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG005/raw_data/PacBio_HiFi/PBmixSequel842_3_C01_PCCD_30hours_15kbV2PD_70pM_HumanHG005_CCS/m64017_200802_073944.fastq.gz
```

Check read coverage for HG005 data
```
{
  "runReadStats.dockerImage": "mobinasri/fai_read_stats:latest",
  "runReadStats.reads": ["/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/HiFi_Sequel2_CCS/m64017_200723_190224.fastq.gz","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/HiFi_Sequel2_CCS/m64017_200730_190124.fastq.gz","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/HiFi_Sequel2_CCS/m64017_200801_011415.fastq.gz","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/HiFi_Sequel2_CCS/m64017_200802_073944.fastq.gz","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/HiFi_Sequel2_CCS/m64109_200304_195708.fastq.gz","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/HiFi_Sequel2_CCS/m64109_200309_192110.fastq.gz","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/HiFi_Sequel2_CCS/m64109_200311_013444.fastq.gz"]
}
```

```
#!/bin/bash
#SBATCH --job-name=read_stats_Hg5
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --exclude=phoenix-06
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

cd /private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/HiFi_Sequel2_CCS

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium"
export TOIL_COORDINATION_DIR=/data/tmp


LOCAL_FOLDER=/data/tmp/$(whoami)/read_stats_HG005
mkdir -p ${LOCAL_FOLDER}

toil-wdl-runner \
    --jobStore ${LOCAL_FOLDER}/jobstore \
    --batchSystem single_machine \
    --maxCores 32 \
    --batchLogsDir /private/groups/patenlab/mira/WashU_pedigree/data/read_stats/HiFi_DCv1.2/toil_logs \
    ~/progs/master_hpp_production_workflows/QC/wdl/tasks/read_stats.wdl \
    /private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/HiFi_Sequel2_CCS/read_stats.inputs.json  \
    --outputDirectory /private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/HiFi_Sequel2_CCS/outputs \
    --outputFile /private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/HiFi_Sequel2_CCS/outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    2>&1 | tee log.txt

toil clean ${LOCAL_FOLDER}/jobstore
```

```
grep "total_reads" sample_all.report.tsv > coverage_per_file.tsv
grep "N50" sample_all.report.tsv | paste -d "\t" coverage_per_file.tsv - > tmp; mv tmp coverage_per_file.tsv
awk '{ print $1"\t"$3*$6 / 3000000000 }' coverage_per_file.tsv


m64017_200723_190224.fastq	6.59322
m64017_200730_190124.fastq	6.69309
m64017_200801_011415.fastq	5.57823
m64017_200802_073944.fastq	4.31361
m64109_200304_195708.fastq	8.82954
m64109_200309_192110.fastq	8.58724
m64109_200311_013444.fastq	8.11014
sample.fastq	48.3507

m64017_200723_190224.fastq	6.59322
m64017_200730_190124.fastq	6.69309
m64017_200802_073944.fastq	4.31361
m64109_200304_195708.fastq	8.82954
m64109_200309_192110.fastq	8.58724
m64109_200311_013444.fastq	8.11014

42 x coverage

```

DeepPolisher submitted in batch: https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/hprc_DeepPolisher/GIAB_samples_manuscript

#### Sequel - DCv1.1

Redownloading original winnowmap alignments Mobin made: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/04ab5ae4-0170-11ee-904c-0a13c5208311--HPRC_Polishing/HPRC_Y1/HG002/hifi/winnowmap/

Original source fastq files:
```
https://storage.googleapis.com/brain-genomics/awcarroll/share/t2t/v1.1RC/m64011_190830_220126.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/awcarroll/share/t2t/v1.1RC/m64011_190901_095311.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/awcarroll/share/t2t/v1.1RC/m64012_190920_173625.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/awcarroll/share/t2t/v1.1RC/m64015_190920_185703.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/awcarroll/share/t2t/v1.1RC/m64012_190921_234837.dc.q20.fastq.gz
https://storage.googleapis.com/brain-genomics/awcarroll/share/t2t/v1.1RC/m64015_190922_010918.dc.q20.fastq.gz
```
```
/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1
```

Check coverage

```
#!/bin/bash
#SBATCH --job-name=HG002_mosdepth
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -n --fast-mode -t 4 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/mosdepth/HG002 \
    /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.bam

python3 ~/progs/mosdepth/scripts/plot-dist.py /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/mosdepth/HG002.mosdepth.global.dist.txt
```

Downsample from 66x to 40x
```
#!/bin/bash
#SBATCH --job-name=downsampling_HG002_DCv1.1
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

samtools view -s 0.65 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Sequel_Revio_DCv1.1/HG002_DC_v1.1.0_Q20.hifi_dc_1.1.winnowmap_v2.03.40x.bam

```

## HG002 HiFi Revio DCv1.1 from pacbio

https://downloads.pacbcloud.com/public/revio/2022Q4/

Check read coverage for HG002 data
```
{
  "runReadStats.dockerImage": "mobinasri/fai_read_stats:latest",
  "runReadStats.reads": ["/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/m84005_220919_232112_s2.hifi_reads.bam","/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/m84011_220902_175841_s1.hifi_reads.bam"]
}
```

```
#!/bin/bash
#SBATCH --job-name=read_stats_Hg2_DCv1.1_revio
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

cd /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=high_priority"
export TOIL_COORDINATION_DIR=/data/tmp


LOCAL_FOLDER=/data/tmp/$(whoami)/read_stats_HG002
mkdir -p ${LOCAL_FOLDER}

toil-wdl-runner \
    --jobStore ${LOCAL_FOLDER}/jobstore \
    --batchSystem single_machine \
    --maxCores 32 \
    --batchLogsDir /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/toil_logs \
    ~/progs/master_hpp_production_workflows/QC/wdl/tasks/read_stats.wdl \
    /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/read_stats.inputs.json  \
    --outputDirectory /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/read_stats_outputs \
    --outputFile /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/read_stats_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    2>&1 | tee log.txt

toil clean ${LOCAL_FOLDER}/jobstore

```
Get coverage
```
/private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/read_stats_outputs

grep "total_reads" sample_all.report.tsv > coverage_per_file.tsv
grep "N50" sample_all.report.tsv | paste -d "\t" coverage_per_file.tsv - > tmp; mv tmp coverage_per_file.tsv
awk '{ print $1"\t"$3*$6 / 3000000000 }' coverage_per_file.tsv

m84005_220919_232112_s2.hifi_reads.fq	33.8921
m84011_220902_175841_s1.hifi_reads.fq	36.9876
sample.fastq	70.8804
```
Downsample coverage
```
#!/bin/bash
#SBATCH --job-name=downsampling_HG002_Revio_DCv1.1
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

# 70x (delete after)
samtools merge -@32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_70x.bam /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/m84011_220902_175841_s1.hifi_reads.bam /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/m84005_220919_232112_s2.hifi_reads.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_70x.bam

samtools view -s 0.86 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_70x.bam > /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_60x.bam

samtools view -s 0.72 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_70x.bam > /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_50x.bam

samtools view -s 0.6 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_70x.bam > /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_40x.bam

samtools view -s 0.43 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_70x.bam > /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_30x.bam

samtools view -s 0.29 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_70x.bam > /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_20x.bam

samtools view -s 0.15 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_70x.bam > /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/HiFi_Revio_DCv1.1/pacbio_2022Q4/HG002_unaligned_HiFi_Revio_DCv1.1_10x.bam

```
