## This document contains notes for polishing the t2t primates assemblies

### Check read coverages for all data

Location of data in this spreadsheet: https://docs.google.com/spreadsheets/d/17HVGUJ7cvf8SvxkLNVg4cGcIvDTaDCJvnGf56CT27d4/edit#gid=437306236

phoenix batch submission for read stats: https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/read_stats/t2t_primates

```sh
for data in hifi ont ilm ; do
    for sample in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 ; do \
      grep "total_reads" ${data}/${sample}/analysis/read_stats_outputs/sample_all.report.tsv > ${data}/${sample}/analysis/read_stats_outputs/coverage_per_file.tsv; \
      grep "N50" ${data}/${sample}/analysis/read_stats_outputs/sample_all.report.tsv | paste -d "\t" ${data}/${sample}/analysis/read_stats_outputs/coverage_per_file.tsv - > ${data}/${sample}/analysis/read_stats_outputs/tmp; mv ${data}/${sample}/analysis/read_stats_outputs/tmp ${data}/${sample}/analysis/read_stats_outputs/coverage_per_file.tsv; \
      cov=`awk '{ print $1"\t"$3*$6 / 3000000000 }' ${data}/${sample}/analysis/read_stats_outputs/coverage_per_file.tsv | cut -f2 | head -n 1`; \
      echo ${sample} ${data} ${cov} >> coverages_listed.csv ;done ;done
```

Added coverages to column in spreadsheet.

Download illumina files and check coverage with mosdepth
```sh
#!/bin/bash

#SBATCH --job-name=t2t_primates_ilm_mosdepth
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=200gb
#SBATCH --threads-per-core=1
#SBATCH --output=slurm_logs/downsample_%x_%j_%A_%a.log
#SBATCH --time=12:00:00
#SBATCH --array=[1-6]%6

set -ex

## Pull samples names from CSV passed to script
sample_file=$1

# Read the CSV file and extract the sample ID for the current job array task
# Skip first row to avoid the header
sample_id=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $1}' "${sample_file}")
ilm_reads=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $6}' "${sample_file}" | sed 's/\"//g' | sed 's/[][]//g')

# Ensure a sample ID is obtained
if [ -z "${sample_id}" ]; then
    echo "Error: Failed to retrieve a valid sample ID. Exiting."
    exit 1
fi

echo "${sample_id}"

OUTDIR=/private/groups/patenlab/mira/t2t_primates_polishing/reads/ilm/${sample_id}
mkdir -p ${OUTDIR}
cd ${OUTDIR}

source /private/home/mmastora/progs/miniconda3/etc/profile.d/conda.sh
conda activate awscli

## make folder on local node for data
aws s3 cp ${ilm_reads} ${OUTDIR}
aws s3 cp ${ilm_reads}.bai ${OUTDIR}

ILM_NAME=`basename ${ilm_reads}`

mkdir -p ${OUTDIR}/mosdepth/

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v ${OUTDIR}/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -n --fast-mode -t 4 ${OUTDIR}/mosdepth/${sample_id} \
    ${OUTDIR}/${ILM_NAME}

source /private/home/mmastora/progs/miniconda3/etc/profile.d/conda.sh
conda activate analysis

python3 ~/progs/mosdepth/scripts/plot-dist.py ${OUTDIR}/mosdepth/${sample_id}.mosdepth.global.dist.txt &> ${OUTDIR}/mosdepth/coverages.txt
```

```
  cd /private/groups/patenlab/mira/t2t_primates_polishing/reads/ilm/

mkdir -p slurm_logs/

sbatch ilm_coverage.slurm.sh /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/hprc_DeepPolisher/T2T_primates/T2T_primates_all_manuscript.csv

ls | grep "m*" | while read line ; do tail -n1 ${line}/mosdepth/coverages.txt; done
```

Only need to downsample mSymSyn1
```sh
#!/bin/bash

#SBATCH --job-name=mSymSyn1_downsample_ilm
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=200gb
#SBATCH --threads-per-core=1
#SBATCH --time=12:00:00

samtools view -s .6 -b -h -@ 32 /private/groups/patenlab/mira/t2t_primates_polishing/reads/ilm/mSymSyn1/bwa.analysis-dip.illumina.dedup.bam > /private/groups/patenlab/mira/t2t_primates_polishing/reads/ilm/mSymSyn1/bwa.analysis-dip.illumina.dedup.30x.bam

samtools index /private/groups/patenlab/mira/t2t_primates_polishing/reads/ilm/mSymSyn1/bwa.analysis-dip.illumina.dedup.30x.bam
```

### Downsampling and downloading hifi data

```sh
#!/bin/bash

#SBATCH --job-name=t2t_primates_hifi_downsample
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=200gb
#SBATCH --threads-per-core=1
#SBATCH --output=slurm_logs/downsample_%x_%j_%A_%a.log
#SBATCH --time=12:00:00
#SBATCH --array=[1-6]%6

set -ex

## Pull samples names from CSV passed to script
sample_file=$1

# Read the CSV file and extract the sample ID for the current job array task
# Skip first row to avoid the header
sample_id=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $1}' "${sample_file}")
hifi_reads=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $4}' "${sample_file}" | sed 's/\"//g' | sed 's/[][]//g')
downsample=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $8}' "${sample_file}")

echo $hifi_reads

# Ensure a sample ID is obtained
if [ -z "${sample_id}" ]; then
    echo "Error: Failed to retrieve a valid sample ID. Exiting."
    exit 1
fi

echo "${sample_id}"

OUTDIR=/private/groups/patenlab/mira/t2t_primates_polishing/reads/hifi/${sample_id}
mkdir -p ${OUTDIR}
cd ${OUTDIR}

source /private/home/mmastora/progs/miniconda3/etc/profile.d/conda.sh
conda activate awscli

## make folder on local node for data
LOCAL_FOLDER=/data/tmp/$(whoami)/${sample_id}
mkdir -p ${LOCAL_FOLDER}

aws s3 cp ${hifi_reads} ${LOCAL_FOLDER}
aws s3 cp ${hifi_reads}.bai ${LOCAL_FOLDER}

HIFI_NAME=`basename ${hifi_reads}`

samtools view -s ${downsample} -b -h -@ 32 ${LOCAL_FOLDER}/${HIFI_NAME} > ${OUTDIR}/${HIFI_NAME}.40x.bam
samtools index ${OUTDIR}/${HIFI_NAME}.40x.bam

rm -rf /data/tmp/$(whoami)/${sample_id}
```

```sh
cd /private/groups/patenlab/mira/t2t_primates_polishing/reads

mkdir -p slurm_logs/
sbatch downsample_hifi.sh T2T_primates_all_manuscript.csv
```

List files for spreasheet:
```
cd /private/groups/patenlab/mira/t2t_primates_polishing/reads/hifi

for sample in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 ; do realpath ${sample}/*.bam ; done
for sample in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 ; do realpath ${sample}/*.bai ; done
```

Download hap1 and hap2 raw fasta files
```
cd /private/groups/patenlab/mira/t2t_primates_polishing/assemblies

conda activate awscli

cut -f 1 -d"," /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/hprc_DeepPolisher/T2T_primates/T2T_primates_all_manuscript.csv | grep -v "sample_id" | while read line ; do
    mkdir -p ${line}
    hap1_fa=`grep ${line} /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/hprc_DeepPolisher/T2T_primates/T2T_primates_all_manuscript.csv | cut -f2 -d"," | sed 's/.fasta.gz/.fasta/g'`
    hap2_fa=`grep ${line} /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/hprc_DeepPolisher/T2T_primates/T2T_primates_all_manuscript.csv | cut -f3 -d"," | sed 's/.fasta.gz/.fasta/g'`
    hap1_fai=`grep ${line} /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/hprc_DeepPolisher/T2T_primates/T2T_primates_all_manuscript.csv | cut -f2 -d"," | sed 's/.fasta.gz/.fasta.fai/g'`
    hap2_fai=`grep ${line} /private/groups/patenlab/mira/phoenix_batch_submissions/polishing/hprc_DeepPolisher/T2T_primates/T2T_primates_all_manuscript.csv | cut -f3 -d"," | sed 's/.fasta.gz/.fasta.fai/g'`

    aws s3 cp ${hap1_fa} ${line}/
    aws s3 cp ${hap2_fa} ${line}/
    aws s3 cp ${hap1_fai} ${line}/
    aws s3 cp ${hap2_fai} ${line}/
  done
```


Submitted to deeppolisher https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/hprc_DeepPolisher/T2T_primates
