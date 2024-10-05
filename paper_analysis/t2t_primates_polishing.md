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

### Optimizing GQ filters for primates

https://github.com/miramastoras/DeepPolisher_manuscript/blob/main/scripts/Optimize_GQ_filters_T2T_primates_verkko_model2.ipynb

### Apply GQ filtered polishing edits to primates

Filter vcf
```
for sample in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1; do \
    UNFILT=/private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/${sample}_verkko_model2/analysis/DeepPolisher_outputs/polisher_output.no_filters.vcf.gz

    bcftools view -Oz -i 'FORMAT/GQ>24 && (ILEN = 1)' ${UNFILT} > /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_INS1_GQ24.vcf.gz
    tabix -p vcf /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_INS1_GQ24.vcf.gz

    bcftools view -Oz -i 'FORMAT/GQ>21 && (ILEN = -1)' ${UNFILT} > /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_DEL1_GQ21.vcf.gz
    tabix -p vcf /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_DEL1_GQ21.vcf.gz

    bcftools view -Oz -e 'FORMAT/GQ<=13 || (ILEN = 1) || (ILEN = -1)' ${UNFILT} > /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_GQ13.vcf.gz
    tabix -p vcf /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_GQ13.vcf.gz

    bcftools concat -a -Oz /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_INS1_GQ24.vcf.gz \
    /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_DEL1_GQ21.vcf.gz \
    /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_GQ13.vcf.gz \
    > /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_DP_primate_filters.vcf.gz
  done

for sample in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1; do \
    tabix -p vcf /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_DP_primate_filters.vcf.gz
  done
```

Apply polishing edits
```
for sample in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1; do \
    echo $sample
    bcftools consensus -f /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/diploid/${sample}.dip.20230906.fasta.gz -H 2 /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/verkko_model2_primate_optimized_filters/${sample}_DP_primate_filters.vcf.gz > /private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/${sample}_verkko_model2/${sample}_verkko_model2.T2T_GQ_DP_polished.dip.fasta
  done

# mPanPan1
# Applied 3841 variants
# mPanTro3
# Applied 2193 variants
# mPonAbe1
# Applied 2205 variants
# mPonPyg2
# Applied 2777 variants
# mSymSyn1
# Applied 3224 variants
```

Count edits with HPRC filters
```
for sample in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1; do \
    HPRCFILT=/private/groups/patenlab/mira/t2t_primates_polishing/DeepPolisher/${sample}_verkko_model2/analysis/DeepPolisher_outputs/polisher_output.vcf.gz
    count=`zcat $HPRCFILT | grep -v "^#" | wc -l`
    echo ${sample},${count}
    done
```
### Run Merqury QV stratifications

Run mosdepth on bam files
```
#!/bin/bash
#SBATCH --job-name=mosdepth_primate
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00

ASM=$1
BAM=$2
SAMPLE=$3

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -Q 1 --threads 4 \
    -f ${ASM} \
    --quantize 0:5:10:150: \
    /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/mosdepth/${SAMPLE}_mosdepth \
    ${BAM}

zcat /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/mosdepth/${SAMPLE}_mosdepth.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/mosdepth/${SAMPLE}_mosdepth_quantized.tsv

awk 'BEGIN {FS=OFS="\t"} {
        if ($4 == "0:5") {
            $4 = "NO_COVERAGE"
        } else if ($4 == "5:10") {
            $4 = "LOW_COVERAGE"
        } else if ($4 == "10:150") {
            $4 = "CALLABLE"
        } else if ($4 == "150:inf") {
            $4 = "HIGH_COVERAGE"
        }
        print }' /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/mosdepth/${SAMPLE}_mosdepth_quantized.tsv > /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/mosdepth/${SAMPLE}_mosdepth_quantized.quant.bed

grep NO_COVERAGE /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/mosdepth/${SAMPLE}_mosdepth_quantized.quant.bed > /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/mosdepth/${SAMPLE}_mosdepth_quantized.quant.lt5x_cov.bed

awk -v OFS="\t" '{if ($3 -$2 >100) print $1,$2,$3 }' /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/mosdepth/${SAMPLE}_mosdepth_quantized.quant.bed  > /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/mosdepth/${SAMPLE}_mosdepth_quantized.quant.lt5x_cov.gt100bp.MAPQ1.bed

```

```
cd /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/mosdepth

for sample in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1; do
  echo sbatch slurm.sh /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/diploid/${sample}.dip.20230906.fasta.gz /private/groups/patenlab/mira/t2t_primates_polishing/hprc_DeepPolisher/${sample}/analysis/hprc_DeepPolisher_outputs/${sample}.hifi.to.diploid.asm.PHARAOH.bam ${sample}
done
```

#### Re-polishing primates with higher coverage

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
downsample=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $9}' "${sample_file}")

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

samtools view -s ${downsample} -b -h -@ 32 ${LOCAL_FOLDER}/${HIFI_NAME} > ${OUTDIR}/${HIFI_NAME}.60x.bam
samtools index ${OUTDIR}/${HIFI_NAME}.60x.bam

rm -rf /data/tmp/$(whoami)/${sample_id}
```

```sh
cd /private/groups/patenlab/mira/t2t_primates_polishing/reads

mkdir -p slurm_logs/
sbatch downsample_hifi_60x.sh T2T_primates_all_manuscript.csv
```

List files for spreasheet:
```
cd /private/groups/patenlab/mira/t2t_primates_polishing/reads/hifi

for sample in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 ; do realpath ${sample}/*.60x.bam ; done
for sample in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 ; do realpath ${sample}/*.60x.bam.bai ; done
```
