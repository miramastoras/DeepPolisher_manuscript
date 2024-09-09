cd /private/groups/migalab/juklucas/polishing

###############################################################################
##                       create preprocess input jsons                       ##
###############################################################################

deactivate
conda activate pandas_env


mkdir verkko_preprocess_input_jsons
cd verkko_preprocess_input_jsons

cp /private/groups/hprc/verkko/batch1/verkko_preprocess_input_jsons/verkko_preprocess_input_mapping.csv ./

## Note that I put in dummy data for HiC since we don't need it
python3 /private/groups/hprc/hprc_intermediate_assembly/hpc/launch_from_table.py \
     --data_table ../giab_polishing.csv \
     --field_mapping verkko_preprocess_input_mapping.csv \
     --workflow_name verkko_trio_plus_hic_preprocess


###############################################################################
##                             launch preprocess                             ##
###############################################################################

## on HPC...
cd /private/groups/migalab/juklucas/polishing

## check that github repos are up to date
git -C /private/groups/hprc/hprc_intermediate_assembly pull
git -C /private/home/juklucas/github/hpp_production_workflows/ pull

mkdir slurm_logs

sbatch \
     --job-name=HPRC-verkko-polishing \
     --array=[1-2] \
     --cpus-per-task=64 \
     --mem=200gb \
     --partition=high_priority \
     /private/groups/hprc/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
     --wdl /private/home/juklucas/github/hpp_production_workflows/assembly/wdl/tasks/verkko_trio_plus_hic_preprocess.wdl \
     --sample_csv giab_polishing.csv \
     --input_json_path '../verkko_preprocess_input_jsons/${SAMPLE_ID}_verkko_trio_plus_hic_preprocess.json' 


###############################################################################
##                       recreate preprocess input jsons for HG002           ##
###############################################################################

## We realized that Mobin ran Hifiasm with 300X for the yak DBs so it's probably
## best to run Verkko with the same data.

## Download the 300X data and put into one file (since each sample has 1800 
## files or so) then rerun preprocessing.


cd /private/groups/migalab/juklucas/polishing/HG002/inputs

mkdir ilmn_HG002/
cd ilmn_HG002/

srun \
  --job-name "HG002" \
  --cpus-per-task 4 \
  --mem 20G \
  --time 24:00:00 \
  --partition high_priority \
  --pty bash \
  -i

wget -i ../HG002_ilmn_ftp.txt
cd ..

mkdir ilmn_HG003/
cd ilmn_HG003/

srun \
  --job-name "HG003" \
  --cpus-per-task 4 \
  --mem 20G \
  --time 24:00:00 \
  --partition high_priority \
  --pty bash \
  -i

wget -i ../HG003_ilmn_ftp.txt
cd ..

mkdir ilmn_HG004/
cd ilmn_HG004/

srun \
  --job-name "HG004" \
  --cpus-per-task 4 \
  --mem 20G \
  --time 24:00:00 \
  --partition high_priority \
  --pty bash \
  -i

wget -i ../HG004_ilmn_ftp.txt
cd .. 

srun \
  --job-name "HG002" \
  --cpus-per-task 24 \
  --mem 60G \
  --time 24:00:00 \
  --partition high_priority \
  --pty bash \
  -i

zcat ilmn_HG002/* | pigz > HG002_300X_ilmn.fastq.gz &
zcat ilmn_HG003/* | pigz > HG003_300X_ilmn.fastq.gz &
zcat ilmn_HG004/* | pigz > HG004_300X_ilmn.fastq.gz &


###############################################################################
##                             Run Preprocessing                             ##
############################################################################### 

## Run WDL that localizes the data and creates meryl DBs
cd /private/groups/migalab/juklucas/polishing/verkko_preprocess_input_jsons

## Note that I put in dummy data for HiC since we don't need it
python3 /private/groups/hprc/hprc_intermediate_assembly/hpc/launch_from_table.py \
     --data_table ../giab_polishing.csv \
     --field_mapping verkko_preprocess_input_mapping.csv \
     --workflow_name verkko_trio_plus_hic_preprocess


cd /private/groups/migalab/juklucas/polishing

sbatch \
     --job-name=HPRC-verkko-polishing \
     --array=[2] \
     --cpus-per-task=64 \
     --mem=200gb \
     --partition=high_priority \
     /private/groups/hprc/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
     --wdl /private/home/juklucas/github/hpp_production_workflows/assembly/wdl/tasks/verkko_trio_plus_hic_preprocess.wdl \
     --sample_csv giab_polishing.csv \
     --input_json_path '../verkko_preprocess_input_jsons/${SAMPLE_ID}_verkko_trio_plus_hic_preprocess.json' 


###############################################################################
##                             Fix Data Problem.                             ##
############################################################################### 

## I got the following error from Toil for HG002
# gzip: ilmn_HG003/3F2_GCCAAT_L001_R1_003.fastq.gz.11: invalid compressed data--format violated

cd /private/groups/migalab/juklucas/polishing/HG002/inputs

grep "3F2_GCCAAT_L001_R1_003.fastq.gz" HG003_ilmn_ftp.txt > HG003_ilmn_failures_ftp.txt

cd ilmn_HG003
rm 3F2_GCCAAT_L001_R1_003.fastq.gz*
wget -i ../HG003_ilmn_failures_ftp.txt

# Function to check if a gzip file is corrupted
check_gzip_file() {
    local file_path="$1"
    if gzip -t "$file_path" 2>/dev/null; then
        echo "$file_path is not corrupted."
    else
        echo "$file_path is corrupted."
    fi
}

# Iterate through each file and check for corruption
for file_path in 3F2_GCCAAT_L001_R1_003.fastq.gz*; do
    #echo $file_path
    check_gzip_file "$file_path"
done

## looks good, run concatenation again...
srun \
  --job-name "HG002" \
  --cpus-per-task 24 \
  --mem 60G \
  --time 24:00:00 \
  --partition high_priority \
  --pty bash \
  -i

cat ilmn_HG003/* > HG003_300X_ilmn.fastq.gz &

###############################################################################
##                      Run Preprocessing For HG002 Again                    ##
############################################################################### 

cd /private/groups/migalab/juklucas/polishing

## get rid of prior results
rm -rf  HG002/analysis/

sbatch \
     --job-name=HPRC-verkko-polishing \
     --array=[1] \
     --cpus-per-task=64 \
     --mem=200gb \
     --partition=high_priority \
     /private/groups/hprc/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
     --wdl /private/home/juklucas/github/hpp_production_workflows/assembly/wdl/tasks/verkko_trio_plus_hic_preprocess.wdl \
     --sample_csv giab_polishing.csv \
     --input_json_path '../verkko_preprocess_input_jsons/${SAMPLE_ID}_verkko_trio_plus_hic_preprocess.json' 

###############################################################################
##                          Launch Verkko (Trio) For HG002                   ##
###############################################################################     

cd /private/groups/migalab/juklucas/polishing/HG002

## Don't launch from the head node...
srun \
  --job-name "HG002" \
  --cpus-per-task 24 \
  --mem 60G \
  --time 96:00:00 \
  --partition high_priority \
  --pty bash \
  -i

## create input directories then pull in data from preprocessing WDL
mkdir -p verkko/hifi 
mkdir -p verkko/ont 
mkdir -p verkko/meryl/mat
mkdir -p verkko/meryl/pat

PREPROCESS_OUTPUT_JSON=HG002_verkko_trio_plus_hic_preprocess_outputs.json

## Soft link in ultralong reads
jq -r '.["verkko_preprocess_wf.ultralong_fq"][]' "$PREPROCESS_OUTPUT_JSON" \
| while read -r file; do
    ln -s "$file" "verkko/ont/$(basename "$file")"
done

## Soft link in hifi reads
jq -r '.["verkko_preprocess_wf.hifi_fq"][]' "$PREPROCESS_OUTPUT_JSON" \
| while read -r file; do
    ln -s "$file" "verkko/hifi/$(basename "$file")"
done

## meryl hapmers are stored as tars, extract them
PAT_HAPMER_TAR=$(jq -r '.["verkko_preprocess_wf.pat_compr_hapmer"]' "$PREPROCESS_OUTPUT_JSON")
tar xvf "${PAT_HAPMER_TAR}" --directory verkko/meryl/pat/

MAT_HAPMER_TAR=$(jq -r '.["verkko_preprocess_wf.mat_compr_hapmer"]' "$PREPROCESS_OUTPUT_JSON")
tar xvf "${MAT_HAPMER_TAR}" --directory verkko/meryl/mat/


## Actually launch...
cd verkko

source ~/.bashrc
deactivate
conda activate verkko_2.0

TRIO_ASM_DIR=trio_asm
MAT_HAPMER_DB=meryl/mat/maternal.hapmers.meryl
PAT_HAPMER_DB=meryl/pat/paternal.hapmers.meryl


verkko \
    --correct-overlap-batches 64 \
    --screen human \
    --slurm \
    --snakeopts '--cluster "/private/home/juklucas/miniconda3/envs/verkko_2.0/lib/verkko/profiles/slurm-sge-submit.sh {threads} {resources.mem_gb} {resources.time_h} {rulename} {resources.job_id}"' \
    -d $TRIO_ASM_DIR \
    --hifi hifi/*.gz \
    --nano ont/*.gz \
    --hap-kmers $MAT_HAPMER_DB $PAT_HAPMER_DB \
    trio

pigz -p 12 -c trio_asm/assembly.haplotype1.fasta > trio_asm/HG002_mat_verkko_2.0_2024_06_05.fa.gz &
pigz -p 12 -c trio_asm/assembly.haplotype2.fasta > trio_asm/HG002_pat_verkko_2.0_2024_06_05.fa.gz &


cd /private/groups/migalab/juklucas/polishing/HG002/verkko

mkdir qc 
cd qc 

###############################################################################
##                                QC HG002 Assembly                          ##
###############################################################################  

cat <<EOF > qc_inputs.json
{
  "comparisonQC.annotationBed": "/private/groups/hprc/ref_files/chm13/chm13v2.0_GenomeFeature_v1.0.bed",
  "comparisonQC.annotationCENSAT": "/private/groups/hprc/ref_files/chm13/chm13v2.0_censat_v2.0.bed",
  "comparisonQC.annotationSD": "/private/groups/hprc/ref_files/chm13/chm13v2.0_SD.flattened.bed",
  "comparisonQC.genesFasta": "/private/groups/hprc/ref_files/grch38/Homo_sapiens.GRCh38.cdna.all.fa",
  "comparisonQC.hs38Paf": "/private/groups/hprc/ref_files/grch38/hs38.paf",
  "comparisonQC.sampleName": "HG002_for_polishing",
  "comparisonQC.breakFasta": false,
  "comparisonQC.patYak": "/private/groups/patenlab/masri/hprc/polishing/HG002/parents/pat.HG003.yak",
  "comparisonQC.matYak": "/private/groups/patenlab/masri/hprc/polishing/HG002/parents/mat.HG004.yak",
  "comparisonQC.hap1Fasta": "/private/groups/migalab/juklucas/polishing/HG002/verkko/trio_asm/HG002_pat_verkko_2.0_2024_06_05.fa.gz",
  "comparisonQC.hap2Fasta": "/private/groups/migalab/juklucas/polishing/HG002/verkko/trio_asm/HG002_mat_verkko_2.0_2024_06_05.fa.gz",
  "comparisonQC.reference": "/private/groups/hprc/ref_files/chm13/chm13v2.0.fa.gz",
  "childReadsILM": [
    "/private/groups/migalab/juklucas/polishing/HG002/inputs/HG002_300X_ilmn.fastq.gz"
  ],
  "extractReadsReferenceFasta": "/private/groups/hprc/ref_files/grch38/hs38DH.fa",
  "comparisonQC.x_hap_compleasm_db": "/private/groups/hprc/ref_files/compleasm/primatesnoY_odb10.tar.gz",
  "comparisonQC.y_hap_compleasm_db": "/private/groups/hprc/ref_files/compleasm/primatesnoX_odb10.tar.gz",
  "comparisonQC.x_hap_compleasm_lineage": "primatesnoY_odb10",
  "comparisonQC.y_hap_compleasm_lineage": "primatesnoX_odb10",
  "comparisonQC.isMaleSample": true
}
EOF



mkdir -p batch_logs 
mkdir -p analysis

export SINGULARITY_CACHEDIR=`pwd`/analysis/cache/.singularity/cache 
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/analysis/cache/.cache/miniwdl 
export TOIL_SLURM_ARGS="--time=24:00:00 --partition=high_priority"
export TOIL_COORDINATION_DIR=/data/tmp


toil-wdl-runner \
    --jobStore ./bigstore \
    --batchSystem slurm \
    --batchLogsDir ./batch_logs \
    /private/home/juklucas/github/hpp_production_workflows/QC/wdl/workflows/comparison_qc.wdl \
    qc_inputs.json \
    --outputDirectory analysis \
    --outputFile qc_outputs.json \
    --runLocalJobsOnWorkers \
    --disableProgress \
    --caching=false \
    2>&1 | tee qc_log.txt  

###############################################################################
##                          Launch Verkko (Trio) For HG005                   ##
###############################################################################     

cd /private/groups/migalab/juklucas/polishing/HG005

## Don't launch from the head node...
srun \
  --job-name "HG005" \
  --cpus-per-task 24 \
  --mem 60G \
  --time 96:00:00 \
  --partition high_priority \
  --pty bash \
  -i

## create input directories then pull in data from preprocessing WDL
mkdir -p verkko/hifi 
mkdir -p verkko/ont 
mkdir -p verkko/meryl/mat
mkdir -p verkko/meryl/pat

PREPROCESS_OUTPUT_JSON=HG005_verkko_trio_plus_hic_preprocess_outputs.json

## Soft link in ultralong reads
jq -r '.["verkko_preprocess_wf.ultralong_fq"][]' "$PREPROCESS_OUTPUT_JSON" \
| while read -r file; do
    ln -s "$file" "verkko/ont/$(basename "$file")"
done

## Soft link in hifi reads
jq -r '.["verkko_preprocess_wf.hifi_fq"][]' "$PREPROCESS_OUTPUT_JSON" \
| while read -r file; do
    ln -s "$file" "verkko/hifi/$(basename "$file")"
done

## meryl hapmers are stored as tars, extract them
PAT_HAPMER_TAR=$(jq -r '.["verkko_preprocess_wf.pat_compr_hapmer"]' "$PREPROCESS_OUTPUT_JSON")
tar xvf "${PAT_HAPMER_TAR}" --directory verkko/meryl/pat/

MAT_HAPMER_TAR=$(jq -r '.["verkko_preprocess_wf.mat_compr_hapmer"]' "$PREPROCESS_OUTPUT_JSON")
tar xvf "${MAT_HAPMER_TAR}" --directory verkko/meryl/mat/


## Actually launch...
cd verkko

source ~/.bashrc
deactivate
conda activate verkko_2.0

TRIO_ASM_DIR=trio_asm
MAT_HAPMER_DB=meryl/mat/maternal.hapmers.meryl
PAT_HAPMER_DB=meryl/pat/paternal.hapmers.meryl


verkko \
    --correct-overlap-batches 64 \
    --screen human \
    --slurm \
    --snakeopts '--cluster "/private/home/juklucas/miniconda3/envs/verkko_2.0/lib/verkko/profiles/slurm-sge-submit.sh {threads} {resources.mem_gb} {resources.time_h} {rulename} {resources.job_id}"' \
    -d $TRIO_ASM_DIR \
    --hifi hifi/*.gz \
    --nano ont/*.gz \
    --hap-kmers $MAT_HAPMER_DB $PAT_HAPMER_DB \
    trio


###############################################################################
##                          Check HG006 Ilmn Coverage                        ##
###############################################################################   

cd /private/groups/migalab/juklucas/polishing/

mkdir check_hg005_ilmn
cd check_hg005_ilmn

nano HG006_inputs.json
{
  "fastqReadCoverage.sumFastqReads.diskSizeGB": 528,
  "fastqReadCoverage.sumFastqReads.out_prefix": "HG006",
  "fastqReadCoverage.sumFastqReads.inputFastq": [
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_001.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_002.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_003.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_004.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_005.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_006.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_007.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_008.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_009.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_010.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_011.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_012.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_013.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_014.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_015.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_016.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_017.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_018.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_019.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_020.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_021.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_022.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_023.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_024.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_025.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_026.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_027.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_028.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_029.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_030.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_031.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_032.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_033.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_034.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R1_035.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_001.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_002.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_003.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_004.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_005.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_006.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_007.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_008.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_009.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_010.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_011.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_012.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_013.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_014.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_015.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_016.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_017.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_018.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_019.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_020.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_021.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_022.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_023.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_024.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_025.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_026.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_027.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_028.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_029.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_030.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_031.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_032.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_033.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_034.fastq.gz",
        "s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/141008_D00360_0058_AHB675ADXX/Project_NA24694/Sample_NA24694/NA24694_GCCAAT_L001_R2_035.fastq.gz"
    ],
  "fastqReadCoverage.sumFastqReads.memSizeGB": 128
}

mkdir -p hg006_logs 
mkdir -p hg006_ilmn

export SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache 
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl 
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=high_priority"
export TOIL_COORDINATION_DIR=/data/tmp


toil-wdl-runner \
    --jobStore ./hg006_bigstore \
    --batchSystem slurm \
    --batchLogsDir ./hg006_logs \
    /private/home/juklucas/github/hpp_production_workflows/data_processing/wdl/tasks/fastqReadCoverage.wdl \
    HG006_inputs.json \
    --outputDirectory hg006_ilmn \
    --outputFile HG006_outputs.json \
    --runLocalJobsOnWorkers \
    --disableProgress \
    --caching=false \
    2>&1 | tee hg006_log.txt  

cat hg006_ilmn/2764ddde-b19e-43c5-aa99-e0bebf3f5a6f/HG006.txt
# 41326288048

## the coverage is only 15X!!! We need higher coverage and Mobin looked into his yak
## files and found that he was actually using the full dataset, so I will do that too.

cd /private/groups/migalab/juklucas/polishing/

# rm -rf check_hg005_ilmn


###############################################################################
##                         Get Full Data Listing For HG005 Ilmn              ##
###############################################################################   

cd /private/groups/migalab/juklucas/polishing/HG005

mkdir inputs 
cd inputs 


aws s3 ls \
    s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG006/ \
    --recursive \
    --human-readable \
    --summarize \
    | grep "fastq.gz" \
    | grep -v "ubuntu" \
    | awk '{print "\"s3://human-pangenomics/"$5 "\","}' \
     > HG006_ilmn_all.txt

aws s3 ls \
    s3://human-pangenomics/working/HPRC_PLUS/HG005/raw_data/Illumina/parents/HG007/ \
    --recursive \
    --human-readable \
    --summarize \
    | grep "fastq.gz" \
    | grep -v "ubuntu" \
    | awk '{print "\"s3://human-pangenomics/"$5 "\","}' \
     > HG007_ilmn_all.txt    

###############################################################################
##                      Run Preprocessing For HG005 Again                    ##
############################################################################### 

cd /private/groups/migalab/juklucas/polishing

## get rid of prior results
rm -rf  HG005/analysis/
rm -rf  HG005/toil_logs/ HG005/HG005_verkko_trio_plus_hic_preprocess_outputs.json HG005/HG005_verkko_trio_plus_hic_preprocess_log.txt

sbatch \
     --job-name=HPRC-verkko-polishing \
     --array=[2] \
     --cpus-per-task=64 \
     --mem=200gb \
     --partition=high_priority \
     /private/groups/hprc/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
     --wdl /private/home/juklucas/github/hpp_production_workflows/assembly/wdl/tasks/verkko_trio_plus_hic_preprocess.wdl \
     --sample_csv giab_polishing.csv \
     --input_json_path '../verkko_preprocess_input_jsons/${SAMPLE_ID}_verkko_trio_plus_hic_preprocess.json' \

## try with Docker (still failed for some reason)
sbatch \
     --job-name=HPRC-verkko-polishing \
     --array=[2] \
     --partition=high_priority \
     /private/groups/hprc/hprc_intermediate_assembly/hpc/toil_sbatch_slurm.sh \
     --wdl /private/home/juklucas/github/hpp_production_workflows/assembly/wdl/tasks/verkko_trio_plus_hic_preprocess.wdl \
     --sample_csv giab_polishing.csv \
     --input_json_path '../verkko_preprocess_input_jsons/${SAMPLE_ID}_verkko_trio_plus_hic_preprocess.json' \
    --toil_args '--container docker'


###############################################################################
##                      Download HG005 Data Again Manually                   ##
############################################################################### 


srun \
  --job-name "HG004" \
  --cpus-per-task 4 \
  --mem 20G \
  --time 24:00:00 \
  --partition high_priority \
  --pty bash \
  -i

wget -i ../HG004_ilmn_ftp.txt
cd .. 

srun \
  --job-name "HG002" \
  --cpus-per-task 24 \
  --mem 60G \
  --time 24:00:00 \
  --partition high_priority \
  --pty bash \
  -i

zcat ilmn_HG002/* | pigz > HG002_300X_ilmn.fastq.gz &
zcat ilmn_HG003/* | pigz > HG003_300X_ilmn.fastq.gz &
zcat ilmn_HG004/* | pigz > HG004_300X_ilmn.fastq.gz &


###############################################################################
##                             Call Meryl Manually                           ##
###############################################################################

## given to Mira (bless her)

# srun \
#   --job-name "meryl" \
#   --cpus-per-task 64 \
#   --partition high_priority \
#   --mem 250G \
#   --time 24:00:00 \
#   --pty bash \
#   -i

# cd /private/groups/migalab/juklucas/polishing/HG005/inputs


# docker run \
#     -it \
#     -v /private/groups/patenlab/mira/hprc_polishing/data/reads/:/reads/ \
#     -v $(pwd):/data/ \
#     -u 30145:620 \
#     juklucas/hpp_merqury@sha256:387bb58c820c5825e8cf29bc0f8678975ae6a9e6d350c56d5eb5865a6d3d5e82 \
#     /bin/bash


# meryl count \
#     compress \
#     k=30 \
#     threads=42 \
#     memory=64 \
#     /reads/HG007/HG007_ilm_100x.fastq.gz \
#     output maternal_compress.k30.meryl \
#     >maternal_merqury.stdout 2>maternal_merqury.stderr

# meryl count \
#     compress \
#     k=30 \
#     threads=42 \
#     memory=64 \
#     /reads/HG006/HG006_ilm_100x.fastq.gz \
#     output paternal_compress.k30.meryl \
#     >paternal_merqury.stdout 2>paternal_merqury.stderr

# meryl count \
#     compress \
#     k=30 \
#     threads=42 \
#     memory=64 \
#     /reads/HG005/illumina/HG005.ilm.fastq.gz \
#     output child_compress.k30.meryl \
#     >child_merqury.stdout 2>child_merqury.stderr

# $MERQURY/trio/hapmers.sh \
#     maternal_compress.k30.meryl \
#     paternal_compress.k30.meryl \
#     child_compress.k30.meryl \
#     >hapmer.stdout 2>hapmer.sterr

###############################################################################
##                        launch preprocess for HG005 again                  ##
###############################################################################

## on HPC...
cd /private/groups/migalab/juklucas/polishing

## I need to download the files for HG005 again...
## make a backup of the old input json (just in case) then remove most illumina
## files because we won't actually use them for hapmer creation.
cp verkko_preprocess_input_jsons/HG005_verkko_trio_plus_hic_preprocess.json \
    verkko_preprocess_input_jsons/HG005_verkko_trio_plus_hic_preprocess_backup.json

mkdir -p slurm_logs

sbatch \
     --job-name=HPRC-verkko-polishing \
     --array=[2] \
     --cpus-per-task=64 \
     --mem=200gb \
     --partition=high_priority \
     /private/groups/hprc/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
     --wdl /private/home/juklucas/github/hpp_production_workflows/assembly/wdl/tasks/verkko_trio_plus_hic_preprocess.wdl \
     --sample_csv giab_polishing.csv \
     --input_json_path '../verkko_preprocess_input_jsons/${SAMPLE_ID}_verkko_trio_plus_hic_preprocess.json' 


###############################################################################
##                             Relaunch HG005 To Rephase...                  ##
###############################################################################

cd /private/groups/migalab/juklucas/polishing/HG005/

## Don't launch from the head node...
srun \
  --job-name "HG005" \
  --cpus-per-task 48 \
  --mem 300G \
  --time 96:00:00 \
  --partition high_priority \
  --pty bash \
  -i


PREPROCESS_OUTPUT_JSON=HG005_verkko_trio_plus_hic_preprocess_outputs.json

## Soft link in ultralong reads
jq -r '.["verkko_preprocess_wf.ultralong_fq"][]' "$PREPROCESS_OUTPUT_JSON" \
| while read -r file; do
    ln -s "$file" "verkko/ont/$(basename "$file")"
done

## Soft link in hifi reads
jq -r '.["verkko_preprocess_wf.hifi_fq"][]' "$PREPROCESS_OUTPUT_JSON" \
| while read -r file; do
    ln -s "$file" "verkko/hifi/$(basename "$file")"
done

## make the files look old so snakemake doesn't restart verkko
touch -a -m -t 202001011205.02 verkko/hifi/*
touch -a -m -t 202001011205.02 verkko/ont/*

find /private/groups/migalab/juklucas/polishing/HG005/analysis/verkko_trio_plus_hic_preprocess_outputs/ -type f -exec touch -a -m -t 202001011205.02 {} +


source ~/.bashrc
deactivate
conda activate verkko_2.0

TRIO_ASM_DIR=trio_asm

rm -rf trio_asm/6-rukki
rm -rf trio_asm/7-consensus/

verkko \
    -d trio_asm \
    --hifi hifi/*.gz \
    --nano ont/*.gz \
    --correct-overlap-batches 64 \
    --local-memory 300 \
    --local-cpus 48 \
    --screen human \
    --snakeopts "-R rukki" \
    --hap-kmers /private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/illumina/maternal_compress.k30.only.meryl /private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/illumina/paternal_compress.k30.only.meryl trio
    


pigz -p 12 -c trio_asm/assembly.haplotype1.fasta > trio_asm/HG005_mat_verkko_2.0_2024_07_29.fa.gz &
pigz -p 12 -c trio_asm/assembly.haplotype2.fasta > trio_asm/HG005_pat_verkko_2.0_2024_07_29.fa.gz &


mkdir qc 
cd qc 


cat <<EOF > qc_inputs.json
{
  "comparisonQC.annotationBed": "/private/groups/hprc/ref_files/chm13/chm13v2.0_GenomeFeature_v1.0.bed",
  "comparisonQC.annotationCENSAT": "/private/groups/hprc/ref_files/chm13/chm13v2.0_censat_v2.0.bed",
  "comparisonQC.annotationSD": "/private/groups/hprc/ref_files/chm13/chm13v2.0_SD.flattened.bed",
  "comparisonQC.genesFasta": "/private/groups/hprc/ref_files/grch38/Homo_sapiens.GRCh38.cdna.all.fa",
  "comparisonQC.hs38Paf": "/private/groups/hprc/ref_files/grch38/hs38.paf",
  "comparisonQC.sampleName": "HG005_for_polishing",
  "comparisonQC.breakFasta": false,
  "comparisonQC.patYak": "/private/groups/patenlab/masri/hprc/polishing/HG005/parents/pat.HG006.yak",
  "comparisonQC.matYak": "/private/groups/patenlab/masri/hprc/polishing/HG005/parents/mat.HG007.yak",
  "comparisonQC.hap1Fasta": "/private/groups/migalab/juklucas/polishing/HG005/verkko/trio_asm/HG005_pat_verkko_2.0_2024_07_29.fa.gz",
  "comparisonQC.hap2Fasta": "/private/groups/migalab/juklucas/polishing/HG005/verkko/trio_asm/HG005_mat_verkko_2.0_2024_07_29.fa.gz",
  "comparisonQC.reference": "/private/groups/hprc/ref_files/chm13/chm13v2.0.fa.gz",
  "childReadsILM": [
    "/private/groups/patenlab/mira/hprc_polishing/data/reads/HG005/illumina/HG005.ilm.fastq.gz"
  ],
  "extractReadsReferenceFasta": "/private/groups/hprc/ref_files/grch38/hs38DH.fa",
  "comparisonQC.x_hap_compleasm_db": "/private/groups/hprc/ref_files/compleasm/primatesnoY_odb10_mb_download.tar.gz",
  "comparisonQC.y_hap_compleasm_db": "/private/groups/hprc/ref_files/compleasm/primatesnoX_odb10_mb_download.tar.gz",
  "comparisonQC.x_hap_compleasm_lineage": "primatesnoY_odb10",
  "comparisonQC.y_hap_compleasm_lineage": "primatesnoX_odb10",
  "comparisonQC.isMaleSample": true
}
EOF



mkdir -p batch_logs 
mkdir -p analysis

export TOIL_SLURM_ARGS="--time=24:00:00 --partition=high_priority"
export TOIL_COORDINATION_DIR=/data/tmp


toil-wdl-runner \
    --jobStore ./bigstore \
    --batchSystem slurm \
    --batchLogsDir ./batch_logs \
    /private/home/juklucas/github/hpp_production_workflows/QC/wdl/workflows/comparison_qc.wdl \
    qc_inputs.json \
    --outputDirectory analysis \
    --outputFile qc_outputs.json \
    --runLocalJobsOnWorkers \
    --disableProgress \
    --caching=false \
    2>&1 | tee qc_log.txt  
