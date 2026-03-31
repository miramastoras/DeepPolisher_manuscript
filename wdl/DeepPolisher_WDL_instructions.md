## Mira's procedure for running DeepPolisher and QV WDLs
3/31/2026

### 1. Download Kishwar's model files to phoenix and put them in a tarball:

The WDL requires them to be tar'd in a folder called "checkpoint" and the name of the tar file needs to match
the prefix for the model files.

```sh
mkdir -p /private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint/
gsutil cp gs://brain-genomics-public/research/polisher-models/hybrid_ont_hifi_coverage_60/* . /private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint

cd /private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/
tar -zcvf checkpoint-460.tar.gz checkpoint/
```

I keep a list of all deeppolisher models here: https://docs.google.com/spreadsheets/d/1RuTgdq1HHeBAAMTr7QTV1m2bqwtenBoMiA5w2jtP8IE/edit?gid=0#gid=0

### 2. Create a github repo for running WDL pipelines in batch

Here is my example: https://github.com/miramastoras/phoenix_batch_submissions/tree/main/workflows

Every folder is a WDL workflow name. Within each folder there are four items:

```sh
ls phoenix_batch_submissions/workflows/DeepPolisher

DeepPolisher.csv
DeepPolisher_input_jsons
DeepPolisher_input_mapping.csv
run_DeepPolisher.sh
```
In this example {WDL_NAME} is DeepPolisher.

`{WDL_NAME}.csv` is a csv file containing the inputs for each pipeline execution, with one line per run.

`{WDL_NAME}_input_jsons` contains input jsons for each pipeline run, autopopulated from the CSV.

`{WDL_NAME}_input_mapping.csv` Gives the mapping of the CSV columns to the input parameters for the WDL. Basically the instructions on how to autopopulate the input files.

`run_{WDL_NAME}.sh` bash code for autopopulating the input jsons, and submitting the workflow on phoenix.

For this example, here is the DeepPolisher WDL we are using: https://github.com/miramastoras/hpp_production_workflows/blob/master/QC/wdl/tasks/DeepPolisher.wdl

### 3. First, set up input_mapping.csv based on the WDL you are going to run.

Here is what DeepPolisher_input_mapping.csv looks like:
```sh
input,type,value
runDeepPolisher.useOptimalGQFilter,scalar,true
runDeepPolisher.ModelFilesTarGZ,scalar,$input.DeepPolisherModelFiles
runDeepPolisher.Bai,scalar,$input.Bai
runDeepPolisher.Bam,scalar,$input.Bam
runDeepPolisher.dockerImage,scalar,$input.DeepPolisherDocker
runDeepPolisher.sampleName,scalar,$input.sample_id
runDeepPolisher.Fasta,scalar,$input.Fasta
runDeepPolisher.diskSize,scalar,600
runDeepPolisher.memSize,scalar,256
```

The first column "input" has the names of the inputs exactly as listed in the WDL. Here's the first few lines:
```sh
workflow runDeepPolisher {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Run DeepPolisher on phased reads aligned to one haplotype"
    }
    input {
        File Bam
        File Bai
        File Fasta
        File ModelFilesTarGZ
        String sampleName
        String dockerImage
        Boolean useOptimalGQFilter=true
        String customGQFilter=""
        Int memSize=128
        Int diskSize=128
        Int threadCount=32
    }
}
```
The second column "type" refers to whether the input parameter is a single value or input file "scalar" or an array of input values/files. scalar and array are the two options.

You can tell by looking at the input definition in the WDL

```sh
# scalar
File ModelFilesTarGZ
String sampleName
Boolean useOptimalGQFilter=true
Int memSize=128
# array
Array inputFastqs
```

The third column "value" tells which column from the input CSV to find that input parameter in.

```sh
# This says the input parameter "runDeepPolisher.ModelFilesTarGZ"
# should be autopopulated by the column named DeepPolisherModelFiles
# from DeepPolisher.csv
runDeepPolisher.ModelFilesTarGZ,scalar,$input.DeepPolisherModelFiles
```

> if you are setting this up for a new WDL, you can run
 ```
 java -jar /private/home/mmastora/progs/womtool-85.jar inputs DeepPolisher.wdl
 ```
 >to print the required inputs to the WDL in json format.


### 4. Now set up a Google Docs csv file to store input parameters

Example: https://docs.google.com/spreadsheets/d/1zQ1ICWvyYBqy_7L7vS6QsvlwGJyMIbwQrcZxH_Cbu1g/edit?gid=629196395#gid=629196395

I use a different tab for each WDL workflow or pipeline, and keep adding to it every time I have a new run.

### 5. Enter the full filepaths on phoenix to the required input data in the Google docs.

For the DeepPolisher model hybrid_ont_hifi_coverage_60_checkpoint_460, I added input data to column 13 here: https://docs.google.com/spreadsheets/d/1zQ1ICWvyYBqy_7L7vS6QsvlwGJyMIbwQrcZxH_Cbu1g/edit?gid=1300297667#gid=1300297667

### 6. Download the DeepPolisher tab and put it into your github repo as "DeepPolisher.csv"

```sh
ls phoenix_batch_submissions/workflows/DeepPolisher

DeepPolisher.csv # now updated
DeepPolisher_input_jsons
DeepPolisher_input_mapping.csv # now updated
run_DeepPolisher.sh
```

### 7. Follow steps in run_DeepPolisher.sh to autopopulate the input jsons, then submit the pipeline on phoenix

Here is run_DeepPolisher.sh:
```sh
###############################################################################
##                             create input jsons                            ##
###############################################################################
## workflow name = DeepPolisher

## I do this on personal computer but you could do it on phoenix as well, just push to github

# Remove top up data from data table

mkdir -p ~/Desktop/github_repos/phoenix_batch_submissions/workflows/DeepPolisher/DeepPolisher_input_jsons
cd ~/Desktop/github_repos/phoenix_batch_submissions/workflows/DeepPolisher/DeepPolisher_input_jsons

python3 /Users/miramastoras/Desktop/Paten_lab/hprc_intermediate_assembly/hpc/launch_from_table.py \
     --data_table ../DeepPolisher.csv \
     --field_mapping ../DeepPolisher_input_mapping.csv \
     --workflow_name DeepPolisher

## Add/commit/push to github (hprc_intermediate_assembly)

###############################################################################
##                             create launch workflow                      ##
###############################################################################

## on HPC...

## check that github repo is up to date
git -C  /private/groups/patenlab/mira/phoenix_batch_submissions pull

# move to working dir
mkdir -p /private/groups/patenlab/mira/phoenix_batch_executions/workflows/DeepPolisher
cd /private/groups/patenlab/mira/phoenix_batch_executions/workflows/DeepPolisher

## get files
cp -r /private/groups/patenlab/mira/phoenix_batch_submissions/workflows/DeepPolisher/* ./

mkdir -p slurm_logs
export PYTHONPATH="/private/home/juklucas/miniconda3/envs/toil/bin/python"

# submit job
sbatch \
     --job-name=DeepPolisher \
     --array=[11]%1 \ # Here is where you give it 1-the row number in your CSV you want to run.
     --partition=long \
     --cpus-per-task=32 \
     --exclude=phoenix-[09,10,22,23,24,18] \
     --mem=400gb \
     --mail-type=FAIL,END \
     --mail-user=mmastora@ucsc.edu \
     /private/groups/hprc/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
     --wdl ~/progs/hpp_production_workflows/QC/wdl/tasks/DeepPolisher.wdl \
     --sample_csv DeepPolisher.csv \
     --input_json_path '../DeepPolisher_input_jsons/${SAMPLE_ID}_DeepPolisher.json'

###############################################################################
##                             write output files to csv                     ##
###############################################################################

# on hpc after entire batch has finished
cd /private/groups/patenlab/mira/phoenix_batch_executions/workflows/DeepPolisher

python3 /private/groups/hprc/polishing/hprc_intermediate_assembly/hpc/update_table_with_outputs.py \
      --input_data_table ./DeepPolisher.csv \
      --output_data_table ./DeepPolisher.results.csv \
      --json_location '{sample_id}_DeepPolisher_outputs.json'
```

Nextime you run DeepPolisher, you just need to make a new column in the csv, update it on github and submit on phoenix.


### 8. Repeat the procedure for the next steps:

Apply polishing edits to assembly
https://github.com/miramastoras/hpp_production_workflows/blob/master/QC/wdl/tasks/applyPolish.wdl

Run QV QC:
https://github.com/miramastoras/hpp_production_workflows/blob/master/QC/wdl/workflows/hprc_polishing_QC_no_meryl.wdl
