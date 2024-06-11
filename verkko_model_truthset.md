## Creating a truthset for training a DeepPolisher hifi verkko model

1. Run dipcall

```
#!/bin/bash
#SBATCH --job-name=dipcall_truthset
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
    humanpangenomics/hpp_dipcall_v0.3:latest /opt/dipcall/dipcall.kit/run-dipcall \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_pat \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_pat_verkko_2.0_2024_06_05.fa \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.pat_X_EBV_MT.fasta \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.pat_X_EBV_MT.fasta > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_pat/HG002_t2t_v0.9_2_hprc_pat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
        humanpangenomics/hpp_dipcall_v0.3:latest make -j2 -f /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_pat/HG002_t2t_v0.9_2_hprc_pat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
    humanpangenomics/hpp_dipcall_v0.3:latest /opt/dipcall/dipcall.kit/run-dipcall \
    /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_mat \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_mat_verkko_2.0_2024_06_05.fa \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.mat_Y_EBV_MT.fasta \
    /private/groups/patenlab/mira/hprc_polishing/HG002_v0.9_truth/hg002v0.9.mat_Y_EBV_MT.fasta > /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_mat/HG002_t2t_v0.9_2_hprc_mat.mak

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
        humanpangenomics/hpp_dipcall_v0.3:latest make -j2 -f /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/dipcall_verkkov2.0_mat/HG002_t2t_v0.9_2_hprc_mat.mak
```

Prepare bed file

Dipcall raw assembly against GRCh38
```

{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_mat_verkko_2.0_2024_06_05.fa",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/HG002_pat_verkko_2.0_2024_06_05.fa",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall

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


```
# 2. removed unmapped records in paf files
cd /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall/dipcall_outfiles/HG002_verkko_2.0_2024_06_05.dipcall

grep "tp:A:P" HG002_verkko_2.0_2024_06_05.hap1.paf > HG002_verkko_2.0_2024_06_05.hap1.AP.paf
grep "tp:A:P" HG002_verkko_2.0_2024_06_05.hap2.paf > HG002_verkko_2.0_2024_06_05.hap2.AP.paf

# 3. Project GIAB confidence bedfile to y2 assembly, each haplotype separately
# Hap1
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall/dipcall_outfiles/HG002_verkko_2.0_2024_06_05.dipcall/HG002_verkko_2.0_2024_06_05.hap1.AP.paf \
  --blocks /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap1_GIAB_conf.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap1_GIAB_conf.projection.bed

# Hap 2
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/data/HG002_verkko_2.0/dipcall/dipcall_outfiles/HG002_verkko_2.0_2024_06_05.dipcall/HG002_verkko_2.0_2024_06_05.hap2.AP.paf \
  --blocks /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap2_GIAB_conf.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/verkko_model_truthset/hifi_model/HG002_verkko_2.0_hap2_GIAB_conf.projection.bed
```

Uploaded for kishwar here:
```
gsutil -o GSUtil:parallel_composite_upload_threshold=50GB -m cp -r verkko_model_truthset_HG002_v0.9 gs://pepper-deepvariant/mira/
```
