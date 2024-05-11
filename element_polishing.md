## Polishing HG002 assembly with a second round of element-deepvariant calls

The purpose of this analysis is to test whether adding an additional round of element-deepvariant polishing can clean up calls missed by DeepPolisher-hifi.

### 0. Align element data provided by Google to polished assembly

Using the HG002 hprc hifiasm v0.19.5 assembly polished by DeepPolisher minimap2 model1 with the HPRC GQ filters

using 1000bp insert sizes reads

```
#!/bin/bash
#SBATCH --job-name=element_all_to_one
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00

MAT_ASM=/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.polished.fasta

PAT_ASM=/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.fasta

# align ins 1000 to polished assembly - mat
docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/long_read_aligner:v0.3.3 minimap2 --cs --eqx -t 32 -ax sr ${MAT_ASM} /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/HG002-20220606.grch38.bam.element.fq.gz | samtools view -b -h > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_mat/HG002.50x_element_google_mat.hprc_polished.bam

# align ins 1000 to polished assembly - pat
docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/long_read_aligner:v0.3.3 minimap2 --cs --eqx -t 32 -ax sr ${PAT_ASM} /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/HG002-20220606.grch38.bam.element.fq.gz | samtools view -b -h > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_pat/HG002.50x_element_google_pat.hprc_polished.bam

# sort files and index
samtools sort -@32 /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_pat/HG002.50x_element_google_pat.hprc_polished.bam > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_pat/HG002.50x_element_google_pat.hprc_polished.srt.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_pat/HG002.50x_element_google_pat.hprc_polished.srt.bam

samtools sort -@32 /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_mat/HG002.50x_element_google_mat.hprc_polished.bam > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_mat/HG002.50x_element_google_mat.hprc_polished.srt.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_mat/HG002.50x_element_google_mat.hprc_polished.srt.bam
```
Remove reads with divergence > 0.04
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl
```

```
{
  "runCorrectBam.correctBam.dockerImage": "mobinasri/secphase:dev-v0.2.0-hom",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_mat/HG002.50x_element_google_mat.hprc_polished.srt.bam",
  "runCorrectBam.correctBam.options": "--maxDiv 0.04 --minReadLen 0 --minAlignmentLen 0",
  "runCorrectBam.correctBam.suffix": "maxDiv.04"
}
```
```
{
  "runCorrectBam.correctBam.dockerImage": "mobinasri/secphase:dev-v0.2.0-hom",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_pat/HG002.50x_element_google_pat.hprc_polished.srt.bam",
  "runCorrectBam.correctBam.options": "--maxDiv 0.04 --minReadLen 0 --minAlignmentLen 0",
  "runCorrectBam.correctBam.suffix": "maxDiv.04"
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_pat/correct_bam
/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_mat/correct_bam


export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=24:00:00 --partition=high_priority"
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
### 2. Run DeepVariant on element alignments

```
#!/bin/bash
#SBATCH --job-name=element_dv_all2pat
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=24:00:00

BIN_VERSION="1.6.0"

time docker run \
  -u `id -u`:`id -g` \
  -v "/private/groups/patenlab/mira":"/private/groups/patenlab/mira" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.fasta \
  --reads=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_pat/correct_bam/correct_bam_outfiles/HG002.50x_element_google_pat.hprc_polished.srt.maxDiv.04.bam \
  --output_vcf=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.pat.vcf \
  --output_gvcf=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.pat.gvcf \
  --num_shards=32
```

```
#!/bin/bash
#SBATCH --job-name=element_dv_all2mat
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=24:00:00

BIN_VERSION="1.6.0"

time docker run \
  -u `id -u`:`id -g` \
  -v "/private/groups/patenlab/mira":"/private/groups/patenlab/mira" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.polished.fasta \
  --reads=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_mat/correct_bam/correct_bam_outfiles/HG002.50x_element_google_mat.hprc_polished.srt.maxDiv.04.bam \
  --output_vcf=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.mat.vcf \
  --output_gvcf=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.mat.gvcf \
  --num_shards=32
```

Filter variants by GQ7 and extract just homalt variants

```
bcftools view -Oz -f "PASS" -e 'FORMAT/GQ<=7' /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.mat.vcf > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.mat.GQ7.vcf.gz

bcftools view -Oz -f "PASS" -e 'FORMAT/GQ<=7' /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.pat.vcf > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.pat.GQ7.vcf.gz

# Maternal
zcat /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.mat.GQ7.vcf.gz | grep "^#" > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.mat.GQ7.homalt.vcf

zcat /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.mat.GQ7.vcf.gz | grep -v "^#" | grep "1/1" >> /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.mat.GQ7.homalt.vcf

bcftools stats /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.mat.GQ7.homalt.vcf | head -n 30
# SN	0	number of records:	3276

# Paternal
zcat /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.pat.GQ7.vcf.gz | grep "^#" > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.pat.GQ7.homalt.vcf

zcat /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.pat.GQ7.vcf.gz | grep -v "^#" | grep "1/1" >> /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.pat.GQ7.homalt.vcf

bcftools stats /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/google_50x/deepvariant_1.6.WGS.pat.GQ7.homalt.vcf | head -n 30
# 2310 records
```
### 3. Polish with element homalt GQ7 variants, run evaluation

located under phoenix batch submissions github: https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/applyPolish_dipcall_happy/GIAB_samples_manuscript

### 4. Investigate missing variants 
