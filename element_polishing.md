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

# align element to polished assembly - mat
docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/long_read_aligner:v0.3.3 minimap2 --cs --eqx -t 32 -ax sr ${MAT_ASM} /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/HG002-20220606.grch38.bam.element.fq.gz | samtools view -b -h > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_mat/HG002.50x_element_google_mat.hprc_polished.bam

# align element to polished assembly - pat
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


## Polishing with heterozygous variants - adapting PHARAOH to phase element reads using hifi information

Align element reads to the diploid polished assembly
```
#!/bin/bash
#SBATCH --job-name=element_all_to_dip
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00

ASM=/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_dip.polished.fasta

# align element to polished assembly
docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/long_read_aligner:v0.3.3 minimap2 --cs --eqx -t 32 -ax sr ${ASM} /private/groups/patenlab/mira/hprc_polishing/data/element_HG002/HG002-20220606.grch38.bam.element.fq.gz | samtools view -b -h > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_dip/HG002.50x_element_google_dip.hprc_polished.bam

ulimit -Sn 5000
mkdir -p /data/tmp/mira/

# sort files and index
samtools sort -T /data/tmp/mira/ -@32 /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_dip/HG002.50x_element_google_dip.hprc_polished.bam > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_dip/HG002.50x_element_google_dip.hprc_polished.srt.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_dip/HG002.50x_element_google_dip.hprc_polished.srt.bam
```

Align HiFi reads all to mat and all to pat
```
{
  "longReadAlignmentScattered.preset": "map-hifi",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.fasta",
  "longReadAlignmentScattered.enableRunningSecphase": false,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -L -Y -I8g",
  "longReadAlignmentScattered.kmerSize": 19,
  "longReadAlignmentScattered.sampleName": "HG002",
  "longReadAlignmentScattered.suffix": "DCv1.2_40x.all_to_pat.minimap2v2.26",
  "longReadAlignmentScattered.readFiles": ["/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190714_120746.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190830_220126.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190901_095311.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64012_190920_173625.dc.q20.fastq.gz"],
  "longReadAlignmentScattered.enableAddingMDTag": false
}
```

Align HiFi reads all to mat and all to pat
```
{
  "longReadAlignmentScattered.preset": "map-hifi",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.polished.fasta",
  "longReadAlignmentScattered.enableRunningSecphase": false,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.alignerOptions": "--cs --eqx -L -Y -I8g",
  "longReadAlignmentScattered.kmerSize": 19,
  "longReadAlignmentScattered.sampleName": "HG002",
  "longReadAlignmentScattered.suffix": "DCv1.2_40x.all_to_mat.minimap2v2.26",
  "longReadAlignmentScattered.readFiles": ["/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190714_120746.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190830_220126.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64011_190901_095311.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/downsampled_40x/m64012_190920_173625.dc.q20.fastq.gz"],
  "longReadAlignmentScattered.enableAddingMDTag": false
}
```

```
#!/bin/bash
#SBATCH --job-name=align_HiFi_all2mat_polished
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --time=7-00:00

cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/HiFi_DCv1.2_all_to_pat
cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/HiFi_DCv1.2_all_to_mat

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

Run PHARAOH using 1000bp as the "homozygous" region length

```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/groups/hprc/polishing/hpp_production_workflows/QC/wdl/workflows/PHARAOH.wdl
```

```
{
  "PHARAOH.allHifiToDiploidBai": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_dip/HG002.50x_element_google_dip.hprc_polished.srt.bam.bai",
  "PHARAOH.minWindowSizeBp": "1000",
  "PHARAOH.extendBp": "50000",
  "PHARAOH.allONTToHap2Bam": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/HiFi_DCv1.2_all_to_mat/long_read_aligner_scattered_outfiles/HG002.DCv1.2_40x.all_to_mat.minimap2v2.26.bam",
  "PHARAOH.allONTToHap1Bam": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/HiFi_DCv1.2_all_to_pat/long_read_aligner_scattered_outfiles/HG002.DCv1.2_40x.all_to_pat.minimap2v2.26.bam",
  "PHARAOH.diploidFaGz": "/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_dip.polished.fasta.gz",
  "PHARAOH.Hap1Fasta": "/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.fasta",
  "PHARAOH.allONTToHap2Bai": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/HiFi_DCv1.2_all_to_mat/long_read_aligner_scattered_outfiles/HG002.DCv1.2_40x.all_to_mat.minimap2v2.26.bam.bai",
  "PHARAOH.Hap2Fasta": "/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.polished.fasta",
  "PHARAOH.Hap1FastaIndex": "/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.fasta.fai",
  "PHARAOH.allHifiToDiploidBam": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/50x_element_google_dip/HG002.50x_element_google_dip.hprc_polished.srt.bam",
  "PHARAOH.allONTToHap1Bai": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/HiFi_DCv1.2_all_to_pat/long_read_aligner_scattered_outfiles/HG002.DCv1.2_40x.all_to_pat.minimap2v2.26.bam.bai",
  "PHARAOH.Hap2FastaIndex": "/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.polished.fasta.fai",
  "PHARAOH.sampleName": "HG002",
  "PHARAOH.PharaohHiFiPreset": "sr",
  "PHARAOH.hifiAlignmentOptions": "--cs --eqx",
  "PHARAOH.hifiAlignmentOptions": "--cs --eqx"
}
```

```
#!/bin/bash
#SBATCH --job-name=PHARAOH_element_test1
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00


export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=7-00:00 --partition=long"
export TOIL_COORDINATION_DIR=/data/tmp


LOCAL_FOLDER=/data/tmp/$(whoami)/PHARAOH_element
mkdir -p ${LOCAL_FOLDER}

mkdir -p /private/groups/patenlab/mira/hprc_polishing/element_polishing/PHARAOH/default_1kb_extend_20kb/toil_logs

toil-wdl-runner \
    --jobStore ${LOCAL_FOLDER}/jobstore \
    --stats \
    --clean=never \
    --batchSystem single_machine \
    --maxCores 32 \
    --batchLogsDir /private/groups/patenlab/mira/hprc_polishing/element_polishing/PHARAOH/default_1kb_extend_20kb/toil_logs \
    /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/workflows/PHARAOH.wdl \
    /private/groups/patenlab/mira/hprc_polishing/element_polishing/PHARAOH/default_1kb_extend_20kb/PHARAOH.inputs.json  \
    --outputDirectory /private/groups/patenlab/mira/hprc_polishing/element_polishing/PHARAOH/default_1kb_extend_20kb/PHARAOH_outputs \
    --outputFile /private/groups/patenlab/mira/hprc_polishing/element_polishing/PHARAOH/default_1kb_extend_20kb/PHARAOH_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    2>&1 | tee log.txt

toil clean ${LOCAL_FOLDER}/jobstore
```

Run deepvariant
```

```
