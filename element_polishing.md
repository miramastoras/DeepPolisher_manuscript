## Polishing HG002 assembly with a second round of element-deepvariant calls

The purpose of this analysis is to test whether adding an additional round of element-deepvariant polishing can clean up calls missed by DeepPolisher-hifi.

### 0. Download element data from the Q100 project

Element data obtained from Q100 project: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/element_trio_all_to_one/

```
#!/bin/bash
cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/element_trio_all_to_one/hg002v1.0.mat.Y.EBV_hg002_element_ins1000.dedup.bam

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/element_trio_all_to_one/hg002v1.0.mat.Y.EBV_hg002_element_ins1000.dedup.bam.bai

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/element_trio_all_to_one/hg002v1.0.mat.Y.EBV_hg002_element_ins500_600.dedup.bam

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/element_trio_all_to_one/hg002v1.0.mat.Y.EBV_hg002_element_ins500_600.dedup.bam.bai
```

### 1. Generating alignments of all element reads to each haplotype

Using the HG002 hprc hifiasm v0.19.5 assembly polished by DeepPolisher minimap2 model1 with the HPRC GQ filters

using 1000bp insert sizes reads

```
#!/bin/bash
#SBATCH --job-name=element_all_to_one
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=20:00:00

MAT_ASM=/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.polished.fasta

PAT_ASM=/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.fasta

OUTDIR=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ

#samtools fastq -@32 /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/hg002v1.0.mat.Y.EBV_hg002_element_ins500_600.dedup.bam > ${OUTDIR}/hg002.element_ins500_600.fastq

#samtools fastq -@32 ${OUTDIR}/hg002v1.0.mat.Y.EBV_hg002_element_ins1000.dedup.bam > ${OUTDIR}/hg002.element_ins1000.fastq

# align ins 1000 to polished assembly - mat
docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/long_read_aligner:v0.3.3 minimap2 --cs --eqx -t 32 -ax sr ${MAT_ASM} ${OUTDIR}/hg002.element_ins1000.fastq | samtools view -b -h > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_mat/HG002.ins1000.mm2.all_to_y2_mat.polished.bam

# align ins 1000 to polished assembly - pat
docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/long_read_aligner:v0.3.3 minimap2 --cs --eqx -t 32 -ax sr ${PAT_ASM} ${OUTDIR}/hg002.element_ins1000.fastq | samtools view -b -h > ${OUTDIR}/element_ins1000_pat/HG002.ins1000.mm2.all_to_y2_pat.polished.bam

# sort files and index
samtools sort -@32 /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_pat/HG002.ins1000.mm2.all_to_y2_pat.polished.bam > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/HG002.ins1000.mm2.all_to_y2_pat.polished.srt.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/HG002.ins1000.mm2.all_to_y2_pat.polished.srt.bam

samtools sort -@32 /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_mat/HG002.ins1000.mm2.all_to_y2_mat.polished.bam > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_mat/HG002.ins1000.mm2.all_to_y2_mat.polished.srt.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_mat/HG002.ins1000.mm2.all_to_y2_mat.polished.srt.bam
```


Remove reads with divergence > 0.4
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl
```

```
{
  "runCorrectBam.correctBam.dockerImage": "mobinasri/secphase:dev-v0.2.0-hom",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_mat/HG002.ins1000.mm2.all_to_y2_mat.polished.srt.bam",
  "runCorrectBam.correctBam.options": "--maxDiv 0.04 --minReadLen 0 --minAlignmentLen 0",
  "runCorrectBam.correctBam.suffix": "maxDiv.04"
}
```
```
{
  "runCorrectBam.correctBam.dockerImage": "mobinasri/secphase:dev-v0.2.0-hom",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_pat/HG002.ins1000.mm2.all_to_y2_pat.polished.srt.bam",
  "runCorrectBam.correctBam.options": "--maxDiv 0.04 --minReadLen 0 --minAlignmentLen 0",
  "runCorrectBam.correctBam.suffix": "maxDiv.04"
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_pat
/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_mat


export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=24:00:00 --partition=high_priority"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats
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

Check coverage
```
#!/bin/bash
#SBATCH --job-name=HG002_mosdepth
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_pat/correct_bam_outfiles/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -n --fast-mode -t 4 /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_pat/correct_bam_outfiles/mosdepth/HG002 \
    /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_pat/correct_bam_outfiles/HG002.ins1000.mm2.all_to_y2_pat.polished.srt.maxDiv.04.bam

python3 ~/progs/mosdepth/scripts/plot-dist.py /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_pat/correct_bam_outfiles/mosdepth/HG002.mosdepth.global.dist.txt
```


### 2. Run DeepVariant on element alignments

```
#!/bin/bash
#SBATCH --job-name=element_dv_all2pat
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
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
  --reads=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_pat/correct_bam_outfiles/HG002.ins1000.mm2.all_to_y2_pat.polished.srt.maxDiv.04.bam \
  --output_vcf=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.pat.vcf \
  --output_gvcf=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.pat.gvcf \
  --num_shards=32
```

```
#!/bin/bash
#SBATCH --job-name=element_dv_all2mat
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
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
  --reads=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_mat/correct_bam_outfiles/HG002.ins1000.mm2.all_to_y2_mat.polished.srt.maxDiv.04.bam \
  --output_vcf=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.mat.vcf \
  --output_gvcf=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.mat.gvcf \
  --num_shards=32
```

Filter variants by GQ7 and extract just homalt variants

```
bcftools view -Oz -f "PASS" -e 'FORMAT/GQ<=7' /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.mat.vcf > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.mat.GQ7.vcf.gz

bcftools view -Oz -f "PASS" -e 'FORMAT/GQ<=7' /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.pat.vcf > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.pat.GQ7.vcf.gz

# Maternal
zcat /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.mat.GQ7.vcf.gz | grep "^#" > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.mat.GQ7.homalt.vcf

zcat /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.mat.GQ7.vcf.gz | grep -v "^#" | grep "1/1" >> /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.mat.GQ7.homalt.vcf

bcftools stats /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.mat.GQ7.homalt.vcf | head -n 30
# SN	0	number of records:	2194

# Paternal
zcat /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.pat.GQ7.vcf.gz | grep "^#" > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.pat.GQ7.homalt.vcf

zcat /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.pat.GQ7.vcf.gz | grep -v "^#" | grep "1/1" >> /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.pat.GQ7.homalt.vcf

bcftools stats /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/deepvariant/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/deepvariant_1.6.WGS.pat.GQ7.homalt.vcf | head -n 30
# 2987 records
```

### 3. Polish with element homalt GQ7 variants, run evaluation

located under phoenix batch submissions github: https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/applyPolish_dipcall_happy/GIAB_samples_manuscript

### 4. Rerun GIAB variant evaluation on just places where GIAB v4.2.1 and Q100 agree.

```
# intersect concordant bed with dipcall bed
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/HG002_GRCh38_T2T_concordant.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent_intersect_HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.dipcall.bed > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/hom_polish_T2T_GIAB_concordant/HG002_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/HG002_GRCh38_T2T_concordant_GQ20_INS1_GQ12_DEL1_GQ5_else.dipcall.bed

# DeepPolisher
bash /private/home/mmastora/progs/scripts/GIAB_happy.sh /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.dipcall.vcf.gz /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/hom_polish_T2T_GIAB_concordant/HG002_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/HG002_GRCh38_T2T_concordant_GQ20_INS1_GQ12_DEL1_GQ5_else.dipcall.bed /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/hom_polish_T2T_GIAB_concordant/HG002_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/happy_out HG002

# DP + element homalt

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/HG002_GRCh38_T2T_concordant.bed -b /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/applyPolish_dipcall_outputs/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent_intersect_HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_hap1.polished.dipcall.bed > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/hom_polish_T2T_GIAB_concordant/HG002_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ_element_homaltGQ7/HG002_T2T_GIAB_concordant_HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7.dipcall.bed

bash /private/home/mmastora/progs/scripts/GIAB_happy.sh /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/applyPolish_dipcall_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_hap1.polished.dipcall.vcf.gz /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/hom_polish_T2T_GIAB_concordant/HG002_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ_element_homaltGQ7/HG002_T2T_GIAB_concordant_HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7.dipcall.bed /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/hom_polish_T2T_GIAB_concordant/HG002_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ_element_homaltGQ7/happy_out HG002

# 2529996443 / 3099922541 = 81.614%
# 2530008724 / 3099922541 = 81.615%
```

### 5. Annotate missing variants by the following categories:

Extract FP/FN variants
```
# before element
zcat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.vcf.gz | grep "^#" >  /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.FPFN.vcf

zcat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.vcf.gz | grep -v "^#" | grep :F >> /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.FPFN.vcf

# 11404 false records

# after element
zcat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.vcf.gz | grep "^#" > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf

zcat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.vcf.gz | grep -v "^#" | grep :F >> /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf

bcftools stats /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf

# 10330 false records
```

#### Records concordant with Q100 regions:

Before element polishing
```
bedtools intersect -header -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.FPFN.vcf -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/HG002_GRCh38_T2T_concordant.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.FPFN.T2T.concordant.vcf
```

After element polishing
```
bedtools intersect -header -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf -b /private/groups/patenlab/mira/hprc_polishing/GIAB_T2T_platinum_truthset/HG002_GRCh38_T2T_concordant.bed > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.T2T.concordant.vcf

bcftools stats
# 3878 records fall in concordant regions
```

#### What is the genotype makeup of the missing variants?

https://github.com/miramastoras/element_polishing/blob/main/scripts/count_missing_vars_happy.sh

All FP/FN before element polishing
```
bcftools view -Oz /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.FPFN.vcf > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.FPFN.vcf.gz

bash ~/progs/element_polishing/scripts/count_missing_vars_happy.sh /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.FPFN.vcf.gz /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/count_missing_vars_happy.txt

TruthHet_Query_WrongHet,TruthHet_QueryHom,TruthHom_QueryHet,TruthHom_QueryWrongHom
461,3596,6724,623
```

All FP/FN after element polishing
```
bcftools view -Oz /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf.gz

bash ~/progs/element_polishing/scripts/count_missing_vars_happy.sh /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf.gz /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/count_missing_vars_happy.txt

# TruthHet_Query_WrongHet,TruthHet_QueryHom,TruthHom_QueryHet,TruthHom_QueryWrongHom
# 459,3700,5612,559
```

Just GIAB concordant FP/FN before element polishing
```
bcftools view -Oz /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.FPFN.T2T.concordant.vcf > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.FPFN.T2T.concordant.vcf.gz

bash ~/progs/element_polishing/scripts/count_missing_vars_happy.sh /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_happy_out.FPFN.T2T.concordant.vcf.gz /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/happy_outputs/count_missing_vars_happy.FPFN.GIAB_T2T_concordant.txt

TruthHet_Query_WrongHet,TruthHet_QueryHom,TruthHom_QueryHet,TruthHom_QueryWrongHom
424,2527,1822,177
```

Just GIAB concordant FP/FN after element polishing
```
bcftools view -Oz /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.T2T.concordant.vcf > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.T2T.concordant.vcf.gz

bash ~/progs/element_polishing/scripts/count_missing_vars_happy.sh /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.T2T.concordant.vcf.gz /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/count_missing_vars_happy.FPFN.GIAB_T2T_concordant.txt

TruthHet_Query_WrongHet,TruthHet_QueryHom,TruthHom_QueryHet,TruthHom_QueryWrongHom
424,2630,712,112
```

#### Check how many missing variants are in callable regions

Pat
```
#!/bin/bash
#SBATCH --job-name=ele_pat_mosdepth_callable
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/shnegi/.conda/envs/bedtools

/private/groups/migalab/shnegi/apps/mosdepth --threads 4 \
    -f /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.fasta \
    --use-median --fast-mode --mapq 10 \
    --quantize 0:1:10:60: \
    HG002.element_pat.quantized \
    /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_pat/correct_bam_outfiles/HG002.ins1000.mm2.all_to_y2_pat.polished.srt.maxDiv.04.bam

zcat /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/HG002.element_pat.quantized.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/HG002.element_pat.quantized.tsv

awk 'BEGIN {FS=OFS="\t"} {
        if ($4 == "0:1") {
            $4 = "NO_COVERAGE"
        } else if ($4 == "1:10") {
            $4 = "LOW_COVERAGE"
        } else if ($4 == "10:60") {
            $4 = "CALLABLE"
        } else if ($4 == "60:inf") {
            $4 = "HIGH_COVERAGE"
        }
        print
}' /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/HG002.element_pat.quantized.tsv > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/HG002.element_pat.quantized.bed
```

Mat
```
#!/bin/bash
#SBATCH --job-name=ele_mat_mosdepth_callable
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/shnegi/.conda/envs/bedtools

/private/groups/migalab/shnegi/apps/mosdepth --threads 4 \
    -f /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.polished.fasta \
    --use-median --fast-mode --mapq 10 \
    --quantize 0:1:10:60: \
    HG002.element_mat.quantized \
    /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ/element_ins1000_mat/correct_bam_outfiles/HG002.ins1000.mm2.all_to_y2_mat.polished.srt.maxDiv.04.bam

zcat /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/HG002.element_mat.quantized.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/HG002.element_mat.quantized.tsv

awk 'BEGIN {FS=OFS="\t"} {
        if ($4 == "0:1") {
            $4 = "NO_COVERAGE"
        } else if ($4 == "1:10") {
            $4 = "LOW_COVERAGE"
        } else if ($4 == "10:60") {
            $4 = "CALLABLE"
        } else if ($4 == "60:inf") {
            $4 = "HIGH_COVERAGE"
        }
        print
}' /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/HG002.element_mat.quantized.tsv > /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/HG002.element_mat.quantized.bed
```
Separate bed file by category
```
grep "NO_COVERAGE" HG002.element_mat.quantized.bed > HG002.element.dip.NO_COVERAGE.quantized.bed
grep "NO_COVERAGE" HG002.element_pat.quantized.bed >> HG002.element.dip.NO_COVERAGE.quantized.bed

grep "LOW_COVERAGE" HG002.element_mat.quantized.bed > HG002.element.dip.LOW_COVERAGE.quantized.bed
grep "LOW_COVERAGE" HG002.element_pat.quantized.bed >> HG002.element.dip.LOW_COVERAGE.quantized.bed

grep "CALLABLE" HG002.element_mat.quantized.bed > HG002.element.dip.CALLABLE.quantized.bed
grep "CALLABLE" HG002.element_pat.quantized.bed >> HG002.element.dip.CALLABLE.quantized.bed

grep "HIGH_COVERAGE" HG002.element_mat.quantized.bed > HG002.element.dip.HIGH_COVERAGE.quantized.bed
grep "HIGH_COVERAGE" HG002.element_pat.quantized.bed >> HG002.element.dip.HIGH_COVERAGE.quantized.bed
```
Project bedfiles to HG38
```
# switch naming on paf files, as asm files were switched originally
grep "tp:A:P" /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/applyPolish_dipcall_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_hap1.polished.dipcall/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_hap1.polished.hap1.paf > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/applyPolish_dipcall_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7.polished.hap2.AP.paf

grep "tp:A:P" /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/applyPolish_dipcall_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_hap1.polished.dipcall/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_hap1.polished.hap2.paf > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/applyPolish_dipcall_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7.polished.hap1.AP.paf
```
Convert vcf files to bed
```
export PATH=$PATH:/private/home/mmastora/progs/bin/

# all GIAB
/private/home/mmastora/progs/bin/vcf2bed --do-not-split < /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf | awk '{print $1"\t"$2-10"\t"$3+10"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12}' > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf.10bp.bed

bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf.10bp.bed | bedtools merge -i - -c 1 -o count > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf.10bp.srt.mrg.bed

# T2T concordant
/private/home/mmastora/progs/bin/vcf2bed --do-not-split < /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.T2T.concordant.vcf | awk '{print $1"\t"$2-10"\t"$3+10"\t"$6"\t"$7"\t"$10"\t"$11"\t"$12}' > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.T2T.concordant.vcf.10bp.bed

bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.T2T.concordant.vcf.10bp.bed | bedtools merge -i - -c 1 -o count > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.T2T.concordant.vcf.10bp.srt.mrg.bed
```

Project GIAB FP/FN to pat
```
docker run --rm -u `id -u`:`id -g` \
-v /private/groups:/private/groups \
mobinasri/flagger:latest python3 \
/home/programs/src/project_blocks_multi_thread.py \
--threads 16 --mode 'ref2asm' \
--paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/applyPolish_dipcall_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7.polished.hap1.AP.paf \
--blocks /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf.10bp.srt.mrg.bed \
--outputProjectable /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/projections/GIAB.FPFN.projectable.pat.bed \
--outputProjection /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/projections/GIAB.FPFN.projection.pat.bed
```
Project GIAB FP/FN to mat
```
docker run --rm -u `id -u`:`id -g` \
-v /private/groups:/private/groups \
mobinasri/flagger:latest python3 \
/home/programs/src/project_blocks_multi_thread.py \
--threads 16 --mode 'ref2asm' \
--paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/applyPolish_dipcall_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7.polished.hap2.AP.paf \
--blocks /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.vcf.10bp.srt.mrg.bed \
--outputProjectable /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/projections/GIAB.FPFN.projectable.mat.bed \
--outputProjection /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/projections/GIAB.FPFN.projection.mat.bed
```
Project GIAB-T2T concordant FP/FN to pat
```
docker run --rm -u `id -u`:`id -g` \
-v /private/groups:/private/groups \
mobinasri/flagger:latest python3 \
/home/programs/src/project_blocks_multi_thread.py \
--threads 16 --mode 'ref2asm' \
--paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/applyPolish_dipcall_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7.polished.hap1.AP.paf \
--blocks /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.T2T.concordant.vcf.10bp.srt.mrg.bed \
--outputProjectable /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/projections/GIAB.T2T_concordant.FPFN.projectable.pat.bed \
--outputProjection /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/projections/GIAB.T2T_concordant.FPFN.projection.pat.bed
```
Project GIAB-T2T concordant FP/FN to mat
```
docker run --rm -u `id -u`:`id -g` \
-v /private/groups:/private/groups \
mobinasri/flagger:latest python3 \
/home/programs/src/project_blocks_multi_thread.py \
--threads 16 --mode 'ref2asm' \
--paf /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/applyPolish_dipcall_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7.polished.hap2.AP.paf \
--blocks /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7/happy_outputs/HG002_y2_DCv1.2_mm2model1_hprcGQ_elementGQ7_happy_out.FPFN.T2T.concordant.vcf.10bp.srt.mrg.bed \
--outputProjectable /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/projections/GIAB.T2T_concordant.FPFN.projectable.mat.bed \
--outputProjection /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/projections/GIAB.T2T_concordant.FPFN.projection.mat.bed
```

combine haplotypes
```
cat GIAB.FPFN.projection.mat.bed GIAB.FPFN.projection.pat.bed > GIAB.FPFN.projection.dip.bed
cat GIAB.T2T_concordant.FPFN.projection.mat.bed GIAB.T2T_concordant.FPFN.projection.pat.bed > GIAB.T2T_concordant.FPFN.projection.dip.bed
```

Intersect with callable regions: GIAB

total projected variants:
```
awk '{sum+=$4;} END{print sum;}' GIAB.FPFN.projection.dip.bed
# 20031
```
```
cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/projections

bedtools intersect -a GIAB.FPFN.projection.dip.bed -b ../HG002.element.dip.CALLABLE.quantized.bed | sort | uniq | awk '{sum+=$4;} END{print sum;}'
# 15970

bedtools intersect -a GIAB.FPFN.projection.dip.bed -b ../HG002.element.dip.HIGH_COVERAGE.quantized.bed | sort | uniq | awk '{sum+=$4;} END{print sum;}'
# 330

bedtools intersect -a GIAB.FPFN.projection.dip.bed -b ../HG002.element.dip.LOW_COVERAGE.quantized.bed | sort | uniq | awk '{sum+=$4;} END{print sum;}'
# 1820

bedtools intersect -a GIAB.FPFN.projection.dip.bed -b ../HG002.element.dip.NO_COVERAGE.quantized.bed | sort | uniq | awk '{sum+=$4;} END{print sum;}'
# 3493
```

Intersect with callable regions: T2T concordant

total projected variants:
```
awk '{sum+=$4;} END{print sum;}' GIAB.T2T_concordant.FPFN.projection.dip.bed
# 7725
```
```
cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/callable_regions/HG002_y2_DCv1.2_PHv6_DPmm2_model1_dockerv0.8_HPRC_GQ_ele_GQ7_hom/projections

bedtools intersect -a GIAB.T2T_concordant.FPFN.projection.dip.bed -b ../HG002.element.dip.CALLABLE.quantized.bed | sort | uniq | awk '{sum+=$4;} END{print sum;}'
# 7591

bedtools intersect -a GIAB.T2T_concordant.FPFN.projection.dip.bed -b ../HG002.element.dip.HIGH_COVERAGE.quantized.bed | sort | uniq | awk '{sum+=$4;} END{print sum;}'
# 59

bedtools intersect -a GIAB.T2T_concordant.FPFN.projection.dip.bed -b ../HG002.element.dip.LOW_COVERAGE.quantized.bed | sort | uniq | awk '{sum+=$4;} END{print sum;}'
# 156

bedtools intersect -a GIAB.T2T_concordant.FPFN.projection.dip.bed -b ../HG002.element.dip.NO_COVERAGE.quantized.bed | sort | uniq | awk '{sum+=$4;} END{print sum;}'
# 94
```

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
