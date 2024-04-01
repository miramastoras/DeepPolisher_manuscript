## Polishing HG002 with a second round of element-deepvariant calls

Element data obtained from Q100 project: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/element_trio_all_to_one/

```
#!/bin/bash
cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/element_trio_all_to_one/hg002v1.0.mat.Y.EBV_hg002_element_ins1000.dedup.bam

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/element_trio_all_to_one/hg002v1.0.mat.Y.EBV_hg002_element_ins1000.dedup.bam.bai

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/element_trio_all_to_one/hg002v1.0.mat.Y.EBV_hg002_element_ins500_600.dedup.bam

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/element_trio_all_to_one/hg002v1.0.mat.Y.EBV_hg002_element_ins500_600.dedup.bam.bai
```

### 1. Generating alignments of all element reads to each haplotype of the HG002 hprc hifiasm v0.19.5 assembly polished by DeepPolisher minimap2 model1

```
#!/bin/bash
#SBATCH --job-name=element_pat_mm2
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=20:00:00

MAT_ASM=/private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/HG002_y2_DCv1.2_PHv6_DPmm2model1_dockerv0.0.8_12122023/applyPolish_dipcall_happy/applyPolish_dipcall_outfiles/HG002_hap2.polished.fasta

PAT_ASM=/private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/HG002_y2_DCv1.2_PHv6_DPmm2model1_dockerv0.0.8_12122023/applyPolish_dipcall_happy/applyPolish_dipcall_outfiles/HG002_hap1.polished.fasta

OUTDIR=/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1

samtools fastq -@32 /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1/hg002v1.0.mat.Y.EBV_hg002_element_ins500_600.dedup.bam > ${OUTDIR}/hg002.element_ins500_600.fastq

samtools fastq -@32 ${OUTDIR}/hg002v1.0.mat.Y.EBV_hg002_element_ins1000.dedup.bam > ${OUTDIR}/hg002.element_ins1000.fastq

# align ins 500 to polished assembly - mat
docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/long_read_aligner:v0.3.3 minimap2 --cs --eqx -t 32 -ax sr ${MAT_ASM} ${OUTDIR}/hg002.element_ins500_600.fastq | samtools view -b -h > ${OUTDIR}/element_ins500_600_mat/HG002.ins500_600.mm2.all_to_y2_mat.polished.bam

# align ins 500 to polished assembly - pat
docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/long_read_aligner:v0.3.3 minimap2 --cs --eqx -t 32 -ax sr ${PAT_ASM} ${OUTDIR}/hg002.element_ins500_600.fastq | samtools view -b -h > ${OUTDIR}/element_ins500_600_pat/HG002.ins500_600.mm2.all_to_y2_pat.polished.bam

# align ins 1000 to polished assembly - mat
docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/long_read_aligner:v0.3.3 minimap2 --cs --eqx -t 32 -ax sr ${MAT_ASM} ${OUTDIR}/hg002.element_ins1000.fastq | samtools view -b -h > ${OUTDIR}/element_ins1000_mat/HG002.ins1000.mm2.all_to_y2_mat.polished.bam

# align ins 1000 to polished assembly - pat
docker run -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira mobinasri/long_read_aligner:v0.3.3 minimap2 --cs --eqx -t 32 -ax sr ${PAT_ASM} ${OUTDIR}/hg002.element_ins1000.fastq | samtools view -b -h > ${OUTDIR}/element_ins1000_pat/HG002.ins1000.mm2.all_to_y2_pat.polished.bam
```

Remove reads with divergence > 0.4
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/correct_bam.wdl
```

```
{
  "runCorrectBam.correctBam.dockerImage": "mobinasri/secphase:dev-v0.2.0-hom",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1/element_ins1000_mat/HG002.ins1000.mm2.all_to_y2_mat.polished.bam",
  "runCorrectBam.correctBam.options": "--maxDiv 0.04 --minReadLen 0 --minAlignmentLen 0",
  "runCorrectBam.correctBam.suffix": "maxDiv.04"
}
```
```
{
  "runCorrectBam.correctBam.dockerImage": "mobinasri/secphase:dev-v0.2.0-hom",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1/element_ins1000_pat/HG002.ins1000.mm2.all_to_y2_pat.polished.bam",
  "runCorrectBam.correctBam.options": "--maxDiv 0.04 --minReadLen 0 --minAlignmentLen 0",
  "runCorrectBam.correctBam.suffix": "maxDiv.04"
}
```
```
{
  "runCorrectBam.correctBam.dockerImage": "mobinasri/secphase:dev-v0.2.0-hom",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1/element_ins500_600_pat/HG002.ins500_600.mm2.all_to_y2_pat.polished.bam",
  "runCorrectBam.correctBam.options": "--maxDiv 0.04 --minReadLen 0 --minAlignmentLen 0",
  "runCorrectBam.correctBam.suffix": "maxDiv.04"
}
```
```
{
  "runCorrectBam.correctBam.dockerImage": "mobinasri/secphase:dev-v0.2.0-hom",
  "runCorrectBam.correctBam.Bam": "/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1/element_ins500_600_mat/HG002.ins500_600.mm2.all_to_y2_mat.polished.bam",
  "runCorrectBam.correctBam.options": "--maxDiv 0.04 --minReadLen 0 --minAlignmentLen 0",
  "runCorrectBam.correctBam.suffix": "maxDiv.04"
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1/element_ins1000_pat
/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1/element_ins1000_mat
/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1/element_ins500_600_pat
/private/groups/patenlab/mira/hprc_polishing/element_polishing/DeepPolisher_assemblies/alignments/HG002_y2_DCv1.2_PHv6_DPmm2_model1/element_ins500_600_mat


export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=4:00:00 --partition=high_priority"
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
    --restart \
    2>&1 | tee log.txt
```
