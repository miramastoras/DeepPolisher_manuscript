### Evaluation of DeepPolisher performance in the centromere

1. Dipcall HG002 Q100 v 1.1 to HG002 raw assembly

`dipcall_inputs.json`:

```json
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/data/hg002v1.1.mat1.fasta.gz",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/data/hg002v1.1.mat.fasta.gz",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}

{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/data/hg002v1.1.pat1.fasta.gz",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/data/hg002v1.1.pat.fasta.gz",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}

{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/hg002v1.1.mat.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/hg002v1.1.mat.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat2.fa",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}

{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/hg002v1.1.pat.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/hg002v1.1.pat.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat2.fa",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}
```
Dipcall polished assembly to Q100
```
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/hg002v1.1.mat.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/hg002v1.1.mat.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.polished.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.copy.polished.fasta",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}

{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/hg002v1.1.pat.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/hg002v1.1.pat.fasta",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.fasta",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.copy.polished.fasta",
  "runDipcall.dipcall.referenceIsHS38": false,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
#!/bin/bash
#SBATCH --job-name=dipcall_truthset
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --exclude=phoenix-[09,10,22,23,24,18]

cd /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_pat
cd /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_mat
cd /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_pat
cd /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_mat
cd /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_pol_to_Q100_pat
cd /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_pol_to_Q100_mat

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24,18]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir ./toil_logs \
    /private/groups/hprc/polishing/hpp_production_workflows/QC/wdl/tasks/dipcall.wdl \
    dipcall_inputs.json \
    --outputDirectory ./dipcall_outfiles \
    --outputFile dipcall_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```



2. Project Q100 centromere annotations to HG002 raw assembly

```
docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode "asm2ref" \
    --paf /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_pat/dipcall_outfiles/hg002v1.1.pat.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/project_censat/Q100_censat.projectable_to_HG002_raw.pat.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/project_censat/Q100_censat.projection_to_HG002_raw.pat.bed \
    --threads 16

# count projectable bases
grep PATERNAL /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed | awk '{sum += $3-$2}END{print sum}'
# 302183823
awk '{sum += $3-$2}END{print sum}'  /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/project_censat/Q100_censat.projection_to_HG002_raw.pat.bed
# 297847515

docker run --rm -it -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode "asm2ref" \
    --paf /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_mat/dipcall_outfiles/hg002v1.1.mat.AP.paf \
    --blocks /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/project_censat/Q100_censat.projectable_to_HG002_raw.mat.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/project_censat/Q100_censat.projection_to_HG002_raw.mat.bed \
    --threads 16

grep MATERNAL /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed | awk '{sum += $3-$2}END{print sum}'
# 273720002
awk '{sum += $3-$2}END{print sum}'  /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/project_censat/Q100_censat.projection_to_HG002_raw.mat.bed
# 251352603
```

define regions
```
cat /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/project_censat/Q100_censat.projection_to_HG002_raw.pat.bed /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/project_censat/Q100_censat.projection_to_HG002_raw.mat.bed > /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/project_censat/Q100_censat.projection_to_HG002_raw.dip.bed

```
Prepare truth vcf
```
cd /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_mat/dipcall_outfiles/hg002v1.1.dipcall

bcftools view -s hg002v1.1.dipcall/hg002v1.1.hap1.bam /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_mat/dipcall_outfiles/hg002v1.1.dipcall/hg002v1.1.pair.vcf.gz | sed 's/hg002v1.1.dipcall\/hg002v1.1.hap1.bam/HG002v1.1/g' > /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_mat/dipcall_outfiles/hg002v1.1.dipcall.mat.vcf

cd /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_pat/dipcall_outfiles/hg002v1.1.dipcall

bcftools view -s hg002v1.1.dipcall/hg002v1.1.hap1.bam /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_pat/dipcall_outfiles/hg002v1.1.dipcall/hg002v1.1.pair.vcf.gz | sed 's/hg002v1.1.dipcall\/hg002v1.1.hap1.bam/HG002v1.1/g' > /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_mat/dipcall_outfiles/hg002v1.1.dipcall.pat.vcf

# combine pat and mat
bcftools concat -Oz /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_mat/dipcall_outfiles/hg002v1.1.dipcall.pat.vcf  /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_Q100_to_HG002_raw_mat/dipcall_outfiles/hg002v1.1.dipcall.mat.vcf > /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.dipcall.HG002_raw.vcf.gz

```

Run hap.py

```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
    /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.dipcall.HG002_raw.vcf.gz /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GQ20_INS1_GQ12_DEL1_GQ5_else/HG002.mm2_model1.polisher_output.GQ20_INS1_GQ12_DEL1_GQ5_else.vcf.gz \
    -r /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    -f /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/project_censat/Q100_censat.projection_to_HG002_raw.dip.bed \
    -o /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/happy_HG002_DP_vs_Q100_all_censat/happy \
    --pass-only --no-roc --no-json --engine=vcfeval --threads=16
```

This didn't work at all, seems that the dipcall is messed up .

### Dipcall variants in centromeres against Q100

#### Raw, mat

subset bam file to regions of mapq 60
```
samtools view -bh -q 60 /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_mat/dipcall_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.bam > /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_mat/dipcall_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.Q60.bam

docker run --rm -it  -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest bedtools bamtobed \
    -i /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_mat/dipcall_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.Q60.bam \
    | bedtools sort -i - | bedtools merge -i - \
    > /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_mat/dipcall_outfiles/dipcall_HG002_raw_to_Q100_mat.Q60.bed
```
Intersect with censat bed
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_mat/dipcall_outfiles/dipcall_HG002_raw_to_Q100_mat.Q60.bed -b hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed | awk '{sum += $3-$2}END{print sum}'
243921318

grep MATERNAL hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed | awk '{sum += $3-$2}END{print sum}'
273720002
```
Intersect with dipcall vcf
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_mat/dipcall_outfiles/dipcall_HG002_raw_to_Q100_mat.Q60.bed -b /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed | bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_mat/dipcall_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pair.vcf.gz -b - | sort | uniq | wc -l

17901

# Count number of SNPs and indels in dipcall vcf

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_mat/dipcall_outfiles/dipcall_HG002_raw_to_Q100_mat.Q60.bed -b /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed | bedtools intersect -header -a /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_mat/dipcall_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pair.vcf.gz -b - | bcftools stats | head -n 30

#
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	2
SN	0	number of records:	17906
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	13219
SN	0	number of MNPs:	0
SN	0	number of indels:	4789
SN	0	number of others:	0
SN	0	number of multiallelic sites:	1077
```

#### Repeat for raw, pat

subset bam file to regions of mapq 60
```
samtools view -bh -q 60 /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_pat/dipcall_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x2.dipcall/HG002.trio_hifiasm_0.19.5.DC_1.2_40x2.hap1.bam > /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_pat/dipcall_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x2.pat.MAPQ60.bam

docker run --rm -it  -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    mobinasri/flagger:latest bedtools bamtobed \
    -i /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_pat/dipcall_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x2.pat.MAPQ60.bam \
    | bedtools sort -i - | bedtools merge -i - \
    > /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_pat/dipcall_outfiles/dipcall_HG002_raw_to_Q100_pat.Q60.bed
```
Intersect with censat bed
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_pat/dipcall_outfiles/dipcall_HG002_raw_to_Q100_pat.Q60.bed -b /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed | awk '{sum += $3-$2}END{print sum}'
270440101

grep PATERNAL /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed| awk '{sum += $3-$2}END{print sum}'
302183823
```
Intersect with dipcall vcf
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_pat/dipcall_outfiles/dipcall_HG002_raw_to_Q100_pat.Q60.bed -b /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed | bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_pat/dipcall_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x2.dipcall/HG002.trio_hifiasm_0.19.5.DC_1.2_40x2.pair.vcf.gz -b - | sort | uniq | wc -l

11500

# Count number of SNPs and indels in dipcall vcf

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_pat/dipcall_outfiles/dipcall_HG002_raw_to_Q100_pat.Q60.bed -b /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/hg002v1.1.cenSatv2.0.noheader.lifted_from_v1.0.1.bed | bedtools intersect -header -a /private/groups/patenlab/mira/hprc_polishing/centromere_evaluation/dipcall_HG002_raw_to_Q100_pat/dipcall_outfiles/HG002.trio_hifiasm_0.19.5.DC_1.2_40x2.dipcall/HG002.trio_hifiasm_0.19.5.DC_1.2_40x2.pair.vcf.gz -b - | bcftools stats | head -n 30

#
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	2
SN	0	number of records:	11501
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	7854
SN	0	number of MNPs:	0
SN	0	number of indels:	3727
SN	0	number of others:	0
SN	0	number of multiallelic sites:	870
```
