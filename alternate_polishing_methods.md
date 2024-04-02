## Benchmarking DeepPolisher against alternate polishing methods

Alternate polishing approaches:

- [T2T polishing pipeline](https://github.com/arangrhie/T2T-Polish/blob/master/doc/T2T_polishing_case_study.md)
- [HiFi-DeepVariant](https://github.com/google/deepvariant)
- [NextPolish2](https://github.com/Nextomics/NextPolish2)

### 1. T2T polishing pipeline

#### 1.1 HG002

**Generate illumina alignments**

BWA index
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa index -p \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa
```

BWA mem alignment
```
#!/bin/bash
#SBATCH --job-name=HG002_ilm
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa mem -t32 \
    /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R1.fastq.gz \
    /private/groups/patenlab/mira/hprc_polishing/data/reads/HG002/illumina/HG002_HiSeq30x_subsampled_R2.fastq.gz \
    | samtools view -b -h \
    > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/illumina_all_to_dip/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.HiSeq30x.bwa.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/illumina_all_to_dip/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.HiSeq30x.bwa.bam
```

#### 1.2 HG005

**Generate illumina alignments**

BWA index
```
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa index -p \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa
```

Convert bam to fastq
```
#!/bin/bash
#SBATCH --job-name=HG005_ilm
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
    quay.io/masri2019/hpp_bwa:latest bwa mem -t32 \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_one/reads/HG005.ilm.fastq.gz \
    | samtools view -b -h \
    > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_dip/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.bwa.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/illumina_all_to_dip/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.ilm.bwa.bam
```
