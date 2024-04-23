# Experimenting with different read coverages for DeepPolisher

## Titrating HiFi coverage for the PHARAOH + DeepPolisher pipeline

### HG002 coverage titrations

Downsample HiFi read alignment file using samtools - using reads already aligned to assembly with winnowmap from previous experiment.
```
#!/bin/bash
#SBATCH --job-name=downsampling_HG002_DCv1.2
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

# 60x
samtools view -s 0.75 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/downsampled_bams/HG002.DCv1.2.60x.winnowmapv2.03.bam

# 50x
samtools view -s 0.625 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/downsampled_bams/HG002.DCv1.2.50x.winnowmapv2.03.bam

# 30 x - using minimap2 bam so it can be used for DP without pharaoh
samtools view -s 0.75 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/toil_hifi_dv1.2_mm2_all2dip_out/HG002.DCv1.2.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/downsampled_bams/HG002.DCv1.2.30x.minimap2v2.26.bam

# 20x
samtools view -s 0.5 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/toil_hifi_dv1.2_mm2_all2dip_out/HG002.DCv1.2.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/downsampled_bams/HG002.DCv1.2.20x.minimap2v2.26.bam

# 10x
samtools view -s 0.25 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/toil_hifi_dv1.2_mm2_all2dip_out/HG002.DCv1.2.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_minimap2_all2diploid/downsampled_bams/HG002.DCv1.2.10x.minimap2v2.26.bam
```

Run HPRC_deepPolisher workflow on different coverage titrations

### HG005 coverage titrations

Downsample HiFi read alignment file using samtools - using reads already aligned to assembly with winnowmap from previous experiment.
```
#!/bin/bash
#SBATCH --job-name=downsampling_HG005_DCv1.2
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

# 50x
samtools view -s 0.833 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_winnowmap/downsample_60x/HG005.DCv1.2.60x.winnowmapv2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_winnowmap/downsampled_bams/HG005.DCv1.2.50x.winnowmapv2.03.bam

# 30 x
samtools view -s 0.75 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/long_read_aligner_scattered_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/downsampled_bams/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.30x.bam

# 20x
samtools view -s 0.5 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/long_read_aligner_scattered_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/downsampled_bams/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.20x.bam

# 10x
samtools view -s 0.25 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/long_read_aligner_scattered_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_minimap2/downsampled_bams/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.minimap2v2.26.10x.bam

```

## Titrating HiFi coverage for DeepPolisher without PHARAOH
