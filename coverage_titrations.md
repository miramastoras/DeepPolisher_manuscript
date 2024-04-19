## Experimenting with different read coveraged for DeepPolisher pipeline

### HG002 coverage titrations

Downsample PHARAOH bam files
```
/private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG002_y2_DCv1.2_mm2_PHARAOH_whatshap/toil_hprc_pipeline_out/HG002_y2_DCv1.2_PHARAOH_minimap_alignments/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.minimap2v2.26.PHARAOH.bam

/private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_PHv6_DPmm2model1/toil_hprc_deepPolisher_out/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.PHARAOHv6.bam
```

```
#!/bin/bash
#SBATCH --job-name=downsample_HG002
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

samtools view -s 0.5 -b -h -@ 32 /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG002_y2_DCv1.2_mm2_PHARAOH_whatshap/toil_hprc_pipeline_out/HG002_y2_DCv1.2_PHARAOH_minimap_alignments/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.minimap2v2.26.PHARAOH.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_PHARAOHv6_downsampled/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.minimap2v2.26.PHARAOH.20x.downsampled.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/HiFi_DCv1.2_PHARAOHv6_downsampled/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.minimap2v2.26.PHARAOH.20x.downsampled.bam

```
