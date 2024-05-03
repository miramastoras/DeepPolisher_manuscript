### Adding a PL tag to distinguish between HiFi and ONT reads

```
cd /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/alignments/add_PL_tag
```
Bam files:
```
# HiFi
/private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/PHARAOH/DCv1.2_DoradoR10/minimap2/toil_pharaoh_out/hifi_ONT_combined_model_minimap2/HG002.trio_hifiasm_0.19.5.DCv1.2.minimap2v2.26.PHARAOH.bam

# ONT
/private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/alignments/Dorado_R10_UL/minimap2/R10_all_to_dip/toil_correctbam_out/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.R1041_Dorado_v0.1.1_43x_alignment.minimap2.26.secphase_v0.4.3_marker_mode.bam
```

Create header mapping file. RG ID must be the prefix of each .bam file. PL tag can only be PACBIO or ONT (https://samtools.github.io/hts-specs/SAMv1.pdf)
```
printf '@RG\tID:HG002.trio_hifiasm_0.19.5.DCv1.2.minimap2v2.26.PHARAOH\tPL:PACBIO\n@RG\tID:HG002.trio_hifiasm_0.19.5.DC_1.2_40x.R1041_Dorado_v0.1.1_43x_alignment.minimap2.26.secphase_v0.4.3_marker_mode\tPL:ONT\n' > rg.txt
```

Merge hifi+ont bam
```
#!/bin/bash
#SBATCH --job-name=merge_hifi_ont
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

cd /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/alignments/add_PL_tag

samtools merge -@64 -rh rg.txt /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/alignments/merged_hybrid_bam/HG002.trio_hifiasm_0.19.5.DCv1.2.PHARAOH.Dorado.R10.secphase.mm2v2.26.merged.PL_tag.bam /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/PHARAOH/DCv1.2_DoradoR10/minimap2/toil_pharaoh_out/hifi_ONT_combined_model_minimap2/HG002.trio_hifiasm_0.19.5.DCv1.2.minimap2v2.26.PHARAOH.bam /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/alignments/Dorado_R10_UL/minimap2/R10_all_to_dip/toil_correctbam_out/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.R1041_Dorado_v0.1.1_43x_alignment.minimap2.26.secphase_v0.4.3_marker_mode.bam

```

View new header tags
```
samtools view /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/alignments/merged_hybrid_bam/HG002.trio_hifiasm_0.19.5.DCv1.2.PHARAOH.Dorado.R10.secphase.mm2v2.26.merged.PL_tag.bam -H | grep @RG
```

Each read will have the following RG tag:
```
RG:Z:HG002.trio_hifiasm_0.19.5.DCv1.2.minimap2v2.26.PHARAOH
RG:Z:HG002.trio_hifiasm_0.19.5.DC_1.2_40x.R1041_Dorado_v0.1.1_43x_alignment.minimap2.26
```

### Testing hybrid model 2, docker polisher_v0.2.0_04172024 whole genome

```
#!/bin/bash
#SBATCH --job-name=DP-hybrid_model2
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --time=7-00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    google/deepconsensus:polisher_v0.2.0_04172024 \
    polisher make_images \
    --bam /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/alignments/add_PL_tag/HG002.trio_hifiasm_0.19.5.DCv1.2.PHARAOH.Dorado.R10.secphase.mm2v2.26.merged.PL_tag.bam \
    --fasta /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --output /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/DeepPolisher/HG002_y2_mm2_DCv1.2_R10_Dorado_hybrid_model2/images/images \
    --cpus 64 \
    --use_pl_tags \
    --pl_tag_height 30 \
    --pl_tag_order pacbio,ont

# Inference on images to generate VCFs
docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    google/deepconsensus:polisher_v0.2.0_04172024 \
    polisher inference \
    --input_dir /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/DeepPolisher/HG002_y2_mm2_DCv1.2_R10_Dorado_hybrid_model2/images \
    --out_dir /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/DeepPolisher/HG002_y2_mm2_DCv1.2_R10_Dorado_hybrid_model2/vcf/ \
    --checkpoint /private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-226/checkpoint-226 \
    --reference_file /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --sample_name HG002 \
    --cpus 64
```

### Testing hybrid model 2, docker polisher_v0.2.0_04172024 chr20

Subset merged bam to chr20
```
#!/bin/bash
#SBATCH --job-name=samtools_view
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

samtools view -bh /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/alignments/add_PL_tag/HG002.trio_hifiasm_0.19.5.DCv1.2.PHARAOH.Dorado.R10.secphase.mm2v2.26.merged.PL_tag.bam h2tg000002l h2tg000005l h2tg000006l h2tg000021l h2tg000022l h2tg000024l h2tg000031l h2tg000035l h2tg000062l h2tg000068l h1tg000004l h1tg000011l h1tg000012l h1tg000016l h1tg000017l h1tg000020l h1tg000026l h1tg000028l h1tg000029l h1tg000032l h1tg000036l h1tg000040l h1tg000044l h1tg000054l h1tg000101l h1tg000116l h1tg000204l h1tg000217l > /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/alignments/add_PL_tag/HG002.trio_hifiasm_0.19.5.DCv1.2.PHARAOH.Dorado.R10.secphase.mm2v2.26.merged.PL_tag.chr20.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/alignments/add_PL_tag/HG002.trio_hifiasm_0.19.5.DCv1.2.PHARAOH.Dorado.R10.secphase.mm2v2.26.merged.PL_tag.chr20.bam
```

```
#!/bin/bash
#SBATCH --job-name=DP-hybrid_model2_chr20
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --exclude=phoenix-[09,10,22,23,24]
#SBATCH --time=7-00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    google/deepconsensus:polisher_v0.2.0_04172024 \
    polisher make_images \
    --bam /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/alignments/add_PL_tag/HG002.trio_hifiasm_0.19.5.DCv1.2.PHARAOH.Dorado.R10.secphase.mm2v2.26.merged.PL_tag.chr20.bam \
    --fasta /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.chr20.full_contigs.fa \
    --output /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/DeepPolisher/HG002_y2_mm2_DCv1.2_R10_Dorado_hybrid_model2_chr20/images/images \
    --cpus 64 \
    --use_pl_tags \
    --pl_tag_height 30 \
    --pl_tag_order pacbio,ont

# Inference on images to generate VCFs
docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    google/deepconsensus:polisher_v0.2.0_04172024 \
    polisher inference \
    --input_dir /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/DeepPolisher/HG002_y2_mm2_DCv1.2_R10_Dorado_hybrid_model2_chr20/images \
    --out_dir /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/DeepPolisher/HG002_y2_mm2_DCv1.2_R10_Dorado_hybrid_model2_chr20/vcf/ \
    --checkpoint /private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-226/checkpoint-226 \
    --reference_file /private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.chr20.full_contigs.fa \
    --sample_name HG002 \
    --cpus 64
```
