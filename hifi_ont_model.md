### Adding a PL tag to distinguish between HiFi and ONT reads

Bam files:
```
# HiFi
HG002.trio_hifiasm_0.19.5.DCv1.2.minimap2v2.26.PHARAOH.bam
# ONT
HG002.trio_hifiasm_0.19.5.DC_1.2_40x.R1041_Dorado_v0.1.1_43x_alignment.minimap2.26.bam
```

Create header mapping file. RG ID must be the prefix of each .bam file 
```
printf '@RG\tID:HG002.trio_hifiasm_0.19.5.DCv1.2.minimap2v2.26.PHARAOH\tPL:PACBIO\n@RG\tID:HG002.trio_hifiasm_0.19.5.DC_1.2_40x.R1041_Dorado_v0.1.1_43x_alignment.minimap2.26\tPL:ONT\n' > rg.txt
```

Merge hifi+ont bam
```
samtools merge -rh rg.txt hybrid.bam HG002.trio_hifiasm_0.19.5.DCv1.2.minimap2v2.26.PHARAOH.bam HG002.trio_hifiasm_0.19.5.DC_1.2_40x.R1041_Dorado_v0.1.1_43x_alignment.minimap2.26.bam
```

View new header tags
```
samtools view merged.bam -H | grep @RG

@RG	ID:HG002.trio_hifiasm_0.19.5.DCv1.2.minimap2v2.26.PHARAOH	PL:PACBIO
@RG	ID:HG002.trio_hifiasm_0.19.5.DC_1.2_40x.R1041_Dorado_v0.1.1_43x_alignment.minimap2.26	PL:ONT
```

Each read will have the following RG tag:
```
RG:Z:HG002.trio_hifiasm_0.19.5.DCv1.2.minimap2v2.26.PHARAOH
RG:Z:HG002.trio_hifiasm_0.19.5.DC_1.2_40x.R1041_Dorado_v0.1.1_43x_alignment.minimap2.26
```
