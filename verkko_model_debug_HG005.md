## Debugging verkko model with HG005

### 1. Rerun polishing QC in only the alignable regions between the assemblies

Align raw hifiasm and verkko assemblies
```
#!/bin/bash
#SBATCH --job-name=align_pat
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/migalab/juklucas/polishing/HG005/verkko/trio_asm/assembly.haplotype1.fasta \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa -o /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/HG005.hifiasm_raw_to_verkko.pat.paf

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/migalab/juklucas/polishing/HG005/verkko/trio_asm/assembly.haplotype2.fasta \
    /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa -o /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/HG005_hifiasm_raw_to_verkko.mat.paf
```

Align polished assemblies

```
#!/bin/bash
#SBATCH --job-name=align_pol_mat
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=high_priority
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_verkko_model1_hprc_filters/applyPolish_dipcall_outputs/HG005_verkko_model1_hprc_filters_hap1.polished.fasta \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_20x_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/applyPolish_dipcall_outputs/HG002_y2_DCv1.2_20x_PHv6_mm2_model1_dockerv0.8_HPRC_GQ_hap1.polished.fasta -o /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/HG005.hifiasm_polished_to_verkko.pat.paf

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_verkko_model1_hprc_filters/applyPolish_dipcall_outputs/HG005_verkko_model1_hprc_filters_hap2.polished.fasta \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG002_y2_DCv1.2_20x_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/applyPolish_dipcall_outputs/HG002_y2_DCv1.2_20x_PHv6_mm2_model1_dockerv0.8_HPRC_GQ_hap2.polished.fasta -o  /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/HG005_hifiasm_polished_to_verkko.mat.paf
```

Convert paf to bed format
```
cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions

# Convert paf file of alignment to bed file

cut -f 1,3,4 /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/HG005_hifiasm_polished_to_verkko.dip.paf > HG005_hifiasm_polished_alignable.bed

cut -f 6,8,9 /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/HG005_hifiasm_polished_to_verkko.dip.paf > HG005_verkko_polished_alignable.bed

cut -f 1,3,4 /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/HG005.hifiasm_raw_to_verkko.dip.paf > HG005_hifiasm_raw_alignable.bed

cut -f 6,8,9 /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/HG005.hifiasm_raw_to_verkko.dip.paf > HG005_verkko_raw_alignable.bed

```

Subset assemblies to alignable regions only
```
cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies

# verkko raw
bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_verkko/HG005_verkko_paternal.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/HG005_verkko_raw_alignable.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Raw.dip.WG.alignable.subBed.pat.fasta

bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_verkko/HG005_verkko_maternal.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/HG005_verkko_raw_alignable.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Raw.dip.WG.alignable.subBed.mat.fasta

# verkko polished
bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_verkko_model1_hprc_filters/applyPolish_dipcall_outputs/HG005_verkko_model1_hprc_filters_hap2.polished.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/HG005_verkko_raw_alignable.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Pol.pat.WG.alignable.subBed.fasta

bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_verkko_model1_hprc_filters/applyPolish_dipcall_outputs/HG005_verkko_model1_hprc_filters_hap1.polished.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/HG005_verkko_raw_alignable.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Pol.mat.WG.alignable.subBed.fasta

# hifiasm raw
bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/HG005_hifiasm_raw_alignable.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.raw.pat.WG.alignable.subBed.fasta

bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/HG005_hifiasm_raw_alignable.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.raw.mat.WG.alignable.subBed.fasta

# hifiasm polished
bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/HG005_hifiasm_raw_alignable.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.pol.pat.WG.alignable.subBed.fasta

bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.polished.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/HG005_hifiasm_raw_alignable.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.pol.mat.WG.alignable.subBed.fasta
```

align subsetted assemblies against grch38
```
srun \
  --job-name "align" \
  --cpus-per-task 32 \
  --partition high_priority \
  --mem 250G \
  --time 1:00:00 \
  --pty bash \
  -i

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/hprc/ref_files/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.pol.mat.WG.alignable.subBed.fasta -o /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/hifiasm.pol.mat.WG.alignable.subBed.grch38.paf

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/hprc/ref_files/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.pol.pat.WG.alignable.subBed.fasta -o /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/hifiasm.pol.pat.WG.alignable.subBed.grch38.paf

s
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/hprc/ref_files/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.raw.mat.WG.alignable.subBed.fasta -o /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/hifiasm.raw.mat.WG.alignable.subBed.grch38.paf

s
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/hprc/ref_files/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.raw.pat.WG.alignable.subBed.fasta -o /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/hifiasm.raw.pat.WG.alignable.subBed.grch38.paf

s
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/hprc/ref_files/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Pol.mat.WG.alignable.subBed.fasta -o /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/Verkko.Pol.mat.WG.alignable.subBed.grch38.paf

s
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/hprc/ref_files/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Pol.pat.WG.alignable.subBed.fasta -o /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/Verkko.Pol.pat.WG.alignable.subBed.grch38.paf

s
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/hprc/ref_files/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Raw.dip.WG.alignable.subBed.mat.fasta -o /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/Verkko.Raw.dip.WG.alignable.subBed.mat.grch38.paf

s
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    mobinasri/long_read_aligner:v0.3.3 \
    minimap2 --eqx -x asm5 -t32 -c --cs --secondary=no \
    /private/groups/hprc/ref_files/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Raw.dip.WG.alignable.subBed.pat.fasta -o /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/Verkko.Raw.dip.WG.alignable.subBed.pat.grch38.paf
```


Project confidence regions to each assembly
```
# hifiasm polished
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/:/private/groups/ \
    mobinasri/flagger:latest \
    python3 /home/programs/src/project_blocks_multi_thread.py --threads 32 --mode 'ref2asm' \
    --paf /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/hifiasm.pol.pat.WG.alignable.subBed.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projectable.polished.hap1.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projection.polished.hap1.bed

s
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/:/private/groups/ \
    mobinasri/flagger:latest \
    python3 /home/programs/src/project_blocks_multi_thread.py --threads 32 --mode 'ref2asm' \
    --paf /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/hifiasm.pol.mat.WG.alignable.subBed.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projectable.polished.hap2.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projection.polished.hap2.bed


# hifiasm raw
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/:/private/groups/ \
    mobinasri/flagger:latest \
    python3 /home/programs/src/project_blocks_multi_thread.py --threads 32 --mode 'ref2asm' \
    --paf /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/hifiasm.raw.pat.WG.alignable.subBed.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projectable.raw.hap1.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projection.raw.hap1.bed

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/:/private/groups/ \
    mobinasri/flagger:latest \
    python3 /home/programs/src/project_blocks_multi_thread.py --threads 32 --mode 'ref2asm' \
    --paf /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/hifiasm.raw.mat.WG.alignable.subBed.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projectable.raw.hap2.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projection.raw.hap2.bed


# verkko polished
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/:/private/groups/ \
    mobinasri/flagger:latest \
    python3 /home/programs/src/project_blocks_multi_thread.py --threads 32 --mode 'ref2asm' \
    --paf /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/Verkko.Pol.mat.WG.alignable.subBed.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projectable.polished.mat.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projection.polished.mat.bed

s
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/:/private/groups/ \
    mobinasri/flagger:latest \
    python3 /home/programs/src/project_blocks_multi_thread.py --threads 32 --mode 'ref2asm' \
    --paf /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/Verkko.Pol.pat.WG.alignable.subBed.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projectable.polished.pat.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projection.polished.pat.bed


# verkko raw
docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/:/private/groups/ \
    mobinasri/flagger:latest \
    python3 /home/programs/src/project_blocks_multi_thread.py --threads 32 --mode 'ref2asm' \
    --paf /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/Verkko.Raw.dip.WG.alignable.subBed.mat.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projectable.raw.mat.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projection.raw.mat.bed

docker run --rm -u `id -u`:`id -g` \
    -v /private/groups/:/private/groups/ \
    mobinasri/flagger:latest \
    python3 /home/programs/src/project_blocks_multi_thread.py --threads 32 --mode 'ref2asm' \
    --paf /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/grch38/Verkko.Raw.dip.WG.alignable.subBed.pat.grch38.paf \
    --blocks /private/groups/hprc/ref_files/giab/HG002_intersect_HG005_GIAB_v4.2.1.bed \
    --outputProjectable /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projectable.raw.pat.bed \
    --outputProjection /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projection.raw.pat.bed
```

Subset assemblies to confidence regions

```
bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Raw.dip.WG.alignable.subBed.mat.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projection.raw.mat.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/Verkko.Raw.dip.WG.alignable.subBed.mat.GIABconf.fasta

bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Raw.dip.WG.alignable.subBed.pat.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projection.raw.pat.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/Verkko.Raw.dip.WG.alignable.subBed.pat.GIABconf.fasta

# verkko polished
bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Pol.mat.WG.alignable.subBed.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projection.polished.mat.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/Verkko.Pol.mat.WG.alignable.subBed.GIABconf.asta

bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Pol.pat.WG.alignable.subBed.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.verkko.GIABconf.projection.polished.pat.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/Verkko.Pol.pat.WG.alignable.subBed.GIABconf.asta


# hifiasm raw
bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.raw.mat.WG.alignable.subBed.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projection.raw.hap2.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/hifiasm.raw.mat.WG.alignable.subBed.GIABconf.fasta

bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.raw.pat.WG.alignable.subBed.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projection.raw.hap1.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/hifiasm.raw.pat.WG.alignable.subBed.GIABconf.fasta


# hifiasm polished
bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.pol.mat.WG.alignable.subBed.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projection.polished.hap2.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/hifiasm.pol.mat.WG.alignable.subBed.GIABconf.fasta

bedtools getfasta -fi /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.pol.pat.WG.alignable.subBed.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/HG005.hifiasm.GIABconf.projection.polished.hap1.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/hifiasm.pol.pat.WG.alignable.subBed.GIABconf.fasta


```

Run merqury on whole genome "alignable" only assemblies
```
srun \
  --job-name "merqury" \
  --cpus-per-task 32 \
  --partition high_priority \
  --mem 250G \
  --time 5:00:00 \
  --pty bash \
  -i


mkdir -p /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_hifiasm_polished_alignable_WG

cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_hifiasm_polished_alignable_WG

docker run \
    -it \
    -v /private/groups/:/private/groups/ \
    -v $(pwd):/data/ \
    -u 30162:620 \
    juklucas/hpp_merqury@sha256:387bb58c820c5825e8cf29bc0f8678975ae6a9e6d350c56d5eb5865a6d3d5e82 \
    /bin/bash

merqury.sh \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k31.meryl \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.pol.mat.WG.alignable.subBed.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.pol.pat.WG.alignable.subBed.fasta \
    HG005_hifiasm_polished_alignable_WG
```

```
srun \
  --job-name "merqury" \
  --cpus-per-task 32 \
  --partition high_priority \
  --mem 250G \
  --time 5:00:00 \
  --pty bash \
  -i


mkdir -p /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_hifiasm_raw_alignable_WG

cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_hifiasm_raw_alignable_WG

docker run \
    -it \
    -v /private/groups/:/private/groups/ \
    -v $(pwd):/data/ \
    -u 30162:620 \
    juklucas/hpp_merqury@sha256:387bb58c820c5825e8cf29bc0f8678975ae6a9e6d350c56d5eb5865a6d3d5e82 \
    /bin/bash

merqury.sh \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k31.meryl \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.raw.mat.WG.alignable.subBed.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/hifiasm.raw.pat.WG.alignable.subBed.fasta \
    HG005_hifiasm_raw_alignable_WG
```

```
srun \
  --job-name "merqury" \
  --cpus-per-task 32 \
  --partition high_priority \
  --mem 250G \
  --time 5:00:00 \
  --pty bash \
  -i


mkdir -p /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_verkko_polished_alignable_WG

cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_verkko_polished_alignable_WG

docker run \
    -it \
    -v /private/groups/:/private/groups/ \
    -v $(pwd):/data/ \
    -u 30162:620 \
    juklucas/hpp_merqury@sha256:387bb58c820c5825e8cf29bc0f8678975ae6a9e6d350c56d5eb5865a6d3d5e82 \
    /bin/bash

merqury.sh \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k31.meryl \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Pol.mat.WG.alignable.subBed.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Pol.pat.WG.alignable.subBed.fasta \
    HG005_verkko_polished_alignable_WG
```

```
srun \
  --job-name "merqury" \
  --cpus-per-task 32 \
  --partition high_priority \
  --mem 250G \
  --time 5:00:00 \
  --pty bash \
  -i

mkdir -p /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_verkko_raw_alignable_WG

cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_verkko_raw_alignable_WG

docker run \
    -it \
    -v /private/groups/:/private/groups/ \
    -v $(pwd):/data/ \
    -u 30162:620 \
    juklucas/hpp_merqury@sha256:387bb58c820c5825e8cf29bc0f8678975ae6a9e6d350c56d5eb5865a6d3d5e82 \
    /bin/bash

merqury.sh \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k31.meryl \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Raw.dip.WG.alignable.subBed.mat.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies/Verkko.Raw.dip.WG.alignable.subBed.pat.fasta \
    HG005_verkko_raw_alignable_WG
```


Run merqury on assemblies subsetted to confidence regions
```
srun \
  --job-name "merqury" \
  --cpus-per-task 32 \
  --partition high_priority \
  --mem 250G \
  --time 5:00:00 \
  --pty bash \
  -i

mkdir -p /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_verkko_raw_alignable_conf

cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_verkko_raw_alignable_conf

docker run \
    -it \
    -v /private/groups/:/private/groups/ \
    -v $(pwd):/data/ \
    -u 30162:620 \
    juklucas/hpp_merqury@sha256:387bb58c820c5825e8cf29bc0f8678975ae6a9e6d350c56d5eb5865a6d3d5e82 \
    /bin/bash

merqury.sh \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k31.meryl \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/Verkko.Raw.dip.WG.alignable.subBed.mat.GIABconf.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/Verkko.Raw.dip.WG.alignable.subBed.mat.GIABconf.fasta \
    HG005_verkko_raw_alignable_conf
```

```
srun \
  --job-name "merqury" \
  --cpus-per-task 32 \
  --partition high_priority \
  --mem 250G \
  --time 5:00:00 \
  --pty bash \
  -i

mkdir -p /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_verkko_pol_alignable_conf

cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_verkko_pol_alignable_conf

docker run \
    -it \
    -v /private/groups/:/private/groups/ \
    -v $(pwd):/data/ \
    -u 30162:620 \
    juklucas/hpp_merqury@sha256:387bb58c820c5825e8cf29bc0f8678975ae6a9e6d350c56d5eb5865a6d3d5e82 \
    /bin/bash

merqury.sh \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k31.meryl \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/Verkko.Pol.mat.WG.alignable.subBed.GIABconf.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/Verkko.Pol.pat.WG.alignable.subBed.GIABconf.fasta \
    HG005_verkko_pol_alignable_conf
```
```
srun \
  --job-name "merqury" \
  --cpus-per-task 32 \
  --partition high_priority \
  --mem 250G \
  --time 5:00:00 \
  --pty bash \
  -i

mkdir -p /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_hifiasm_pol_alignable_conf

cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_hifiasm_pol_alignable_conf

docker run \
    -it \
    -v /private/groups/:/private/groups/ \
    -v $(pwd):/data/ \
    -u 30162:620 \
    juklucas/hpp_merqury@sha256:387bb58c820c5825e8cf29bc0f8678975ae6a9e6d350c56d5eb5865a6d3d5e82 \
    /bin/bash

merqury.sh \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k31.meryl \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/hifiasm.pol.mat.WG.alignable.subBed.GIABconf.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/hifiasm.pol.pat.WG.alignable.subBed.GIABconf.fasta \
    HG005_hifiasm_pol_alignable_conf
```

```
srun \
  --job-name "merqury" \
  --cpus-per-task 32 \
  --partition high_priority \
  --mem 250G \
  --time 5:00:00 \
  --pty bash \
  -i

mkdir -p /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_hifiasm_raw_alignable_conf

cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/merqury_alignable_wg/HG005_hifiasm_raw_alignable_conf

docker run \
    -it \
    -v /private/groups/:/private/groups/ \
    -v $(pwd):/data/ \
    -u 30162:620 \
    juklucas/hpp_merqury@sha256:387bb58c820c5825e8cf29bc0f8678975ae6a9e6d350c56d5eb5865a6d3d5e82 \
    /bin/bash

merqury.sh \
    /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/meryl_dbs/HG005.ilm.k31.meryl \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/hifiasm.raw.mat.WG.alignable.subBed.GIABconf.fasta \
    /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/align_hifiasm_verkko/qv_only_alignable_regions/alignable_assemblies_conf/hifiasm.raw.pat.WG.alignable.subBed.GIABconf.fasta \
    HG005_hifiasm_raw_alignable_conf
```

### Re-run merqury removing any contigs < 100kb

```
cd /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb
```


```
# verkko polished

awk -v OFS="\t" '{if($2>100000)print $1,"0",$2,$1}' /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_verkko_model1_hprc_filters/applyPolish_dipcall_outputs/HG005_verkko_model1_hprc_filters_hap1.polished.fasta.fai > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_model1_hprc_filters_hap1.polished.fasta.100kb.bed

bedtools getfasta -name -fi /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_verkko_model1_hprc_filters/applyPolish_dipcall_outputs/HG005_verkko_model1_hprc_filters_hap1.polished.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_model1_hprc_filters_hap1.polished.fasta.100kb.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_model1_hprc_filters_hap1.polished.100kb.fasta


awk -v OFS="\t" '{if($2>100000)print $1,"0",$2,$1}' /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_verkko_model1_hprc_filters/applyPolish_dipcall_outputs/HG005_verkko_model1_hprc_filters_hap2.polished.fasta.fai > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_model1_hprc_filters_hap2.polished.fasta.100kb.bed

bedtools getfasta -name -fi /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/applyPolish_dipcall_happy/HG005_verkko_model1_hprc_filters/applyPolish_dipcall_outputs/HG005_verkko_model1_hprc_filters_hap2.polished.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_model1_hprc_filters_hap2.polished.fasta.100kb.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_model1_hprc_filters_hap2.polished.100kb.fasta

# verkko raw
awk -v OFS="\t" '{if($2>100000)print $1,"0",$2,$1}' /private/groups/patenlab/mira/hprc_polishing/data/HG005_verkko/HG005_verkko_paternal.fasta.fai > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_paternal.fasta.100kb.bed

bedtools getfasta -name -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_verkko/HG005_verkko_paternal.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_paternal.fasta.100kb.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_paternal.100kb.fasta

awk -v OFS="\t" '{if($2>100000)print $1,"0",$2,$1}' /private/groups/patenlab/mira/hprc_polishing/data/HG005_verkko/HG005_verkko_maternal.fasta.fai > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_maternal.fasta.100kb.bed

bedtools getfasta -name -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_verkko/HG005_verkko_maternal.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_maternal.fasta.100kb.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_verkko_maternal.100kb.fasta

# hifiasm polished
awk -v OFS="\t" '{if($2>100000)print $1,"0",$2,$1}' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.fasta.fai > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_hifiasm_mm2_model1_hap1.fasta.100kb.bed

bedtools getfasta -name -fi /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_hifiasm_mm2_model1_hap1.fasta.100kb.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_hifiasm_mm2_model1_hap1.100kb.fasta

awk -v OFS="\t" '{if($2>100000)print $1,"0",$2,$1}' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.polished.fasta.fai > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_hifiasm_mm2_model1_hap2.fasta.100kb.bed

bedtools getfasta -name -fi /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.polished.fasta -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_hifiasm_mm2_model1_hap2.fasta.100kb.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005_hifiasm_mm2_model1_hap2.100kb.fasta

# hifiasm raw
awk -v OFS="\t" '{if($2>100000)print $1,"0",$2,$1}' /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa.fai > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.100kb.bed

bedtools getfasta -name -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.100kb.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.100kb.fasta

awk -v OFS="\t" '{if($2>100000)print $1,"0",$2,$1}' /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa.fai > /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.100kb.bed

bedtools getfasta -name -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa -bed /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.100kb.bed -fo /private/groups/patenlab/mira/hprc_polishing/investigate_variants/HG005_verkko_model1/merqury_gt_100kb/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.100kb.fasta
```
