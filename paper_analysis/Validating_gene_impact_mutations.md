## Validating impacted genes with GIAB benchmark variant calls

### Liftover locations of mutations from assemblies to GRCh38

1. Download mutation loci provided by Prajna

https://public.gi.ucsc.edu/~pnhebbar/assemblyPolishing_mutationLoci/

Location:
```
cd /private/groups/patenlab/mira/hprc_polishing/gene_impact
```

Count number of nonsynonymous mutations
```
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do num=`cat ${sample}_mutation_loci_${hap}_${asm}/genomic_loci_nonsynonymous.txt | wc -l`
            echo ${sample}_${hap}_${asm},${num}
            done
          done
        done

HG002_hap1_pol,68904
HG002_hap2_pol,39288
HG002_hap1_raw,68906
HG002_hap2_raw,39288
HG005_hap1_pol,67139
HG005_hap2_pol,21324
HG005_hap1_raw,67147
HG005_hap2_raw,21327
```

Convert mutation loci to bed format
```
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do
                awk -v FS=: '{print $1"\t"$2-5"\t"$2+5}' ${sample}_mutation_loci_${hap}_${asm}/genomic_loci_nonsynonymous.txt | bedtools sort -i - | bedtools merge -c 1 -o count -i - > ${sample}_mutation_loci_${hap}_${asm}/genomic_loci_nonsynonymous.bed
                awk -v FS=: '{print $1"\t"$2-5"\t"$2+5}' ${sample}_mutation_loci_${hap}_${asm}/genomic_loci_frameshift.txt | bedtools sort -i - | bedtools merge -c 1 -o count -i - > ${sample}_mutation_loci_${hap}_${asm}/genomic_loci_frameshift.bed
                awk -v FS=: '{print $1"\t"$2-5"\t"$2+5}' ${sample}_mutation_loci_${hap}_${asm}/genomic_loci_early_stop.txt | bedtools sort -i - | bedtools merge -c 1 -o count -i - > ${sample}_mutation_loci_${hap}_${asm}/genomic_loci_early_stop.bed
            done
          done
        done
```

Copy paf files from dipcall to working directory
```
cp /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf /private/groups/patenlab/mira/hprc_polishing/gene_impact/paf_files/HG005_hap1_raw_grch38.paf

cp /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG005_y2_raw/dipCallTar/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf /private/groups/patenlab/mira/hprc_polishing/gene_impact/paf_files/HG005_hap2_raw_grch38.paf

cp /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG002_y2_raw/dipCallTar/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf /private/groups/patenlab/mira/hprc_polishing/gene_impact/paf_files/HG002_hap1_raw_grch38.paf

cp /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/y2_terra_tables/y2_polisher_evaluation/HG002_y2_raw/dipCallTar/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf /private/groups/patenlab/mira/hprc_polishing/gene_impact/paf_files/HG002_hap2_raw_grch38.paf

cp /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.dipcall/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_polished.hap1.AP.paf /private/groups/patenlab/mira/hprc_polishing/gene_impact/paf_files/HG005_hap1_pol_grch38.paf

cp /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.dipcall/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_polished.hap2.AP.paf /private/groups/patenlab/mira/hprc_polishing/gene_impact/paf_files/HG005_hap2_pol_grch38.paf

cp /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.dipcall/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_polished.hap1.AP.paf /private/groups/patenlab/mira/hprc_polishing/gene_impact/paf_files/HG002_hap1_pol_grch38.paf

cp /private/groups/patenlab/mira/hprc_polishing/qv_problems/HPRC_intermediate_asm/GQ_filters/GIAB/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else/applyPolish_dipcall_outputs/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.polished.dipcall/HG002_GQ20_INS1_GQ12_DEL1_GQ5_else_polished.hap2.AP.paf /private/groups/patenlab/mira/hprc_polishing/gene_impact/paf_files/HG002_hap2_pol_grch38.paf
```

Use project blocks script to liftover mutations loci
```
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do
            docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
            mobinasri/flagger:latest \
            python3 /home/programs/src/project_blocks_multi_thread.py \
            --threads 10 \
            --mode 'asm2ref' \
            --paf /private/groups/patenlab/mira/hprc_polishing/gene_impact/paf_files/${sample}_${hap}_${asm}_grch38.paf \
            --blocks /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_early_stop.bed \
            --outputProjectable /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_early_stop.projectable.bed \
            --outputProjection /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_early_stop.projection.bed
            done
          done
        done

# count number of projectable mutations
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do num=`cat /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_early_stop.projectable.bed | awk '{sum += $4}END{print sum}'`
            echo ${sample}_${hap}_${asm},$num
            done
          done
        done

HG002_hap1_pol,40
HG002_hap2_pol,23
HG002_hap1_raw,40
HG002_hap2_raw,23
HG005_hap1_pol,44
HG005_hap2_pol,17
HG005_hap1_raw,44
HG005_hap2_raw,17
```

```
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do
            docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
            mobinasri/flagger:latest \
            python3 /home/programs/src/project_blocks_multi_thread.py \
            --threads 10 \
            --mode 'asm2ref' \
            --paf /private/groups/patenlab/mira/hprc_polishing/gene_impact/paf_files/${sample}_${hap}_${asm}_grch38.paf \
            --blocks /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_frameshift.bed \
            --outputProjectable /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_frameshift.projectable.bed \
            --outputProjection /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_frameshift.projection.bed
            done
          done
        done

# count number of projectable mutations
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do num=`cat /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_frameshift.projectable.bed | awk '{sum += $4}END{print sum}'`
            echo ${sample}_${hap}_${asm},$num
            done
          done
        done
```

```
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do
            docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
            mobinasri/flagger:latest \
            python3 /home/programs/src/project_blocks_multi_thread.py \
            --threads 10 \
            --mode 'asm2ref' \
            --paf /private/groups/patenlab/mira/hprc_polishing/gene_impact/paf_files/${sample}_${hap}_${asm}_grch38.paf \
            --blocks /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_nonsynonymous.bed \
            --outputProjectable /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_nonsynonymous.projectable.bed \
            --outputProjection /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_nonsynonymous.projection.bed
            done
          done
        done

# count number of projectable mutations
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do num=`cat /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_nonsynonymous.projectable.bed | awk '{sum += $4}END{print sum}'`
            echo ${sample}_${hap}_${asm},$num
            done
          done
        done
```

Intersect projection files with GIAB benchmark vcf
```
# nonsynonymous
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do num=`bedtools intersect -wo -a /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_nonsynonymous.projection.bed -b /private/groups/patenlab/mira/data/${sample}_GRCh38_1_22_v4.2.1_benchmark.vcf | sort | uniq | awk '{sum += $4}END{print sum}'`
            echo ${sample}_${hap}_${asm},$num
            done
          done
        done

# frameshift
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do num=`bedtools intersect -wo -a /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_frameshift.projection.bed -b /private/groups/patenlab/mira/data/${sample}_GRCh38_1_22_v4.2.1_benchmark.vcf | sort | uniq | awk '{sum += $4}END{print sum}'`
            echo ${sample}_${hap}_${asm},$num
            done
          done
        done

# early stop
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do num=`bedtools intersect -wo -a /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_early_stop.projection.bed -b /private/groups/patenlab/mira/data/${sample}_GRCh38_1_22_v4.2.1_benchmark.vcf | sort | uniq | awk '{sum += $4}END{print sum}'`
            echo ${sample}_${hap}_${asm},$num
            done
          done
        done
```

Save validated coordinates
```
# nonsynonymous
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do bedtools intersect -wo -a /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_nonsynonymous.projection.bed -b /private/groups/patenlab/mira/data/${sample}_GRCh38_1_22_v4.2.1_benchmark.vcf | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_nonsynonymous.projection.VALIDATED.bed
            done
          done
        done

# frameshift
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do bedtools intersect -wo -a /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_frameshift.projection.bed -b /private/groups/patenlab/mira/data/${sample}_GRCh38_1_22_v4.2.1_benchmark.vcf | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_frameshift.projection.VALIDATED.bed
            done
          done
        done

# early stop
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
            do bedtools intersect -wo -a /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_early_stop.projection.bed -b /private/groups/patenlab/mira/data/${sample}_GRCh38_1_22_v4.2.1_benchmark.vcf | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_early_stop.projection.VALIDATED.bed
            done
          done
        done
```

```
for sample in HG002 HG005
    do for asm in pol raw
        do for hap in hap1 hap2
           do ls /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_early_stop.projection.VALIDATED.bed
           ls /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_frameshift.projection.VALIDATED.bed
           ls /private/groups/patenlab/mira/hprc_polishing/gene_impact/${sample}_mutation_loci_${hap}_${asm}/genomic_loci_nonsynonymous.projection.VALIDATED.bed
           done
         done
       done
```
