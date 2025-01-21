## Running DeepPolisher for the T2T marmoset projects

### 1. Polishing initial hifiasm assembly for Baguette

Assemblies:
```
/private/groups/cgl/pnhebbar/marmoset/Baguette/hifiasm_primary_deeppolisher/baguette_hifiasm_hic_phased.asm.hic.hap1.p_ctg.formatted.fa
/private/groups/cgl/pnhebbar/marmoset/Baguette/hifiasm_primary_deeppolisher/baguette_hifiasm_hic_phased.asm.hic.hap2.p_ctg.formatted.fa
```

HiFi (Revio) reads:
```
/private/nanopore/basecalled/marmoset/Baguette/genomic_data/PacBio_HiFi/quick_stats/n50_Baguette_CJA.tsv

# 61.33 X coverage
```
Using all 61 x coverage of HiFi

ONT:
```

ls | grep "summary_stats" | grep "dorado0.5.3" | while read line ;do cat $line | cut -f1-5 ;done

Sample	read_N50	Gb	coverage	100kb+
Baguette_1 	 65206 	 79.52 	 24.1 	 7.11
Sample	read_N50	Gb	coverage	100kb+
Baguette_1 	 94570 	 51.65 	 15.65 	 7.38
Sample	read_N50	Gb	coverage	100kb+
Baguette_2 	 92579 	 53.14 	 16.1 	 7.42
Sample	read_N50	Gb	coverage	100kb+
Baguette_3 	 93075 	 56.85 	 17.23 	 7.97
```
~ 28x coverage of 100kb+

```
ls | grep "dorado0.5.3" | grep ".bam" | while read line ;do realpath $line ;done
```

### 2. Creating confidence regions for marmosets

According to GIAB, the following criteria were used for creating confidence regions:
```
These excluded regions fall in several categories:
(1) the modeled centromere and heterochromatin in GRCh38 because these are highly repetitive regions and generally differ in structure and copy number between any individual and the reference;
(2) the VDJ, which encodes immune system components and undergoes somatic recombination in B cells;
(3) in GRCh37, regions that are either expanded or collapsed relative to GRCh38;
(4) segmental duplications with greater than 5 copies longer than 10 kb and identity greater than 99 %, where errors are likely in mapping and variant calling, e.g., due to structural or copy number variation resulting in calling paralogous sequence variants;26,27
(5) potential large duplications that are in HG002 relative to GRCh37 or GRCh38;
(6) putative insertions, deletions, and inversions >49bp in size and flanking sequence;
(7) tandem repeats larger than 10,000 bp where variants can be difficult to detect accurately given the length of PacBio HiFi reads.
```

1),2),4),and 7) can be annotated directly to the marmoset genomes.

#### 2.1 Generating annotations (Prajna):

#### 2.2 Combining annotations

```
# fix tabs
sed 's/ /\t/g'  /private/groups/cgl/pnhebbar/assemblyPolish/marmoset_confidence/baguette_hifiasm_hic_hap2_centr.bed > /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/Baguette/confidence_regions/baguette_hifiasm_hic_hap2_centr.bed
sed 's/ /\t/g'  /private/groups/cgl/pnhebbar/assemblyPolish/marmoset_confidence/baguette_hifiasm_hic_hap1_centr.bed > /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/Baguette/confidence_regions/baguette_hifiasm_hic_hap1_centr.bed


cat /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/Baguette/confidence_regions/baguette_hifiasm_hic_hap1_centr.bed /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/Baguette/confidence_regions/baguette_hifiasm_hic_hap2_centr.bed /private/groups/cgl/pnhebbar/assemblyPolish/marmoset_confidence/baguette_hifiasm_hic_hap1_vdj.bed /private/groups/cgl/pnhebbar/assemblyPolish/marmoset_confidence/baguette_hifiasm_hic_hap2_vdj.bed /private/groups/cgl/pnhebbar/assemblyPolish/marmoset_confidence/baguette_hifiasm_hic_hap1_SD_length10k_identity99_copies5.bed /private/groups/cgl/pnhebbar/assemblyPolish/marmoset_confidence/baguette_hifiasm_hic_hap2_SD_length10k_identity99_copies5.bed /private/groups/cgl/pnhebbar/assemblyPolish/marmoset_confidence/baguette_hifiasm_hic_hap2_tandem_repeats_length10k.bed /private/groups/cgl/pnhebbar/assemblyPolish/marmoset_confidence/baguette_hifiasm_hic_hap1_tandem_repeats_length10k.bed | bedtools sort -i - | bedtools merge -i - > /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/Baguette/confidence_regions/baguette_hifiasm_confidence.dip.bed
```

Count bases in bed file
```
awk '{sum += $3-$2}END{print sum}' /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/Baguette/confidence_regions/baguette_hifiasm_confidence.dip.bed

168767160

awk '{sum+=$2;} END{print sum;}' /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/Baguette/Baguette_hic_final.dip.fa.fai

5697787823

# only removing like 2.9% of the genome
```

Run mosdepth to get coverage dropout regions
```
#!/bin/bash
#SBATCH --job-name=mosdepth
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups \
    -v /private/groups/patenlab/mira/phoenix_batch_executions/workflows/hprc_DeepPolisher/Baguette_hifiasm_hic/analysis/hprc_DeepPolisher_outputs/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth -Q1 --threads 4 \
    -f /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/Baguette/Baguette_hic_final.dip.fa \
    --quantize 0:5:10:150: \
    /private/groups/patenlab/mira/phoenix_batch_executions/workflows/hprc_DeepPolisher/Baguette_hifiasm_hic/analysis/hprc_DeepPolisher_outputs/mosdepth/Baguette_hic_hifiasm_dip \
    /private/groups/patenlab/mira/phoenix_batch_executions/workflows/hprc_DeepPolisher/Baguette_hifiasm_hic/analysis/hprc_DeepPolisher_outputs/Baguette_hifiasm_hic.hifi.to.diploid.asm.PHARAOH.bam

cd /private/groups/patenlab/mira/phoenix_batch_executions/workflows/hprc_DeepPolisher/Baguette_hifiasm_hic/analysis/hprc_DeepPolisher_outputs/mosdepth/

zcat Baguette_hic_hifiasm_dip.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > Baguette_hic_hifiasm_dip.quantized.tsv

awk 'BEGIN {FS=OFS="\t"} {
        if ($4 == "0:5") {
            $4 = "NO_COVERAGE"
        } else if ($4 == "5:10") {
            $4 = "LOW_COVERAGE"
        } else if ($4 == "10:150") {
            $4 = "CALLABLE"
        } else if ($4 == "150:inf") {
            $4 = "HIGH_COVERAGE"
        }
        print }' Baguette_hic_hifiasm_dip.quantized.tsv > Baguette_hic_hifiasm_dip.quantized_mapQ1.quant.bed

grep NO_COVERAGE Baguette_hic_hifiasm_dip.quantized_mapQ1.quant.bed > Baguette_hic_hifiasm_dip.quantized_mapQ1.quant.lt5x_cov.bed

awk -v OFS="\t" '{if ($3 -$2 >100) print $1,$2,$3 }' Baguette_hic_hifiasm_dip.quantized_mapQ1.quant.lt5x_cov.bed > Baguette_hic_hifiasm_dip.quantized_mapQ1.quant.lt5x_cov.gt100bp.bed
```

Merge confidence regions with mosdepth regions
```
cat /private/groups/patenlab/mira/phoenix_batch_executions/workflows/hprc_DeepPolisher/Baguette_hifiasm_hic/analysis/hprc_DeepPolisher_outputs/mosdepth/Baguette_hic_hifiasm_dip.quantized_mapQ1.quant.lt5x_cov.gt100bp.bed /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/Baguette/confidence_regions/baguette_hifiasm_confidence.dip.bed | bedtools sort -i - | bedtools merge -i - > /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/Baguette/confidence_regions/baguette_hifiasm_confidence.mosdepth_mapQ1_lt5x_cov.gt100bp.dip.bed

awk '{sum += $3-$2}END{print sum}' /private/groups/patenlab/mira/t2t_primates_polishing/assemblies/Baguette/confidence_regions/baguette_hifiasm_confidence.mosdepth_mapQ1_lt5x_cov.gt100bp.dip.bed

500714241 / 5697787823 = 8.7% of genome
```
