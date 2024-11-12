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

Run illumina meryl
```

```
