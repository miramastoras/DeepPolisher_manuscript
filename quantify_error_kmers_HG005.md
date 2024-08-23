# Quantifying the illumina error kmers that are not being changed by polishing

## 1. Define bed file containing the error kmers unchanged by polishing

the automated wdl workflow [hprc_polishing_QC_no_meryl.wdl](https://github.com/miramastoras/hpp_production_workflows/blob/master/QC/wdl/workflows/hprc_polishing_QC_no_meryl.wdl) was used to run merqury, and to project error kmers from the polished assembly to the raw assembly using [project_blocks_multi_thread.py](https://github.com/mobinasri/flagger/blob/main/programs/src/project_blocks_multi_thread.py)

Combine polished projected error kmers for each haplotype:
```
cd /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs

cat HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.hap1PolToRaw_asm_only.projection.bed HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.hap2PolToRaw_altHap_only.projection.bed > HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ.error_kmers_projected.dip.bed
```

We use bedtools intersect to extract error kmers in the raw assembly that are also in the polished assembly (projected to raw assembly coordinates)
```
bedtools intersect -f 1 -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.dip_only.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ.error_kmers_projected.dip.bed \
    | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed
```

Counting error kmers:
```
# unchanged by polishing
wc -l /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed
# 1574857

# raw assembly
wc -l /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_manuscript/hprc_polishing_QC_no_meryl_k31/HG005_deepvariant/analysis/hprc_polishing_QC_no_meryl_outputs/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.dip_only.bed
# 1617645

# 1574857 / 1617645 = 97.3% of error kmers are unchanged by polishing
```

## 2. Define annotation bed files

Low coverage (<5x)
```
/private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed
```
CHM13 simple repeats and low complexity regions
```
/private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed
```
GC content
```
/private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed
```
Homopolymers > 10bp
```
/private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_homopolymers.gt10.bed
```

### 2.1 Regions of low coverage in the PHARAOH alignments

Defining "no coverage" as <5x coverage to account for the edges of coverage dropout regions, and making the assumption that <5x coverage would reduce DeepPolisher's ability to call edits in those regions. We also restrict to regions of <5x coverage for >100 bp, to avoid picking up long indels in the alignments.

```
#!/bin/bash
#SBATCH --job-name=HG005_mosdepth
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
    -v /private/groups/patenlab/mira/hprc_polishing/hifi_ONT_combined_model/polishing_dropouts/mosdepth/:/opt/mount \
    quay.io/biocontainers/mosdepth:0.2.4--he527e40_0 \
    mosdepth --threads 4 \
    -f /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa \
    --quantize 0:5:10:150: \
    /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized \
    /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_PHv6_DPmm2model1/toil_hprc_deepPolisher_out/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.PHARAOHv6.bam

zcat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.tsv

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
        print }' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.tsv > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.bed

grep NO_COVERAGE /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed

awk -v OFS="\t" '{if ($3 -$2 >100) print $1,$2,$3 }' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.bed  > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed
```
### 2.2 CHM13 low complexity and simple repeats projected to HG005 assembly

Downloaded from:
https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/

```
cd /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation
grep Low_complexity chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed > chm13_simple_low_complexity.bed
grep Simple_repeat chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed >>chm13_simple_low_complexity.bed
```

Dipcall HG005 raw assembly to chm13

```

{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/chm13v2.0.fa.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/chm13v2.0.fa",
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa",
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa",
  "runDipcall.dipcall.referenceIsHS38": true,
  "dipcall.isMaleSample":true
}
```
Run dipcall

```sh
cd /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/dipcall_chm13

export SINGULARITY_CACHEDIR=`pwd`/../cache/.singularity/cache
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/../cache/.cache/miniwdl
export TOIL_SLURM_ARGS="--time=12:00:00 --partition=medium --exclude=phoenix-[09,10,22,23,24,18]"
export TOIL_COORDINATION_DIR=/data/tmp

mkdir -p toil_logs

time toil-wdl-runner \
    --jobStore ./jobstore \
    --stats \
    --clean=never \
    --batchSystem single_machine \
    --maxCores 32 \
    --batchLogsDir ./toil_logs \
    /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/dipcall.wdl \
    dipcall_inputs.json \
    --outputDirectory ./dipcall_outfiles \
    --outputFile dipcall_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --logDebug \
    2>&1 | tee log.txt
```

Project repeats to HG005 with project_blocks_multi_thread.py

https://github.com/mobinasri/flagger/blob/main/programs/src/project_blocks_multi_thread.py.

```
# Removed unmapped records in paf files
cd /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/dipcall_chm13/dipcall_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dipcall/

grep "tp:A:P" HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.paf > ../HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf
grep "tp:A:P" HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.paf > ../HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf

bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/chm13_simple_low_complexity.srt.bed | bedtools merge -i - > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/chm13_simple_low_complexity.srt.mrg.bed

# Project repeat annotation bedfile to y2 assembly, each haplotype separately
# Hap1
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/dipcall_chm13/dipcall_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap1.AP.paf \
  --blocks /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/chm13_simple_low_complexity.srt.mrg.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap1.chm13_simple_low_complexity.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap1.chm13_simple_low_complexity.projection.bed

# Hap 2
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira \
  mobinasri/flagger:latest \
  python3 /home/programs/src/project_blocks_multi_thread.py \
  --threads 10 \
  --mode 'ref2asm' \
  --paf /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/dipcall_chm13/dipcall_outfiles/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.hap2.AP.paf \
  --blocks /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/chm13_simple_low_complexity.srt.mrg.bed \
  --outputProjectable /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap2.chm13_simple_low_complexity.projectable.bed \
  --outputProjection /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap2.chm13_simple_low_complexity.projection.bed

cat /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap2.chm13_simple_low_complexity.projection.bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_hap1.chm13_simple_low_complexity.projection.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed

# Number of bases in repeat annotation
bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed | bedtools merge -i - | awk '{sum += $3-$2}END{print sum}'

# 143921510 / 5972405005 = 2.4% of genome
```

### 2.3 Bed regions of high GC content

```
# calculate GC content for all error kmers unchanged by polishing
bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed
```

### 2.4 Homopolymers > 10 bp

Python script for identifying homopolymers > 10bp
```py
from Bio import SeqIO

def find_homopolymers(sequence, min_length=10):
    """Find homopolymer regions in a sequence."""
    homopolymers = []
    n = len(sequence)
    i = 0
    while i < n:
        start = i
        while i < n and sequence[i] == sequence[start]:
            i += 1
        length = i - start
        if length >= min_length:
            homopolymers.append((start, i, sequence[start]))
    return homopolymers

def fasta_to_bed(fasta_file, bed_file, min_length=10):
    """Convert FASTA file to BED file with homopolymers."""
    with open(bed_file, 'w') as out_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequence = str(record.seq)
            homopolymers = find_homopolymers(sequence, min_length)
            for start, end, base in homopolymers:
                # Output BED format: chromosome, start, end, name, score, strand
                out_file.write(f"{record.id}\t{start}\t{end}\t{base}_homopolymer\t0\t+\n")

fasta_file = '/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa'  # Replace with your FASTA file path
bed_file = '/private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_homopolymers.gt10.bed'     # Replace with desired output BED file path
fasta_to_bed(fasta_file, bed_file)
```

## 3. Count number of error kmers unchanged by polishing in each bed file / annotation

Coverage <5x for >100bp
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed -wa | sort | uniq | wc -l

359392 / 1574857 = 22.8%
```
CHM13 simple and low complexity repeats
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed -wa | sort | uniq | wc -l

569002 / 1574857 = 36.1%
```
GC content
```
# count error kmers with > n % GC content

awk '$5 > 0.9' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 23215 / 1574857 = 1.4%

awk '$5 > 0.8' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 156986 / 1574857 = 10%

awk '$5 > 0.7' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 651829 / 1574857 = 41%

awk '$5 > 0.6' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 1132050 / 1574857 = 71%

# count error kmers with > n % AT content

awk '$4 > 0.9' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 63799 / 1574857 = 4%

awk '$4 > 0.8' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 109984 / 1574857 = 6.9%

awk '$4 > 0.7' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 154277 / 1574857 = 9.7%

awk '$4 > 0.6' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | wc -l
# 204499 / 1574857 = 12.9%
```

Homopolymers > 10bp
```
# homopolymer > 10 bp
bedtools intersect -wa -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_homopolymers.gt10.bed | sort | uniq | wc -l

166080 / 1574857 = 10 %
```

These categories will have overlapping kmers. To obtain the total number of error kmers explained by these categories (and those unexplained) we do the following:

```
# combine all annotated kmers into one file, reporting the entire original kmer in the -a bed file with -wa

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed -wa | sort | uniq >> /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_error_kmers_unchanged_by_polishing.annotated.bed

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed -wa | sort | uniq >> /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_error_kmers_unchanged_by_polishing.annotated.bed

awk '$5 > 0.7' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_ilm_k31_error_kmers_unchanged_by_polishing.nuc.bed | cut -f 1-3 >> /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_error_kmers_unchanged_by_polishing.annotated.bed

bedtools intersect -wa -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_homopolymers.gt10.bed | sort | uniq >> /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_error_kmers_unchanged_by_polishing.annotated.bed

# discard repeats and count all annotated kmers
sort /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/HG005_error_kmers_unchanged_by_polishing.annotated.bed | uniq | wc -l

1221596 / 1574857 = 77%
```

We can account for 77% of the unchanged error kmers with the annotations described above.

## 4. How many error kmers in each annotation category are validated by element data?

hprc polishing QC wdl location and input files : https://github.com/miramastoras/phoenix_batch_submissions/tree/main/polishing/hprc_polishing_QC_no_meryl/GIAB_samples_DP_manuscript_element_k31

Get element FP kmers unchanged by polishing
```
# combine projection bed files
cat /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_DP_manuscript_element_k31/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap1.hap1PolToRaw_asm_only.projection.bed /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_DP_manuscript_element_k31/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/HG005_GQ20_INS1_GQ12_DEL1_GQ5_else_hap2.hap2PolToRaw_altHap_only.projection.bed > /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_DP_manuscript_element_k31/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/HG005_polished_dip_FP_kmers_elementk31_projection.bed

bedtools intersect -f 1 -wa \
    -a /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_DP_manuscript_element_k31/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/raw_results/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.merqury.dip_only.bed \
    -b /private/groups/patenlab/mira/hprc_polishing/polisher_evaluation/GIAB_samples_DP_manuscript_element_k31/hprc_polishing_QC_no_meryl/HG005_y2_DCv1.2_PHv6_mm2_model1_dockerv0.8_HPRC_GQ/analysis/hprc_polishing_QC_no_meryl_outputs/HG005_polished_dip_FP_kmers_elementk31_projection.bed \
    | sort | uniq > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed

# 3097256 / 3142727 = 98.5% unchanged
```

Illumina error kmers unchanged which are also unchanged in element

```
# in illumina, in element
bedtools intersect -f 1 -r \
    -wa -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed \
    -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed

wc -l /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed

# 948307 / 1574857 = 60%
```
How many overlap with each annotation?

Coverage <5x for >100bp
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed -wa | sort | uniq | wc -l

294555 / 948307 = 31 %
```
CHM13 simple and low complexity repeats
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed -wa | sort | uniq | wc -l

315053 / 948307 = 33.2%
```
GC content
```
bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed

# count error kmers with > n % GC content

awk '$5 > 0.7' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | wc -l
# 419247 / 948307 = 44.2%
```

Homopolymers > 10bp
```
# homopolymer > 10 bp
bedtools intersect -wa -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_homopolymers.gt10.bed | sort | uniq | wc -l

34263 / 948307 = 10 %
```

How much of the error kmers validated by element do these annotations account for? :
```
bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth_5/HG005_y2_DCv1.2_40x.PHARAOHv6.quantized.quant.NO_COVERAGE.5x.gt100bp.bed -wa | sort | uniq >> /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.ANNOTATED.bed

bedtools intersect -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/repeat_annotation/HG005_dip.chm13_simple_low_complexity.projection.bed -wa | sort | uniq >> /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.ANNOTATED.bed

awk '$5 > 0.7' /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.nuc.bed | cut -f1-3 >> /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.ANNOTATED.bed

bedtools intersect -wa -a /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed -b /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/homopolymers/HG005_homopolymers.gt10.bed | sort | uniq >> /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.ANNOTATED.bed

sort /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.ANNOTATED.bed | uniq | wc -l

# 749030 / 948307 = 78%
```

## 5. Is there a GC bias in both illumina and element error kmers?
Merge error kmers of the three categories:
```
# illumina unchanged by polishing
bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/mosdepth/HG005_error_kmers_unchanged_by_polishing.bed | bedtools merge -i - > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.bed

# illumina and element
bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.bed | bedtools merge -i - > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.bed

# just element unchanged
bedtools sort -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_element/HG005_error_kmers_unchanged_by_polishing.elementk31.bed | bedtools merge -i - > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.bed
```

Calculate GC and AT content for each
```
bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.nuc.bed

bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.nuc.bed

bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.nuc.bed  
```

Get txt file of counts
```
ls | grep .nuc.bed | while read line
    do cut -f4 $line | grep -v "_pct" > ${line}.AT.txt
    cut -f5 $line | grep -v "_pct" > ${line}.GC.txt
    done  
```

randomly permute same size regions for each file across genome
```
bedtools shuffle -seed 2400 -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.bed -g /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.fai.genome > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.RANDOM.bed

bedtools shuffle -seed 2400 -i /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.bed -g /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.fai.genome > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.RANDOM.bed

bedtools shuffle -seed 2400 -i  /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.bed -g /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.fai.genome > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.RANDOM.bed
```

Calculate GC and AC content for the random shuffled regions
```
bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.RANDOM.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.mrg.srt.RANDOM.nuc.bed

bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.RANDOM.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_k31_illumina_kmers_unchanged_by_pol_in_elementk31.srt.mrg.RANDOM.nuc.bed

bedtools nuc -bed /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.RANDOM.bed -fi /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa > /private/groups/patenlab/mira/hprc_polishing/qv_problems/HG005_coverage/GC_content/HG005_error_kmers_unchanged_by_polishing.elementk31.srt.mrg.RANDOM.nuc.bed
```
Plot histograms: https://colab.research.google.com/drive/18kY97xhAcgFaaGU6tH2L-IYEnmbbvC0E?authuser=1#scrollTo=dtQs42ojHGsY
