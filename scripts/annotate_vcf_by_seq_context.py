'''
Purpose: Annotate variants in a vcf file by sequence context
Author: Mira Mastoras, mmastora@ucsc.edu
Usage: python3 annotate_vcf_by_seq_context.py \
            -v /Users/miramastoras/Desktop/test_files_annotate_vcf/HG04115.polisher_output.h2tg000055l.vcf.gz \
            -f /Users/miramastoras/Desktop/test_files_annotate_vcf/HG04115.mat.h2tg000055l.fa \
            -o /Users/miramastoras/Desktop/test_files_annotate_vcf/output.vcf
'''
import argparse
import pysam
from pysam import VariantFile, FastaFile


def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='annotate_vcf_by_seq_context.py',
        description="""Annotate variants in vcf by sequence context (dimer or homopolymer)""")

    parser.add_argument("-f", "--fasta",
                        required=True,
                        help="input fasta file")
    parser.add_argument("-o", "--out",
                        required=True,
                        help="output file")
    parser.add_argument("-v", "--vcf",
                        required=False,
                        help="input vcf file")

    return parser.parse_args()

def extend_dimers(reference_sequence, pattern, start_index):
    end_index = start_index + 2
    while end_index < len(reference_sequence) - 1:
        di_mer = reference_sequence[end_index] + reference_sequence[end_index + 1]
        if di_mer == pattern:
            end_index += 2
        else:
            break

    return end_index

def extend_homopolymers(reference_sequence, pattern, start_index):
    end_index = start_index + 1
    while end_index < len(reference_sequence) - 1:
        mono_mer = reference_sequence[end_index]
        if mono_mer == pattern:
            end_index += 1
        else:
            break

    return end_index

def annotate_repeats(contig,start,end,FastaFileObject):
    '''
    :param contig: contig of region to annotate
    :param start: start coord of region to annotate
    :param end: end coord of region to annotate
    :param FastaFileObject: FastaFile object
    :return: homopolymer_coords: list of regions and their annotation, in the format [[contig,start,end,pattern]]
    '''
    seq = FastaFileObject.fetch(reference=contig, start=start, end=end)
    repeat_coords=[]

    end_index = 0
    dimer_base = '**'
    i = 0
    while i < len(seq) - 1 :
        # if current base is not the same as next base, check for dimer repeats
        if seq[i] != seq[i + 1]:
            dimer_base = seq[i] + seq[i + 1]
            end_index = extend_dimers(seq, dimer_base, i)
            dimer_length=int((end_index - i)/2)

            if dimer_length > 1 and end_index <= len(seq) - 1 :
                repeat_coords.append([contig,i+start,end_index+start,dimer_base])
                i = end_index
            else:
                # if dimer length is only 2 bp, we want to start next check on second base in the pair
                i = end_index - 1

        # if current base is same as next base, check for homopolymer repeats
        else:
            hom_base = seq[i]
            end_index = extend_homopolymers(seq, hom_base, i)
            hom_length = int((end_index - i)/2)

            if hom_length >= 1 and end_index <= len(seq) - 1 :
                repeat_coords.append([contig, i+start,end_index+start, hom_base])

            i = end_index

    return repeat_coords


def main():
    # parse command line arguments
    args = arg_parser()

    fasta_start=12590
    fasta_contig="h2tg000055l"
    fasta_end=12659

    #small_variant_vcf = VariantFile(vcf_file)
    assembly_fasta_file = FastaFile(args.fasta)

    data=annotate_repeats(contig=fasta_contig,start=fasta_start,end=fasta_end,FastaFileObject=assembly_fasta_file)

    with open(args.out, 'w') as outfile:
        for item in data:
            print(*item, sep="\t", file=outfile)


if __name__ == '__main__':
    main()