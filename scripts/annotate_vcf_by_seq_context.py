'''
Purpose: Annotate variants in a vcf file by sequence context
Author: Mira Mastoras, mmastora@ucsc.edu
Usage: python3 annotate_vcf_by_seq_context.py \
            -v /Users/miramastoras/Desktop/test_files_annotate_vcf/HG04115.polisher_output.h2tg000055l.vcf \
            -f /Users/miramastoras/Desktop/test_files_annotate_vcf/HG04115.mat.h2tg000055l.fa \
            -o /Users/miramastoras/Desktop/test_files_annotate_vcf
'''
import argparse
import pysam
from pysam import VariantFile, FastaFile
import os


def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='annotate_vcf_by_seq_context.py',
        description="Annotate variants in vcf by sequence context (dimer or homopolymer).\
        Repeat annotations will be stored in the vcf INFO column with the format REP_ANN=[HOMOPOLYMER/DIMER]:[repeat pattern]:[repeat length]")

    parser.add_argument("-f", "--fasta",
                        required=True,
                        help="input fasta file")
    parser.add_argument("-o", "--out_dir",
                        required=True,
                        help="directory to save output files")
    parser.add_argument("-v", "--vcf",
                        required=True,
                        help="input vcf file")

    return parser.parse_args()

def extend_dimers(reference_sequence, pattern, start_index):
    '''
    :param reference_sequence: string containing genome sequence
    :param pattern: dimer pattern to search for "at gc ct ag"
    :param start_index: start index in sequence
    :return: end index in sequence for end of repeat pattern
    '''
    end_index = start_index + 2
    while end_index < len(reference_sequence) - 1:
        di_mer = reference_sequence[end_index] + reference_sequence[end_index + 1]
        if di_mer == pattern:
            end_index += 2
        else:
            break

    return end_index

def extend_homopolymers(reference_sequence, pattern, start_index):
    '''
    :param reference_sequence: string containing genome sequence
    :param pattern: dimer pattern to search for "at gc ct ag"
    :param start_index: start index in sequence
    :return: end index in sequence for end of repeat pattern
    '''
    end_index = start_index + 1
    while end_index < len(reference_sequence) :
        mono_mer = reference_sequence[end_index]
        if mono_mer == pattern:
            end_index += 1
        else:
            break

    return end_index

def fetch_repeat_annotation(contig,start,end,FastaFileObject):
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
    hom_base='*'
    i = 0

    while i < len(seq) - 1 :

        # if current base is not the same as next base, check for dimer repeats
        if seq[i] != seq[i + 1]:
            dimer_base = seq[i] + seq[i + 1]
            end_index = extend_dimers(seq, dimer_base, i)
            dimer_length=int((end_index - i)/2)

            if dimer_length > 1 and end_index <= len(seq) - 1 :
                repeat_coords.append([contig,i+start,end_index-1+start,dimer_base,"DIMER"])
                #print(contig,start + i,start + (end_index-1),dimer_base)
                i = end_index - 1 # allows for edge case where one dimer goes into a different dimer or a homopolymer (allows overlaps)
            else:
                # if dimer length is only 2 bp, we want to start next check on second base in the pair
                i = end_index - 1

        # if current base is same as next base, check for homopolymer repeats
        else:
            hom_base = seq[i]
            end_index = extend_homopolymers(seq, hom_base, i)
            hom_length = int((end_index - i)/2)

            if hom_length >= 1 :
                repeat_coords.append([contig, i+start,end_index-1+start, hom_base, "HOMOPOLYMER"])
                #print(contig, start+i, start+(end_index-1), hom_base)
            i = end_index - 1 # allows for edge case of end of homopolymer being start of a dimer

    return repeat_coords

def annotate_variants(vcf_file,fasta_file,vcf_out,bed_out):
    # for each record in vcf, fetch annotations in 200 bp windows on each side
    # this is to make sure we can capture lengths of very long repeats around the variants
    for rec in vcf_file.fetch():
        # remove any SVs
        alternate_allele = rec.alleles[1]

        if len(alternate_allele) > 50:
            continue

        # fetch annotations
        seq_start = rec.start - 200
        seq_end = rec.stop + 200
        annotations = fetch_repeat_annotation(contig=rec.contig, start=seq_start, end=seq_end,
                                   FastaFileObject=fasta_file)
        #print("rec start ",rec.start, " rec stop ", rec.stop)
        variant_annotation = ""

        if len(annotations)==0:
            continue
        else:
            for block in annotations:
                print(block[0],block[1],block[2]+1,block[3],block[4],sep="\t",file=bed_out)
                block_start=block[1]
                block_end=block[2]

                block_len=block_end - block_start +1

                # check if block is inside annotation
                if (rec.start in range(block_start,block_end)) or (rec.stop in range(block_start,block_end)):
                    #print("block inside annotation")
                    #print(block)
                    variant_annotation = variant_annotation + ":" + str(block[4]) + ":" + str(block[3]) + ":" + str(block_len)
                # in the case of large deletion, check if there are homopolymers or dimers inside the deleted segment
                elif (block_start in range(rec.start,rec.stop)) or (block_end in range(rec.start,rec.stop)):
                    #print("annotation inside variant")
                    #print(block)
                    variant_annotation = variant_annotation + ":" + str(block[4]) + ":" + str(block[3]) + ":" + str(block_len)
                # check if block is book-ended by an annotation
                elif abs(rec.start - block_end) == 1 or abs(rec.start - block_start) == 1 or abs(rec.stop - block_end) == 1 or abs(rec.stop - block_start) == 1:
                    variant_annotation = variant_annotation + ":" + str(block[4]) + ":" + str(block[3]) + ":" + str(block_len)
                    #print("block bookended by annotation")
                    #print(block)

        rec.info.__setitem__('REP_ANN',variant_annotation[1:len(variant_annotation)])
        vcf_out.write(rec)

    return


def main():
    # parse command line arguments
    args = arg_parser()

    # read in data files
    small_variant_vcf = VariantFile(args.vcf)
    assembly_fasta_file = FastaFile(args.fasta)

    file_basename=os.path.basename(args.vcf).split(".vcf")[0]
    vcf_out_filename=args.out_dir + "/" + file_basename + ".rep_annotated.vcf"
    bed_out_filename=args.out_dir + "/" + file_basename + ".rep_annotations.bed"

    small_variant_vcf.header.add_meta('INFO', items=[('ID', 'REP_ANN'), ('Number', '.'), ('Type', 'String'),
                                          ('Description', "repeat annotation for homopolymer and dimer context of variant")])

    vcf_out = VariantFile(vcf_out_filename, 'w', header=small_variant_vcf.header)
    bed_out = open(bed_out_filename,'w')

    annotate_variants(small_variant_vcf,assembly_fasta_file,vcf_out,bed_out)

    bed_out.close()

if __name__ == '__main__':
    main()