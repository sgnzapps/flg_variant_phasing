#!/usr/bin/env python3

# detect_copies.py

'''
A python script to detect number of repeat copies in FLG exon 3
'''

# Import packages
import os
import sys
from sys import argv
import pysam
from Bio.Seq import Seq
from simplesam import Reader

__version__ = '0.0.5'

def detect_copies(bam):
    # Set default variables
    min_size = 800  # Minimum insertion size
    start = 152303959  # Start coordinate for insertion region
    stop = 152307308  # End coordinate for insertion region
    base_count = 10  # Number of repeats in reference contig
    allele_freq_threshold = 0.20 # Minimum frequency for allele to be considered
    rep10_rep8_separator = 152304650
    # Get path of input file
    path_out = os.path.abspath(bam).strip('.bam')
    read_reps, fasta_dict, qual_dict, read_rep_coords = parse_sam(bam, min_size, start, stop)
    ordered_counts, new_read_reps = calculate_counts(read_reps, read_rep_coords, rep10_rep8_separator, base_count)
    top_alleles = define_main_alleles(ordered_counts, allele_freq_threshold)
    top_alleles_reads = write_data(new_read_reps, top_alleles, fasta_dict, path_out, qual_dict)
    write_summary(top_alleles, ordered_counts, path_out)
    subset_bam(top_alleles_reads, bam, path_out)

def subset_bam(read_list, bam, path_out):
    infile = pysam.AlignmentFile(bam, "r")
    outfile = pysam.AlignmentFile(path_out + '.filt.bam', "wb", template=infile)
    for seg in infile:
        if seg.query_name in read_list:
            outfile.write(seg)
    outfile.close()
    infile.close()

def write_data(read_reps, top_alleles, fasta_dict, path_out, qual_dict):
    top_alleles_reads = {}
    for k, j in top_alleles:
        out = open(path_out + '.' + str(k[0:2]) + 'rep.fastq', 'w')
        for i in read_reps:
            if read_reps[i] == k:
                out.write('@' + i + '\n')
                out.write(fasta_dict[i] + '\n')
                out.write('+\n')
                out.write(qual_dict[i] + '\n')
                top_alleles_reads[i] = 1
        out.close()
    out = open(path_out + '.index.tsv', 'w')
    for i in read_reps:
        out.write(i + '\t' + str(read_reps[i]) + '\n')
    out.close()
    return top_alleles_reads

def write_summary(top_alleles, ordered_counts, path_out):
    out = open(path_out + '.summmary.tsv', 'w')
    out.write('CNV_alleles\tNo_of_reads\tSelected_alleles\n')
    for i in ordered_counts:
        out.write(i[0] + '\t' + str(i[1]) + '\t')
        if i in top_alleles:
            out.write('yes\n')
        else:
            out.write('no\n')
    out.close()

def define_main_alleles(ordered_counts, allele_freq_threshold):
    first = ordered_counts[0][1]
    if len(ordered_counts) > 1:
        second = ordered_counts[1][1]
        if second/(first + second) >= allele_freq_threshold:  # If second allele is at least x%
            return ordered_counts[:2]  # Print both alleles
        else:
            return [ordered_counts[0]]  # Only print top allele
    else:
        return [ordered_counts[0]]  # Only print top allele

def calculate_counts(read_reps, read_rep_coords, rep10_rep8_separator, base_count):
    counts_dict = {}
    new_read_reps = {}
    for i in read_reps:
        reps = add_pos_info(read_reps[i], read_rep_coords[i], rep10_rep8_separator, base_count)
        new_read_reps[i] = reps
        if reps not in counts_dict:
            counts_dict[reps] = 1
        else:
            counts_dict[reps] += 1
    ordered_counts = sorted(counts_dict.items(), key=lambda item: item[1], reverse=True)
    return ordered_counts, new_read_reps

def add_pos_info(ins_count, coord, rep10_rep8_separator, base_count):
    if ins_count == 1:
        if coord[0] <= rep10_rep8_separator:
            return str(ins_count + base_count) + '-rpt10'
        elif coord[0] > rep10_rep8_separator:
            return str(ins_count + base_count) + '-rpt8'
    else:
        return str(ins_count + base_count)

def parse_sam(bam, min_size, start, stop):
    read_reps = {}
    read_rep_coords = {}
    sam = pysam.AlignmentFile(bam, 'rb')
    for seg in sam:
        rcoord = seg.reference_start
        cigar = seg.cigartuples
        if cigar:
            read_reps[seg.query_name], read_rep_coords[seg.query_name] = parse_cigar(cigar, rcoord, min_size, start, stop)
        else:
            raise Exception("Error: %s has no cigar" % seg.query_name)
    fasta_dict, qual_dict = get_fastq(bam)
    return read_reps, fasta_dict, qual_dict, read_rep_coords

def parse_cigar(cigar, rcoord, min_size, start, stop):
    counts = 0
    rcoord_list = []
    for i in cigar:
        if i[0] == 0 or i[0] == 2:  # If Match or Deletion
            rcoord = rcoord + i[1]  # Update reference coordinate
        if start < rcoord < stop:  # If reference coordinate within expected insertion region
            if i[0] == 1 and i[1] > min_size:  # If insertion and larger than min_size
                counts += 1
                rcoord_list.append(rcoord)
    return counts, rcoord_list

def get_fastq(bam):
    fasta_dict = {}
    qual_dict = {}
    with Reader(open(bam)) as sam:
        for read in sam:
            if read.flag == 0:
                fasta_dict[read.qname] = read.seq
                qual_dict[read.qname] = read.qual
            elif read.flag == 16:
                seq = Seq(read.seq)
                fasta_dict[read.qname] = str(seq.reverse_complement())
                qual_dict[read.qname] = read.qual[::-1]
    return fasta_dict, qual_dict

if __name__ == '__main__':
    if len(argv) != 2:
        sys.exit("Input Error - usage: python detect_copies.py data.bam")

    bam = argv[1]
    detect_copies(bam)
