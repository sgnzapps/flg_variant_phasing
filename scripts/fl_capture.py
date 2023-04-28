#!/usr/bin/env python3

# fl_capture.py

'''
A python script to capture reads that span across a specified genomic region
'''

# Import packages
import os
import sys
from sys import argv
import pysam
import random

__version__ = '0.0.4'

# Main function to capture full-length reads
def fl_capture(bed, bam, loc, lcov=0.98, depth=500, f_strandness=0.5, max_sclip_size=300, seed=7):
    pass_reads = []
    forward = 0
    reverse = 0
    f_max = depth * f_strandness
    r_max = depth - f_max
    f = open(bed).read().splitlines()
    random.seed(seed)
    random.shuffle(f)
    sam = pysam.AlignmentFile(bam, 'rb')
    cigar_sclip_dict = soft_clip_check(sam, max_sclip_size)
    for i in f:
        strand = i.split('\t')[5].strip()
        if strand == '+':
            if forward < f_max:
                rname = i.split('\t')[3].strip()
                seg = (
                    i.split('\t')[0].strip(), 
                    int(i.split('\t')[1].strip()), 
                    int(i.split('\t')[2].strip())
                    )
                if seg_within(seg, loc, lcov):
                    if cigar_sclip_dict[rname]:
                        pass_reads.append(rname)
                        forward += 1
        elif strand == '-':
            if reverse < r_max:
                rname = i.split('\t')[3].strip()
                seg = (
                    i.split('\t')[0].strip(), 
                    int(i.split('\t')[1].strip()), 
                    int(i.split('\t')[2].strip())
                    )
                if seg_within(seg, loc, lcov):
                    if cigar_sclip_dict[rname]:
                        pass_reads.append(rname)
                        reverse += 1
        else:
            raise Exception('Error: Unknown strandess %s in BED file.' % strand)
        if forward + reverse == depth:
            break
    path_out = os.path.abspath(bam).strip('.bam')
    read_extract(bam, pass_reads, depth, path_out)
    print('Captured %i full-length reads' % len(pass_reads))

# Evaluate if segment is spanning within region
def seg_within(seg, loc, lcov):
    ref_c, ref_s, ref_e = loc
    seg_c, seg_s, seg_e = seg
    if (
        seg_c == ref_c and 
        seg_s >= ref_s and 
        seg_e <= ref_e and
        (seg_e-seg_s)/(ref_e-ref_s) >= lcov
        ):
        return True

# Evaluate if CIGAR has large soft clipping
def soft_clip_check(sam, sclip_size):
    cigar_dict = {}
    for seg in sam:
        rname = seg.query_name
        cigar = seg.cigartuples
        if cigar[0][0] == 4 and cigar[0][1] > sclip_size:  # if start soft-clip more than threshold
            cigar_dict[rname] = False
        elif cigar[-1][0] == 4 and cigar[-1][1] > sclip_size:  # if end soft-clip more than threshold
            cigar_dict[rname] = False
        else:
            cigar_dict[rname] = True
    return cigar_dict

# Extract reads from BAM and save into a new BAM
def read_extract(bam, read_names, depth, path_out):
    sam = pysam.AlignmentFile(bam, 'rb')
    sam_index = pysam.IndexedReads(sam)
    sam_index.build()
    header = sam.header.copy()
    out = pysam.Samfile(path_out + '.fl' + str(depth) + 'sr.unsort.bam', 'wb', header=header)
    for n in read_names:
        try:
            sam_index.find(n)
        except KeyError:
            raise Exception('Error: Read %s not found in BAM file.' % n)
        else:
            iter = sam_index.find(n)
            for i in iter:
                out.write(i)
    out.close()
    sam.close()
    pysam.sort("-o", path_out + '.fl' + str(depth) + 'sr.bam', path_out + '.fl' + str(depth) + 'sr.unsort.bam')
    os.remove(path_out + '.fl' + str(depth) + 'sr.unsort.bam')

if __name__ == '__main__':
    if len(argv) != 6:
        sys.exit("Input Error - usage: python fl_capture.py data.bed data.bam region(chrN:10000-20000) length-wise-coverage(0.98) depth-of-coverage(500)")

    bed = argv[1]
    bam = argv[2]
    loc = argv[3]
    loc = (
        loc.split(':')[0], 
        int(loc.split(':')[1].split('-')[0]), 
        int(loc.split(':')[1].split('-')[1])
        )
    lcov = float(argv[4])
    depth = int(argv[5])
    fl_capture(bed, bam, loc, lcov, depth)
