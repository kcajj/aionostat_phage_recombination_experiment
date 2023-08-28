import pysam
from Bio import SeqIO
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

def decode_flag(flag):
    if flag==0 or flag==16: return 'primary'
    if flag==2048 or flag==2064: return 'supplementary'
    if flag==256 or flag==272: return 'secondary'

def analyse_bam(bam_file, ref_file):
    ref = SeqIO.read(ref_file, "fasta")
    ref_name=ref.name
    ref_seq=ref.seq
    l=len(ref_seq)
    with pysam.AlignmentFile(bam_file, "rb") as bam:

        #mapping and mismatches arrays for each read_name
        mismatches=defaultdict(lambda: np.zeros(l))
        mapping=defaultdict(lambda: {'primary':np.zeros(l,dtype=bool),'supplementary':np.zeros(l,dtype=bool),'secondary':np.zeros(l,dtype=bool)})
        coverage=defaultdict(lambda: np.zeros(l))

        for read in bam.fetch():

            # do not consider the data coming from the reference mapped on itself
            if read.query_name==ref_name:
                continue

            # mapping array
            start=read.reference_start
            end=read.reference_end
            for type,map in mapping[read.query_name].items():
                for pos,value in enumerate(map):
                    if pos>=start:
                        if pos>=end: break
                        alignment_type = decode_flag(read.flag)
                        mapping[read.query_name][alignment_type][pos]=True
            
            # mismatches array
            for (read_pos,reference_pos) in read.get_aligned_pairs():
                if reference_pos!=None and read_pos!=None:
                    if ref_seq[reference_pos]!=read.query_sequence[read_pos]:
                        mismatches[read.query_name][reference_pos]+=1

                    #coverage array
                    coverage[read.query_name][reference_pos]+=1
        
        #normalise for coverage
        for ref,array in mismatches.items():
            for pos, val in enumerate(array):
                if val>0:
                    mismatches[ref][pos]=val/coverage[ref][pos]

    return mismatches,mapping