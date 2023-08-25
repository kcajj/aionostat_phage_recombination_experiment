import pysam
from Bio import SeqIO
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

def analyse_bam(bam_file, ref, l):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        #mapping and mismatches arrays for each read_name
        mismatches=defaultdict(lambda: np.zeros(l))
        mapping=defaultdict(lambda: np.zeros(l,dtype=bool))
        for read in bam.fetch():
            start=read.reference_start
            end=read.reference_end
            for pos,val in enumerate(mapping):
                if pos>=start and pos<end:
                    mapping[read.query_name][pos]=True
            for (read_pos,reference_pos) in read.get_aligned_pairs():
                if reference_pos!=None and read_pos!=None:
                    if ref[reference_pos]!=read.query_sequence[read_pos]:
                        mismatches[read.query_name][reference_pos]+=1 
    return mismatches,mapping

if __name__ == "__main__":
    
    #populations=['P2','P3']
    #isolates=['C1','C2','C3','C4']
    
    populations=['P2']
    isolates=['C1']

    for population in populations:
        for isolate in isolates:
            
            #get the length of the isolate genome
            assembly_file=f'/home/giacomocastagnetti/code/rec_genome_analysis/results/assemblies/{population}/{isolate}.fasta'
            ref = SeqIO.read(assembly_file, "fasta").seq
            l=len(ref)

            bam_file = f'/home/giacomocastagnetti/code/rec_genome_analysis/results/mappings/{population}/{isolate}.bam'
            mismatch_distribution, mapping = analyse_bam(bam_file, ref, l)
            #np.set_printoptions(threshold=99999999)
            #print(mismatch_distribution)

            k=10000

            for reference, distribution in mismatch_distribution.items():
                distribution=np.convolve(distribution,np.ones(k),'valid')/k
                l=len(distribution)
                x=np.linspace(0,l,l)
                plt.plot(distribution)
            plt.legend(mismatch_distribution.keys())
            plt.title('mutation density between assembly and references')
            plt.ylabel('mutation density')
            plt.xlabel('bp')
            plt.show()

            references=[]
            for reference in mismatch_distribution.keys():
                references.append(reference)

            for reference in references:
                ref_file=f'data/references/'
                ref = SeqIO.read(ref_file, "fasta").seq
                l=len(ref)

                sam_file=f'/home/giacomocastagnetti/code/rec_genome_analysis/results/mappings/references/{reference}.sam'

            #plot the alignments between references