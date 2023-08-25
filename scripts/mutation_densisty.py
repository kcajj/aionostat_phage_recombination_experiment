import pysam
from Bio import SeqIO
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

def analyse_bam(bam_file, ref_file):
    ref = SeqIO.read(ref_file, "fasta")
    ref_name=ref.name
    ref_seq=ref.seq
    l=len(ref_seq)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        #mapping and mismatches arrays for each read_name
        mismatches=defaultdict(lambda: np.zeros(l))
        mapping=defaultdict(lambda: np.zeros(l,dtype=bool))
        for read in bam.fetch():
            if read.query_name==ref_name:
                continue
            start=read.reference_start
            end=read.reference_end
            for pos,val in enumerate(mapping):
                if pos>=start and pos<end:
                    mapping[read.query_name][pos]=True
            for (read_pos,reference_pos) in read.get_aligned_pairs():
                if reference_pos!=None and read_pos!=None:
                    if ref_seq[reference_pos]!=read.query_sequence[read_pos]:
                        mismatches[read.query_name][reference_pos]+=1 
    return mismatches,mapping

if __name__ == "__main__":
    
    populations=['P2','P3']
    isolates=['C1','C2','C3','C4']
    k=10000

    for population in populations:
        for isolate in isolates:
            
            #get the length of the isolate genome
            assembly_file=f'/home/giacomocastagnetti/code/rec_genome_analysis/results/assemblies/{population}/{isolate}.fasta'
            bam_file = f'/home/giacomocastagnetti/code/rec_genome_analysis/results/mappings/{population}/{isolate}.bam'

            mismatch_distribution, mapping = analyse_bam(bam_file, assembly_file)

            figure=plt.figure()
            for reference, distribution in mismatch_distribution.items():
                distribution=np.convolve(distribution,np.ones(k),'valid')/k
                l=len(distribution)
                x=np.linspace(0,l,l)
                plt.plot(distribution)
            plt.legend(mismatch_distribution.keys())
            plt.title(f'mutation density between isolate assembly and references (on clone {isolate}), with convolution of {k}')
            plt.ylabel('mutation density')
            plt.xlabel('bp')
            figure.savefig(f'/home/giacomocastagnetti/code/rec_genome_analysis/results/plots/{population}/{isolate}.png')
            plt.close()
            for reference in mismatch_distribution.keys():
                
                ref_file=f'/home/giacomocastagnetti/code/rec_genome_analysis/data/references/{reference}_reference.fa'
                sam_file=f'/home/giacomocastagnetti/code/rec_genome_analysis/results/mappings/references/{reference}.sam'

                mismatch_distribution, mapping = analyse_bam(sam_file, ref_file)

                figure=plt.figure()
                for name, distribution in mismatch_distribution.items():
                    distribution=np.convolve(distribution,np.ones(k),'valid')/k
                    l=len(distribution)
                    x=np.linspace(0,l,l)
                    plt.plot(distribution)
                plt.legend(mismatch_distribution.keys())
                plt.title(f'mutation density between references (on {reference}), with convolution of {k}')
                plt.ylabel('mutation density')
                plt.xlabel('bp')
                figure.savefig(f'/home/giacomocastagnetti/code/rec_genome_analysis/results/plots/references/{reference}.png')
                plt.close()