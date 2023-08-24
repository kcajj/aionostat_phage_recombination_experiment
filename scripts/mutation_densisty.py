import pysam
from Bio import SeqIO
import numpy as np

def analyse_bam(bam_file, l):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        mismatches=np.zeros(l)
        mapping=np.zeros(l,dtype=bool)
        for read in bam.fetch():
            print(read.reference_start)
            for (read_pos,reference_pos) in read.get_aligned_pairs():
                if reference_pos!=None and read_pos!=None:
                    print(read_pos,reference_pos)
                    print(read.query_sequence[read_pos])
                    #if ref_seq==read.query_sequence[read_pos]:

                    #add mismatch count
                    #add mapping positions
                        
    return mismatches,mapping

if __name__ == "__main__":
    
    #populations=['P2','P3']
    #isolates=['C1','C2','C3','C4']
    
    populations=['P2']
    isolates=['C1']

    for population in populations:
        for isolate in isolates:
            
            #get the length of the isolate genome
            assembly=f'/home/giacomocastagnetti/code/rec_genome_analysis/results/assemblies/{population}/{isolate}.fasta'
            seq = SeqIO.read(assembly, "fasta").seq
            l=len(seq)

            bam_file_path = f'/home/giacomocastagnetti/code/rec_genome_analysis/results/mappings/{population}/{isolate}.bam'
            mismatch_distribution, mapping = analyse_bam(bam_file_path, l)
            print(mismatch_distribution)
            print(mapping)