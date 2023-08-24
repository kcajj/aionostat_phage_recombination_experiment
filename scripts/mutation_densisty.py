import pysam
import numpy as np

def analyse_bam(bam_file, l):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        mismatches=np.zeros(l)
        mapping=np.zeros(l,dtype=bool)
        for read in bam.fetch():
            print(read.reference_start)
    return mismatches,mapping

if __name__ == "__main__":
    
    #populations=['P2','P3']
    #isolates=['C1','C2','C3','C4']
    
    populations=['P2']
    isolates=['C1']

    for population in populations:
        for isolate in isolates:
            
            #get the length of the isolate genome

            l=120000
            bam_file_path = f'/home/giacomocastagnetti/code/rec_genome_analysis/results/mappings/{population}/{isolate}.bam'
            mismatch_distribution, mapping = analyse_bam(bam_file_path, l)
            print(mismatch_distribution)
            print(mapping)