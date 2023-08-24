import pysam

def analyse_bam(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        print(bam.get_reference_length())
        for read in bam.fetch():
            print(read)
    return 0

if __name__ == "__main__":
    
    #populations=['P2','P3']
    #isolates=['C1','C2','C3','C4']
    
    populations=['P2']
    isolates=['C1']

    for population in populations:
        for isolate in isolates:
            
            bam_file_path = f'results/mappings/{population}/new_chemistry/{isolate}.bam'
            mismatch_distribution, mapping = analyse_bam(bam_file_path)