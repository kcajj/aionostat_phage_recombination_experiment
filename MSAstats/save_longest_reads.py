from Bio import SeqIO
import gzip
import pandas as pd

def get_longest_seq(path, longest_reads, n):
    reads={}
    for i_read,row in enumerate(longest_reads.iterrows()):
        reads[row[1].read_name]=[i_read,row[1].is_reverse_primary]
        if i_read==n-1: break

    output=[0] * n

    c=0
    with gzip.open(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if record.id in reads.keys():
                c+=1
                if bool(reads[record.id][1]): #take the flag information from the bam file
                    output[reads[record.id][0]]=[record.id, str(record.seq.reverse_complement())]
                else:
                    output[reads[record.id][0]]=[record.id, str(record.seq)]
                if c==n: break #to speed up the code
    return output

populations=['P2','P3']
timepoints=['1','3','5','7']
for population in populations:
    for timepoint in timepoints:
        file=f'MSAstats/data/population_reads/{population}_{timepoint}.fastq.gz'

        longest_reads=pd.read_csv(f'/home/giacomocastagnetti/code/rec_genome_analysis/chimeric_reads/longest_matching_reads/{population}/{population}_{timepoint}.csv')

        n=1000
        longest=get_longest_seq(file, longest_reads, n)

        ref1_file='MSAstats/data/references/EM11_assembly.fasta'
        ref2_file='MSAstats/data/references/EM60_assembly.fasta'

        for i_r,read in enumerate(longest):
            out_file=open(f'MSAstats/results/seq_for_msa/{population}/{timepoint}/{population}_{timepoint}_{i_r}.fasta','w')
            out_file.write('>'+read[0]+'\n'+read[1]+'\n')
        
        print(f'saved reads for {population}, {timepoint}')
