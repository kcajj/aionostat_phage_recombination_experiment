from Bio import SeqIO
import gzip
import pandas as pd

'''
takes the list of longest reads, searches them in the fastq file and stores their
sequences in individual files
'''
def get_longest_seq(path, longest_reads, n):
    '''
    reads the csv file to take the names of the longest reads.
    searches the read names in the fastq file
    '''
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
                if c==n:
                    break #to speed up the code

    return output

populations=['P2','P3']
timepoints=['1','3','5','7']
for population in populations:
    for timepoint in timepoints:
        file=f'data/population_reads/{population}_{timepoint}.fastq.gz'#takes file with all reads
        longest_reads=pd.read_csv(f'results/longest_matching_reads/{population}/{population}_{timepoint}.csv')#takes names of longest reads
        out_fasta=f'results/seq_for_msa/{population}/{timepoint}/{population}_{timepoint}'

        n=1000#number of reads to save
        longest=get_longest_seq(file, longest_reads, n)

        for i_r,read in enumerate(longest):
            out_file=open(f'{out_fasta}_{i_r}.fasta','w')
            out_file.write('>'+read[0]+'\n'+read[1]+'\n')
        
        print(f'saved {n} reads for {population}, {timepoint}')