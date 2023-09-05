from Bio import SeqIO
import gzip
import pandas as pd

def get_longest_seq(path, longest_reads, n):
    l=[]
    c=0
    for row in longest_reads.iterrows():
        with gzip.open(path, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                if record.id==row[1].read_name:
                    if bool(row[1].is_reverse_primary): #take the flag information from the bam file
                        l.append([record.id, str(record.seq.reverse_complement())])
                    else:
                        l.append([record.id, str(record.seq)])
                    break
        c+=1
        if c==n: break
        
    return l

populations=['P2','P3']
timepoints=['1','3','5','7']
for population in populations:
    for timepoint in timepoints:
        file=f'MSAstats/data/population_reads/{population}_{timepoint}.fastq.gz'
        n=2

        longest_reads=pd.read_csv(f'/home/giacomocastagnetti/code/rec_genome_analysis/chimeric_reads/longest_matching_reads/{population}/{population}_{timepoint}.csv')

        n=5
        longest=get_longest_seq(file, longest_reads, n)

        ref1_file='MSAstats/data/references/EM11_assembly.fasta'
        ref2_file='MSAstats/data/references/EM60_assembly.fasta'

        for i_r,read in enumerate(longest):
            out_file=open(f'MSAstats/results/seq_for_msa/{population}/{timepoint}/{population}_{timepoint}_{i_r}.fasta','w')
            out_file.write('>'+read[0]+'\n'+read[1]+'\n')
            ref1=SeqIO.read(ref1_file, 'fasta')
            out_file.write('>'+ref1.id+'\n'+str(ref1.seq)+'\n')
            ref2=SeqIO.read(ref2_file, 'fasta')
            out_file.write('>'+ref2.id+'\n'+str(ref2.seq)+'\n')