from Bio import SeqIO
import gzip

def select_for_longest(path,n):
    l=[]
    with gzip.open(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            l.append([record.id, str(record.seq), len(record)])
    l.sort(key=lambda x: x[2])
    return l[-n:]

#populations=['P2','P3']
#timepoints=['1','3','5','7']
populations=['P2']
timepoints=['5','7']
for population in populations:
    for timepoint in timepoints:
        file=f'MSAstats/data/population_reads/{population}_{timepoint}.fastq.gz'
        n=2
        longest=select_for_longest(file, 10)
        ref1_file='MSAstats/data/references/EM11_assembly.fasta'
        ref2_file='MSAstats/data/references/EM60_assembly.fasta'
        for i_r,read in enumerate(longest):
            out_file=open(f'MSAstats/results/reads/{population}/{timepoint}/{population}_{timepoint}_{i_r}.fasta','w')
            out_file.write('>'+read[0]+'\n'+read[1]+'\n')
            ref1=SeqIO.read(ref1_file, 'fasta')
            out_file.write('>'+ref1.id+'\n'+str(ref1.seq)+'\n')
            ref2=SeqIO.read(ref2_file, 'fasta')
            out_file.write('>'+ref2.id+'\n'+str(ref2.seq)+'\n')