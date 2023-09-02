from Bio import AlignIO
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import gzip

def select_for_longest(path,n):
    l=[]
    with gzip.open(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            l.append([record.id, str(record.seq), len(record)])
    l.sort(key=lambda x: x[2])
    return l[-n:]

def read_msa(path):
    alignment = AlignIO.read(open(path), "fasta")
    print(alignment)
    l=alignment.get_alignment_length()
    msa_matrix=np.zeros([3,l],dtype=str)
    for i,record in enumerate(alignment):
        for pos,nuc in enumerate(record.seq):
            msa_matrix[i][pos]=nuc
    return msa_matrix

#populations=['P2','P3']
#timepoints=['1','3','5','7']
populations=['P2']
timepoints=['1']
for population in populations:
    for timepoint in timepoints:
        file=f'MSAstats/data/population_reads/{population}_{timepoint}.fastq.gz'
        longest=select_for_longest(file, 1)
        ref1_file='MSAstats/data/references/EM11_assembly.fasta'
        ref2_file='MSAstats/data/references/EM60_assembly.fasta'
        for i_r,read in enumerate(longest):
            out_file=open(f'MSAstats/results/{population}_{timepoint}_{i_r}','w')
            out_file.write('>'+read[0]+'\n'+read[1]+'\n')
            ref1=SeqIO.read(ref1_file, 'fasta')
            out_file.write('>'+ref1.id+'\n'+str(ref1.seq)+'\n')
            ref2=SeqIO.read(ref2_file, 'fasta')
            out_file.write('>'+ref2.id+'\n'+str(ref2.seq)+'\n')

        '''
        msa_matrix=read_msa(file)
        l=len(msa_matrix[0])
        first_e_distribution=np.zeros(l)
        second_e_distribution=np.zeros(l)
        for pos,array in enumerate(msa_matrix[0]):
            nuc_assembly=msa_matrix[0,pos]
            nuc_first_ref=msa_matrix[1,pos]
            nuc_second_ref=msa_matrix[2,pos]
            #print(nuc_assembly, nuc_first_ref, nuc_second_ref)
            if nuc_assembly!='-' and nuc_first_ref!='-' and nuc_second_ref!='-':
                if (nuc_assembly!=nuc_first_ref and nuc_assembly!=nuc_second_ref) or (nuc_assembly==nuc_first_ref and nuc_assembly==nuc_second_ref):
                    continue
                elif nuc_assembly==nuc_first_ref and nuc_assembly!=nuc_second_ref:
                    first_e_distribution[pos]=1
                elif nuc_assembly!=nuc_first_ref and nuc_assembly==nuc_second_ref:
                    second_e_distribution[pos]=1

        k=1000
        first_to_plot=np.convolve(first_e_distribution, np.ones(k))/k
        second_to_plot=np.convolve(second_e_distribution, np.ones(k))/k
        x=np.linspace(0,len(first_to_plot),len(first_to_plot))
        plt.plot(x, first_to_plot)
        plt.plot(x, second_to_plot)
        plt.show()
        '''