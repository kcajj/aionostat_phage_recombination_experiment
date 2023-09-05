from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt

def read_msa(path):
    alignment = AlignIO.read(open(path), "fasta")
    print(alignment)
    l=alignment.get_alignment_length()
    msa_matrix=np.zeros([3,l],dtype=str)
    for i,record in enumerate(alignment):
        for pos,nuc in enumerate(record.seq):
            msa_matrix[i][pos]=nuc
    return msa_matrix

def get_evidences_distributions(msa_matrix):
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

    return first_e_distribution, second_e_distribution

#populations=['P2','P3']
#timepoints=['1','3','5','7']
populations=['P2']
timepoints=['5']
for population in populations:
    for timepoint in timepoints:
        file=f'MSAstats/results/msa/{population}/{timepoint}/{population}_{timepoint}_3_msa.fasta'
        #file='/home/giacomocastagnetti/code/rec_genome_analysis/MSAstats/results/msa/clones/P3_C3_msa.fasta'

        msa_matrix=read_msa(file)
        
        first_e_distribution, second_e_distribution = get_evidences_distributions(msa_matrix)

        k=1000
        first_to_plot=np.convolve(first_e_distribution, np.ones(k))/k
        second_to_plot=np.convolve(second_e_distribution, np.ones(k))/k
        x=np.linspace(0,len(first_to_plot),len(first_to_plot))
        plt.title('distribution of evidences')
        plt.plot(x, first_to_plot,label='EM11')
        plt.plot(x, second_to_plot, label='EM60')
        plt.xlabel('bp')
        plt.ylabel('evidence score')
        plt.legend()
        plt.show()