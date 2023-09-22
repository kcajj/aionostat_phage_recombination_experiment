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
        if nuc_assembly!='-' and nuc_first_ref!='-' and nuc_second_ref!='-':
            if (nuc_assembly!=nuc_first_ref and nuc_assembly!=nuc_second_ref) or (nuc_assembly==nuc_first_ref and nuc_assembly==nuc_second_ref):
                continue
            elif nuc_assembly==nuc_first_ref and nuc_assembly!=nuc_second_ref:
                first_e_distribution[pos]=1
            elif nuc_assembly!=nuc_first_ref and nuc_assembly==nuc_second_ref:
                second_e_distribution[pos]=1

    return first_e_distribution, second_e_distribution

populations=['P2','P3']
clones=['C1','C2','C3','C4']
references=['EM11','EM60']

for population in populations:
    for clone in clones:
        file=f'MSAstats/results/msa/clones/{population}_{clone}_msa.fasta'
        out_folder=f'results/plots/recombination_evidences/clones/{population}_{clone}_msa.png'

        msa_matrix=read_msa(file)

        first_e_distribution, second_e_distribution = get_evidences_distributions(msa_matrix)
        distribution_evidence_sites=first_e_distribution+second_e_distribution

        k=10
        first_convoluted=np.convolve(first_e_distribution, np.ones(k), mode='same')
        second_convoluted=np.convolve(second_e_distribution, np.ones(k), mode='same')
        normaliser_convoluted=np.convolve(distribution_evidence_sites, np.ones(k), mode='same')

        first_normalised=np.divide(first_convoluted,normaliser_convoluted,out=np.zeros_like(first_convoluted), where=normaliser_convoluted!=0)
        second_normalised=np.divide(second_convoluted,normaliser_convoluted,out=np.zeros_like(second_convoluted), where=normaliser_convoluted!=0)

        x=np.linspace(0,len(first_normalised),len(second_normalised))
        plt.title('distribution of evidences')
        plt.plot(x, first_normalised,label=references[0])
        plt.plot(x, second_normalised, label=references[1])
        plt.xlabel('bp')
        plt.ylabel('evidence score')
        plt.xlim([28100,28500])
        plt.legend()
        plt.savefig(out_folder, bbox_inches='tight')
        plt.close()

        out_folder=f'results/plots/recombination_evidences/clones/{population}_{clone}_non_normalised_msa.png'
        first_non_normalised=first_convoluted/k
        second_non_normalised=second_convoluted/k

        x=np.linspace(0,len(first_normalised),len(second_normalised))
        plt.title('distribution of evidences')
        plt.plot(x, first_non_normalised,label=references[0])
        plt.plot(x, second_non_normalised, label=references[1])
        plt.xlabel('bp')
        plt.ylabel('evidence score')
        plt.xlim([28100,28500])
        plt.legend()
        plt.savefig(out_folder, bbox_inches='tight')
        plt.close()
