from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
from handle_msa import read_msa, get_evidences_distributions

'''
    takes a MSA, plots the distribution of evidences (mismatches)
'''

if __name__ == "__main__":

    populations=['P2','P3']
    clones=['C1','C2','C3','C4']
    references=['EM11','EM60']

    for population in populations:
        for clone in clones:

            file=f'results/msa/clones/{population}_{clone}_msa.fasta'
            out_folder=f'results/plots/recombination_evidences/clones/{population}_{clone}_msa.png'

            k=1000 #convolution window

            msa_matrix=read_msa(file)

            first_e_distribution, second_e_distribution = get_evidences_distributions(msa_matrix)
            distribution_evidence_sites=first_e_distribution+second_e_distribution

            first_convoluted=np.convolve(first_e_distribution, np.ones(k), mode='same')
            second_convoluted=np.convolve(second_e_distribution, np.ones(k), mode='same')
            normaliser_convoluted=np.convolve(distribution_evidence_sites, np.ones(k), mode='same')

            ### normalised plot ###
            first_normalised=np.divide(first_convoluted,normaliser_convoluted,out=np.zeros_like(first_convoluted), where=normaliser_convoluted!=0)
            second_normalised=np.divide(second_convoluted,normaliser_convoluted,out=np.zeros_like(second_convoluted), where=normaliser_convoluted!=0)

            x=np.linspace(0,len(first_normalised),len(second_normalised))
            plt.title('distribution of evidences')
            plt.plot(x, first_normalised,label=references[0])
            plt.plot(x, second_normalised, label=references[1])
            plt.xlabel('bp')
            plt.ylabel('evidence score')
            plt.legend()
            plt.savefig(out_folder, bbox_inches='tight')
            plt.close()

            ### non-normalised plot ###
            out_folder=f'results/plots/recombination_evidences/clones/{population}_{clone}_non_normalised_msa.png'
            first_non_normalised=first_convoluted/k
            second_non_normalised=second_convoluted/k

            x=np.linspace(0,len(first_normalised),len(second_normalised))
            plt.title('distribution of evidences')
            plt.plot(x, first_non_normalised,label=references[0])
            plt.plot(x, second_non_normalised, label=references[1])
            plt.xlabel('bp')
            plt.ylabel('evidence score')
            plt.legend()
            plt.savefig(out_folder, bbox_inches='tight')
            plt.close()
