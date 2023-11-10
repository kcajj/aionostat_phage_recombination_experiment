from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
from handle_msa import read_msa, get_evidences_distributions

'''
takes a lot of MSA, plots the total distribution of evidences (mismatches) summing all MSAs
'''

if __name__ == "__main__":

    #populations=['P2','P3']
    #timepoints=['1','3','5','7']
    #reads=['0','1','2','3','4']
    populations=['P2','P3']
    timepoints=['1','3','5','7']

    n=1000 #number of reads
    reads=[]
    for i in range(n):
        reads.append(str(i))

    references=[]
    references_msa='results/msa/refs_msa.fasta'
    references_alignment = AlignIO.read(open(references_msa), "fasta")
    l=references_alignment.get_alignment_length()
    for i,record in enumerate(references_alignment):
        references.append(record.id)

    for population in populations:
        for timepoint in timepoints:

            out_folder=f'results/plots/recombination_evidences/reads/{population}/{population}_{timepoint}.png'

            fig=plt.figure()
            total_first_evidences=np.zeros(l)
            total_second_evidences=np.zeros(l)
            total_evidences=np.zeros(l)

            for read in reads:
                file=f'results/msa/{population}/{timepoint}/{population}_{timepoint}_{read}_msa.fasta'
                
                msa_matrix=read_msa(file)
                
                first_e_distribution, second_e_distribution = get_evidences_distributions(msa_matrix,l)

                total_first_evidences=total_first_evidences+first_e_distribution
                total_second_evidences=total_second_evidences+second_e_distribution
                total_evidences=total_evidences+first_e_distribution+second_e_distribution

            print(f'processed {n} reads from {population} {timepoint}')
            
            k=1000 #convolution window

            ## normalised plot ##

            first_convoluted=np.convolve(total_first_evidences, np.ones(k), mode='same')
            second_convoluted=np.convolve(total_second_evidences, np.ones(k), mode='same')
            total_convoluted=np.convolve(total_evidences, np.ones(k), mode='same')
            
            first_normalised=np.divide(first_convoluted,total_convoluted,out=np.zeros_like(first_convoluted), where=total_convoluted!=0)
            second_normalised=np.divide(second_convoluted,total_convoluted,out=np.zeros_like(second_convoluted), where=total_convoluted!=0)

            x=np.linspace(0,len(first_normalised),len(first_normalised))
            plt.title('distribution of evidences')
            plt.plot(x, first_normalised, label=references[0])
            plt.plot(x, second_normalised, label=references[1])
            plt.xlabel('bp')
            plt.ylabel('evidence score')
            plt.legend()
            fig.savefig(out_folder, bbox_inches='tight')
            plt.close(fig)
            
            ## non-normalised plot ##

            out_folder=f'results/plots/recombination_evidences/reads/{population}/{population}_{timepoint}_non_normalised.png'

            fig=plt.figure()

            first_to_plot=np.convolve(total_first_evidences, np.ones(k), mode='same')/k
            second_to_plot=np.convolve(total_second_evidences, np.ones(k), mode='same')/k

            x=np.linspace(0,len(first_to_plot),len(first_to_plot))
            plt.title('distribution of evidences')
            plt.plot(x, first_to_plot,label=references[0])
            plt.plot(x, second_to_plot, label=references[1])
            plt.xlabel('bp')
            plt.ylabel('evidence score')
            plt.legend()
            fig.savefig(out_folder, bbox_inches='tight')
            plt.close(fig)