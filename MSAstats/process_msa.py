from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt

def read_msa(path):
    alignment = AlignIO.read(open(path), "fasta")
    l=alignment.get_alignment_length()
    msa_matrix=np.zeros([3,l],dtype=str)
    for i,record in enumerate(alignment):
        for pos,nuc in enumerate(record.seq):
            msa_matrix[i][pos]=nuc
    return msa_matrix

def get_evidences_distributions(msa_matrix,l):
    first_e_distribution=np.zeros(l)
    second_e_distribution=np.zeros(l)
    all_changes=np.zeros(l)

    for pos,array in enumerate(msa_matrix[0]):
        nuc_first_ref=msa_matrix[0,pos]
        nuc_second_ref=msa_matrix[1,pos]
        nuc_assembly=msa_matrix[2,pos]
        #print(nuc_assembly, nuc_first_ref, nuc_second_ref)
        if nuc_assembly!='-' and nuc_first_ref!='-' and nuc_second_ref!='-':
            if (nuc_assembly!=nuc_first_ref and nuc_assembly!=nuc_second_ref) or (nuc_assembly==nuc_first_ref and nuc_assembly==nuc_second_ref):
                continue
            elif nuc_assembly==nuc_first_ref and nuc_assembly!=nuc_second_ref:
                first_e_distribution[pos]=1
                all_changes[pos]=1
            elif nuc_assembly!=nuc_first_ref and nuc_assembly==nuc_second_ref:
                second_e_distribution[pos]=1
                all_changes[pos]=1

    return first_e_distribution, second_e_distribution, all_changes

#populations=['P2','P3']
#timepoints=['1','3','5','7']
#reads=['0','1','2','3','4']
populations=['P2','P3']
timepoints=['1','3','5','7']

reads=[]
for i in range(100):
    reads.append(str(i))

references=[]
references_msa='MSAstats/results/msa/refs_msa.fasta'
references_alignment = AlignIO.read(open(references_msa), "fasta")
l=references_alignment.get_alignment_length()
for i,record in enumerate(references_alignment):
    references.append(record.id)

for population in populations:
    for timepoint in timepoints:
        
        print(f'processing 100 reads from {population} {timepoint}')

        out_folder=f'results/plots/recombination_evidences/reads/{population}/{population}_{timepoint}.png'
        #out_folder=f'results/plots/recombination_evidences/clones/C4.png'

        fig=plt.figure()
        total_first_evidences=np.zeros(l)
        total_second_evidences=np.zeros(l)
        total_evidences=np.zeros(l)
        total_changes=np.zeros(l)

        for read in reads:
            file=f'MSAstats/results/msa/{population}/{timepoint}/{population}_{timepoint}_{read}_msa.fasta'
            #file='/home/giacomocastagnetti/code/rec_genome_analysis/MSAstats/results/msa/clones/P2_C4_msa.fasta'
            
            msa_matrix=read_msa(file)
            
            first_e_distribution, second_e_distribution, all_changes = get_evidences_distributions(msa_matrix,l)

            total_first_evidences=total_first_evidences+first_e_distribution
            total_second_evidences=total_second_evidences+second_e_distribution
            total_evidences=total_evidences+first_e_distribution+second_e_distribution
            total_changes=total_changes+all_changes

        k=1000
        first_convoluted=np.convolve(total_first_evidences, np.ones(k), mode='same')
        second_convoluted=np.convolve(total_second_evidences, np.ones(k), mode='same')
        total_convoluted=np.convolve(total_evidences, np.ones(k), mode='same')
        normaliser_convoluted=np.convolve(total_changes, np.ones(k), mode='same')
        
        first_normalised=np.divide(first_convoluted,normaliser_convoluted,out=np.zeros_like(first_convoluted), where=normaliser_convoluted!=0)
        second_normalised=np.divide(second_convoluted,normaliser_convoluted,out=np.zeros_like(second_convoluted), where=normaliser_convoluted!=0)
        total_normalised=np.divide(total_convoluted,normaliser_convoluted,out=np.zeros_like(total_convoluted), where=normaliser_convoluted!=0)

        first_to_plot=np.divide(first_normalised,total_normalised,out=np.zeros_like(first_normalised), where=total_normalised!=0)
        second_to_plot=np.divide(second_normalised,total_normalised,out=np.zeros_like(second_normalised), where=total_normalised!=0)

        x=np.linspace(0,len(first_to_plot),len(first_to_plot))
        plt.title('distribution of evidences')
        plt.plot(x, first_to_plot,label=references[0])
        plt.plot(x, second_to_plot, label=references[1])
        plt.xlabel('bp')
        plt.ylabel('evidence score')
        plt.legend()
        fig.savefig(out_folder, bbox_inches='tight')
        plt.close(fig)