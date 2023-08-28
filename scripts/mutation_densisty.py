import numpy as np
import matplotlib.pyplot as plt

from analyse_bam import analyse_bam

if __name__ == "__main__":
    
    populations=['P2','P3']
    isolates=['C1','C2','C3','C4']
    k=1000

    for population in populations:
        for isolate in isolates:
            
            assembly_file=f'/home/giacomocastagnetti/code/rec_genome_analysis/results/assemblies/{population}/{isolate}.fasta'
            bam_file = f'/home/giacomocastagnetti/code/rec_genome_analysis/results/mappings/{population}/{isolate}.bam'

            mismatch_distribution, mapping = analyse_bam(bam_file, assembly_file)

            fig, axs =plt.subplots(2,sharex=True,constrained_layout = True)
            for reference, distribution in mismatch_distribution.items():
                distribution=np.convolve(distribution,np.ones(k),'valid')/k
                l=len(distribution)
                x=np.linspace(0,l,l)
                axs[0].plot(distribution)
            axs[0].legend(mismatch_distribution.keys())
            fig.suptitle('population: '+population)
            axs[0].set_title(f'mutation density distribution between {isolate} and references, with convolution window of {k}')
            axs[0].set_ylabel('mutation density')
            axs[0].set_xlabel('bp')

            for reference, all_alignments in mapping.items():
                for alignment_type, map in all_alignments.items():
                    y=[]
                    x=[]
                    for pos, val in enumerate(map):
                        if val==True:
                            x.append(pos)
                            y.append(reference+' '+alignment_type)
                    colors={'primary':'red','supplementary':'purple','secondary':'pink'}
                    axs[1].scatter(x,y,
                                   s=5,
                                   alpha=0.2,
                                   c=colors[alignment_type],
                                   label=alignment_type)

            #axs[1].legend()
            axs[1].set_title('mapping information of the references on the assembly')
            fig.savefig(f'/home/giacomocastagnetti/code/rec_genome_analysis/results/plots/{population}/{isolate}.png',
                        bbox_inches='tight')

            for reference in mismatch_distribution.keys():
                
                ref_file=f'/home/giacomocastagnetti/code/rec_genome_analysis/data/references/{reference}_assembly.fasta'
                sam_file=f'/home/giacomocastagnetti/code/rec_genome_analysis/results/mappings/references/{reference}.sam'

                mismatch_distribution, mapping = analyse_bam(sam_file, ref_file)

                figure=plt.figure()
                for name, distribution in mismatch_distribution.items():
                    distribution=np.convolve(distribution,np.ones(k),'valid')/k
                    l=len(distribution)
                    x=np.linspace(0,l,l)
                    plt.plot(distribution)
                plt.legend(mismatch_distribution.keys())
                plt.title(f'mutation density between references (on {reference}), with convolution of {k}')
                plt.ylabel('mutation density')
                plt.xlabel('bp')
                figure.savefig(f'/home/giacomocastagnetti/code/rec_genome_analysis/results/plots/references/{reference}.png')
                plt.close()