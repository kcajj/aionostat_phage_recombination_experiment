import numpy as np
import matplotlib.pyplot as plt

from analyse_bam import analyse_bam

def plot_mappings(mismatch_distribution, mapping, population, isolate, k, out_folder):
    number_of_plots=len(mismatch_distribution.keys())+1
    fig, axs =plt.subplots(number_of_plots,sharex=True,constrained_layout = True, figsize=(8,10))

    fig.suptitle(f'mutation density distribution between {population}-{isolate} and references, with convolution window of {k}')

    for reference, distribution in mismatch_distribution.items():
        distribution=np.convolve(distribution,np.ones(k),'valid')/k
        l=len(distribution)
        x=np.linspace(0,l,l)
        axs[0].plot(x,distribution,color=phage_colors[reference])
    axs[0].legend(mismatch_distribution.keys())
    axs[0].set_title(reference)
    axs[0].set_ylabel('mutation density')
    axs[0].set_xlabel('bp')

    c=1
    for reference, distribution in mismatch_distribution.items():
        axs[c].set_title(reference)
        axs[c].set_ylabel('mutation density')
        axs[c].set_xlabel('bp')

        distribution=np.convolve(distribution,np.ones(k),'valid')/k
        l=len(distribution)
        x=np.linspace(0,l,l)
        axs[c].plot(x,distribution,c=phage_colors[reference])

        for alignment_type, maps in mapping[reference].items():
            for map in maps:
                axs[c].axvspan(map[0],map[1],facecolor=mapping_colors[alignment_type],alpha=0.2)

        c+=1
        
    fig.savefig(out_folder, bbox_inches='tight')
    plt.close()
    
if __name__ == "__main__":
    
    populations=['P2','P3']
    isolates=['C1','C2','C3','C4']
    k=1000

    phage_colors={'EC2D2':'C2','EM11':'C0','EM60':'C1'}
    mapping_colors={'primary':'g','supplementary':'b','secondary':'r'}

    for population in populations:
        for isolate in isolates:
            
            assembly_file=f'/home/giacomocastagnetti/code/rec_genome_analysis/results/assemblies/{population}/{isolate}.fasta'
            bam_file = f'/home/giacomocastagnetti/code/rec_genome_analysis/results/mappings/{population}/{isolate}.bam'
            out_folder = f'/home/giacomocastagnetti/code/rec_genome_analysis/results/plots/{population}/{isolate}.png'

            mismatch_distribution, mapping = analyse_bam(bam_file, assembly_file)

            plot_mappings(mismatch_distribution, mapping, population, isolate, k, out_folder)

    for reference in mismatch_distribution.keys():
        
        ref_file=f'/home/giacomocastagnetti/code/rec_genome_analysis/data/references/{reference}_assembly.fasta'
        sam_file=f'/home/giacomocastagnetti/code/rec_genome_analysis/results/mappings/references/{reference}.sam'
        out_folder=f'/home/giacomocastagnetti/code/rec_genome_analysis/results/plots/references/{reference}.png'

        mismatch_distribution, mapping = analyse_bam(sam_file, ref_file)

        plot_mappings(mismatch_distribution, mapping, reference, '', k, out_folder)