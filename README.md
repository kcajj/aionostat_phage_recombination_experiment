# detect genome recombination

the objective of this repository is to carry out and track the exploratory work done on the data from an Aionostat experiment that aims to find recombination between phages.

this repository allows to:

1. create the plots of mutation density distribution for phage clones assembly with respect to multiple reference genomes.
    - [analyse_bam.py](scripts/analyse_bam.py) defines a function that takes in input a bam file (alignment between assembly and reference) and a reference. It returnsthe distribution of mismatches between the assembly and the reference in an array.
    - [mutation_density.py](scripts/mutation_densisty.py) uses the function defined in analyse_bam.py, it takes the arrays of the distribution of mutation density of an assembly against multiple references and it creates a plot out of it.

2. create plots of recombination evidences
    1. in single phage clones
        - create an MSA between the assembly and the references with MAFFT, for now you have to do it manually
        - [process_clone_msa.py](scripts/process_clone_msa.py) takes in input the msa of the alignment of clone assemblies with the references, computes the evidence distribution and builds its plot.
    2. in multiple reads (coming from the population alignment data) at the same time.
        - [longest_matching_reads.py](scripts/longest_matching_reads.py) takes in input the bam file of the mapping of a population of reads to multiple references. Among the reads that have a secondary mapping it selects the ones that have the longest mapping length, considering the average between primary and secondary mapping. It saves the names and other useful data of these reads in a csv file.
        - [save_longest_reads.py](scripts/save_longest_reads.py) takes in input a fastq file containing all the reads and the list of the best reads. it saves the sequence of each long matching read into an individual fasta file, that will be used by MAFFT
        - build MSA with MAFFT: for now you have to do it manually, i created a for loop in the shell
        - [process_msa.py](scripts/process_msa.py) takes in input multiple MSAs from the best matching reads and extracts the recombination evidences from each of them. then plots all the recombination evidences together in a single plot.