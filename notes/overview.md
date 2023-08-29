# using the pipeline

I'm planning how to use which pipeline.

I should start from the data i need to use

# data

the Aionostat run that generated the data had 2 vials with 3 phages in each of them. A third vial was set up as a negative control without phages.
phage samples:
1. phage populations: 2 populations, 4 timepoints each. samples from each vial during the experiment.
2. phage isolates: lest timepoint, 2 populations, 4 isolates each. single phages isolated from the populations of the last day, 4 phages for each population.

we want to look at population alignment data and look at the assemblies of the isolates.

Valentin ha anche sequenziato vari batteri dell'esperimento:
1. initial bacteria (wbbl+)
2. bacterial culture samples: samples from bacterial vials, at the end of the experiment. these samples were collected because bacteria started replicating faster and we want to understand why.
3. resistant bacteria: 3 clones. some bacteria survived in the phage vial, 3 colonies were isolated and sequenced.

we want to assemble the initial bacterial genome (1), then use it as a reference to map on it the data of the bacteria from the culture vials of the experiment. the resistant bacteria should be assembled and compared with initial bacteria.

the idea is to keep bacteria and phage data anaysis separated.

## folder organisation

we want to have a data folder shared on the cluster, then we create links to the data in my folder, this is called symlink, this is the command: ln -s source destination(file_name)

when i create folders in the group folder, to put there the data and the results i have to give group permissions to my folder. you can do it with chmod g+rwx filename

## recombination genome analysis

1. assembly: we want to assemble genomes sequenced with ont, they are divided in populations, each population has a different reference.

2. mappings:
    1. references: each reference mapped against the others
        ref1
    2. assemblies: map each reference onto the assembly
        p2-c1: ref1, ref2