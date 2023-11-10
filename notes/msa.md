# create MSA that are used to produce the plots of recombination evidences 

## clone assemblies

we want to create a msa between assembly and references

mafft --auto results/seq_for_msa/P2C1_refs.fasta > results/msa/clones/P2_C1_msa.fasta

this method works, now we want to load this file in biopython and put the sequences in a numpy matrix to analyse it with a script.

## reads

### select for the best reads

#### longest

we assume that the best reads are the longest ones.
- we select the reads for the longest read.
- we have the longest read, now we want to put it in a file with the references.

now we produce the msa of this file.

mafft --auto results/reads/P2_1_0.fasta > results/msa/P2_1_0_msa.fasta

then we process the msa with the script produced before.

we need ot select for the reads in the last populations, they are the ones with highest probability of being recombinants. we will execute the whole thing on 10 reads at timepoints 5 and 7 of population two with the following command:

mafft --retree 2 --maxiterate 2 results/reads/P2/5/P2_5_0.fasta > results/msa/P2/5/P2_5_0_msa.fasta

no results, it is too slow. the problem is that the longest reads are not always the best aligned ones (by using minimap we get a shitty alignment).

#### best matching reads on the basis of read population alignment

we need another way to select reads and to align them, i think i will use the data of population alignment (from evo genome analysis). we select reads on the basis of the length of the matching on the reference, we consider the average between primary and secondary mapping.

with the new reads, we run this script:

populations=("P2" "P3")
timepoints=("1" "3" "5" "7")
reads=("0" "1" "2" "3" "4")
for population in "${populations[@]}"
do
    for timepoint in "${timepoints[@]}"
    do
        for read in "${reads[@]}"
        do
        mafft --retree 2 --maxiterate 2 results/seq_for_msa/$population/$timepoint/${population}_${timepoint}_${read}.fasta \
        > results/msa/$population/$timepoint/${population}_${timepoint}_${read}_msa.fasta
        done
    done
done

the procedure is still super slow. the problem is that some reads are mapping in the reverse, mafft is having hard time with that. by making the reverse complement of the read we can get  the msa of the read.

minimap2 -a data/P2C1_refs.fasta results/longest_reads/P2/5/P2_5_0.fasta > test_alignment_longest_reads_references/5_0_on3refs.sam

we will know if reads are reverse or not from the initial population alignment, if they are reverse mapped we build the msa files with the reverse complement of these sequences

the procedure is still super slow but let's go deeper in last timepoint

reads=("5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19")
for read in "${reads[@]}"
do
mafft --retree 2 --maxiterate 2 results/seq_for_msa/P2/7/P2_7_${read}.fasta > results/msa/P2/7/P2_7_${read}_msa.fasta
done

# the approach followed so far was completely wrong, we don't have to create a new msa for each read, we can use MAFFT --addfragments

## creating the MSA

creating msa of the two references:

mafft --auto results/seq_for_msa/refs.fasta > results/msa/refs_msa.fasta

add each read to it:

mafft --auto --addfragments fragments --reorder --thread -1 existing_alignment > output

    Sequences in fragments are ungapped and then aligned to existing_alignment. 
    fragments is a single multi-FASTA format file.
    existing_alignment is a single multi-FASTA format file. 
    Gaps in existing_alignment are preserved, but the alignment length may be changed in the default setting (see example below).
    If the --keeplength option is given, then the alignment length is unchanged.  Insertions at the fragmentary sequences are deleted. 
    Add --mapout to see a correspondence table of positions, fragments.map, between before and after the calculation.  The --mapout option automatically turns on the --keeplength option, to keep the numbering of sites in the reference alignment (explanation added, 2016/Aug). 
    --auto automatically switches algorithm according to data size.  Safer to always use this flag.  (added 2020/Sep)
    --multipair uses a high-cost (in time and memory usage) option.  Same as default.  Applicable to ∼<30,000 sites × ∼<1,000 sequences.
    --6merpair uses a low-cost option.
    Omit --reorder to preserve the original sequence order. 
    Described in Katoh & Frith 2012
    Can be used off-label to align closely-related sequences to a reference to build an MSA.

mafft --auto --addfragments results/longest_reads_seq/P2/7/P2_7_2.fasta --keeplength results/msa/refs_msa.fasta > results/test.fasta

the execution is super fast, we can go further and test multiple reads of the last timepoint:

reads=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19")
for read in "${reads[@]}"
do
mafft --auto --addfragments results/seq_for_msa/P2/7/P2_7_${read}.fasta --keeplength results/msa/refs_msa.fasta > results/msa/P2/7/P2_7_${read}_msa.fasta
done

## total msa production

populations=("P2" "P3")
timepoints=("1" "3" "5" "7")
for population in "${populations[@]}"
do
    for timepoint in "${timepoints[@]}"
    do
        for read in {0..999}
        do
        mafft --auto --addfragments results/seq_for_msa/${population}/${timepoint}/${population}_${timepoint}_${read}.fasta --keeplength results/msa/refs_msa.fasta \
        > results/msa/${population}/${timepoint}/${population}_${timepoint}_${read}_msa.fasta
        done
    done
done
