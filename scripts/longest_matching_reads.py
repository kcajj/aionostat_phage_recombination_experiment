import pysam
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd

'''
finds the reads in a bam file that have the longest matching length
(average between primary and secondary mapping)
'''
def find_reads_with_secondary_mapping(bam_file):
    '''
    looks at all the secondary reads and stores read name and an array of useful information about the mapping in a dictionary
    (not all the info in the array are useful actually, it's just a long array with random information).
    then it looks at all the reads again to find the primary mappings of the secondary reads. it stores them in another
    dictionary together with their mapping information
    '''
    primary_positions = {}
    secondary_positions = defaultdict(list)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            # Check if the read has a secondary alignment
            if read.is_secondary:
                matching_length=len(read.get_aligned_pairs(matches_only=True))
                secondary_positions[read.query_name].append([read.reference_name,read.reference_start,matching_length,read.is_reverse])
        for read in bam.fetch():
            if read.query_name in secondary_positions.keys():
                if not(read.is_secondary) and not(read.is_supplementary):
                    matching_length=len(read.get_aligned_pairs(matches_only=True))
                    primary_positions[read.query_name]=[read.reference_name,read.reference_start,matching_length, read.is_reverse]

    return primary_positions, secondary_positions

def clean_secondary(primary, secondary):
    '''
    cleans the two dictionary produced by the previous function.
    it removes the read if the primary and secondary mapping are on the same reference.
    '''
    to_delete=[]

    for primary_query_name, [primary_reference_name, primary_reference_start,l1,rev] in primary.items():
        for secondary_query_name, secondary_mappings in secondary.items():
            if primary_query_name==secondary_query_name:
                for [secondary_reference_name, secondary_reference_start,l2,rev] in secondary_mappings:
                    if secondary_reference_name==primary_reference_name:
                        secondary[secondary_query_name].remove([secondary_reference_name, secondary_reference_start,l2,rev])
                        if secondary[secondary_query_name]==[]:
                            to_delete.append(primary_query_name)
    
    for name in to_delete:
        primary.pop(name)
        secondary.pop(name)

    return primary, secondary

def get_longest_reads(primary, secondary, n):
    '''
    saves only the longest secondary mapping for each read among the remaining ones in the dictionary
    '''
    longest_matching_reads={}
    for primary_query_name, [ref_name1, ref_start1, primary_alignment_length, is_reverse] in primary.items():
        for secondary_query_name, secondary_mappings in secondary.items():
            if primary_query_name==secondary_query_name:
                sorted_secondary_alignments=sorted(secondary_mappings, key=lambda x: x[2], reverse=True)
                secondary[secondary_query_name]=sorted_secondary_alignments[0]
                longest_matching_reads[primary_query_name]=[(primary_alignment_length+sorted_secondary_alignments[0][2])/2,ref_start1,sorted_secondary_alignments[0][1],ref_name1,sorted_secondary_alignments[0][0],is_reverse,sorted_secondary_alignments[0][3]]

    sorted_dict=sorted(longest_matching_reads.items(), key=lambda x: x[1][0], reverse=True)
    sorted_dict=sorted_dict[:n]

    return secondary, sorted_dict

if __name__ == "__main__":
    
    populations=['P2','P3']
    timepoints=['1','3','5','7']

    for population in populations:
        for timepoint in timepoints:
            
            print(f'selecting the best reads from {population} {timepoint}')

            bam_file_path = f'data/population_alignments/{population}/{population}_{timepoint}.bam'

            primary_positions, secondary_positions = find_reads_with_secondary_mapping(bam_file_path)

            primary_positions, secondary_positions = clean_secondary(primary_positions, secondary_positions)

            n=1000 #number of best reads to select
            secondary_positions, longest_matching_reads = get_longest_reads(primary_positions, secondary_positions, n)

            out_csv=f'results/longest_matching_reads/{population}/{population}_{timepoint}.csv'

            d=[]
            for read_name, data in longest_matching_reads:
                d.append({})
                d[-1]['read_name']=read_name
                d[-1]['avg_match']=data[0]
                d[-1]['ref_start_primary']=data[1]
                d[-1]['ref_start_secondary']=data[2]
                d[-1]['ref_name_primary']=data[3]
                d[-1]['ref_name_secondary']=data[4]
                d[-1]['is_reverse_primary']=data[5]
                d[-1]['is_reverse_secondary']=data[6]
            
            df=pd.DataFrame.from_dict(d)

            df.to_csv(out_csv)