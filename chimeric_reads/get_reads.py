import pysam
import matplotlib.pyplot as plt
from collections import defaultdict
from Bio import SeqIO

def find_reads_with_secondary_mapping(bam_file):
    primary_positions = {}
    secondary_positions = defaultdict(list)
    c=0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            # Check if the read has a secondary alignment
            if read.is_secondary:
                matching_length=len(read.get_aligned_pairs(matches_only=True))
                secondary_positions[read.query_name].append([read.reference_name,read.reference_start,matching_length])
                c+=1
                #if c==10: break
        for read in bam.fetch():
            if read.query_name in secondary_positions.keys():
                if not(read.is_secondary) and not(read.is_supplementary):
                    matching_length=len(read.get_aligned_pairs(matches_only=True))
                    primary_positions[read.query_name]=[read.reference_name,read.reference_start,matching_length]

    return primary_positions, secondary_positions

def clean_secondary(primary, secondary):
    to_delete=[]

    for primary_query_name, [primary_reference_name, primary_reference_start,l1] in primary.items():
        for secondary_query_name, secondary_mappings in secondary.items():
            if primary_query_name==secondary_query_name:
                for [secondary_reference_name, secondary_reference_start,l2] in secondary_mappings:
                    if secondary_reference_name==primary_reference_name:
                        secondary[secondary_query_name].remove([secondary_reference_name, secondary_reference_start,l2])
                        if secondary[secondary_query_name]==[]:
                            to_delete.append(primary_query_name)
    
    for name in to_delete:
        primary.pop(name)
        secondary.pop(name)

    return primary, secondary

def get_longest_reads(primary, secondary):
    longest_matching_reads={}
    for primary_query_name, [ref_name1, ref_start1, primary_alignment_length] in primary.items():
        for secondary_query_name, secondary_mappings in secondary.items():
            if primary_query_name==secondary_query_name:
                sorted_secondary_alignments=sorted(secondary_mappings, key=lambda x: x[2], reverse=True)
                secondary[secondary_query_name]=sorted_secondary_alignments[0]
                longest_matching_reads[primary_query_name]=[(primary_alignment_length+sorted_secondary_alignments[0][2])/2,ref_start1,sorted_secondary_alignments[0][1],ref_name1,sorted_secondary_alignments[0][0]]

    sorted_dict=sorted(longest_matching_reads.items(), key=lambda x: x[1][0], reverse=True)
    sorted_dict=sorted_dict[:5] #the third longest read matching in the file has some common matches
    return secondary, sorted_dict

def find_chimeric_reads(bam_file, longest_matching_reads, ref1_seq, ref2_seq):
    for [read_name,[mean_length,ref_start1,ref_start2,ref1,ref2]] in longest_matching_reads:
        unique_primary=0
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam.fetch():
                # Check if the read has a secondary alignment
                if read.query_name==read_name and read.reference_start==ref_start1 and not(read.is_secondary) and not(read.is_supplementary):
                    alignment1=read.get_aligned_pairs()
                    read_seq=read.query_sequence
                    if read.reference_name=='EM60':
                        ref1_seq, ref2_seq = ref2_seq, ref1_seq
                    unique_primary+=1
        if unique_primary!=1: print('error primary')
        unique_secondary=0
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam.fetch():
                # Check if the read has a secondary alignment
                if read.is_secondary and read.query_name==read_name and read.reference_name==ref2 and read.reference_start==ref_start2:
                    alignment2=read.get_aligned_pairs()
                    unique_secondary+=1
                    for ic, (block_type, block_len) in enumerate(read.cigar):
                        hard_clip2=0
                        if block_type==5:
                            hard_clip2=block_len
                            print('clip secondary', block_type, block_len)
                        break
        if unique_secondary!=1: print('error secondary')
        
        iref_1=0
        iref_2=0
        c=0

        for read_pos, read_nuc in enumerate(read_seq):

            #get the tuple in which the position of the read is pairing
            while alignment1[iref_1][0]!=read_pos:
                iref_1+=1
            #check that the position of the read is not pairing to a gap
            if alignment1[iref_1][1]==None:
                ref1_nuc=None
            else:
                ref1_nuc=ref1_seq[alignment1[iref_1][1]]

            while alignment2[iref_2+hard_clip2][0]!=read_pos:
                iref_2+=1
            if alignment2[iref_2+hard_clip2][1]==None:
                ref2_nuc=None
            else:
                ref2_nuc=ref2_seq[alignment2[iref_2+hard_clip2][1]]

            if ref1_nuc!=None and ref2_nuc!=None:
                c+=1
                #print(read_pos, read_nuc)
                #print(iref_1, ref1_nuc)
                #print(iref_2, ref2_nuc)
                #print()
        print(c)
        #for read_pos, read_nuc_2 in enumerate(read_seq2):
        #    while alignment2[i_2][0]!=read_pos:
        #        i_2+=1
        #    ref2_nuc=ref2_seq[alignment2[i_2][1]]
        #    print(read_nuc_1, ref1_nuc)
        #    print(read_nuc_2, ref2_nuc)'''

def extract_seq(path):
    return str(SeqIO.read(path, 'fasta').seq)

def get_points(pp,sp):
    points=[]
    for name,secondary_mappings in sp.items():
        for secondary_mapping in secondary_mappings:
            points.append((pp[name],secondary_mapping))
    return points


if __name__ == "__main__":
    
    #populations=['P2','P3']
    #timepoints=['1','3','5','7']
    populations=['P2']
    timepoints=['7']

    for population in populations:
        for timepoint in timepoints:
            bam_file_path = f'chimeric_reads/data/{population}_{timepoint}.bam'
            ref1_seq=extract_seq('/home/giacomocastagnetti/code/rec_genome_analysis/data/references/EM11_assembly.fasta')
            ref2_seq=extract_seq('/home/giacomocastagnetti/code/rec_genome_analysis/data/references/EM60_assembly.fasta')

            primary_positions, secondary_positions = find_reads_with_secondary_mapping(bam_file_path)

            primary_positions, secondary_positions = clean_secondary(primary_positions, secondary_positions)

            n=1
            secondary_positions, longest_matching_reads = get_longest_reads(primary_positions, secondary_positions)

            find_chimeric_reads(bam_file_path, longest_matching_reads,ref1_seq,ref2_seq)



            '''
            points = get_points(primary_positions,secondary_positions)
            
            figure=plt.figure()
            x=[]
            y=[]
            for p in points:
                x.append(p[0])
                y.append(p[1])

            plt.scatter(x,y,alpha=0.2)
            plt.ylabel('secondary mapping')
            plt.xlabel('primary mapping')
            figure.savefig(f'plots/secondary_mapping/{phage}/{time}.png')
            plt.close()
            '''