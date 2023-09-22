import gzip
from Bio import SeqIO,bgzf
import numpy as np
import matplotlib.pyplot as plt

'''

this script gives some statistics on the fastq files of the nanopore run,
it is not important for the recombination analysis

'''
len_dict={}

lens=[]

qual_dict={}
for n in range(100):
    qual_dict[n]=0

populations=['P2','P3']
clones=['C1','C2','C3','C4']

for population in populations:
    for clone in clones:

        '''
        path_in = f"./data/nanopore/{population}_{clone}.fastq.gz"
        path_out = f"./data/nanopore_length_threshold/{population}_{clone}.fastq.gz"
        handle_in = gzip.open(path_in, "rt")
        handle_out = gzip.open(path_out, "wt")

        fq = SeqIO.parse(handle_in, "fastq")
        for read in fq:
            # Only export reads that have a G in positions 7, 8,
            # and 9
            if len(read.seq)>500:
                handle_out.write(read.format("fastq"))

        handle_in.close()
        handle_out.close()

        '''
        
        with gzip.open(f"data/nanopore/{population}_{clone}.fastq.gz", "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                l=len(record.seq)
                if l not in len_dict.keys():
                    len_dict[l]=1
                else:
                    len_dict[l]+=1
                lens.append(l)

                #q=record.format('qual')
                #q_string=q.replace('\n',' ')
                #q_list=q_string.split(' ')[8:-1]
                #for nuc_qscore in q_list:
                #    qual_dict[int(nuc_qscore)]+=1

        print('the total number of reads is:', sum(len_dict.values()))
        max_len=max(len_dict.keys())

        b=np.linspace(0, max_len, max_len)
        len_distr=np.zeros(max_len)
        for pos,val in len_dict.items():
            len_distr[pos-1]+=val

        len_distr=np.convolve(len_distr, np.ones(100), mode='same')

        plt.plot(b,len_distr)
        plt.title(f'read length distribution {population} {clone}')
        plt.ylabel('number of reads')
        plt.xlabel('length of the read')
        plt.show()

        plt.hist(lens, len(lens), density=True, histtype='step', cumulative=True)
        plt.xlim([0,4000])
        plt.title(f'cumulative read length distibution {population} {clone}')
        plt.xlabel('length of the read')
        plt.ylabel('fraction of reads')
        plt.show()

        '''
        x=np.linspace(0,100,100)
        qual_distr=[]
        for q in qual_dict.values():
            qual_distr.append(q)

        plt.plot(x,qual_distr)
        plt.ylabel('number of nucleotides')
        plt.xlabel('quality score')
        plt.show()
        '''