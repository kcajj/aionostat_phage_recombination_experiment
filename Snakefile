configfile: "config.yml"

nanopore_reads = 'data/nanopore/{phage}_{tag}.fastq.gz'
reference = 'data/references/{phage}_reference.fasta'

rule flye:
    input:
        reads = nanopore_reads
    output:
        flye_folder = directory('results/{phage}/assemblies/{tag}_flye'),
        assembly = 'results/{phage}/assemblies/{tag}.fasta'
    params:
        genome_size = lambda w : config["genome-size"][w.phage],
        cores = 4,
        coverage = 40
    conda:
        'conda_envs/genome_assembly.yml'
    shell:
        """
        flye --nano-hq {input.reads} \
            --out-dir {output.flye_folder} \
            --threads {params.cores} \
            --genome-size {params.genome_size} \
            --asm-coverage {params.coverage}
        
        cp {output.flye_folder}/assembly.fasta {output.assembly}
        """

rule minimap:
    input:
        reads = lambda w : expand(rules.flye.input.reads, tag=w.qry_tag, phage=w.phage),
        reference = lambda w : expand(rules.flye.output.assembly, tag=w.ref_tag, phage=w.phage)
    output:
        alignment = 'results/{phage}/mapping/{ref_tag}/{qry_tag}.sam'
    conda:
        'conda_envs/read_mapping.yml'
    shell:
        """
        minimap2 -ax map-ont \
            {input.reference} \
            {input.reads} \
            > {output.alignment}
        """

rule bam:
    input:
        sam = rules.minimap.output.alignment
    output:
        bam = 'results/{phage}/mapping/{ref_tag}/{qry_tag}.bam',
        bai = 'results/{phage}/mapping/{ref_tag}/{qry_tag}.bam.bai'
    conda:
        'conda_envs/read_mapping.yml'
    params:
        cores = 4
    shell:
        """
        samtools sort -@ {params.cores} \
            -o {output.bam} \
            {input.sam}
        samtools index {output.bam} \
            {output.bai}
        """

rule alginment_to_ref:
    input:
        reference = reference,
        assembly = rules.flye.output.assembly
    output:
        alignment = 'results/{phage}/mapping/reference/alignment_with_reference_{tag}.bam'
    conda:
        'conda_envs/read_mapping.yml'
    shell:
        """
        minimap2 -a \
            {input.reference} \
            {input.assembly} \
            > {output.alignment}
        """

rule all:
    input:
        new_chemistry_assemblies = expand(rules.plot_pileup.output.plot_folder,ref_tag='new_chemistry',qry_tag=['new_chemistry','1','3','5'],phage=['EC2D2','EM11','EM60']),
        old_chemistry_assemblies = expand(rules.plot_pileup.output.plot_folder,ref_tag='old_chemistry',qry_tag='old_chemistry',phage=['EC2D2','EM11','EM60']),
        reference_alignments = expand(rules.alginment_to_ref.output.alignment,phage=['EC2D2','EM11','EM60'],tag=['new_chemistry','old_chemistry']),