nanopore_reads = 'data/nanopore/{population}_{isolate}.fastq.gz'
reference = 'data/references/{phage}_assembly.fasta'

rule flye:
    input:
        reads = nanopore_reads
    output:
        flye_folder = directory('results/assemblies/{population}/{isolate}_flye'),
        assembly = 'results/assemblies/{population}/{isolate}.fasta'
    params:
        genome_size = '0.150m',
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

rule alginment_to_assembly:
    input:
        assembly = rules.flye.output.assembly,
        references = expand(reference, phage=['EC2D2','EM11','EM60'])
    output:
        alignment = 'results/mappings/{population}/{isolate}.sam',
        bam = 'results/mappings/{population}/{isolate}.bam',
        bai = 'results/mappings/{population}/{isolate}.bam.bai'
    conda:
        'conda_envs/read_mapping.yml'
    params:
        cores = 4
    shell:
        """
        minimap2 -a \
            {input.assembly} \
            {input.references} \
            > {output.alignment}

        samtools sort -@ {params.cores} \
            -o {output.bam} \
            {output.alignment}
        samtools index {output.bam} \
            {output.bai}
        """

rule alignment_references:
    input:
        reference = lambda w: expand(reference, phage=w.ref),
        query = lambda w: expand(reference, phage=['EC2D2','EM11','EM60'])
    output:
        alignment = 'results/mappings/references/{ref}.sam'
    conda:
        'conda_envs/read_mapping.yml'
    params:
        cores = 4
    shell:
        """
        minimap2 -a \
            {input.reference} \
            {input.query} \
            > {output.alignment}
        """

rule all:
    input:
        assemblies = expand(rules.alginment_to_assembly.output.alignment,population=['P2','P3'],isolate=['C1','C2','C3','C4']),
        references_alignments = expand(rules.alignment_references.output.alignment,ref=['EC2D2','EM11','EM60'])