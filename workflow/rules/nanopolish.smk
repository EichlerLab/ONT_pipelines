

rule concat_fastq:
    input:
        fastq=concatenate_fastq,
    output:
        fastq=temp("tmp/reads/{sample}.{bc_vers}.{seq}.fastq.gz"),
    threads: 1
    conda:
        "../envs/nanopolish.yaml"
    log:
        "log/{sample}_{bc_vers}_{seq}.concat_fastq.log",
    resources:
        mem=8,
        hrs=12,
    shell:
        """
        cat {input.fastq} > {output.fastq}
        """


rule index_fastq:
    input:
        fastq=rules.concat_fastq.output.fastq,
    output:
        index=temp("tmp/reads/{sample}.{bc_vers}.{seq}.fastq.gz.index"),
        fai=temp("tmp/reads/{sample}.{bc_vers}.{seq}.fastq.gz.index.fai"),
        gzi=temp("tmp/reads/{sample}.{bc_vers}.{seq}.fastq.gz.index.gzi"),
        readdb=temp("tmp/reads/{sample}.{bc_vers}.{seq}.fastq.gz.index.readdb"),
    params:
        directory=find_fast5_dir,
        summary=find_summary_fofn,
    log:
        "log/{sample}_{bc_vers}_{seq}.index_fastq.log",
    threads: 1
    conda:
        "../envs/nanopolish.yaml"
    resources:
        mem=8,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "nanopolish/0.13.3",
    shell:
        """
        nanopolish index {params.directory} {params.summary} {input.fastq}
        """


rule nanopolish:
    input:
        fastq=rules.concat_fastq.output.fastq,
        sorted_bam=rules.merge_run_aln.output.merged_bam,
        bam_index=rules.index_aln.output.merged_bai,
        reference=REF,
        index=rules.index_fastq.output.index,
        fai=rules.index_fastq.output.fai,
        gzi=rules.index_fastq.output.gzi,
        readdb=rules.index_fastq.output.readdb,
    output:
        methyl=temp("tmp/{sample}/{window}.{bc_vers}.nanopolish.{seq}.tsv"),
    threads: 8
    conda:
        "../envs/nanopolish.yaml"
    log:
        "log/{sample}_{bc_vers}_{seq}_{window}.nanopolish.log",
    resources:
        mem=2,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "nanopolish/0.13.3",
    shell:
        """
        nanopolish call-methylation -t {threads} -r {input.fastq} -b {input.sorted_bam} -g {input.reference} -w {wildcards.window} > {output.methyl}
        """


rule gather_nanopolish:
    input:
        tsv=find_all_windows,
    log:
        "log/{sample}_{bc_vers}_{seq}.nanopolish.log",
    output:
        tsv_all="alignments/{sample}/{sample}.{bc_vers}.nanopolish.{seq}.tsv.gz",
    resources:
        mem=2,
        hrs=72,
    threads: 1
    run:
        df = pd.concat([pd.read_csv(file, sep="\t") for file in input.tsv])
        df.to_csv(output.tsv_all, sep="\t", index=False)


rule h5_nanopolish:
    input:
        tsv=rules.gather_nanopolish.output.tsv_all
    log:
        "log/{sample}_{bc_vers}_{seq}.nanopolish_h5.log",
    output:
        h5="alignments/{sample}/{sample}.{bc_vers}.nanopolish.{seq}.m5",
    resources:
        mem=2,
        hrs=72,
    conda:
        "../envs/nanopolish_h5.yaml"   
    threads: 1
    shell:
        '''
        meth5 create_m5 --input_paths {input.tsv} --output_file {output.h5}
        '''
