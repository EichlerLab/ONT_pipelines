rule link_bam:
    input:
        aln_bam=find_aln_bam,
        aln_bai=find_aln_bai,
        methyl_bam=find_aln_list,
    output:
        linked_bam="methyl_aln/{sample}/{sample}.{phase}.methyl.bam",
    resources:
        mem=4,
        hrs=24,
        disk_free=250,
    threads: 16
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "methylink/0.4.0",
    conda:
        "../envs/methylink.yaml"
    log:
        "log/{sample}.{phase}.methlylink.log",
    shell:
        """
        methylink \
            --threads {threads} \
            --aln {input.aln_bam} \
            --sample {wildcards.sample}_{wildcards.phase} \
            --methyl_bams "$(echo {input.methyl_bam})" \
            --output {output.linked_bam} 2>&1 | tee {log}
        """
