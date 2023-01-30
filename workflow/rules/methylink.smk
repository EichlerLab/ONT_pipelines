rule link_bam:
    input:
        bam=find_aln_bam,
        methyl_bam=find_methyl_list,
    output:
        linked_bam="methyl_aln/{sample}/{sample}.{phase}.methyl.bam",
    resources:
        mem=4,
        hrs=24,
        disk_free=1,
    threads: 12
    log:
        "log/{sample}.methlylink.log",
    conda:
        "../envs/align.yaml"
    script:
        """
        scripts/append_mod_tags.py
        """
