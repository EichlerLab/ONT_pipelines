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
    conda:
        "../envs/methylink.yaml"
    log:
        "log/{sample}.{phase}.methlylink.log",
    script:
        "../scripts/append_mod_tags.py"
