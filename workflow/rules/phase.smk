rule long_phase:
    input:
        snv_vcf="alignments/{sample}/{sample}.clair3.vcf.gz",
        bam=rules.merge_run_aln.output.merged_bam,
        bai=rules.index_aln.output.merged_bai,
        sv_vcf="alignments/{sample}/{sample}.sniffles.vcf.gz",
    output:
        bam="phased_aln/{sample}/{sample}.minimap2.longphase.bam",
    resources:
        mem=4,
        hrs=24,
        disk_free=1,
    threads: 12
    log:
        "log/{sample}.merge_all.log",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "samtools/1.12",
    shell:
        """
        {SDIR}/scripts/longphase haplotag —snp-file={input.snv_vcf} —bam-file={input.bam} —qualityThreshold=1 -t {threads} —sv-file={input.sv_vcf} -o {output.bam}
        """
