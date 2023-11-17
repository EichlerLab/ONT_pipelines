
# rule clair_chr:
#     input:
#         merged_bam=rules.merge_run_aln.output.merged_bam,
#         index=rules.index_aln.output.merged_bai,
#         ref=REF,
#     output:
#         vcf="tmp/variants/{sample}/{chrom}/phased_merge_output.vcf.gz",
#         gvcf="tmp/variants/{sample}/{chrom}/merge_output.gvcf.gz",
#     log:
#         "log/{sample}_{chrom}.clair.log",
#     conda:
#         "../envs/clair3.yaml"
#     threads: 1
#     resources: 
#         mem=lambda wildcards, input, attempt: max(input.size//1e9//8, 16),
#         hrs=24,
#         disk_free=1,
#     shell:
#         """
#         run_clair3.sh --bam_fn={input.merged_bam} --sample_name={wildcards.sample} --ref_fn={input.ref} --threads={threads} --platform=ont --model_path=$(dirname $( which run_clair3.sh  ) )/models/ont_guppy5 --output=$(dirname {output.vcf}) --ctg_name={wildcards.chrom} --enable_phasing --gvcf
#         if [[ ( -f $( echo {output.vcf} | sed 's/phased_//' ) ) && ( ! -f {output.vcf} ) ]]; then
#             cp $( echo {output.vcf} | sed 's/phased_//' ) {output.vcf}
#         elif [[ ( ! -f $( echo {output.vcf} | sed 's/phased_//' ) ) && ( ! -f {output.vcf} ) ]]; then
#             zcat $(dirname {output.vcf})/pileup.vcf.gz | grep ^# | bgzip -c > {output.vcf}
#         fi
#         if [[ ( -f $( echo {output.vcf} | sed 's/phased_//' ) ) && ( ! -f {output.gvcf} ) ]]; then
#             cp $( echo {output.vcf} | sed 's/phased_//' ) {output.gvcf}
#         elif [[ ( ! -f $( echo {output.vcf} | sed 's/phased_//' ) ) && ( ! -f {output.gvcf} ) ]]; then
#             zcat $(dirname {output.vcf})/pileup.vcf.gz | grep ^# | bgzip -c > {output.gvcf}
#         fi
#         """


# rule concat_clair:
#     input:
#         vcf=find_clair_chrs,
#         gvcf=find_clair_gvcf,
#     output:
#         vcf=temp("variants/{sample}/{sample}.clair3.vcf"),
#         gvcf="variants/{sample}/{sample}.clair3.gvcf.gz",
#     envmodules:
#         "modules",
#         "modules-init",
#         "modules-gs/prod",
#         "modules-eichler/prod",
#     log:
#         "log/{sample}.clair.log",
#     conda:
#         "../envs/vcf.yaml"
#     resources:
#         mem=10,
#         hrs=24,
#         disk_free=1,
#     threads: 1
#     shell:
#         """
#         bcftools concat -O v -o {output.vcf} {input.vcf}
#         bcftools concat -O v {input.gvcf} | bcftools sort -O v | bgzip -c > {output.gvcf} 
#         """


rule sniffles_somatic:
    input:
        merged_bam=rules.merge_run_aln.output.merged_bam,
        index=rules.index_aln.output.merged_bai,
        ref=REF,
    output:
        vcf=temp("variants/{sample}/{sample}.sniffles-somatic.vcf"),
    log:
        "log/{sample}.sniffles.log",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "sniffles/2.2",
    conda:
        "../envs/sniffles.yaml"
    resources:
        mem=24,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        sniffles --somatic -i {input.merged_bam} --reference {input.ref} --output-rnames -v {output.vcf}
        """


# rule cuteSV:
#     input:
#         merged_bam=rules.merge_run_aln.output.merged_bam,
#         index=rules.index_aln.output.merged_bai,
#         ref=REF,
#     output:
#         cuteSV_vcf=temp("variants/{sample}/{sample}.cuteSV.vcf"),
#         cuteSV_ref=temp("tmp/cuteSV/{sample}/ref.out"),
#     conda:
#         "../envs/cutesv.yaml"
#     log:
#         "log/{sample}.cutesv.log",
#     envmodules:
#         "modules",
#         "modules-init",
#         "modules-gs/prod",
#         "modules-eichler/prod",
#         "cuteSV/1.0.11",
#     resources:
#         mem=10,
#         hrs=24,
#         disk_free=1,
#     threads: 8
#     shell:
#         """
#         zcat {input.ref} > {output.cuteSV_ref} || cat {input.ref} > {output.cuteSV_ref}
#         cuteSV -t {threads} --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 {input.merged_bam} {output.cuteSV_ref} {output.cuteSV_vcf} $( dirname {output.cuteSV_vcf} )
#         """


# rule svim:
#     input:
#         merged_bam=rules.merge_run_aln.output.merged_bam,
#         index=rules.index_aln.output.merged_bai,
#         ref=REF,
#     output:
#         vcf_tmp=temp("tmp/variants/{sample}/variants.vcf"),
#         vcf="variants/{sample}/{sample}.svim.vcf",
#     envmodules:
#         "modules",
#         "modules-init",
#         "modules-gs/prod",
#         "modules-eichler/prod",
#         "svim/2.0.0",
#     log:
#         "log/{sample}.svim.log",
#     conda:
#         "../envs/svim.yaml"
#     resources:
#         mem=16,
#         hrs=24,
#         disk_free=1,
#     threads: 1
#     shell:
#         """
#         svim alignment --sample {wildcards.sample} $( dirname {output.vcf_tmp} ) {input.merged_bam} {input.ref}
#         cp -l {output.vcf_tmp} {output.vcf}
#         """


# rule bgzip_vcf:
#     input:
#         vcf="variants/{sample}/{sample}.{var_caller}.vcf",
#     output:
#         zipped=("variants/{sample}/{sample}.{var_caller}.vcf.gz"),
#     log:
#         "log/{sample}.{var_caller}_zip.log",
#     envmodules:
#         "modules",
#         "modules-init",
#         "modules-gs/prod",
#         "modules-eichler/prod",
#         "tabix/0.2.6",
#     conda:
#         "../envs/vcf.yaml"
#     resources:
#         mem=10,
#         hrs=24,
#         disk_free=1,
#     threads: 1
#     shell:
#         """
#         bcftools sort -o /dev/stdout -O v {input.vcf} | bgzip -c > {output.zipped}
#         tabix {output.zipped}
#         """
