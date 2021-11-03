

rule clair:
	input:
		merged_bam = rules.merge_run_aln.output.merged_bam,
		index = rules.index_aln.output.merged_bai,
		ref = REF
	output:
		vcf = 'alignments/{sample}/{sample}.{bc_vers}.minimap2.{seq}.clair3.vcf'
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
	conda:
		"../envs/env.yaml"
	resources:
		mem=10,
		hrs=24
	threads: 1 
	shell:
		"""
		{CLAIR_DIR}/clair3.sh -s {wildcards.sample} -r {input.ref} -b {input.merged_bam} -o $( dirname {output.vcf} )
		"""


rule sniffles:
	input:
		merged_bam = rules.merge_run_aln.output.merged_bam,
		index = rules.index_aln.output.merged_bai,
		ref = REF
	output:
		vcf = 'alignments/{sample}/{sample}.{bc_vers}.minimap2.{seq}.sniffles.vcf'
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
		'sniffles/202109'
	conda:
		"../envs/env.yaml"
	resources:
		mem=10,
		hrs=24
	threads: 1 
	shell:
		"""
		sniffles -m {input.merged_bam} -v {output.vcf}
		"""


rule cuteSV:
	input:
		merged_bam = rules.merge_run_aln.output.merged_bam,
		index = rules.index_aln.output.merged_bai,
		ref = REF
	output:
		cuteSV_vcf = 'alignments/{sample}/{sample}.{bc_vers}.minimap2.{seq}.cuteSV.vcf'
	conda:
		"../envs/env.yaml"
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
		'cuteSV/1.0.11'
	resources:
		mem=10,
		hrs=24
	threads: 8
	shell:
		"""
		cuteSV -t {threads} --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 {input.merged_bam} {input.ref} {output.cuteSV_vcf} $( dirname {output.cuteSV_vcf} )
		"""

rule svim:
	input:
		merged_bam = rules.merge_run_aln.output.merged_bam,
		index = rules.index_aln.output.merged_bai,
		ref = REF
	output:
		vcf = temp('alignments/{sample}/{sample}.{bc_vers}.minimap2.{seq}.svim.vcf')
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
		'svim/1.4.2'
	conda:
		"../envs/env.yaml"
	resources:
		mem=16,
		hrs=24
	threads: 1
	shell:
		"""
		svim --sample {wildcards.sample} $( dirname {output.vcf} ) {input.merged_bam} {input.ref}
		"""


rule bgzip_vcf:
	input:
		vcf = 'alignments/{sample}/{sample}.{bc_vers}.minimap2.{seq}.{var_caller}.vcf'
	output:
		zipped = 'alignments/{sample}/{sample}.{bc_vers}.minimap2.{seq}.{var_caller}.vcf.gz'
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
		'tabix/0.2.6'
	conda:
		"../envs/env.yaml"
	resources:
		mem=10,
		hrs=24
	threads: 1
	shell:
		'''
		bgzip -c {input.vcf} > {output.zipped}
		tabix {output.zipped}
		'''	