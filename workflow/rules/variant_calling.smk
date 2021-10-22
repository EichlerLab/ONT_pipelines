import os
import pandas as pd
import numpy as np
from pathlib import Path

import pdb;
tr=pdb.set_trace

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
shell.prefix("source %s/env.cfg; " % SNAKEMAKE_DIR)
configfile: "svcalling.snake_config.yaml"

CLAIR_DIR = config['clair_dir']


MANIFEST = pd.read_csv(config['manifest'], sep='\t')



def find_input_file(wildcards):
	return MANIFEST.loc[(MANIFEST['SAMPLE'] == wildcards.sample) & (MANIFEST['ALN_TYPE'] == wildcards.aln_type), 'BAM'].values[0]


localrules: all, clair


rule all:
	input:
		# expand('svcalling/clair/{sample}_{aln_type}/run_clair3.log', zip, aln_type=MANIFEST['ALN_TYPE'].values, sample=MANIFEST['SAMPLE'].values),
		expand('svcalling/sniffles/{sample}_{aln_type}_sniffles.vcf', zip, aln_type=MANIFEST['ALN_TYPE'].values, sample=MANIFEST['SAMPLE'].values),
		expand('svcalling/cuteSV/{sample}_{aln_type}_cuteSV.vcf', zip, aln_type=MANIFEST['ALN_TYPE'].values, sample=MANIFEST['SAMPLE'].values)



rule clair:
	input:
		merged_bam = find_input_file,
		ref = config['ref_hg38']
	output:
		clair_vcf = 'svcalling/clair/{sample}_{aln_type}/run_clair3.log'
	params:
		out_dir = 'svcalling/clair/{sample}_{aln_type}/'
	resources:
		mem=10,
		hrs=24
	threads: 1 
	shell:
		"""
		{CLAIR_DIR}/clair3.sh -s {wildcards.sample} -r {input.ref} -b {input.merged_bam} -o {params.out_dir}
		"""


rule sniffles:
	input:
		merged_bam = find_input_file,
		ref = config['ref_hg38']
	output:
		sniffles_vcf = 'svcalling/sniffles/{sample}_{aln_type}_sniffles.vcf'
	resources:
		mem=10,
		hrs=24
	threads: 1 
	shell:
		"""
		sniffles -m {input.merged_bam} -v {output.sniffles_vcf}
		"""


rule cuteSV:
	input:
		merged_bam = find_input_file,
		ref = config['ref_hg38']
	output:
		cuteSV_vcf = 'svcalling/cuteSV/{sample}_{aln_type}_cuteSV.vcf'
	params:
		wd = 'svcalling/cuteSV/{sample}_{aln_type}/'
	resources:
		mem=10,
		hrs=24
	threads: 1 
	shell:
		"""
		mkdir {params.wd}

		cuteSV {input.merged_bam} {input.ref} {output.cuteSV_vcf} {params.wd} --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3
		"""

