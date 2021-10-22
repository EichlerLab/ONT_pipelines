import os
import pandas as pd
import numpy as np
from pathlib import Path

import pdb;
tr=pdb.set_trace

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
shell.prefix("source %s/env.cfg; " % SNAKEMAKE_DIR)
configfile: "alignment.snake_config.yaml"


manifest_df = pd.read_csv(config['manifest'], sep='\t', index_col='SAMPLE')

ALN_TYPES=[config['aln_types']][0].split(",")


def find_input_file(wildcards):
	FOFN = manifest_df.at[wildcards.sample, 'FOFN']
	fofn_df = pd.read_csv(FOFN, header=None, names=['FILE'], sep='\t')
	return fofn_df.at[int(wildcards.run), 'FILE']

def find_parts(wildcards):
	fofn_df = pd.read_csv(manifest_df.at[wildcards.sample, 'FOFN'], header=None, names=['FILE'])
	return expand('alignments/{aln_type}/{run}_{sample}_{aln_type}_alignment.bam', aln_type=ALN_TYPES, run=fofn_df.index, sample=wildcards.sample)


wildcard_constraints:
    sample='|'.join(manifest_df.index) 


scattergather:
    split=100

localrules: all


rule all:
	input:
		expand('alignments/{aln_type}/{sample}_{aln_type}_alignment.bam.bai', aln_type=ALN_TYPES, sample=manifest_df.index.values)
		
		


rule fastq_unzip:
	input:
		fastq_gz = find_input_file
	output:
		fastq = 'fastqfiles/{run}_{sample}.fastq'
	resources:
		mem=2,
		hrs=5,
		disk_free = 1
	threads: 1 
	run:
		if input.fastq_gz.endswith("gz"):
			shell("gunzip < {input.fastq_gz} > {output.fastq}")
		else:
			shell("ln -s $(readlink -f {input.fastq_gz}) {output.fastq}")



rule fastq_index:
	input:
		fastq = rules.fastq_unzip.output.fastq
	output:
		fai = "fastqfiles/{run}_{sample}.fastq.fai"
	resources:
		mem=2,
		hrs=5,
		disk_free = 1
	threads: 1 
	shell:
		"""
		samtools fqidx {input.fastq} --output {output.fai}
		"""


rule get_batchesIDs:
	input:
		fai = "fastqfiles/{run}_{sample}.fastq.fai"
	output:
		batches_IDsfiles = temp(scatter.split('tempFiles/splitBatchID/{{run}}_{{sample}}/{{run}}_{{sample}}_{scatteritem}.txt'))
	params:
		nbatches = config['nbatches']
	resources:
		mem=4,
		hrs=5,
		disk_free = 5
	threads: 1 
	shell:
		"""
		Rscript {SNAKEMAKE_DIR}/scripts/indexSplit.R {input.fai} {wildcards.run} {wildcards.sample} {params.nbatches}
		"""



rule minimap_aln:
	input:
		fastq = 'fastqfiles/{run}_{sample}.fastq', 
		batches_IDsfiles = 'tempFiles/splitBatchID/{run}_{sample}/{run}_{sample}_{scatteritem}.txt',
		ref = config['ref_hg38']
	output:
		sorted_bam = temp('alignments/minimap/{run}_{sample}_batches/{run}_{sample}_{scatteritem}_minimap_alignment.sorted.bam')
	resources:
		mem=4,
		hrs=96,
		disk_free = 5
	threads: 4 
	shell:
		"""
		samtools fqidx {input.fastq} -r {input.batches_IDsfiles} >> $TMPDIR/scatteritem.bam
		
		minimap2 -t {threads} --MD --secondary=no --eqx -x map-ont -a {input.ref} $TMPDIR/scatteritem.bam | samtools view -S -b | samtools sort -T $TMPDIR > {output.sorted_bam}
		"""

rule nglmr_aln:
	input:
		fastq = 'fastqfiles/{run}_{sample}.fastq', 
		batches_IDsfiles = 'tempFiles/splitBatchID/{run}_{sample}/{run}_{sample}_{scatteritem}.txt',
		ref = config['ref_hg38']
	output:
		sorted_bam = temp('alignments/ngmlr/{run}_{sample}_batches/{run}_{sample}_{scatteritem}_ngmlr_alignment.sorted.bam')
	resources:
		mem=4,
		hrs=96,
		disk_free = 5
	threads: 4 
	shell:
		"""
		samtools fqidx {input.fastq} -r {input.batches_IDsfiles} >> $TMPDIR/scatteritem.bam

		ngmlr -t {threads} -r {input.ref} -q $TMPDIR/scatteritem.bam -x ont | samtools view -S -b | samtools sort -T $TMPDIR > {output.sorted_bam}
		"""




rule merge_scatter_aln:
	input:
		sorted_bams = gather.split('alignments/{{aln_type}}/{{run}}_{{sample}}_batches/{{run}}_{{sample}}_{scatteritem}_{{aln_type}}_alignment.sorted.bam')
	output:
		scatter_merged_bam = temp('alignments/{aln_type}/{run}_{sample}_{aln_type}_alignment.bam')
	resources:
		mem=4,
		hrs=24,
		disk_free = 1
	threads: 12 
	shell:
		"""
		samtools merge -@{threads} {output.scatter_merged_bam} {input.sorted_bams}
		"""


rule merge_run_aln:
	input:
		scatter_merged_bam = find_parts
	output:
		merged_bam = 'alignments/{aln_type}/{sample}_{aln_type}_alignment.bam'
	resources:
		mem=4,
		hrs=24,
		disk_free = 1
	threads: 12
	shell:
		"""
		samtools merge -@{threads} {output.merged_bam} {input.scatter_merged_bam}
		"""


rule index_aln:
	input:
		merged_bam = 'alignments/{aln_type}/{sample}_{aln_type}_alignment.bam'
	output:
		merged_bai = 'alignments/{aln_type}/{sample}_{aln_type}_alignment.bam.bai'
	resources:
		mem=4,
		hrs=24,
		disk_free = 1
	threads: 1
	shell:
		"""
		samtools index -b {input.merged_bam}
		"""

