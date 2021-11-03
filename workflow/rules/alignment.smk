
wildcard_constraints:
    sample='|'.join(manifest_df.index) 


scattergather:
    split=NBATCHES

		
rule get_batch_ids:
	input:
		fai = find_input_index
	output:
		batches = temp(scatter.split('tmp/splitBatchID/{{run}}_{{sample}}/{{run}}_{{sample}}_{scatteritem}.txt'))
	resources:
		mem=4,
		hrs=5,
		disk_free = 5
	log: 'log/{run}_{sample}_batching.log'
	threads: 1 
	run:
		fai_df = pd.read_csv(input.fai, sep='\t', header=None, names=['contig', 'len', 'byte_start', 'byte', 'offset'])
		fai_df['batch'] = fai_df.index % NBATCHES
		for i in range(1,NBATCHES+1):
			with open(f'tmp/splitBatchID/{wildcards.run}_{wildcards.sample}/{wildcards.run}_{wildcards.sample}_{i}-of-{NBATCHES}.txt') as outfile:
				outfile.write('\n'.join(fai_df.loc[fai_df['batch'] == i-1]['contig']))
				outfile.write('\n')



rule minimap_aln:
	input:
		fastq = find_input_file, 
		batch_file = 'tmp/splitBatchID/{run}_{sample}/{run}_{sample}_{scatteritem}.txt',
		ref = REF
	output:
		sorted_bam = temp('tmp/alignments/minimap2/{run}_{sample}_batches/{run}_{sample}_{scatteritem}.sorted.bam')
	resources:
		mem=4,
		hrs=96,
		disk_free = 5
	log: 'log/{run}_{sample}_{scatteritem}_aln.log'
	conda:
		'../envs/align.yaml'
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
		'minimap2/2.21'
	threads: 4 
	shell:
		'''
		samtools fqidx {input.fastq} -r {input.batch_file} > $TMPDIR/scatteritem.fastq
		minimap2 -t {threads} --MD --secondary=no --eqx -x map-ont -a {input.ref} $TMPDIR/scatteritem.fastq | samtools view -S -b | samtools sort -T $TMPDIR > {output.sorted_bam}
		'''

rule merge_scatter_aln:
	input:
		sorted_bams = gather.split('tmp/alignments/minimap2/{{run}}_{{sample}}_batches/{{run}}_{{sample}}_{scatteritem}.sorted.bam')
	output:
		scatter_merged_bam = temp('tmp/alignments/minimap2/{run}_{sample}_minimap2_alignment.bam')
	resources:
		mem=4,
		hrs=24,
		disk_free = 1
	threads: 12
	log: 'log/{run}_{sample}_merging.log'
	conda:
		'../envs/align.yaml'
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
		'samtools/1.12'
	shell:
		'''
		samtools merge -@{threads} {output.scatter_merged_bam} {input.sorted_bams}
		'''


rule merge_run_aln:
	input:
		scatter_merged_bam = find_parts
	output:
		merged_bam = 'alignments/{sample}/{sample}.{bc_vers}.minimap2.bam'
	resources:
		mem=4,
		hrs=24,
		disk_free = 1
	threads: 12
	log: 'log/{sample}_{bc_vers}.merge_all.log'
	conda:
		'../envs/align.yaml'
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
		'samtools/1.12'
	shell:
		'''
		samtools merge -@{threads} {output.merged_bam} {input.scatter_merged_bam}
		'''


rule index_aln:
	input:
		merged_bam = rules.merge_run_aln.output.merged_bam
	output:
		merged_bai = 'alignments/{sample}/{sample}.{bc_vers}.minimap2.bam.bai'
	resources:
		mem=4,
		hrs=24,
		disk_free = 1
	threads: 1
	log: 'log/{sample}_{bc_vers}.merge_all.log'	
	conda:
		'../envs/align.yaml'
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
		'samtools/1.12'
	shell:
		'''
		samtools index {input.merged_bam}
		'''

