import pandas as pd
import os

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

configfile : "nanopolish.json"

REFERENCE = config['REFERENCE']
MANIFEST = config['MANIFEST']
REGIONS = config['REGIONS']

manifest_df = pd.read_csv(MANIFEST, sep='\t', header=0, index_col=0)

WINDOWS = []

with open(REGIONS, 'r') as bedFile:
	for line in bedFile:
		line = line.rstrip().split('\t')
		WINDOWS.append(line[0]+":"+line[1]+"-"+line[2])


def find_fastq(wildcards):
	df = pd.read_csv(manifest_df.at[wildcards.sample, 'FASTQ_FOFN'], sep='\t', header=None, names=['files'])
	return df['files'].values()

def find_fast5_dir(wildcards):
	df = pd.read_csv(manifest_df.at[wildcards.sample, 'FAST5_FOFN'], sep='\t', header=None, names=['files'])
	return '-d '+' -d '.join(df['files'].values())

def findBam(wildcards):
	return manifest_df.at[wildcards.sample, 'BAM']

def findHapOne(wildcards):
	return manifest_df.at[wildcards.sample, 'HAP_ONE']

def findHapTwo(wildcards):
	return manifest_df.at[wildcards.sample, 'HAP_TWO']


def find_summary_fofn(wildcards):
	if os.path.isfile(manifest_df.at[wildcards.sample, 'SUMMARY']):
		return '-f '+manifest_df.at[wildcards.sample, 'SUMMARY']
	else:
		return ''

rule all:
	input: 
		# expand('final/{sample}.share.tsv', sample=manifest_df.index), 
		expand('calls/{sample}/{window}.hap.methyl.compiled.tsv', sample=manifest_df.index, window=WINDOWS)


rule pre_process:
	input:
		fastq = find_fastq
	output:
		flag = temp(touch('{sample}.prep'))
	params:
		directory = find_fast5_dir,
		summary = find_summary_fofn
	envmodules:
		'nanopolish/0.13.3'
	shell:
		'''
		nanopolish index {params.directory} {params.summary} {input.fastq}
		'''	



# rule merge:
# 	input:
# 		shared = expand('calls/{{sample}}/{window}.share.tsv', window=WINDOWS),
# 		uniq = expand('calls/{{sample}}/{window}.uniq.tsv', window=WINDOWS)
# 	output:
# 		shared = 'final/{sample}.share.tsv',
# 		uniq = 'final/{sample}.uniq.tsv'
# 	params:
# 		sge_opts = '-l mfree=4G'
# 	run:
# 		for i, file in enumerate(input.shared):
# 			if i == 0:
# 				share_df = pd.read_csv(file, sep='\t', header=0, dtype=str)
# 			else:
# 				share_df = share_df.append(pd.read_csv(file, sep='\t', header=0, dtype=str))

# 		share_df.to_csv(output.shared, sep='\t', header=True, index=False)

# 		for i, file in enumerate(input.uniq):
# 			if i == 0:
# 				uniq_df = pd.read_csv(file, sep='\t', header=0, dtype=str)
# 			else:
# 				uniq_df = uniq_df.append(pd.read_csv(file, sep='\t', header=0, dtype=str))

# 		uniq_df.to_csv(output.uniq, sep='\t', header=True, index=False)

# rule compare:
# 	input: 
# 		father = expand('calls/{sample}/{{window}}.methyl.compiled.tsv', sample='14455.fa'),
# 		mother = expand('calls/{sample}/{{window}}.methyl.compiled.tsv', sample='14455.mo'),
# 		sib = expand('calls/{sample}/{{window}}.methyl.compiled.tsv', sample='14455.s1'),
# 		pro = expand('calls/{sample}/{{window}}.methyl.compiled.tsv', sample='14455.p1')
# 	output:
# 		father = expand('calls/{sample}/{{window}}.uniq.tsv', sample='14455.fa'),
# 		mother = expand('calls/{sample}/{{window}}.uniq.tsv', sample='14455.mo'),
# 		sib = expand('calls/{sample}/{{window}}.uniq.tsv', sample='14455.s1'),
# 		pro = expand('calls/{sample}/{{window}}.uniq.tsv', sample='14455.p1'),
# 		father_shared = expand('calls/{sample}/{{window}}.share.tsv', sample='14455.fa'),
# 		mother_shared = expand('calls/{sample}/{{window}}.share.tsv', sample='14455.mo'),
# 		sib_shared = expand('calls/{sample}/{{window}}.share.tsv', sample='14455.s1'),
# 		pro_shared = expand('calls/{sample}/{{window}}.share.tsv', sample='14455.p1')
# 	params:
# 		sge_opts = '-l mfree=16G'
# 	run:
# 		fa_df = pd.read_csv(input.father[0], sep='\t', header=0, index_col=[0,1,2])
# 		mo_df = pd.read_csv(input.mother[0], sep='\t', header=0, index_col=[0,1,2])
# 		sib_df = pd.read_csv(input.sib[0], sep='\t', header=0, index_col=[0,1,2])
# 		pro_df = pd.read_csv(input.pro[0], sep='\t', header=0, index_col=[0,1,2])

# 		fa_df['shared_mo'] = 'False'
# 		fa_df['shared_sib'] = 'False'
# 		fa_df['shared_pro'] = 'False'

# 		mo_df['shared_fa'] = 'False'
# 		mo_df['shared_sib'] = 'False'
# 		mo_df['shared_pro'] = 'False'		

# 		sib_df['shared_fa'] = 'False'
# 		sib_df['shared_mo'] = 'False'
# 		sib_df['shared_pro'] = 'False'		

# 		pro_df['shared_mo'] = 'False'
# 		pro_df['shared_sib'] = 'False'
# 		pro_df['shared_fa'] = 'False'

# 		for index in fa_df.index:
# 			if index in mo_df.index:
# 				fa_df.at[index, 'shared_mo'] = 'True'
# 			if index in sib_df.index:
# 				fa_df.at[index, 'shared_sib'] = 'True'
# 			if index in pro_df.index:
# 				fa_df.at[index, 'shared_pro'] = 'True'

# 		fa_df.to_csv(output.father_shared[0], sep='\t', header=True, index=True)
# 		fa_uniq = fa_df[(fa_df['shared_mo'] == 'False') & (fa_df['shared_sib'] == 'False') & (fa_df['shared_pro'] == 'False')]
# 		fa_uniq.to_csv(output.father[0], sep='\t', header=True, index=True)

# 		for index in mo_df.index:
# 			if index in fa_df.index:
# 				mo_df.at[index, 'shared_fa'] = 'True'
# 			if index in sib_df.index:
# 				mo_df.at[index, 'shared_sib'] = 'True'
# 			if index in pro_df.index:
# 				mo_df.at[index, 'shared_pro'] = 'True'

# 		mo_df.to_csv(output.mother_shared[0], sep='\t', header=True, index=True)
# 		mo_uniq = mo_df[(mo_df['shared_fa'] == 'False') & (mo_df['shared_sib'] == 'False') & (mo_df['shared_pro'] == 'False')]
# 		mo_uniq.to_csv(output.mother[0], sep='\t', header=True, index=True)
		
# 		for index in sib_df.index:
# 			if index in fa_df.index:
# 				sib_df.at[index, 'shared_fa'] = 'True'
# 			if index in mo_df.index:
# 				sib_df.at[index, 'shared_mo'] = 'True'
# 			if index in pro_df.index:
# 				sib_df.at[index, 'shared_pro'] = 'True'

# 		sib_df.to_csv(output.sib_shared[0], sep='\t', header=True, index=True)
# 		sib_uniq = sib_df[(sib_df['shared_mo'] == 'False') & (sib_df['shared_fa'] == 'False') & (sib_df['shared_pro'] == 'False')]
# 		sib_uniq.to_csv(output.sib[0], sep='\t', header=True, index=True)
		
# 		for index in pro_df.index:
# 			if index in fa_df.index:
# 				pro_df.at[index, 'shared_fa'] = 'True'
# 			if index in mo_df.index:
# 				pro_df.at[index, 'shared_mo'] = 'True'
# 			if index in sib_df.index:
# 				pro_df.at[index, 'shared_sib'] = 'True'

# 		pro_df.to_csv(output.pro_shared[0], sep='\t', header=True, index=True)
# 		pro_uniq = pro_df[(pro_df['shared_mo'] == 'False') & (pro_df['shared_fa'] == 'False') & (pro_df['shared_sib'] == 'False')]
# 		pro_uniq.to_csv(output.pro[0], sep='\t', header=True, index=True)




rule calc_frequency:
	input:
		callset = 'calls/{sample}/{window}.hap.methyl.tsv'
	output:
		compiled = 'calls/{sample}/{window}.hap.methyl.compiled.tsv'
	params:
		sge_opts = '-l mfree=8G'
	shell:
		'''
		module load miniconda/4.5.12
		/net/eichler/vol26/7200/software/modules-sw/nanopolish/0.11.1/Linux/RHEL6/x86_64/scripts/calculate_methylation_frequency.py {input.callset} > {output.compiled}
		'''

rule resolve_haplotypes:
	input:
		hapOne = findHapOne,
		hapTwo = findHapTwo,
		callset = 'calls/{sample}/{window}.methylation.tsv'
	output:
		haplotypes = 'calls/{sample}/{window}.hap.methyl.tsv'
	params:
		sge_opts = '-l mfree=16G'
	run:
		with open(input.hapOne, 'r') as oneFile:
			h_one = [ line.split('\t')[0] for line in oneFile ]

		with open(input.hapTwo, 'r') as twoFile:
			h_two = [ line.split('\t')[0] for line in twoFile ]

		call_df = pd.read_csv(input.callset, header=0, sep='\t')

		call_df['haplotype'] = 'NA'

		hap_one_df = call_df[call_df['read_name'].isin(h_one)]
		hap_one_df['haplotype'] = 'H1'
		hap_two_df = call_df[call_df['read_name'].isin(h_two)]
		hap_two_df['haplotype'] = 'H2'

		hap_combined = hap_one_df.append(hap_two_df)
		hap_combined.to_csv(output.haplotypes, sep='\t', header=True, index=False)



rule nanopolish:
	input:
		fastq = findFastq,
		sorted_bam = findBam,
		reference = REFERENCE,
		flag = '{sample}.prep'
	output:
		methyl = 'calls/{sample}/{window}.methylation.tsv'
	params:
		sge_opts = '-l mfree=4G -pe serial 8',
	priority:20
	shell:
		'''
		module load htslib/latest gcc/latest nanopolish/0.11.1
		nanopolish call-methylation -t 8 -r {input.fastq} -b {input.sorted_bam} -g {REFERENCE} -w {wildcards.window} > {output.methyl}
		'''

# rule map_sort:
# 	input: 
# 		fastq = findFastq,
# 		reference = REFERENCE,
# 		flag = '{sample}.prep'
# 	output:
# 		sorted_bam = 'bam/{sample}.sorted.bam',
# 		bai = 'bam/{sample}.sorted.bam.bai'
# 	params:
# 		sge_opts = '-pe serial 8 -l mfree=8G'
# 	shell:
# 		'''
# 		module load htslib/latest gcc/latest nanopolish/0.11.1 minimap2/2.6 samtools/latest
# 		minimap2 -a -x map-ont -t 8 {REFERENCE} {input.fastq} | samtools sort -T {wildcards.sample} -o {output.sorted_bam}
# 		samtools index {output.sorted_bam}
# 		'''
