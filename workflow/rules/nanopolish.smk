

rule concat_fastq:
    input:
        fastq=concatenate_fastq,
    output:
        fastq=temp("tmp/reads/{sample}.{bc_vers}.{seq}.fastq.gz"),
    threads: 1
    conda:
        "../envs/nanopolish.yaml"
    log:
        "log/{sample}_{bc_vers}_{seq}.concat_fastq.log"
    resources:
        mem=8,
        hrs=12,
    shell:
        """
        cat {input.fastq} > {output.fastq}
        """


rule index_fastq:
    input:
        fastq=rules.concat_fastq.output.fastq,
    output:
        index=temp("tmp/reads/{sample}.{bc_vers}.{seq}.fastq.gz.index"),
        fai=temp("tmp/reads/{sample}.{bc_vers}.{seq}.fastq.gz.index.fai"),
        gzi=temp("tmp/reads/{sample}.{bc_vers}.{seq}.fastq.gz.index.gzi"),
        readdb=temp("tmp/reads/{sample}.{bc_vers}.{seq}.fastq.gz.index.readdb"),
    params:
        directory=find_fast5_dir,
        summary=find_summary_fofn
    log:
        "log/{sample}_{bc_vers}_{seq}.index_fastq.log"
    threads: 1
    conda:
        "../envs/nanopolish.yaml"
    resources:
        mem=8,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "nanopolish/0.13.3",
    shell:
        """
        nanopolish index {params.directory} {params.summary} {input.fastq}
        """


rule nanopolish:
    input:
        fastq=rules.concat_fastq.output.fastq,
        sorted_bam=rules.merge_run_aln.output.merged_bam,
        bam_index=rules.index_aln.output.merged_bai,
        reference=REF,
        index=rules.index_fastq.output.index,
        fai=rules.index_fastq.output.fai,
        gzi=rules.index_fastq.output.gzi,
        readdb=rules.index_fastq.output.readdb,
    output:
        methyl="calls/{sample}/{window}.{bc_vers}.nanopolish.{seq}.tsv",
    threads: 8
    conda:
        "../envs/nanopolish.yaml"
    log:
        "log/{sample}_{bc_vers}_{seq}.nanopolish.log"
    resources:
        mem=2,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "nanopolish/0.13.3",
    shell:
        """
        nanopolish call-methylation -t {threads} -r {input.fastq} -b {input.sorted_bam} -g {input.reference} -w {wildcards.window} > {output.methyl}
        """


rule gather_nanopolish:
    input:
        find_all_windows,
    log:
        "log/{sample}_{bc_vers}_{seq}.nanopolish.log"
    output:
        touch(".{sample}.{bc_vers}.{seq}_nanopolish.done"),


# rule merge:
# 	input:
# 		shared = expand("calls/{{sample}}/{window}.share.tsv", window=WINDOWS),
# 		uniq = expand("calls/{{sample}}/{window}.uniq.tsv", window=WINDOWS)
# 	output:
# 		shared = "final/{sample}.share.tsv",
# 		uniq = "final/{sample}.uniq.tsv"
# 	params:
# 		sge_opts = "-l mfree=4G"
# 	run:
# 		for i, file in enumerate(input.shared):
# 			if i == 0:
# 				share_df = pd.read_csv(file, sep="\t", header=0, dtype=str)
# 			else:
# 				share_df = share_df.append(pd.read_csv(file, sep="\t", header=0, dtype=str))

# 		share_df.to_csv(output.shared, sep="\t", header=True, index=False)

# 		for i, file in enumerate(input.uniq):
# 			if i == 0:
# 				uniq_df = pd.read_csv(file, sep="\t", header=0, dtype=str)
# 			else:
# 				uniq_df = uniq_df.append(pd.read_csv(file, sep="\t", header=0, dtype=str))

# 		uniq_df.to_csv(output.uniq, sep="\t", header=True, index=False)

# rule compare:
# 	input:
# 		father = expand("calls/{sample}/{{window}}.methyl.compiled.tsv", sample="14455.fa"),
# 		mother = expand("calls/{sample}/{{window}}.methyl.compiled.tsv", sample="14455.mo"),
# 		sib = expand("calls/{sample}/{{window}}.methyl.compiled.tsv", sample="14455.s1"),
# 		pro = expand("calls/{sample}/{{window}}.methyl.compiled.tsv", sample="14455.p1")
# 	output:
# 		father = expand("calls/{sample}/{{window}}.uniq.tsv", sample="14455.fa"),
# 		mother = expand("calls/{sample}/{{window}}.uniq.tsv", sample="14455.mo"),
# 		sib = expand("calls/{sample}/{{window}}.uniq.tsv", sample="14455.s1"),
# 		pro = expand("calls/{sample}/{{window}}.uniq.tsv", sample="14455.p1"),
# 		father_shared = expand("calls/{sample}/{{window}}.share.tsv", sample="14455.fa"),
# 		mother_shared = expand("calls/{sample}/{{window}}.share.tsv", sample="14455.mo"),
# 		sib_shared = expand("calls/{sample}/{{window}}.share.tsv", sample="14455.s1"),
# 		pro_shared = expand("calls/{sample}/{{window}}.share.tsv", sample="14455.p1")
# 	params:
# 		sge_opts = "-l mfree=16G"
# 	run:
# 		fa_df = pd.read_csv(input.father[0], sep="\t", header=0, index_col=[0,1,2])
# 		mo_df = pd.read_csv(input.mother[0], sep="\t", header=0, index_col=[0,1,2])
# 		sib_df = pd.read_csv(input.sib[0], sep="\t", header=0, index_col=[0,1,2])
# 		pro_df = pd.read_csv(input.pro[0], sep="\t", header=0, index_col=[0,1,2])

# 		fa_df["shared_mo"] = "False"
# 		fa_df["shared_sib"] = "False"
# 		fa_df["shared_pro"] = "False"

# 		mo_df["shared_fa"] = "False"
# 		mo_df["shared_sib"] = "False"
# 		mo_df["shared_pro"] = "False"

# 		sib_df["shared_fa"] = "False"
# 		sib_df["shared_mo"] = "False"
# 		sib_df["shared_pro"] = "False"

# 		pro_df["shared_mo"] = "False"
# 		pro_df["shared_sib"] = "False"
# 		pro_df["shared_fa"] = "False"

# 		for index in fa_df.index:
# 			if index in mo_df.index:
# 				fa_df.at[index, "shared_mo"] = "True"
# 			if index in sib_df.index:
# 				fa_df.at[index, "shared_sib"] = "True"
# 			if index in pro_df.index:
# 				fa_df.at[index, "shared_pro"] = "True"

# 		fa_df.to_csv(output.father_shared[0], sep="\t", header=True, index=True)
# 		fa_uniq = fa_df[(fa_df["shared_mo"] == "False") & (fa_df["shared_sib"] == "False") & (fa_df["shared_pro"] == "False")]
# 		fa_uniq.to_csv(output.father[0], sep="\t", header=True, index=True)

# 		for index in mo_df.index:
# 			if index in fa_df.index:
# 				mo_df.at[index, "shared_fa"] = "True"
# 			if index in sib_df.index:
# 				mo_df.at[index, "shared_sib"] = "True"
# 			if index in pro_df.index:
# 				mo_df.at[index, "shared_pro"] = "True"

# 		mo_df.to_csv(output.mother_shared[0], sep="\t", header=True, index=True)
# 		mo_uniq = mo_df[(mo_df["shared_fa"] == "False") & (mo_df["shared_sib"] == "False") & (mo_df["shared_pro"] == "False")]
# 		mo_uniq.to_csv(output.mother[0], sep="\t", header=True, index=True)

# 		for index in sib_df.index:
# 			if index in fa_df.index:
# 				sib_df.at[index, "shared_fa"] = "True"
# 			if index in mo_df.index:
# 				sib_df.at[index, "shared_mo"] = "True"
# 			if index in pro_df.index:
# 				sib_df.at[index, "shared_pro"] = "True"

# 		sib_df.to_csv(output.sib_shared[0], sep="\t", header=True, index=True)
# 		sib_uniq = sib_df[(sib_df["shared_mo"] == "False") & (sib_df["shared_fa"] == "False") & (sib_df["shared_pro"] == "False")]
# 		sib_uniq.to_csv(output.sib[0], sep="\t", header=True, index=True)

# 		for index in pro_df.index:
# 			if index in fa_df.index:
# 				pro_df.at[index, "shared_fa"] = "True"
# 			if index in mo_df.index:
# 				pro_df.at[index, "shared_mo"] = "True"
# 			if index in sib_df.index:
# 				pro_df.at[index, "shared_sib"] = "True"

# 		pro_df.to_csv(output.pro_shared[0], sep="\t", header=True, index=True)
# 		pro_uniq = pro_df[(pro_df["shared_mo"] == "False") & (pro_df["shared_fa"] == "False") & (pro_df["shared_sib"] == "False")]
# 		pro_uniq.to_csv(output.pro[0], sep="\t", header=True, index=True)


# rule calc_frequency:
# 	input:
# 		callset = "calls/{sample}/{window}.hap.methyl.tsv"
# 	output:
# 		compiled = "calls/{sample}/{window}.hap.methyl.compiled.tsv"
# 	params:
# 		sge_opts = "-l mfree=8G"
# 	shell:
# 		"""
# 		module load miniconda/4.5.12
# 		/net/eichler/vol26/7200/software/modules-sw/nanopolish/0.11.1/Linux/RHEL6/x86_64/scripts/calculate_methylation_frequency.py {input.callset} > {output.compiled}
# 		"""

# rule resolve_haplotypes:
# 	input:
# 		hapOne = findHapOne,
# 		hapTwo = findHapTwo,
# 		callset = "calls/{sample}/{window}.methylation.tsv"
# 	output:
# 		haplotypes = "calls/{sample}/{window}.hap.methyl.tsv"
# 	params:
# 		sge_opts = "-l mfree=16G"
# 	run:
# 		with open(input.hapOne, "r") as oneFile:
# 			h_one = [ line.split("\t")[0] for line in oneFile ]

# 		with open(input.hapTwo, "r") as twoFile:
# 			h_two = [ line.split("\t")[0] for line in twoFile ]

# 		call_df = pd.read_csv(input.callset, header=0, sep="\t")

# 		call_df["haplotype"] = "NA"

# 		hap_one_df = call_df[call_df["read_name"].isin(h_one)]
# 		hap_one_df["haplotype"] = "H1"
# 		hap_two_df = call_df[call_df["read_name"].isin(h_two)]
# 		hap_two_df["haplotype"] = "H2"
# 		hap_combined = hap_one_df.append(hap_two_df)
# 		hap_combined.to_csv(output.haplotypes, sep="\t", header=True, index=False)
# rule map_sort:
# 	input:
# 		fastq = findFastq,
# 		reference = REF,
# 		flag = "{sample}.prep"
# 	output:
# 		sorted_bam = "bam/{sample}.sorted.bam",
# 		bai = "bam/{sample}.sorted.bam.bai"
# 	params:
# 		sge_opts = "-pe serial 8 -l mfree=8G"
# 	shell:
# 		"""
# 		module load htslib/latest gcc/latest nanopolish/0.11.1 minimap2/2.6 samtools/latest
# 		minimap2 -a -x map-ont -t 8 {REF} {input.fastq} | samtools sort -T {wildcards.sample} -o {output.sorted_bam}
# 		samtools index {output.sorted_bam}
# 		"""
