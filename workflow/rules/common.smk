def find_input_index(wildcards):
    FOFN = manifest_df.at[wildcards.sample, "FOFN"]
    fofn_df = pd.read_csv(FOFN, header=None, names=["FILE"], sep="\t")
    file = fofn_df.at[int(wildcards.run), "FILE"]
    return f"{file}.fai"


def find_input_file(wildcards):
    FOFN = manifest_df.at[wildcards.sample, "FOFN"]
    fofn_df = pd.read_csv(FOFN, header=None, names=["FILE"], sep="\t")
    return fofn_df.at[int(wildcards.run), "FILE"]


def find_parts(wildcards):
    fofn_df = pd.read_csv(
        manifest_df.at[wildcards.sample, "FOFN"], header=None, names=["FILE"]
    )
    return expand(
        rules.merge_scatter_aln.output.scatter_merged_bam,
        aln_type="minimap2",
        run=fofn_df.index,
        sample=wildcards.sample,
    )


def find_clair_chrs(wildcards):
    if config.get("CHRS") == None:
        with open(f"{REF}.fai", "r") as infile:
            chroms = [line.split("\t")[0] for line in infile]
    else:
        chroms = config.get("CHRS")
    return expand(
        "tmp/variants/{{sample}}/{chrom}/phased_merge_output.vcf.gz",
        chrom=chroms,
    )


def concatenate_fastq(wildcards):
    FOFN = manifest_df.at[wildcards.sample, "FOFN"]
    fofn_df = pd.read_csv(FOFN, header=None, names=["FILE"], sep="\t")
    return fofn_df["FILE"]


def find_aln_bam(wildcards):
    if wildcards.phase == "longphase":
        return rules.long_phase.output.bam
    elif wildcards.phase == "minimap2":
        return rules.merge_run_aln.output.merged_bam
    else:
        print("phase state must be either minimap2 (unphased) or longphase (phased)")
        exit(1)


def find_aln_list(wildcards):
    FOFN = manifest_df.at[wildcards.sample, "METHYL_FOFN"]
    with open(FOFN, "r") as infile:
        return "\n".join([line.rstrip() for line in infile])
