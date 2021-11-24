def find_input_index(wildcards):
    FOFN = manifest_df.at[wildcards.sample, "FOFN"]
    fofn_df = pd.read_csv(FOFN, header=None, names=["FILE"], sep="\t")
    return fofn_df.at[int(wildcards.run), "FILE"] + ".fai"


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
        "tmp/alignments/{{sample}}/{chrom}/{{bc_vers}}/{{seq}}/merge_output.vcf.gz",
        chrom=chroms,
    )


def concatenate_fastq(wildcards):
    FOFN = manifest_df.at[wildcards.sample, "FOFN"]
    fofn_df = pd.read_csv(FOFN, header=None, names=["FILE"], sep="\t")
    return fofn_df["FILE"]


def find_fast5_dir(wildcards):
    df = pd.read_csv(
        manifest_df.at[wildcards.sample, "FAST5_FOFN"],
        sep="\t",
        header=None,
        names=["files"],
    )
    return "-d " + " -d ".join(df["files"])


def find_summary_fofn(wildcards):
    if os.path.isfile(manifest_df.at[wildcards.sample, "SUMMARY"]):
        return "-f " + manifest_df.at[wildcards.sample, "SUMMARY"]
    else:
        return ""


def find_all_windows(wildcards):
    if config.get("METHYL_WINDOW_FILE") != None:
        with open(config.get("WINDOW_FILE"), "r") as infile:
            windows = [x.rstrip() for x in infile]
        return expand(
            rules.nanopolish.output.methyl,
            window=windows,
            sample=wildcards.sample,
            bc_vers=wildcards.bc_vers,
            seq=wildcards.seq,
        )
    elif config.get("METHYL_CLAIR_CHRS") != None:
        windows = config.get("CHRS")
        return expand(
            rules.nanopolish.output.methyl,
            window=windows,
            sample=wildcards.sample,
            bc_vers=wildcards.bc_vers,
            seq=wildcards.seq,
        )
    else:
        with open(f"{REF}.fai", "r") as infile:
            windows = [x.split("\t")[0] for x in infile]
        return expand(
            rules.nanopolish.output.methyl,
            window=windows,
            sample=wildcards.sample,
            bc_vers=wildcards.bc_vers,
            seq=wildcards.seq,
        )
