#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: William T. Harvey
import pandas as pd

if __name__ == "__main__":  
    NIDS = len(snakemake.output.batches)

    fai_df = pd.read_csv(
        snakemake.input.fai,
        sep="\t",
        header=None,
        names=["contig", "len", "byte_start", "byte", "byte_len", "offset"]
    )

    fai_df["batch"] = fai_df.index % NIDS

    outs = [open(f, "w+") for f in snakemake.output.batches]

    for i in range(len(outs)):
        outs[i].write("\n".join(fai_df.loc[fai_df["batch"] == i]["contig"]) + "\n")
        outs[i].close()
