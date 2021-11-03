#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: William T. Harvey
import argparse
import pysam
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="input fai file")
    parser.add_argument(
        "--outputs", nargs="+", help="list of output files", required=True
    )

    args = parser.parse_args()
    NIDS = len(args.outputs)

    fai_df = pd.read_csv(
        args.infile,
        sep="\t",
        header=None,
        names=["contig", "len", "byte_start", "byte", "offset"],
    )

    fai_df["batch"] = fai_df.index % NBATCHES

    outs = [open(f, "w+") for f in args.outputs]

    for i in range(len(outs)):
        outs[i].write("\n".join(fai_df.loc[fai_df["batch"] == i]["contig"]) + "\n")
        outs[i].close()
