#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: William T. Harvey

if __name__ == "__main__":
    NIDS = len(snakemake.output.batches)

    batch_dict = {}

    for i in range(NIDS):
        batch_dict[i] = []

    with open(snakemake.input.fai, "r") as infile:
        fai_list = [line.split("\t")[0] for line in infile]

    for j in len(fai_list):
        batch_dict[j % NIDS].append(fai_list[j])

    outs = [open(f, "w+") for f in snakemake.output.batches]

    for i in range(NIDS):
        outs[i].write("\n".join(batch_dict[i]) + "\n")
        outs[i].close()
