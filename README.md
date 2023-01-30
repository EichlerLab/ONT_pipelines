# ONT Pipelines

[![Actions Status](https://github.com/EichlerLab/ONT_pipelines/workflows/CI/badge.svg)](https://github.com/mrvollger/EichlerLab/ONT_pipelines/actions)
[![Actions Status](https://github.com/EichlerLab/ONT_pipelines/workflows/Linting/badge.svg)](https://github.com/EichlerLab/ONT_pipelines/actions)
[![Actions Status](https://github.com/EichlerLab/ONT_pipelines/workflows/black/badge.svg)](https://github.com/EichlerLab/ONT_pipelines/actions)

This is a Snakemake pipeline designed to take Oxford Nanopore Technologies data from fastq's to variant calls. In additions to traditional SNVs and indels, this pipeline will also call methylation using nanopolish. 

The `Snakefile` is under `workflow`.


The functions performed by this tool are as follows:

## Alignment ##
 - minimap2

## SV Calling ##
 - Sniffles
 - CuteSV
 - SVIM

## SNV/indel Calling ##
 - Clair3

## Read-phasing ##
 - longphase

## Methylation tag linking ##
 - In-house pipeline


There are dummy rules which can be added to the snakemake command which run specific tests. They are listed below:

 - all_align
   - Runs minimap2 for all samples
 - all_sv
   - Runs all SV callers for all samples
 - all_snv
   - Runs Clair3 for all samples
 - all_haplotag
   - Runs haplotype phasing using longphase for all samples
 - all_haplotag_methyl
   - Haplotags bams and links methylation bams for all samples
 - all_minimap_methyl
   - Links methylation to unphased bams for all samples
