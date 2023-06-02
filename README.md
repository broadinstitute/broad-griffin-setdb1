# Transposable Elements Pipeline Analysis

## Introduction

This repository contains a collection of pipelines and scripts to analyze RNA-seq and ChIP-seq data taking into account multi-mapping reads for accurate transposable elements quantification.

### RNA
[TEtranscript toolkit](https://github.com/mhammell-laboratory/TEtranscripts) is used for multi-mappers-awar family-level quantifications and for differential analyses.

[TElocal](https://github.com/mhammell-laboratory/TElocal) is used for multi-mappers-awar loci-level quantifications.

### Histone modifcations
[SmartMap](https://github.com/shah-rohan/SmartMap) is used for multi-mappers-awar quantifications.


## Getting started

The first step is to download this repository to your computer using the following commands:

```bash
$ git clone git@github.com:broadinstitute/broad-epi-repeats-analysis.git
$ cd broad-epi-repeats-analysis
```

## Genome, genes, and repeats annotations

The script in `src/bash/create-mus-musculus-annotations-mm10.bash` will download and create all the necessary annotations.
The annotations will be place in a new folder `mm10` with three sub-folders called `genome`, `genes`, and `repeats`.

## Genome index

Using the genome FASTA file insided `mm10/genome` you can build the index for *STAR* and *bowtie2* using the `star-build-index` wdl  or the `bowtie2-build-index` wdl respectively. Both are located in the `workflows` folder.

## Quantifications

`star-repeats-quantification` workflow will align the FASTQs to the genome using *STAR* and compute genes and TEs quantifications. Four count files will be reported:
 - **family-level-unique-counts**
 - **family-level-multimappers-counts**
 - **loci-level-unique-counts**
 - **loci-level-multimappers-counts**.

`bowtie2-repeats-quantification` workflow will align the FASTQs to the genome using *bowtie2* and a weighted *bedgraph* accounting for multi-mappers will be generated.


TODO:
- Report unique-mappers track for ChIP
- Create tracks for RNA using unique-mappers and apportioning multi-mappers. It will not be the same counts as in the count files but it will give us an idea.
- Annotate each loci with the percentage of gained counts when including multi-mappers. If unique counts were 3 and with multi-mappers is 100 it is suspicious compared to something going from 50 to 100.

