[![Build Status](https://travis-ci.com/supermaxiste/ARPEGGIO.svg?token=auqzHDuyLxkuTyyxwdvA&branch=master)](https://travis-ci.com/supermaxiste/ARPEGGIO) \
<img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48"><img src="images/harp.png" height="48">

## ARPEGgIO: Automated Reproducible Polyploid EpiGenomIc wOrkflow


ARPEGGIO is a snakemake workflow that analyzes whole genome bisulfite sequencing (WGBS) data coming from (allo)polyploid species. The workflow includes all basic steps in WGBS data analysis (trimming, quality check and alignment), a read sorting tool specific for allopolyploids, the most comprehensive statistical tool for Differential Methylation (DM) analysis and a set of downstream analyses to obtain a list of genes showing differential methylation.

## Motivation

In the last decade, the use of Next-Generation Sequencing technologies has become widespread across life sciences. With technology not being a bottleneck anymore, the new challenge with NGS data has shifted towards data analysis.
To process and analyze WGBS data, many tools exist, but most of them were developed and/or tested with a focus on model species. For non-model species, especially polyploids, there are complexities that are often not taken into account and workflows are almost non-existent.
To fill in this gap we developed ARPEGGIO: an automated and reproducible workflow for polyploid species.

## Why ARPEGGIO?

ARPEGGIO is based on Snakemake: a human readable, Python based language. Snakemake tries to offer the best trade-off between expert and less-experienced users. At its heart, Snakemake aims at being easily interpretable and adaptable to create workflows that are simple to re-run with new data. ARPEGGIO follows this philosophy and provides a workflow that is easy to adapt for all WGBS coming from allopolyploids.

## What's new in ARPEGGIO?

Besides the workflow itself (which is already quite a lot of new), ARPEGGIO includes an allopolyploid specific read-sorting algorithm that has been adapted to deal with BS-seq data: [EAGLE-RC](https://github.com/tony-kuo/eagle). Check out the papers ["Homeolog expression quantification methods for allopolyploids"](https://doi.org/10.1093/bib/bby121) and ["EAGLE: Explicit Alternative Genome Likelihood Evaluator"](https://doi.org/10.1186/s12920-018-0342-1) by Kuo _et al._ for more details. Together with EAGLE-RC, there's also `dmrseq`: an R package for differential methylation analysis. This package has one of the most comprehensive approaches to deal with WGBS data problems: mainly statistical and computational. If you're curious check out ["Detection and accurate false discovery rate control of differentially methylated regions from whole genome bisulfite sequencing"](https://doi.org/10.1093/biostatistics/kxy007) by Korthauer _et al_.

## Workflow overview

![WorkflowOverview](images/WorkflowV2.png)

## Installation

To install this workflow you first need to [install Snakemake via Conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). Once everything is set up, run the following commands to clone the ARPEGGIO repository to your computer and run the workflow through conda.

```
git clone https://github.com/supermaxiste/ARPEGGIO
cd ARPEGGIO
snakemake --use-conda
```
## Setup and run

Check out the [Wiki](https://github.com/supermaxiste/ARPEGGIO/wiki) to set up and run ARPEGGIO. The Wiki includes also a map of the output with explanations about the content of each folder.

## Troubleshooting and support

Google doesn't help? Are you stuck on an error that no one else seems to be having? Have you checked all the pages mentioning your problem but the solutions are not suitable? On the Wiki there's a list of common problems together with some general solutions. If that didn't work either, feel free to open an issue. Please make sure to describe your problem/errors and your trials in detail so that you can get the best help possible.

## Credits

This project was inspired by ![ARMOR](https://github.com/csoneson/ARMOR), if you work with RNA-seq data check it out!
