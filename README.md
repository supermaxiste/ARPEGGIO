## PlantEpiSnake

This repository contains a SnakeMake workflow to analyze BS-seq data coming from plants. The workflow goes through several steps: adapter trimming, quality check, mapping, deduplication, methylation calling and statistics for different methylation contexts.

The following tools are used in each step:  

- Adapter trimming: `TrimGalore`
- Quality check: `FastQC`
- Mapping and deduplication: `Bismark`
- Methylation calling: `Bismark`
- Statistics: for CG and CHG context `dmrseq`, for CHH context `DMRcaller`

A final report with read statistics is also included.


