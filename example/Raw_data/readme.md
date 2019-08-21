This folder includes reads extracted from alignments coming from real Whole-Genome BS-Seq (WGBS) data. This WGBS data was obtained from 8 individuals: 4 individuals from two different (parental) diploid species and four synthetic allopolyploids. From this data we generated paired-end reads for each individual. For each alignment file, the following commands were used:

```
samtools sort input.bam > input_sorted.bam
samtools index input_sorted.bam
samtools view input_sorted.bam -bh "scaffold_X" > output.bam
samtools fastq -1 output_R1.fastq -2 output_R2.fastq -s single.fastq output.bam
```

For the polyploid alignments we picked scaffold_1 and scaffold_2240 (the first, longest scaffold from each one of the assemblies) and for each one of the parent the corresponding scaffold was picked.
