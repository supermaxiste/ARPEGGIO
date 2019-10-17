# This file includes all rules related to checking the quality of the data
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: FastQC,

## Run FastQC on raw untrimmed reads for quality control

rule quality_control:
	input:
		fastq = RAW_DATA_DIR + "{sample}." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		OUTPUT_DIR + "FastQC/{sample}_fastqc.zip"
	params:
		FastQC = OUTPUT_DIR + "FastQC"
	log:
		OUTPUT_DIR + "logs/fastqc_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/qc_{sample}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		" fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## Run FastQC on SE trimmed reads resulting from trim_galore

rule quality_control_trimmed_SE:
	input:
		fastq = OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed.fq.gz"
	output:
		OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed_fastqc.zip"
	params:
		FastQC = OUTPUT_DIR + "FASTQtrimmed/"
	log:
		OUTPUT_DIR + "logs/fastqc_trimmed_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/qc_trim_se_{sample}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"


## Run FastQC on PE trimmed reads resulting from trim_galore

rule quality_control_trimmed_PE:
	input:
		OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz",
		OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz"
	output:
		OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1_fastqc.zip",
		OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2_fastqc.zip"
	params:
		FastQC = OUTPUT_DIR + "FASTQtrimmed/"
	log:
		OUTPUT_DIR + "logs/fastqc_trimmed_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/qc_trim_pe_{sample}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input}"

## Run Bismark to prepare synthetic bisulfite converted control genome to check conversion efficiency

rule bismark_prepare_control:
	input:
		control = config["CONTROL_GENOME"]
	output:
		CONTROL_GENOME + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"
	benchmark:
		OUTPUT_DIR + "benchmark/prepare_control_genome.txt"
	conda:
		"../envs/environment.yaml"
	shell:
		"bismark_genome_preparation {input.control}"

## Run Bismark to perform alignment to the control genome to check conversion efficiency if reads are single-end.

rule bismark_alignment_SE_control:
	input:
		CONTROL_GENOME + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		fastq = OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		sample = OUTPUT_DIR + "Conversion_efficiency/{sample}/{sample}_trimmed_bismark_bt2.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Conversion_efficiency/{sample}/{sample}_bismark_bt2.bam",
		report = OUTPUT_DIR + "Conversion_efficiency/{sample}/{sample}_trimmed_bismark_bt2_SE_report.txt" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Conversion_efficiency/{sample}/{sample}_bismark_bt2_SE_report.txt"
	params:
		output = OUTPUT_DIR + "Conversion_efficiency/{sample}/",
		control = CONTROL_GENOME
	log:
		OUTPUT_DIR + "logs/Conversion_efficiency_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/Conversion_efficiency_{sample}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
		"bismark --multicore {BISMARK_CORES}  -o {params.output} --temp_dir {params.output} --genome {params.control} {input.fastq}"

rule bismark_alignment_PE_control:
	input:
		CONTROL_GENOME + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		fastq1 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_1"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz",
		fastq2 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_2"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		sample = OUTPUT_DIR + "Conversion_efficiency/{sample}/{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Conversion_efficiency/{sample}/{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.bam",
		report = OUTPUT_DIR + "Conversion_efficiency/{sample}/{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_PE_report.txt" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Conversion_efficiency/{sample}/{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_PE_report.txt"
	params:
		output = OUTPUT_DIR + "Conversion_efficiency/{sample}/",
		control = CONTROL_GENOME
	log:
		OUTPUT_DIR + "logs/Conversion_efficiency_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/Conversion_efficiency_{sample}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
		"bismark --multicore {BISMARK_CORES} --genome {params.control} -1 {input.fastq1} -2 {input.fastq2} -o {params.output} --temp_dir {params.output}"

# sort bam files for multiqc (parent1)

rule bam_sorting_p1:
	input:
		p1 = OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bam"
	output:
		o1 = OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated_sorted.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated_sorted.bam"
	conda:
		"../envs/environment.yaml"
	shell:
		"samtools sort {input.p1} > {output.o1}"

# sort bam files for multiqc (parent2)

rule bam_sorting_p2:
	input:
		p2 = OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bam"
	output:
		o2 = OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated_sorted.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated_sorted.bam"
	conda:
		"../envs/environment.yaml"
	shell:
		"samtools sort {input.p2} > {output.o2}"

# sort bam files for multiqc (allopolyploid)

rule bam_sorting_allo:
	input:
		allo1 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] else OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified1.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bam",
		allo2 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] else OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified2.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bam"
	output:
		o3 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1_sorted.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated_sorted_allo.bam" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated_sorted_allo.bam" if config["IS_PAIRED"] else OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified1_sorted.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated_sorted_allo.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated_sorted_allo.bam",
		o4 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2_sorted.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated_sorted_allo.bam" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated_sorted_allo.bam" if config["IS_PAIRED"] else OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified2_sorted.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated_sorted_allo.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated_sorted_allo.bam"
	conda:
		"../envs/environment.yaml"
	shell:
		"samtools sort {input.allo1} > {output.o3};"
		"samtools sort {input.allo2} > {output.o4}"

## Qualimap rule to get statistics about bam files for parent1

rule qualimap_p1:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated_sorted.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated_sorted.bam")
	output:
		OUTPUT_DIR + "qualimap/{sample}_p1/qualimapReport.html"
	benchmark:
		OUTPUT_DIR + "benchmark/qualimap_p1_{sample}.txt"
	params:
		output = OUTPUT_DIR + "qualimap/{sample}_p1"
	conda:
		"../envs/environment2.yaml"
	threads:
		CORES
	shell:
		"qualimap bamqc -bam {input} -outdir {params.output} -nt {CORES} --java-mem-size=4G"

## Qualimap rule to get statistics about bam files for parent2

rule qualimap_p2:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] (OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated_sorted.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated_sorted.bam")
	output:
		OUTPUT_DIR + "qualimap/{sample}_p2/qualimapReport.html"
	benchmark:
		OUTPUT_DIR + "benchmark/qualimap_p2_{sample}.txt"
	params:
		output = OUTPUT_DIR + "qualimap/{sample}_p2"
	conda:
		"../envs/environment2.yaml"
	threads:
		CORES
	shell:
		"qualimap bamqc -bam {input} -outdir {params.output} -nt {CORES} --java-mem-size=4G"

## Qualimap rule to get statistics about bam files for allopolyploid SE reads

rule qualimap_allo_se:
	input:
		genome1 = OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified1_sorted.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated_sorted_allo.bam",
		genome2 = OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified2_sorted.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated_sorted_allo.bam"
	output:
		OUTPUT_DIR + "qualimap/{sample}_allo_se_1/qualimapReport.html",
		OUTPUT_DIR + "qualimap/{sample}_allo_se_2/qualimapReport.html"
	benchmark:
		OUTPUT_DIR + "benchmark/qualimap_allo_se_{sample}.txt"
	params:
		output1 = OUTPUT_DIR + "qualimap/{sample}_allo_se_1/",
		output2 = OUTPUT_DIR + "qualimap/{sample}_allo_se_2/"
	conda:
		"../envs/environment2.yaml"
	threads:
		CORES
	shell:
		"qualimap bamqc -bam {input.genome1} -outdir {params.output1} -nt {CORES} --java-mem-size=4G;"
		"qualimap bamqc -bam {input.genome2} -outdir {params.output2} -nt {CORES} --java-mem-size=4G"

## Qualimap rule to get statistics about bam files for allopolyploid PE reads

rule qualimap_allo_pe:
	input:
		genome1 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1_sorted.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam",
		genome2 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2_sorted.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "qualimap/{sample}_allo_pe_1/qualimapReport.html",
		OUTPUT_DIR + "qualimap/{sample}_allo_pe_2/qualimapReport.html"
	benchmark:
		OUTPUT_DIR + "benchmark/qualimap_allo_pe_{sample}.txt"
	params:
		output1 = OUTPUT_DIR + "qualimap/{sample}_allo_pe_1/",
		output2 = OUTPUT_DIR + "qualimap/{sample}_allo_pe_2/"
	conda:
		"../envs/environment2.yaml"
	threads:
		CORES
	shell:
		"qualimap bamqc -bam {input.genome1} -outdir {params.output1} -nt {CORES} --java-mem-size=4G;"
		"qualimap bamqc -bam {input.genome2} -outdir {params.output2} -nt {CORES} --java-mem-size=4G"

## Run MultiQC to combine all the outputs from QC, trimming and alignment in a single nice report

rule multiqc:
	input:
		multiqc_input
	output:
		OUTPUT_DIR + "MultiQC/multiqc_report.html"
	params:
		inputdir = multiqc_params,
		multiqcdir = OUTPUT_DIR + "MultiQC"
	log:
		OUTPUT_DIR + "logs/multiqc.log"
	benchmark:
		OUTPUT_DIR + "benchmark/multiqc.txt"
	conda:
		"../envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.inputdir} -f -o {params.multiqcdir}"
