# This file includes all rules related to read alignment
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: Bismark

## Run Bismark to prepare synthetic bisulfite converted genomes

rule bismark_prepare_genome:
	input:
		genome1 = config["GENOME_PARENT_1"],
		genome2 = config["GENOME_PARENT_2"]
	output:
		GENOME_1 + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		GENOME_2 + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"
	benchmark:
		OUTPUT_DIR + "benchmark/prepare_genome.txt"
	conda:
		"../envs/environment.yaml"
	shell:
		"bismark_genome_preparation {input.genome1}; "
		"bismark_genome_preparation {input.genome2}"

## Run Bismark to perform alignment to the first parental genome (GENOME_PARENT_1) if reads are single-end.

rule bismark_alignment_SE_1:
	input:
		GENOME_1 + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		fastq = OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		sample = OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_trimmed_bismark_bt2.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_bismark_bt2.bam",
		report = OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_trimmed_bismark_bt2_SE_report.txt" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_bismark_bt2_SE_report.txt"
	params:
		output = OUTPUT_DIR + "Bismark/{sample}_1/",
		genome1 = config["GENOME_PARENT_1"],
		prefix = "1"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/bismark_se1_{sample}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
		"bismark --prefix {params.prefix} --multicore {BISMARK_CORES}  -o {params.output} --temp_dir {params.output} --genome {params.genome1} {input.fastq}"

## Run Bismark to perform alignment to the second parental genome (GENOME_PARENT_2) if reads are single-end.

rule bismark_alignment_SE_2:
	input:
		GENOME_2 + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		fastq = OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		sample = OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_trimmed_bismark_bt2.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_bismark_bt2.bam",
		report = OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_trimmed_bismark_bt2_SE_report.txt" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_bismark_bt2_SE_report.txt"
	params:
		output = OUTPUT_DIR + "Bismark/{sample}_2/",
		genome2 = config["GENOME_PARENT_2"],
		prefix = "2"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/bismark_se2_{sample}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
		"bismark --prefix {params.prefix} --multicore {BISMARK_CORES}  -o {params.output} --temp_dir {params.output} --genome {params.genome2} {input.fastq}"

## Run Bismark to perform alignment to the first parental genome (GENOME_PARENT_1) if reads are paired-end.

rule bismark_alignment_PE_1:
	input:
		GENOME_1 + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		fastq1 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_1"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz",
		fastq2 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_2"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		sample = OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.bam",
		report = OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_PE_report.txt"
	params:
		output = OUTPUT_DIR + "Bismark/{sample}_1/",
		genome1 = config["GENOME_PARENT_1"],
		prefix = "1"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/bismark_pe1_{sample}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
		"bismark --prefix {params.prefix} --multicore {BISMARK_CORES} --genome {params.genome1} -1 {input.fastq1} -2 {input.fastq2} -o {params.output} --temp_dir {params.output}"

## Run Bismark to perform alignment to the second parental genome (GENOME_PARENT_2) if reads are paired-end.

rule bismark_alignment_PE_2:
	input:
		GENOME_2 + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		fastq1 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_1"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz",
		fastq2 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_2"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		sample = OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.bam",
		report = OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_PE_report.txt"
	params:
		output = OUTPUT_DIR + "Bismark/{sample}_2/",
		genome2 = config["GENOME_PARENT_2"],
		prefix = "2"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/bismark_pe2_{sample}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
		"bismark --prefix {params.prefix} --multicore {BISMARK_CORES} --genome {params.genome2} -1 {input.fastq1} -2 {input.fastq2} -o {params.output} --temp_dir {params.output}"
