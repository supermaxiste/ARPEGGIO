import os.path

# This file includes all rules related to read alignment
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: Bismark

## Run Bismark to prepare synthetic bisulfite converted genomes

rule bismark_prepare_genome_1:
	output:
		control = f"{GENOME_1}Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"
	log:
		f"logs/bismark_prepare_genome_1.log"
	params:
		genome1 = lambda w, output: os.path.split(os.path.split(os.path.split(output.control)[0])[0])[0]
	benchmark:
		f"{OUTPUT_DIR}benchmark/prepare_genome.txt"
	conda:
		"../envs/environment.yaml"
	shell:
		"bismark_genome_preparation {params.genome1}"

## Run Bismark to prepare synthetic bisulfite converted genomes

rule bismark_prepare_genome_2:
	output:
		control = f"{GENOME_2}Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"
	log:
		f"logs/bismark_prepare_genome_2.log"
	params:
		genome2 = lambda w, output: os.path.split(os.path.split(os.path.split(output.control)[0])[0])[0]
	benchmark:
		f"{OUTPUT_DIR}benchmark/prepare_genome.txt"
	conda:
		"../envs/environment.yaml"
	shell:
		"bismark_genome_preparation {params.genome2}"

## Run Bismark to perform alignment to the first parental genome (GENOME_PARENT_1) if reads are single-end.

rule bismark_alignment_SE_1:
	input:
		control = f"{GENOME_1}Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		fastq = f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_trimmed.fq.gz" if config["RUN_TRIMMING"] else f"{RAW_DATA_DIR}{{sample}}.{str(config['RAW_DATA_EXTENSION'])}.gz"
	output:
		sample = f"{OUTPUT_DIR}Bismark/{{sample}}_1/1.{{sample}}_trimmed_bismark_bt2.bam" if config["RUN_TRIMMING"] else f"{OUTPUT_DIR}Bismark/{{sample}}_1/1.{{sample}}_bismark_bt2.bam",
		report = f"{OUTPUT_DIR}Bismark/{{sample}}_1/1.{{sample}}_trimmed_bismark_bt2_SE_report.txt" if config["RUN_TRIMMING"] else f"{OUTPUT_DIR}Bismark/{{sample}}_1/1.{{sample}}_bismark_bt2_SE_report.txt"
	log:
		f"logs/bismark_alignment_{{sample}}_SE_1.log"
	params:
		output = lambda w, output: os.path.split(output.sample)[0],
		genome1 = lambda w, input: os.path.split(os.path.split(os.path.split(input.control)[0])[0])[0],
		prefix = "1",
		bismark_cores = BISMARK_CORES
	benchmark:
		f"{OUTPUT_DIR}benchmark/bismark_se1_{{sample}}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark --prefix {params.prefix} --multicore {params.bismark_cores}  -o {params.output} --temp_dir {params.output} --genome {params.genome1} {input.fastq}"

## Run Bismark to perform alignment to the second parental genome (GENOME_PARENT_2) if reads are single-end.

rule bismark_alignment_SE_2:
	input:
		control = f"{GENOME_2}Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		fastq = f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_trimmed.fq.gz" if config["RUN_TRIMMING"] else f"{RAW_DATA_DIR}{{sample}}.{str(config['RAW_DATA_EXTENSION'])}.gz"
	output:
		sample = f"{OUTPUT_DIR}Bismark/{{sample}}_2/2.{{sample}}_trimmed_bismark_bt2.bam" if config["RUN_TRIMMING"] else f"{OUTPUT_DIR}Bismark/{{sample}}_2/2.{{sample}}_bismark_bt2.bam",
		report = f"{OUTPUT_DIR}Bismark/{{sample}}_2/2.{{sample}}_trimmed_bismark_bt2_SE_report.txt" if config["RUN_TRIMMING"] else f"{OUTPUT_DIR}Bismark/{{sample}}_2/2.{{sample}}_bismark_bt2_SE_report.txt"
	log:
		f"logs/bismark_alignment_{{sample}}_SE_2.log"
	params:
		output = lambda w, output: os.path.split(output.sample)[0],
		genome2 = lambda w, input: os.path.split(os.path.split(os.path.split(input.control)[0])[0])[0],
		prefix = "2",
		bismark_cores = BISMARK_CORES
	benchmark:
		f"{OUTPUT_DIR}benchmark/bismark_se2_{{sample}}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark --prefix {params.prefix} --multicore {params.bismark_cores}  -o {params.output} --temp_dir {params.output} --genome {params.genome2} {input.fastq}"

## Run Bismark to perform alignment to the first parental genome (GENOME_PARENT_1) if reads are paired-end.

rule bismark_alignment_PE_1:
	input:
		control = f"{GENOME_1}Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		fastq1 = f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_1'])}_val_1.fq.gz" if config["RUN_TRIMMING"] else f"{RAW_DATA_DIR}{{sample}}_{str(config['PAIR_1'])}.{str(config['RAW_DATA_EXTENSION'])}.gz",
		fastq2 = f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_2'])}_val_2.fq.gz" if config["RUN_TRIMMING"] else f"{RAW_DATA_DIR}{{sample}}_{str(config['PAIR_2'])}.{str(config['RAW_DATA_EXTENSION'])}.gz"
	output:
		sample = f"{OUTPUT_DIR}Bismark/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.bam" if config["RUN_TRIMMING"] else f"{OUTPUT_DIR}Bismark/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.bam",
		report = f"{OUTPUT_DIR}Bismark/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_PE_report.txt" if config["RUN_TRIMMING"] else f"{OUTPUT_DIR}Bismark/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_PE_report.txt"
	log:
		f"logs/bismark_alignment_{{sample}}_PE_1.log"
	params:
		output = lambda w, output: os.path.split(output.sample)[0],
		genome1 = lambda w, input: os.path.split(os.path.split(os.path.split(input.control)[0])[0])[0],
		prefix = "1",
		bismark_cores = BISMARK_CORES
	benchmark:
		f"{OUTPUT_DIR}benchmark/bismark_pe1_{{sample}}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark --prefix {params.prefix} --multicore {params.bismark_cores} --genome {params.genome1} -1 {input.fastq1} -2 {input.fastq2} -o {params.output} --temp_dir {params.output}"

## Run Bismark to perform alignment to the second parental genome (GENOME_PARENT_2) if reads are paired-end.

rule bismark_alignment_PE_2:
	input:
		control = f"{GENOME_2}Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		fastq1 = f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_1'])}_val_1.fq.gz" if config["RUN_TRIMMING"] else f"{RAW_DATA_DIR}{{sample}}_{str(config['PAIR_1'])}.{str(config['RAW_DATA_EXTENSION'])}.gz",
		fastq2 = f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_2'])}_val_2.fq.gz" if config["RUN_TRIMMING"] else f"{RAW_DATA_DIR}{{sample}}_{str(config['PAIR_2'])}.{str(config['RAW_DATA_EXTENSION'])}.gz"
	output:
		sample = f"{OUTPUT_DIR}Bismark/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.bam" if config["RUN_TRIMMING"] else f"{OUTPUT_DIR}Bismark/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.bam",
		report = f"{OUTPUT_DIR}Bismark/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_PE_report.txt" if config["RUN_TRIMMING"] else f"{OUTPUT_DIR}Bismark/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_PE_report.txt"
	log:
		f"logs/bismark_alignment_{{sample}}_PE_2.log"
	params:
		output = lambda w, output: os.path.split(output.sample)[0],
		genome2 = lambda w, input: os.path.split(os.path.split(os.path.split(input.control)[0])[0])[0],
		prefix = "2",
		bismark_cores = BISMARK_CORES
	benchmark:
		f"{OUTPUT_DIR}benchmark/bismark_pe2_{{sample}}.txt"
	conda:
		"../envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark --prefix {params.prefix} --multicore {params.bismark_cores} --genome {params.genome2} -1 {input.fastq1} -2 {input.fastq2} -o {params.output} --temp_dir {params.output}"
