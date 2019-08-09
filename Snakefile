## Configuration file check
import os
if len(config) == 0:
	if os.path.isfile("./config.yaml"):
		configfile: "./config.yaml"
	else:
		sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")

## Import metadata file
import pandas as pd
samples = pd.read_csv(config["METADATA"], sep='\t')

## Clean path to provided input and output directories
import re
def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

# Save output directory, raw data directory and FastQC directory
OUTPUT_DIR = getpath(config["OUTPUT"])
RAW_DATA_DIR = getpath(config["RAW_DATA"])
EAGLE = config["EAGLE"]
CORES = config["CORES_NUMBER"]
BISMARK_CORES = round(config["CORES_NUMBER"]/3) if config["CORES_NUMBER"]>2 else 1

# Count number of samples

n_samples_p1 = len(samples.name[samples.origin == "parent1"])
n_samples_p2 = len(samples.name[samples.origin == "parent2"])
n_samples_allo = len(samples.name[samples.origin == "allopolyploid"])

# General rule to run all analyses, this rule is needed as Snakemake executes only the first rules

## Run all analyses
rule all:
	input:
		OUTPUT_DIR + "MultiQC/multiqc_report.html",
		OUTPUT_DIR + "DMR_analysis/dmrseq/CG_context/parent1_v_allo.txt",
		OUTPUT_DIR + "DMR_analysis/dmrseq/CG_context/parent2_v_allo.txt"

##### The following pseudo-rules generate output files for the main rules #####
# Pseudo-Rule for Quality Control with FastQC on raw read files
rule pseudo_quality_control:
	input:
		expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_1"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()),
		expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_2"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()),
		expand(OUTPUT_DIR + "FastQC/{sample}_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist())

# Pseudo-Rule for Trimming with TrimGalore on raw read files
rule pseudo_trimming:
	input:
		expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz", sample = samples.name[samples.type == 'PE'].values.tolist()),
		expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz", sample = samples.name[samples.type == 'PE'].values.tolist()),
		expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_fq.gz", sample = samples.name[samples.type == 'SE'].values.tolist())


# Pseudo-Rule for Quality Control with FastQC on trimmed read files
rule pseudo_quality_control_trimmed:
	input:
		expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()),
		expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()),
		expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist())


rule pseudo_preparation:
	input:
		OUTPUT_DIR + "Bismark/genome_preparation_done"


# Pseudo-Rule for deduplication with Bismark on bam files
rule pseudo_deduplication:
	input:
		expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_bismark_bt2.deduplicated.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent2')].values.tolist()),
		expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_bismark_bt2.deduplicated.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent1')].values.tolist()),
		expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent2')].values.tolist()),
		expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent1')].values.tolist())

rule pseudo_read_sorting:
	input:
		expand(OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1.ref.bam", sample = samples.name[samples.origin == 'allopolyploid']),
		expand(OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2.ref.bam", sample = samples.name[samples.origin == 'allopolyploid'])

rule pseudo_methylation_extraction:
	input:
		expand(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/{sample}_bismark_bt2.deduplicated.bismark.cov.gz", sample = samples.name[samples.origin == 'parent1']),
		expand(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/{sample}_bismark_bt2.deduplicated.bismark.cov.gz", sample = samples.name[samples.origin == 'parent2']),
		expand(OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_bismark_bt2.deduplicated.bismark.cov.gz", sample = samples.name[samples.origin == 'allopolyploid']),
		expand(OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_bismark_bt2.deduplicated.bismark.cov.gz", sample = samples.name[samples.origin == 'allopolyploid']),
		expand(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz", sample = samples.name[samples.origin == 'parent1']),
		expand(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz", sample = samples.name[samples.origin == 'parent2']),
		expand(OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz", sample = samples.name[samples.origin == 'allopolyploid']),
		expand(OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz", sample = samples.name[samples.origin == 'allopolyploid'])


######################### Main rules of ARPEGGIO #############################

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
	conda:
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		" fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## Run TrimGalore to trim reads, there are two rules depending on paired or single-end status
## For single-end reads

rule trim_galore:
	input:
		fastq1 = RAW_DATA_DIR + "{sample}_" + str(config["PAIR_1"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz" if config["IS_PAIRED"] else RAW_DATA_DIR + "{sample}." + str(config["RAW_DATA_EXTENSION"]) + ".gz",
		fastq2 = RAW_DATA_DIR + "{sample}_" + str(config["PAIR_2"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz" if config["IS_PAIRED"] else ""
	output:
		OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz" if config["IS_PAIRED"] else OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed.fq.gz",
		OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz" if config["IS_PAIRED"] else ""
	params:
		FASTQtrimmeddir = OUTPUT_DIR + "FASTQtrimmed"
	log:
		OUTPUT_DIR + "logs/trimgalore_{sample}.log"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2}" if config["IS_PAIRED"] else "trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq1}"

## Run FastQC on SE trimmed reads resulting from trim_galore

rule quality_control_trimmed_SE:
	input:
		fastq = OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed.fq.gz"
	output:
		OUTPUT_DIR + "FASTQtrimmed/{sample}_fastqc.zip"
	params:
		FastQC = OUTPUT_DIR + "FASTQtrimmed/"
	log:
		OUTPUT_DIR + "logs/fastqc_trimmed_{sample}.log"
	conda:
		"envs/environment.yaml"
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
	conda:
		"envs/environment.yaml"
	threads:
		config["CORES_NUMBER"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input}"


## Run Bismark to prepare synthetic converted genomes

rule bismark_prepare_genome:
	input:
		genome1 = config["GENOME_PARENT_1"],
		genome2 = config["GENOME_PARENT_2"]
	output:
		OUTPUT_DIR + "Bisulfite_Genome"
	conda:
		"envs/environment.yaml"
	shell:
		"bismark_genome_preparation {input.genome1} > {output}; "
		"bismark_genome_preparation {input.genome2} > {output}"

## Run Bismark to perform alignment to the first parental genome (GENOME_PARENT_1) if reads are single-end.

rule bismark_alignment_SE_1:
	input:
		OUTPUT_DIR + "Bisulfite_Genome" if config["GENOME_PREPARATION"] else "",
		fastq = OUTPUT_DIR + "FASTQtrimmed/SE/{sample}_trimmed.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}." + str(config["RAW_DATA_EXTENSION"]) + ".gz"

	output:
		sample = OUTPUT_DIR + "Bismark/{sample}_1/{sample}_bismark_bt2.bam",
		report = OUTPUT_DIR + "Bismark/{sample}_1/{sample}_bismark_SE_report.txt"
	params:
		output = OUTPUT_DIR + "Bismark/{sample}_1/",
		genome1 = config["GENOME_PARENT_1"]
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
		"bismark --multicore {threads}  -o {params.output} --temp_dir {params.output} --genome {params.genome1} {input.fastq}"

## Run Bismark to perform alignment to the second parental genome (GENOME_PARENT_2) if reads are single-end.

rule bismark_alignment_SE_2:
	input:
		OUTPUT_DIR + "Bisulfite_Genome" if config["GENOME_PREPARATION"] else "",
		fastq = OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		sample = OUTPUT_DIR + "Bismark/{sample}_2/{sample}_bismark_bt2.bam",
		report = OUTPUT_DIR + "Bismark/{sample}_2/{sample}_bismark_SE_report.txt"
	params:
		output = OUTPUT_DIR + "Bismark/{sample}_2/",
		genome2 = config["GENOME_PARENT_2"]
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
		"bismark --multicore {threads}  -o {params.output} --temp_dir {params.output} --genome {params.genome2} {input.fastq}"

## Run Bismark to perform alignment to the first parental genome (GENOME_PARENT_1) if reads are paired-end.

rule bismark_alignment_PE_1:
	input:
		OUTPUT_DIR + "Bisulfite_Genome" if config["GENOME_PREPARATION"] else "",
		fastq1 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_1"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz",
		fastq2 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_2"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		sample = OUTPUT_DIR + "Bismark/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.bam",
		report = OUTPUT_DIR + "Bismark/{sample}_1/{sample}_R1_val_1_bismark_bt2_PE_report.txt"
	params:
		output = OUTPUT_DIR + "Bismark/{sample}_1/",
		genome1 = config["GENOME_PARENT_1"]
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
		"bismark --multicore {threads} --genome {params.genome1} -1 {input.fastq1} -2 {input.fastq2} -o {params.output} --temp_dir {params.output}"

## Run Bismark to perform alignment to the second parental genome (GENOME_PARENT_2) if reads are paired-end.

rule bismark_alignment_PE_2:
	input:
		OUTPUT_DIR + "Bisulfite_Genome" if config["GENOME_PREPARATION"] else "",
		fastq1 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_1"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz",
		fastq2 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_2"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		sample = OUTPUT_DIR + "Bismark/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.bam",
		report = OUTPUT_DIR + "Bismark/{sample}_2/{sample}_R1_val_1_bismark_bt2_PE_report.txt"
	params:
		output = OUTPUT_DIR + "Bismark/{sample}_2/",
		genome2 = config["GENOME_PARENT_2"]
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
		"bismark --multicore {threads} --genome {params.genome2} -1 {input.fastq1} -2 {input.fastq2} -o {params.output} --temp_dir {params.output}"

## Run deduplication of the alignments to remove duplicated reads for SE reads (GENOME_PARENT_1)

rule deduplication_SE_1:
	input:
		OUTPUT_DIR + "Bismark/{sample}_1/{sample}_bismark_bt2.bam"
	output:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_bismark_bt2.deduplicated.bam"
	params:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	conda:
		"envs/environment.yaml"
	shell:
		"deduplicate_bismark -s --output_dir {params} --bam {input}"

## Run deduplication of the alignments to remove duplicated reads for SE reads (GENOME_PARENT_2)

rule deduplication_SE_2:
	input:
		OUTPUT_DIR + "Bismark/{sample}_2/{sample}_bismark_bt2.bam"
	output:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_bismark_bt2.deduplicated.bam"
	params:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	conda:
		"envs/environment.yaml"
	shell:
		"deduplicate_bismark -s --output_dir {params} --bam {input}"

## Run deduplication of the alignments to remove duplicated reads for PE reads (GENOME_PARENT_1)

rule deduplication_PE_1:
	input:
		OUTPUT_DIR + "Bismark/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.bam"
	output:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	params:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	conda:
		"envs/environment.yaml"
	shell:
		"deduplicate_bismark -p --output_dir {params} --bam {input}"

## Run deduplication of the alignments to remove duplicated reads for PE reads (GENOME_PARENT_2)

rule deduplication_PE_2:
	input:
		OUTPUT_DIR + "Bismark/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.bam"
	output:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	params:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	conda:
		"envs/environment.yaml"
	shell:
		"deduplicate_bismark -p --output_dir {params} --bam {input}"

## Run EAGLE-RC to classify reads to the most probable genome for SE reads

rule read_sorting_SE:
	input:
		reads1 = OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_bismark_bt2.deduplicated.bam",
		reads2 = OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_bismark_bt2.deduplicated.bam"
	output:
		OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1.ref.bam",
		OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2.ref.bam"
	params:
		list = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified_reads.list",
		phred = "--phred64" if config["PHRED_SCORE_64"] else "",
		genome1 = config["ASSEMBLY_PARENT_1"],
		genome2 = config["ASSEMBLY_PARENT_2"],
		output = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified"
	shell:
		"{EAGLE} --ngi {params.phred} --ref1={params.genome1} --bam1={input.reads1} --ref2={params.genome2} --bam2={input.reads2} -o {params.output} --bs=3 > {params.list}"

## Run EAGLE-RC to classify reads to the most probable genome for PE reads

rule read_sorting_PE:
	input:
		reads1 = OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam",
		reads2 = OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1.ref.bam",
		OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2.ref.bam"
	params:
		list = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified_reads.list",
		phred = "--phred64" if config["PHRED_SCORE_64"] else "",
		genome1 = config["ASSEMBLY_PARENT_1"],
		genome2 = config["ASSEMBLY_PARENT_2"],
		output = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified"
	shell:
		"{EAGLE} --ngi --paired {params.phred} --ref1={params.genome1} --bam1={input.reads1} --ref2={params.genome2} --bam2={input.reads2} -o {params.output} --bs=3 > {params.list}"

## Run Bismark methylation extraction on SE bam files for parent species 1

rule methylation_extraction_SE_parent_1:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_bismark_bt2.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p1/{sample}_bismark_bt2.deduplicated.bismark.cov.gz"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_p1/",
		genome = config["GENOME_PARENT_1"]
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {threads} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on SE bam files for parent species 2

rule methylation_extraction_SE_parent_2:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_bismark_bt2.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p2/{sample}_bismark_bt2.deduplicated.bismark.cov.gz"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_p2/",
		genome = config["GENOME_PARENT_2"]
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {threads} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on SE bam files for allopolyploid species (GENOME_PARENT_1)

rule methylation_extraction_SE_allo_1:
	input:
		OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_bismark_bt2.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_bismark_bt2.deduplicated.bismark.cov.gz"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_1/",
		genome = config["GENOME_PARENT_1"]
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {threads} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on SE bam files for allopolyploid species (GENOME_PARENT_2)

rule methylation_extraction_SE_allo_2:
	input:
		OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_bismark_bt2.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_bismark_bt2.deduplicated.bismark.cov.gz"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_2/",
		genome = config["GENOME_PARENT_2"]
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {threads} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for parent species 1

rule methylation_extraction_PE_parent_1:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		OUTPUT_DIR + "Bismark/extraction/{sample}_p1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_p1/",
		genome = config["GENOME_PARENT_1"]
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {threads} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for parent species 2

rule methylation_extraction_PE_parent_2:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		OUTPUT_DIR + "Bismark/extraction/{sample}_p2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_p2/",
		genome = config["GENOME_PARENT_2"]
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {threads} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for allopolyploid species (GENOME_PARENT_1)

rule methylation_extraction_PE_allo_1:
	input:
		OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.bam"
	output:
		 OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_classified1.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		 OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_classified1.ref.bedGraph.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_1/",
		genome = config["GENOME_PARENT_1"]
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {threads} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for allopolyploid species (GENOME_PARENT_2)

rule methylation_extraction_PE_allo_2:
	input:
		 OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_classified2.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_classified2.ref.bedGraph.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_2/",
		genome = config["GENOME_PARENT_2"]
	conda:
		"envs/environment.yaml"
	threads:
		BISMARK_CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {threads} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark coverage2cytosine on extraction output to obtain a single file with information about all cytosines (parent 1)

rule coverage2cytosine_1:
	input:
		f1 = OUTPUT_DIR + "Bismark/extraction/{sample}_p1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
	output:
		o1 = OUTPUT_DIR + "Bismark/extraction/{sample}_p1/{sample}.CX_report.txt"
	params:
		genome1 = config["GENOME_PARENT_1"],
		filename1 = OUTPUT_DIR + "Bismark/extraction/{sample}_p1/{sample}"
	conda:
		"envs/environment.yaml"
	shell:
		"coverage2cytosine -CX --genome_folder {params.genome1} -o {params.filename1} {input.f1}"

## Run Bismark coverage2cytosine on extraction output to obtain a single file with information about all cytosines (parent 2)

rule coverage2cytosine_2:
	input:
		f2 = OUTPUT_DIR + "Bismark/extraction/{sample}_p2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
	output:
		o2 = OUTPUT_DIR + "Bismark/extraction/{sample}_p2/{sample}.CX_report.txt"
	params:
		genome2 = config["GENOME_PARENT_2"],
		filename2 = OUTPUT_DIR + "Bismark/extraction/{sample}_p2/{sample}"
	conda:
		"envs/environment.yaml"
	shell:
		"coverage2cytosine -CX --genome_folder {params.genome2} -o {params.filename2} {input.f2}"

## Run Bismark coverage2cytosine on extraction output to obtain a single file with information about all cytosines (allopolyploid)

rule coverage2cytosine_allo:
	input:
		f1 = OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_classified1.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		f2 = OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_classified2.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
	output:
		o1 = OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}.CX_report.txt",
		o2 = OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}.CX_report.txt"
	params:
		genome1 = config["GENOME_PARENT_1"],
		genome2 = config["GENOME_PARENT_2"],
		filename1 = OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}",
		filename2 = OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}"
	conda:
		"envs/environment.yaml"
	shell:
		"""
		coverage2cytosine -CX --genome_folder {params.genome1} -o {params.filename1} {input.f1}
		coverage2cytosine -CX --genome_folder {params.genome2} -o {params.filename2} {input.f2}
		"""

## Run R script to remove all rows with no cytosines and generate three files for each cytosine context for parent 1

rule context_separation_parent_1:
	input:
		p1 = OUTPUT_DIR + "Bismark/extraction/{sample}_p1/{sample}.CX_report.txt"
	output:
		OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CG.cov",
		OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CHG.cov",
		OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CHH.cov"
	params:
		sample_name = "{sample}",
		output = OUTPUT_DIR + "DMR_analysis/context_separation/parent1/"
	conda:
		"envs/environment.yaml"
	shell:
		"Rscript scripts/CoverageFileGeneratorComplete.R {input.p1} {params.output} {params.sample_name}"

## Run R script to remove all rows with no cytosines and generate three files for each cytosine context for parent 2

rule context_separation_parent_2:
	input:
		p1 = OUTPUT_DIR + "Bismark/extraction/{sample}_p2/{sample}.CX_report.txt"
	output:
		OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CG.cov",
		OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CHG.cov",
		OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CHH.cov"
	params:
		sample_name = "{sample}",
		output = OUTPUT_DIR + "DMR_analysis/context_separation/parent2/"
	conda:
		"envs/environment.yaml"
	shell:
		"Rscript scripts/CoverageFileGeneratorComplete.R {input.p1} {params.output} {params.sample_name}"

## Run R script to remove all rows with no cytosines and generate three files for each cytosine context for parent 1

rule combine_context_allo:
	input:
		p1 = OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}.CX_report.txt",
		p2 = OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}.CX_report.txt"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_12/{sample}.CX_report.txt"
	shell:
		"cat {input.p1} {input.p2} > {output}"

rule context_separation_allo:
	input:
		OUTPUT_DIR + "Bismark/extraction/{sample}_12/{sample}.CX_report.txt"
	output:
		OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CG.cov",
		OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CHG.cov",
		OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CHH.cov"
	params:
		sample_name = "{sample}",
		output = OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/"
	conda:
		"envs/environment.yaml"
	shell:
		"Rscript scripts/CoverageFileGeneratorComplete.R {input} {params.output} {params.sample_name};"
		"Rscript scripts/CoverageFileGeneratorComplete.R {input} {params.output} {params.sample_name}"

# Define a function to create an input for dmrseq_CG to include all the samples from the previous steps

def dmrseq_CG_input_p1(wildcards):
	input = []
	input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CG.cov", sample = samples.name[samples.origin == 'parent1'].values.tolist()))
	return input

def dmrseq_CG_input_p2(wildcards):
	input = []
	input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CG.cov", sample = samples.name[samples.origin == 'parent2'].values.tolist()))
	return input

def dmrseq_CG_input_allo(wildcards):
	input = []
	input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CG.cov", sample = samples.name[samples.origin == 'allopolyploid'].values.tolist()))
	return input

# Run dmrseq for CG context

rule dmrseq_CG:
	input:
		p1 = dmrseq_CG_input_p1,
		p2 = dmrseq_CG_input_p2,
		allo = dmrseq_CG_input_allo
	output:
		comparison1 = OUTPUT_DIR + "DMR_analysis/dmrseq/CG_context/parent1_v_allo.txt",
		comparison2 = OUTPUT_DIR + "DMR_analysis/dmrseq/CG_context/parent2_v_allo.txt"
	params:
		n_samples_p1 = n_samples_p1,
		n_samples_p2 = n_samples_p2,
		n_samples_allo = n_samples_allo,
		script = "scripts/dmrseq.R"
	threads:
		config["CORES_NUMBER"]
	conda:
		"envs/environment_R.yaml"
	shell:
		"Rscript scripts/dmrseq.R {params.n_samples_p1} {params.n_samples_allo} {output.comparison1} {threads} {input.p1} {input.allo};"
		"Rscript scripts/dmrseq.R {params.n_samples_p2} {params.n_samples_allo} {output.comparison2} {threads} {input.p2} {input.allo}"


## Define a function to create an input for MultiQC to include all the settings specified in the config file.

def multiqc_input(wildcards):
	input = []
	input.extend(expand(OUTPUT_DIR + "FastQC/{sample}_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist()))
	input.extend(expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_1"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
	input.extend(expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_2"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
	if config["GENOME_PREPARATION"]:
		input.extend(expand(OUTPUT_DIR + "Bisulfite_Genome"))
	if config["RUN_TRIMMING"]:
		input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed.fq.gz", sample = samples.name[samples.type == 'SE'].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz", sample = samples.name[samples.type == 'PE'].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz", sample = samples.name[samples.type == 'PE'].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
	if config["RUN_BISMARK"]:
		## alignment
		input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_1/{sample}_bismark_bt2.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent2')].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_2/{sample}_bismark_bt2.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent1')].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent2')].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent1')].values.tolist()))
		## deduplication
		input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_bismark_bt2.deduplicated.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent2')].values.tolist())),
		input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_bismark_bt2.deduplicated.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent1')].values.tolist())),
		input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent2')].values.tolist())),
		input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent1')].values.tolist()))
	return input

## Define a function to create input directories based on the settings in the config file
def multiqc_params(wildcards):
	param = [OUTPUT_DIR + "FastQC", OUTPUT_DIR + "Bismark"]
	if config["GENOME_PREPARATION"]:
		param.append(OUTPUT_DIR + "Bisulfite_Genome")
	if config["RUN_TRIMMING"]:
		param.append(OUTPUT_DIR + "FASTQtrimmed")
	if config["RUN_BISMARK"]:
		param.append(OUTPUT_DIR + "Bismark")
	return param

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
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.inputdir} -f -o {params.multiqcdir}"
