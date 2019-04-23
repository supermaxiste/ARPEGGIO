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

# General rule to run all analyses, this rule is needed as Snakemake executes only the first rules

## Run all analyses
rule all:
	input:
		OUTPUT_DIR + "Bismark/allopolyploid1/allopolyploid1_R1_val_1_bismark_bt2_pe.bam"

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
        expand(OUTPUT_DIR + "FastQtrimmed/{sample}_fq.gz", sample = samples.name[samples.type == 'SE'].values.tolist())

# Pseudo-Rule for Quality Control with FastQC on trimmed read files
rule pseudo_quality_control_trimmed:
    input:
        expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_1"]) + "_val_1.fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()),
        expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_2"]) + "_val_2.fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()),
        expand(OUTPUT_DIR + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist())

rule pseudo_preparation:
    input:
        OUTPUT_DIR + "Bismark/genome_preparation_done"

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
		config["CORES_NUMBER"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
        " fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## Run TrimGalore to trim reads, there are two rules depending on paired or single-end status
## For single-end reads

rule trim_galore_SE:
	input:
		fastq = RAW_DATA_DIR + "{sample}." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed.fq.gz"
	params:
		FASTQtrimmeddir = OUTPUT_DIR + "FASTQtrimmed"
	log:
		OUTPUT_DIR + "logs/trimgalore_{sample}.log"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq}"

## For paired-end reads
rule trim_galore_PE:
    input:
        fastq1 = RAW_DATA_DIR + "{sample}_" + str(config["PAIR_1"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz",
        fastq2 = RAW_DATA_DIR + "{sample}_" + str(config["PAIR_2"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
    output:
        OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz",
        OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz"
    params:
        FASTQtrimmeddir = OUTPUT_DIR + "FASTQtrimmed"
    log:
        OUTPUT_DIR + "logs/trimgalore_{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        "echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
        "trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2}"

## Run FastQC on trimmed reads resulting from trim_galore

rule quality_control_trimmed:
    input:
        fastq = OUTPUT_DIR + "FASTQtrimmed/{sample}.fq.gz"
    output:
        OUTPUT_DIR + "FastQC/{sample}_fastqc.zip"
    params:
        FastQC = OUTPUT_DIR + "FastQC"
    log:
        OUTPUT_DIR + "logs/fastqc_trimmed_{sample}.log"
    conda:
        "envs/environment.yaml"
    threads:
        config["CORES_NUMBER"]
    shell:
        "echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
        "fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## Run Bismark to prepare synthetic converted genomes

rule bismark_prepare_genome:
    input:
        genome1 = config["GENOME_PARENT_1"],
        genome2 = config["GENOME_PARENT_2"]
    output:
        OUTPUT_DIR + "Bisulfite_Genome"
    shell:
        "bismark_genome_preparation {input.genome1} > {output}; "
        "bismark_genome_preparation {input.genome2} > {output}"

## Run Bismark to perform alignment to the first parental genome (GENOME_PARENT_1) if reads are paired-end.

rule bismark_alignment_PE_1:
    input:
        OUTPUT_DIR + "Bisulfite_Genome" if config["GENOME_PREPARATION"] else "",
        fastq1 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_1"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz",
        fastq2 = OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz" if config["RUN_TRIMMING"] else RAW_DATA_DIR + "{sample}_" + str(config["PAIR_2"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
    output:
        sample = OUTPUT_DIR + "Bismark/{sample}/{sample}_R1_val_1_bismark_bt2_pe.bam",
        report = OUTPUT_DIR + "Bismark/{sample}/{sample}_R1_val_1_bismark_bt2_PE_report.txt"
    params:
        output = OUTPUT_DIR + "Bismark/{sample}/",
        genome1 = config["GENOME_PARENT_1"]
    log:
        OUTPUT_DIR + "logs/bismark_{sample}.log"
    conda:
        "envs/environment.yaml"
    threads:
        config["CORES_NUMBER"]
    shell:
        "echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
        "bismark --multicore {threads} --genome {params.genome1} -1 {input.fastq1} -2 {input.fastq2} -o {params} --temp_dir {params.output}"
