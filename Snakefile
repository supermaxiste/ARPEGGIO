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
		OUTPUT_DIR + "FastQC/parent1_R1_val_1_fastqc.zip"#,
        #OUTPUT_DIR + "FASTQtrimmed/parent1_R1_val_1.fq.gz"

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
