## Configuration file check
import os
if len(config) == 0:
  if os.path.isfile("./config.yaml"):
    configfile: "./config.yaml"
  else:
    sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")

## Import metadata
import pandas as pd
samples = pd.read_table(config["metatxt"])

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
OUTPUT_DIR = getpath(config["output"])
RAW_DATA_DIR = getpath(config["RAW_DATA"])

##### The following pseudo-rules generate output files for the main rules #####
# Pseudo-Rule for Quality Control with FastQC on raw read files
rule pseudo_quality_control:
    input:
        expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_1"]) + "_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
        expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_2"]) + "_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
        expand(OUTPUT_DIR + "FastQC/{sample}_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

# Pseudo-Rule for Trimming with TrimGalore on raw read files
rule pseudo_trimming:
    input:
		expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["fqext1"]) + "_trimed_1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["fqext2"]) + "_trimed_2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(OUTPUT_DIR + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

# Pseudo-Rule for Quality Control with FastQC on trimmed read files
rule pseudo_quality_control_trimmed:
    input:
        expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_1"]) + "_trimed_1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
        expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_2"]) + "_trimed_2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
        expand(OUTPUT_DIR + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

######################### Main rules of ARPEGGIO #############################

## Run FastQC on raw untrimmed reads for quality control
rule quality_control:
	input:
		fastq = RAW_DATA_DIR + "{sample}." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "FastQC/{sample}_fastqc.zip"
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
"fastqc -o {params.FastQC} -t {threads} {input.fastq}"
