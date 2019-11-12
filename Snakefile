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

# Save settings from config file and sanitizing paths

OUTPUT_DIR = os.path.normpath(config["OUTPUT"]) + "/"
RAW_DATA_DIR = os.path.normpath(config["RAW_DATA"]) + "/"
GENOME_1 = os.path.normpath(config["GENOME_PARENT_1"]) + "/"
GENOME_2 = os.path.normpath(config["GENOME_PARENT_2"]) + "/"
CONTROL_GENOME = os.path.normpath(config["CONTROL_GENOME"]) + "/"
CORES = config["CORES_NUMBER"]
BISMARK_CORES = round(config["CORES_NUMBER"]/3) if config["CORES_NUMBER"]>2 else 1

# Count number of samples for parental species and allopolyploid

n_samples_p1 = len(samples.name[samples.origin == "parent1"])
n_samples_p2 = len(samples.name[samples.origin == "parent2"])
n_samples_allo = len(samples.name[samples.origin == "allopolyploid"])
if config["POLYPLOID_ONLY"]:
	n_samples_A = len(samples.name[(samples.origin == "allopolyploid") & (samples.condition == "A")])
	n_samples_B = len(samples.name[(samples.origin == "allopolyploid") & (samples.condition == "B")])
if config["DIPLOID_ONLY"]:
	n_samples_A = len(samples.name[((samples.origin == "parent1") | (samples.origin == "parent2")) & (samples.condition == "A")])
	n_samples_B = len(samples.name[((samples.origin == "parent1") | (samples.origin == "parent2")) & (samples.condition == "B")])

####################### Functions to define inputs ############################

include: "rules/input_functions.smk"

## Run all analyses
rule all:
	input:
		OUTPUT_DIR + "MultiQC/multiqc_report.html",
		dmr_input

######################### Main rules of ARPEGGIO #############################

include: "rules/build_eagle.smk"
include: "rules/quality_check.smk"
include: "rules/trimming.smk"
include: "rules/alignment.smk"
include: "rules/read_sorting.smk"
include: "rules/deduplication.smk"
include: "rules/methylation_extraction.smk"
include: "rules/differential_methylation_analysis.smk"
include: "rules/downstream.smk"
