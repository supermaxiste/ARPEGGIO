## Check minimum Snakemake version

from snakemake.utils import min_version

min_version("5.20.1")

## Configuration file check
import os
if len(config) == 0:
	if os.path.isfile("./config.yaml"):
		configfile: "./config.yaml"
	else:
		sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")

## Conditional parameters check

if config["RUN_BISMARK"]==False:
	if config["RUN_READ_SORTING"] or config["RUN_DMR_ANALYSIS"] or config["RUN_DOWNSTREAM"]:
		sys.exit("Your conditional parameters are set incorrectly: one step is stopping the workflow and one downstream step is set to True. Please fix this problem before running ARPEGGIO")
if config["RUN_READ_SORTING"]==False:
	if config["RUN_DMR_ANALYSIS"] or config["RUN_DOWNSTREAM"]:
		sys.exit("Your conditional parameters are set incorrectly: one step is stopping the workflow and one downstream step is set to True. Please fix this problem before running ARPEGGIO")
if config["RUN_DMR_ANALYSIS"]==False:
	if config["RUN_DOWNSTREAM"]:
		sys.exit("Your conditional parameters are set incorrectly: one step is stopping the workflow and one downstream step is set to True. Please fix this problem before running ARPEGGIO")

## Import metadata file
import pandas as pd
samples = pd.read_csv(config["METADATA"], sep='\t')

# Save settings from config file and sanitizing paths
import math

OUTPUT_DIR = os.path.normpath(config["OUTPUT"]) + "/"
RAW_DATA_DIR = os.path.normpath(config["RAW_DATA"]) + "/"
GENOME_1 = os.path.normpath(config["GENOME_PARENT_1"]) + "/"
GENOME_2 = os.path.normpath(config["GENOME_PARENT_2"]) + "/"
CONTROL_GENOME = os.path.normpath(config["CONTROL_GENOME"]) + "/"
CORES = workflow.cores
BISMARK_CORES = max(math.trunc(CORES/3), 1)
TRIM_CORES = max(math.trunc(CORES/6 - 1), 1) if config["IS_PAIRED"] else max(math.trunc(CORES/3 - 1), 1)

# Count number of samples for parental species and allopolyploid

n_samples_p1 = len(samples.name[samples.origin == "parent1"])
n_samples_p2 = len(samples.name[samples.origin == "parent2"])
n_samples_allo = len(samples.name[samples.origin == "allopolyploid"])

# Check if metadata file has been setup and read correctly

if ((not config["POLYPLOID_ONLY"]) and (not config["DIPLOID_ONLY"])):
	if n_samples_p1 < 2:
		sys.exit("There seems to be fewer than two samples for your parent1 species, please have at least two samples in your metadata file or make sure the metadata file has the correct formatting (tab-separated columns)")
	elif n_samples_p2 < 2:
		sys.exit("There seems to be fewer than two samples for your parent2 species, please have at least two samples in your metadata file or make sure the metadata file has the correct formatting (tab-separated columns)")
	elif n_samples_allo < 2:
		sys.exit("There seems to be fewer than two samples for your allopolyploid species, please have at least two samples in your metadata file or make sure the metadata file has the correct formatting (tab-separated columns)")
if config["POLYPLOID_ONLY"]:
	n_samples_A = len(samples.name[(samples.origin == "allopolyploid") & (samples.condition == "A")])
	n_samples_B = len(samples.name[(samples.origin == "allopolyploid") & (samples.condition == "B")])
	if n_samples_A < 2 or n_samples_B < 2:
		sys.exit("There seems to be fewer than two samples for either of your two conditions, please have at least two samples per condition in your metadata file or make sure the metadata file has the correct formatting (tab-separated columns)")
if config["DIPLOID_ONLY"]:
	n_samples_A = len(samples.name[((samples.origin == "parent1") | (samples.origin == "parent2")) & (samples.condition == "A")])
	n_samples_B = len(samples.name[((samples.origin == "parent1") | (samples.origin == "parent2")) & (samples.condition == "B")])
	if n_samples_A < 2 or n_samples_B < 2:
		sys.exit("There seems to be fewer than two samples for either of your two conditions, please have at least two samples per condition in your metadata file or make sure the metadata file has the correct formatting (tab-separated columns)")

####################### Docker setup ##########################################

singularity: config["DOCKER_IMAGE"]

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
