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

# Save settings from config file

OUTPUT_DIR = getpath(config["OUTPUT"])
RAW_DATA_DIR = getpath(config["RAW_DATA"])
GENOME_1 = getpath(config["GENOME_PARENT_1"])
GENOME_2 = getpath(config["GENOME_PARENT_2"])
EAGLE = config["EAGLE"]
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

# General rule to run all analyses depending on the config settings

def dmr_input(wildcards):
	input = []
	if config["RUN_DOWNSTREAM"]:
		if config["ONLY_CG_CONTEXT"]:
			if config["POLYPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid.txt", context = ["CG_context"]))
			elif config["DIPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid.txt", context = ["CG_context"]))
			else:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent1_v_allo_{context}.txt", context = ["CG_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent2_v_allo_{context}.txt", context = ["CG_context"]))
		else:
			if config["POLYPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
			elif config["DIPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
			else:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent1_v_allo_{context}.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent2_v_allo_{context}.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
	if config["RUN_DMR_ANALYSIS"]:
		if config["ONLY_CG_CONTEXT"]:
			if config["POLYPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_1.txt", context = ["CG_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_2.txt", context = ["CG_context"]))
			elif config["DIPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_diploid_{context}.txt", context = ["CG_context"]))
			else:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo.txt", context = ["CG_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo.txt", context = ["CG_context"]))
		else:
			if config["POLYPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_1.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_2.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
			elif config["DIPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_diploid_{context}.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
			else:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
	return input

## Run all analyses
rule all:
	input:
		OUTPUT_DIR + "MultiQC/multiqc_report.html",
		dmr_input

######################### Main rules of ARPEGGIO #############################

include: "rules/quality_check.smk"
include: "rules/trimming.smk"
include: "rules/alignment.smk"
include: "rules/read_sorting.smk"
include: "rules/deduplication.smk"
include: "rules/methylation_extraction.smk"
include: "rules/differential_methylation_analysis.smk"


## Define a function to create an input for MultiQC to include all the settings specified in the config file.

def multiqc_input(wildcards):
	input = []
	if config["IS_PAIRED"]:
		input.extend(expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_1"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_2"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
	else:
		input.extend(expand(OUTPUT_DIR + "FastQC/{sample}_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist()))
	if config["RUN_TRIMMING"]:
		if config["IS_PAIRED"]:
			input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz", sample = samples.name[samples.type == 'PE'].values.tolist()))
			input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz", sample = samples.name[samples.type == 'PE'].values.tolist()))
			input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
			input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
		else:
			input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed.fq.gz", sample = samples.name[samples.type == 'SE'].values.tolist()))
			input.extend(expand(OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist()))
	if config["RUN_BISMARK"]:
		if config["IS_PAIRED"]:
			## alignment
			input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent2')].values.tolist()))
			input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent1')].values.tolist()))
			## deduplication
			input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent2')].values.tolist()))
			input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent1')].values.tolist()))
			## qualimap
			input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_p1/qualimapReport.html", sample = samples.name[samples.origin == 'parent1'].values.tolist()))
			input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_p2/qualimapReport.html", sample = samples.name[samples.origin == 'parent2'].values.tolist()))
			input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_allo_pe_1/qualimapReport.html", sample = samples.name[(samples.type == 'PE') & (samples.origin == 'allopolyploid')].values.tolist()))
			input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_allo_pe_2/qualimapReport.html", sample = samples.name[(samples.type == 'PE') & (samples.origin == 'allopolyploid')].values.tolist()))
		else:
			if config["RUN_TRIMMING"]:
				## not paired, trimmed
				## alignment
				input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_trimmed_bismark_bt2.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent2')].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_trimmed_bismark_bt2.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent1')].values.tolist()))
				## deduplication
				input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent2')].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent1')].values.tolist()))
				## qualimap
				input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_p1/qualimapReport.html", sample = samples.name[samples.origin == 'parent1'].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_p2/qualimapReport.html", sample = samples.name[samples.origin == 'parent2'].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_allo_se_1/qualimapReport.html", sample = samples.name[(samples.type == 'SE') & (samples.origin == 'allopolyploid')].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_allo_se_2/qualimapReport.html", sample = samples.name[(samples.type == 'SE') & (samples.origin == 'allopolyploid')].values.tolist()))
			else:
				## not paired, not trimmed
				## alignment
				input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_bismark_bt2.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent2')].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_bismark_bt2.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent1')].values.tolist()))
				## deduplication
				input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent2')].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bam", sample = samples.name[(samples.type == 'SE') & (samples.origin != 'parent1')].values.tolist()))
				## qualimap
				input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_p1/qualimapReport.html", sample = samples.name[samples.origin == 'parent1'].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_p2/qualimapReport.html", sample = samples.name[samples.origin == 'parent2'].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_allo_se_1/qualimapReport.html", sample = samples.name[(samples.type == 'SE') & (samples.origin == 'allopolyploid')].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "qualimap/{sample}_allo_se_2/qualimapReport.html", sample = samples.name[(samples.type == 'SE') & (samples.origin == 'allopolyploid')].values.tolist()))
	return input

## Define a function to create parameter directories for multiqc based on the settings in the config file

def multiqc_params(wildcards):
	param = [OUTPUT_DIR + "FastQC"]
	if config["RUN_TRIMMING"]:
		param.append(OUTPUT_DIR + "FASTQtrimmed")
	if config["RUN_BISMARK"]:
		param.append(OUTPUT_DIR + "Bismark")
		param.append(OUTPUT_DIR + "qualimap")
	return param
