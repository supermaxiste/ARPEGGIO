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
GENOME_1 = getpath(config["GENOME_PARENT_1"])
GENOME_2 = getpath(config["GENOME_PARENT_2"])
EAGLE = config["EAGLE"]
CORES = config["CORES_NUMBER"]
BISMARK_CORES = round(config["CORES_NUMBER"]/3) if config["CORES_NUMBER"]>2 else 1

# Count number of samples

n_samples_p1 = len(samples.name[samples.origin == "parent1"])
n_samples_p2 = len(samples.name[samples.origin == "parent2"])
n_samples_allo = len(samples.name[samples.origin == "allopolyploid"])
if config["POLYPLOID_ONLY"]:
	n_samples_A = len(samples.name[(samples.origin == "allopolyploid") & (samples.condition == "A")])
	n_samples_B = len(samples.name[(samples.origin == "allopolyploid") & (samples.condition == "B")])
if config["DIPLOID_ONLY"]:
	n_samples_A = len(samples.name[((samples.origin == "parent1") | (samples.origin == "parent2")) & (samples.condition == "A")])
	n_samples_B = len(samples.name[((samples.origin == "parent1") | (samples.origin == "parent2")) & (samples.condition == "B")])

# General rule to run all analyses, this rule is needed as Snakemake executes only the first rules

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
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		" fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## Run TrimGalore to trim reads, there are two rules depending on paired or single-end status
## For single-end reads

rule trim_galore_se:
	input:
		fastq1 = RAW_DATA_DIR + "{sample}." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		OUTPUT_DIR + "FASTQtrimmed/{sample}_trimmed.fq.gz"
	params:
		FASTQtrimmeddir = OUTPUT_DIR + "FASTQtrimmed",
		trim_5_r1 = config["CLIP_5_R1"],
		trim_3_r1 = config["CLIP_3_R1"]
	log:
		OUTPUT_DIR + "logs/trimgalore_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/trim_se_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"trim_galore -q 20 --clip_R1 {params.trim_5_r1} --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq1}" if config["RUN_TRIMMING"] and config["TRIM_5_ONLY"] else ("trim_galore -q 20 --clip_R1 {params.trim_5_r1}  --three_prime_clip_R1 {params.trim_3_r1} --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq1}" if config["RUN_TRIMMING"] and config["TRIM_3_ONLY"] else ("trim_galore -q 20 --clip_R1 {params.trim_5_r1}  --three_prime_clip_R1 {params.trim_3_r1} --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq1}" if config["RUN_TRIMMING"] else "trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq1}"))

## For paired-end reads

rule trim_galore_pe:
	input:
		fastq1 = RAW_DATA_DIR + "{sample}_" + str(config["PAIR_1"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz",
		fastq2 = RAW_DATA_DIR + "{sample}_" + str(config["PAIR_2"]) + "." + str(config["RAW_DATA_EXTENSION"]) + ".gz"
	output:
		OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_1"]) + "_val_1.fq.gz",
		OUTPUT_DIR + "FASTQtrimmed/{sample}_" + str(config["PAIR_2"]) + "_val_2.fq.gz"
	params:
		FASTQtrimmeddir = OUTPUT_DIR + "FASTQtrimmed",
		trim_5_r1 = config["CLIP_5_R1"],
		trim_5_r2 = config["CLIP_5_R2"],
		trim_3_r1 = config["CLIP_3_R1"],
		trim_3_r2 = config["CLIP_3_R2"]
	log:
		OUTPUT_DIR + "logs/trimgalore_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/trim_pe_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log};"
		"trim_galore -q 20 --clip_R1 {params.trim_5_r1} --clip_R2 {params.trim_5_r2} --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2}" if config["RUN_TRIMMING"] and config["TRIM_5_ONLY"] else ("trim_galore -q 20 --three_prime_clip_R1 {params.trim_3_r1} --three_prime_clip_R2 {params.trim_3_r2} --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2}" if config["RUN_TRIMMING"] and config["TRIM_3_ONLY"] else ("trim_galore -q 20 --clip_R1 {params.trim_5_r1} --clip_R2 {params.trim_5_r2} --three_prime_clip_R1 {params.trim_3_r1} --three_prime_clip_R2 {params.trim_3_r2} --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2}" if config["RUN_TRIMMING"] else "trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2}"))

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
	benchmark:
		OUTPUT_DIR + "benchmark/qc_trim_pe_{sample}.txt"
	conda:
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input}"


## Run Bismark to prepare synthetic converted genomes

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
		"envs/environment.yaml"
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
		"envs/environment.yaml"
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
		"envs/environment.yaml"
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
		sample = OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.bam",
		report = OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_PE_report.txt"
	params:
		output = OUTPUT_DIR + "Bismark/{sample}_1/",
		genome1 = config["GENOME_PARENT_1"],
		prefix = "1"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/bismark_pe1_{sample}.txt"
	conda:
		"envs/environment.yaml"
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
		sample = OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.bam",
		report = OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_PE_report.txt"
	params:
		output = OUTPUT_DIR + "Bismark/{sample}_2/",
		genome2 = config["GENOME_PARENT_2"],
		prefix = "2"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/bismark_pe2_{sample}.txt"
	conda:
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"echo 'Bismark version:\n' > {log}; bismark --version >> {log}; "
		"bismark --prefix {params.prefix} --multicore {BISMARK_CORES} --genome {params.genome2} -1 {input.fastq1} -2 {input.fastq2} -o {params.output} --temp_dir {params.output}"

## Run deduplication of the alignments to remove duplicated reads for SE reads (GENOME_PARENT_1)

rule deduplication_SE_1:
	input:
		OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_trimmed_bismark_bt2.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_bismark_bt2.bam"
	output:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bam"
	params:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/dedup_se1_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"deduplicate_bismark -s --output_dir {params} --bam {input}"

## Run deduplication of the alignments to remove duplicated reads for SE reads (GENOME_PARENT_2)

rule deduplication_SE_2:
	input:
		OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_trimmed_bismark_bt2.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_bismark_bt2.bam"
	output:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bam"
	params:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/dedup_se2_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"deduplicate_bismark -s --output_dir {params} --bam {input}"

## Run deduplication of the alignments to remove duplicated reads for PE reads (GENOME_PARENT_1)

rule deduplication_PE_1:
	input:
		OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.bam"
	output:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	params:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/dedup_pe1_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"deduplicate_bismark -p --output_dir {params} --bam {input}"

## Run deduplication of the alignments to remove duplicated reads for PE reads (GENOME_PARENT_2)

rule deduplication_PE_2:
	input:
		OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.bam"
	output:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	params:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/"
	log:
		OUTPUT_DIR + "logs/bismark_{sample}.log"
	benchmark:
		OUTPUT_DIR + "benchmark/dedup_pe2_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"deduplicate_bismark -p --output_dir {params} --bam {input}"

## Run EAGLE-RC to classify reads to the most probable genome for SE reads

rule read_sorting_SE:
	input:
		reads1 = OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bam",
		reads2 = OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bam"
	output:
		OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified1.ref.bam",
		OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified2.ref.bam"
	benchmark:
		OUTPUT_DIR + "benchmark/readsorting_se_{sample}.txt"
	params:
		list = OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified_reads.list",
		phred = "--phred64" if config["PHRED_SCORE_64"] else "",
		genome1 = config["ASSEMBLY_PARENT_1"],
		genome2 = config["ASSEMBLY_PARENT_2"],
		output = OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified"
	shell:
		"{EAGLE} --ngi {params.phred} --ref1={params.genome1} --bam1={input.reads1} --ref2={params.genome2} --bam2={input.reads2} -o {params.output} --bs=3 > {params.list}"

## Run EAGLE-RC to classify reads to the most probable genome for PE reads

rule read_sorting_PE:
	input:
		reads1 = OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam",
		reads2 = OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1.ref.bam",
		OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2.ref.bam"
	benchmark:
		OUTPUT_DIR + "benchmark/readsorting_pe_{sample}.txt"
	params:
		list = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified_reads.list",
		phred = "--phred64" if config["PHRED_SCORE_64"] else "",
		genome1 = config["ASSEMBLY_PARENT_1"],
		genome2 = config["ASSEMBLY_PARENT_2"],
		output = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified"
	shell:
		"{EAGLE} --ngi --paired {params.phred} --ref1={params.genome1} --bam1={input.reads1} --ref2={params.genome2} --bam2={input.reads2} -o {params.output} --bs=3 > {params.list}"

# sort bam files for multiqc

rule bam_sorting_p1:
	input:
		p1 = OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bam")
	output:
		o1 = OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated_sorted.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated_sorted.bam")
	conda:
		"envs/environment.yaml"
	shell:
		"samtools sort {input.p1} > {output.o1}"

rule bam_sorting_p2:
	input:
		p2 = OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bam"),
	output:
		o2 = OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated_sorted.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated_sorted.bam")
	conda:
		"envs/environment.yaml"
	shell:
		"samtools sort {input.p2} > {output.o2}"

rule bam_sorting_allo:
	input:
		allo1 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified1.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bam")),
		allo2 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified2.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bam"))
	output:
		o3 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1_sorted.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated_sorted_allo.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified1_sorted.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated_sorted_allo.bam")),
		o4 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2_sorted.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated_sorted_allo.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified2_sorted.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated_sorted_allo.bam"))
	conda:
		"envs/environment.yaml"
	shell:
		"samtools sort {input.allo1} > {output.o3};"
		"samtools sort {input.allo2} > {output.o4}"

## Qualimap rule to get statistics about bam files for parent1

rule qualimap_p1:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated_sorted.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated_sorted.bam")
	output:
		OUTPUT_DIR + "qualimap/{sample}_p1/qualimapReport.html"
	benchmark:
		OUTPUT_DIR + "benchmark/qualimap_p1_{sample}.txt"
	params:
		output = OUTPUT_DIR + "qualimap/{sample}_p1"
	conda:
		"envs/environment2.yaml"
	shell:
		"qualimap bamqc -bam {input} -outdir {params.output}"

## Qualimap rule to get statistics about bam files for parent2

rule qualimap_p2:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated_sorted.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated_sorted.bam")
	output:
		OUTPUT_DIR + "qualimap/{sample}_p2/qualimapReport.html"
	benchmark:
		OUTPUT_DIR + "benchmark/qualimap_p2_{sample}.txt"
	params:
		output = OUTPUT_DIR + "qualimap/{sample}_p2"
	conda:
		"envs/environment2.yaml"
	shell:
		"qualimap bamqc -bam {input} -outdir {params.output}"

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
		"envs/environment2.yaml"
	shell:
		"qualimap bamqc -bam {input.genome1} -outdir {params.output1};"
		"qualimap bamqc -bam {input.genome2} -outdir {params.output2}"

## Qualimap rule to get statistics about bam files for allopolyploid PE reads

rule qualimap_allo_pe:
	input:
		genome1 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1_sorted.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_allo.bam",
		genome2 = OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2_sorted.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated_sorted_allo.bam"
	output:
		OUTPUT_DIR + "qualimap/{sample}_allo_pe_1/qualimapReport.html",
		OUTPUT_DIR + "qualimap/{sample}_allo_pe_2/qualimapReport.html"
	benchmark:
		OUTPUT_DIR + "benchmark/qualimap_allo_pe_{sample}.txt"
	params:
		output1 = OUTPUT_DIR + "qualimap/{sample}_allo_pe_1/",
		output2 = OUTPUT_DIR + "qualimap/{sample}_allo_pe_2/"
	conda:
		"envs/environment2.yaml"
	shell:
		"qualimap bamqc -bam {input.genome1} -outdir {params.output1};"
		"qualimap bamqc -bam {input.genome2} -outdir {params.output2}"

## Run Bismark methylation extraction on SE bam files for parent species 1

rule methylation_extraction_SE_parent_1:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/{sample}_bismark_bt2.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_bismark_bt2.deduplicated.bismark.cov.gz"
	benchmark:
		OUTPUT_DIR + "benchmark/extraction_se_p1_{sample}.txt"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_p1/",
		genome = config["GENOME_PARENT_1"]
	conda:
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"


## Run Bismark methylation extraction on SE bam files for parent species 2

rule methylation_extraction_SE_parent_2:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/{sample}_bismark_bt2.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_bismark_bt2.deduplicated.bismark.cov.gz"
	benchmark:
		OUTPUT_DIR + "benchmark/extraction_se_p2_{sample}.txt"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_p2/",
		genome = config["GENOME_PARENT_2"]
	conda:
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on SE bam files for allopolyploid species (GENOME_PARENT_1)

rule methylation_extraction_SE_allo_1:
	input:
		OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified1.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_se/{sample}_classified1.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bismark.cov.gz"
	benchmark:
		OUTPUT_DIR + "benchmark/extraction_se_allo1_{sample}.txt"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_se/" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/" ,
		genome = config["GENOME_PARENT_1"]
	conda:
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on SE bam files for allopolyploid species (GENOME_PARENT_2)

rule methylation_extraction_SE_allo_2:
	input:
		OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified2.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_se/{sample}_classified2.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bismark.cov.gz"
	benchmark:
		OUTPUT_DIR + "benchmark/extraction_se_allo2_{sample}.txt"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_se/" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/",
		genome = config["GENOME_PARENT_2"]
	conda:
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for parent species 1

rule methylation_extraction_PE_parent_1:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
	benchmark:
		OUTPUT_DIR + "benchmark/extraction_pe_p1_{sample}.txt"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_p1/",
		genome = config["GENOME_PARENT_1"]
	conda:
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for parent species 2

rule methylation_extraction_PE_parent_2:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
	benchmark:
		OUTPUT_DIR + "benchmark/extraction_pe_p2_{sample}.txt"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_p2/",
		genome = config["GENOME_PARENT_2"]
	conda:
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for allopolyploid species (GENOME_PARENT_1)

rule methylation_extraction_PE_allo_1:
	input:
		OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_classified1.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_classified1.ref.bedGraph.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
	benchmark:
 		OUTPUT_DIR + "benchmark/extraction_pe_allo1_{sample}.txt"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_1/",
		genome = config["GENOME_PARENT_1"]
	conda:
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for allopolyploid species (GENOME_PARENT_2)

rule methylation_extraction_PE_allo_2:
	input:
		 OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_classified2.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_classified2.ref.bedGraph.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
	benchmark:
		OUTPUT_DIR + "benchmark/extraction_pe_allo2_{sample}.txt"
	params:
		output = OUTPUT_DIR + "Bismark/extraction/{sample}_2/",
		genome = config["GENOME_PARENT_2"]
	conda:
		"envs/environment.yaml"
	threads:
		CORES
	shell:
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark coverage2cytosine on extraction output to obtain a single file with information about all cytosines (parent 1)

rule coverage2cytosine_1:
	input:
		f1 = OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_bismark_bt2.deduplicated.bismark.cov.gz")
	output:
		o1 = OUTPUT_DIR + "Bismark/extraction/{sample}_p1/{sample}.CX_report.txt"
	benchmark:
		OUTPUT_DIR + "benchmark/c2c_1_{sample}.txt"
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
		f2 = OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_bismark_bt2.deduplicated.bismark.cov.gz")
	output:
		o2 = OUTPUT_DIR + "Bismark/extraction/{sample}_p2/{sample}.CX_report.txt"
	benchmark:
		OUTPUT_DIR + "benchmark/c2c_2_{sample}.txt"
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
		f1 = OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_classified1.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/extraction/{sample}_se/{sample}_classified1.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bismark.cov.gz")),
		f2 = OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_classified2.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/extraction/{sample}_se/{sample}_classified2.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bismark.cov.gz"))
	output:
		o1 = OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}.CX_report.txt",
		o2 = OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}.CX_report.txt"
	benchmark:
		OUTPUT_DIR + "benchmark/c2c_allo_{sample}.txt"
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
	benchmark:
		OUTPUT_DIR + "benchmark/dmrseq_CG.txt"
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

# Define a function to create an input for dmrseq_CHG to include all the samples from the previous steps

def dmrseq_CHG_input_p1(wildcards):
	input = []
	input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CHG.cov", sample = samples.name[samples.origin == 'parent1'].values.tolist()))
	return input

def dmrseq_CHG_input_p2(wildcards):
	input = []
	input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CHG.cov", sample = samples.name[samples.origin == 'parent2'].values.tolist()))
	return input

def dmrseq_CHG_input_allo(wildcards):
	input = []
	input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CHG.cov", sample = samples.name[samples.origin == 'allopolyploid'].values.tolist()))
	return input

# Run dmrseq for CHG context

rule dmrseq_CHG:
	input:
		p1 = dmrseq_CHG_input_p1,
		p2 = dmrseq_CHG_input_p2,
		allo = dmrseq_CHG_input_allo
	output:
		comparison1 = OUTPUT_DIR + "DMR_analysis/dmrseq/CHG_context/parent1_v_allo.txt",
		comparison2 = OUTPUT_DIR + "DMR_analysis/dmrseq/CHG_context/parent2_v_allo.txt"
	benchmark:
		OUTPUT_DIR + "benchmark/dmrseq_CHG.txt"
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

# Define a function to create an input for dmrseq_CHH to include all the samples from the previous steps

def dmrseq_CHH_input_p1(wildcards):
	input = []
	input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CHH.cov", sample = samples.name[samples.origin == 'parent1'].values.tolist()))
	return input

def dmrseq_CHH_input_p2(wildcards):
	input = []
	input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CHH.cov", sample = samples.name[samples.origin == 'parent2'].values.tolist()))
	return input

def dmrseq_CHH_input_allo(wildcards):
	input = []
	input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CHH.cov", sample = samples.name[samples.origin == 'allopolyploid'].values.tolist()))
	return input

# Run dmrseq for CHH context

rule dmrseq_CHH:
	input:
		p1 = dmrseq_CHH_input_p1,
		p2 = dmrseq_CHH_input_p2,
		allo = dmrseq_CHH_input_allo
	output:
		comparison1 = OUTPUT_DIR + "DMR_analysis/dmrseq/CHH_context/parent1_v_allo.txt",
		comparison2 = OUTPUT_DIR + "DMR_analysis/dmrseq/CHH_context/parent2_v_allo.txt"
	benchmark:
		OUTPUT_DIR + "benchmark/dmrseq_CHH.txt"
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

## Rules for dmrseq comparing polyploids to polyploids or diploids to diploids based on conditions given in the metadata file

def dmrseq_CG_special_input_A(wildcards):
	input = []
	if config["POLYPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CG.cov", sample = samples.name[(samples.origin == 'allopolyploid') & (samples.condition == 'A')].values.tolist()))
	if config["DIPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CG.cov", sample = samples.name[(samples.origin == 'parent1') & (samples.condition == 'A')].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CG.cov", sample = samples.name[(samples.origin == 'parent2') & (samples.condition == 'A')].values.tolist()))
	return input

def dmrseq_CG_special_input_B(wildcards):
	input = []
	if config["POLYPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CG.cov", sample = samples.name[(samples.origin == 'allopolyploid') & (samples.condition == 'B')].values.tolist()))
	if config["DIPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CG.cov", sample = samples.name[(samples.origin == 'parent1') & (samples.condition == 'B')].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CG.cov", sample = samples.name[(samples.origin == 'parent2') & (samples.condition == 'B')].values.tolist()))
	return input

rule dmrseq_CG_special:
	input:
		condA = dmrseq_CG_special_input_A,
		condB = dmrseq_CG_special_input_B
	output:
		comparison = OUTPUT_DIR + "DMR_analysis/dmrseq/CG_context/A_v_B_diploid.txt" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/CG_context/A_v_B_polyploid.txt"
	benchmark:
		OUTPUT_DIR + "benchmark/dmrseq_CG_special.txt"
	threads:
		config["CORES_NUMBER"]
	conda:
		"envs/environment_R.yaml"
	shell:
		"Rscript scripts/dmrseq.R {n_samples_B} {n_samples_A} {output.comparison} {threads} {input.condB} {input.condA}"

def dmrseq_CHG_special_input_A(wildcards):
	input = []
	if config["POLYPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CHG.cov", sample = samples.name[(samples.origin == 'allopolyploid') & (samples.condition == 'A')].values.tolist()))
	if config["DIPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CHG.cov", sample = samples.name[(samples.origin == 'parent1') & (samples.condition == 'A')].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CHG.cov", sample = samples.name[(samples.origin == 'parent2') & (samples.condition == 'A')].values.tolist()))
	return input

def dmrseq_CHG_special_input_B(wildcards):
	input = []
	if config["POLYPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CHG.cov", sample = samples.name[(samples.origin == 'allopolyploid') & (samples.condition == 'B')].values.tolist()))
	if config["DIPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CHG.cov", sample = samples.name[(samples.origin == 'parent1') & (samples.condition == 'B')].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CHG.cov", sample = samples.name[(samples.origin == 'parent2') & (samples.condition == 'B')].values.tolist()))
	return input

rule dmrseq_CHG_special:
	input:
		condA = dmrseq_CHG_special_input_A,
		condB = dmrseq_CHG_special_input_B
	output:
		comparison = OUTPUT_DIR + "DMR_analysis/dmrseq/CHG_context/A_v_B_diploid.txt" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/CHG_context/A_v_B_polyploid.txt"
	threads:
		config["CORES_NUMBER"]
	conda:
		"envs/environment_R.yaml"
	shell:
		"Rscript scripts/dmrseq.R {n_samples_B} {n_samples_A} {output.comparison} {threads} {input.condB} {input.condA}"

def dmrseq_CHH_special_input_A(wildcards):
	input = []
	if config["POLYPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CHH.cov", sample = samples.name[(samples.origin == 'allopolyploid') & (samples.condition == 'A')].values.tolist()))
	if config["DIPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CHH.cov", sample = samples.name[(samples.origin == 'parent1') & (samples.condition == 'A')].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CHH.cov", sample = samples.name[(samples.origin == 'parent2') & (samples.condition == 'A')].values.tolist()))
	return input

def dmrseq_CHH_special_input_B(wildcards):
	input = []
	if config["POLYPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/allopolyploid/{sample}_CHH.cov", sample = samples.name[(samples.origin == 'allopolyploid') & (samples.condition == 'B')].values.tolist()))
	if config["DIPLOID_ONLY"]:
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent1/{sample}_CHH.cov", sample = samples.name[(samples.origin == 'parent1') & (samples.condition == 'B')].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "DMR_analysis/context_separation/parent2/{sample}_CHH.cov", sample = samples.name[(samples.origin == 'parent2') & (samples.condition == 'B')].values.tolist()))
	return input

rule dmrseq_CHH_special:
	input:
		condA = dmrseq_CHH_special_input_A,
		condB = dmrseq_CHH_special_input_B
	output:
		comparison = OUTPUT_DIR + "DMR_analysis/dmrseq/CHH_context/A_v_B_diploid.txt" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/CHH_context/A_v_B_polyploid.txt"
	benchmark:
		OUTPUT_DIR + "benchmark/dmrseq_CHG_special.txt"
	threads:
		config["CORES_NUMBER"]
	conda:
		"envs/environment_R.yaml"
	shell:
		"Rscript scripts/dmrseq.R {n_samples_B} {n_samples_A} {output.comparison} {threads} {input.condB} {input.condA}"


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

## Define a function to create input directories based on the settings in the config file

def multiqc_params(wildcards):
	param = [OUTPUT_DIR + "FastQC"]
	if config["RUN_TRIMMING"]:
		param.append(OUTPUT_DIR + "FASTQtrimmed")
	if config["RUN_BISMARK"]:
		param.append(OUTPUT_DIR + "Bismark")
		param.append(OUTPUT_DIR + "qualimap")
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
	benchmark:
		OUTPUT_DIR + "benchmark/multiqc.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.inputdir} -f -o {params.multiqcdir}"

# The following rules are for the dowstream analyses of the workflow

# The first downstream rule takes the dmrseq output and creates a bed file with all significant DMRs (specifically DMRs with q-values < 0.05)

rule dm_regions_bed:
	input:
		i1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo.txt",
		i2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo.txt"
	output:
		o1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo_sig_sorted.bed",
		o2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo_sig_sorted.bed"
	params:
		p1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo_sig",
		p2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo_sig"
	conda:
		"envs/environment_downstream.yaml"
	shell:
		"Rscript scripts/significantGenesToBed.R {input.i1} {params.p1};"
		"Rscript scripts/significantGenesToBed.R {input.i2} {params.p2}"

# The first downstream rule for special mode: diploid vs diploid or polyploid vs polyploid. Read above for a short explanation of the rule.

rule dm_regions_bed_special:
	input:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid.txt" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid.txt"
	output:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid_sig_sorted.txt" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid_sig_sorted.txt"
	params:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid_sig" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid_sig"
	conda:
		"envs/environment_downstream.yaml"
	shell:
		"Rscript scripts/significantGenesToBed.R {input} {params}"

# The second downsteam rule checks the gene regions provided in the annotation file and finds any overlaps with DMRs. The output includes all genes showing an overlap with the overlapping DMR.

rule bedtools_intersect:
	input:
		i1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo_sig_sorted.bed",
		i2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo_sig_sorted.bed"
	output:
		o1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo_genes_overlap.txt",
		o2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo_genes_overlap.txt"
	params:
		anno1 = config["ANNOTATION_PARENT_1"],
		anno2 = config["ANNOTATION_PARENT_2"]
	conda:
		"envs/environment_downstream.yaml"
	shell:
		"bedtools intersect -a {params.anno1} -b {input.i1} -wo > {output.o1};"
		"bedtools intersect -a {params.anno2} -b {input.i2} -wo > {output.o2}"

# The second downstream rule for special mode: diploid vs diploid or polyploid vs polyploid. Read above for a short explanation of the rule.

rule bedtools_intersect_special:
	input:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid_sig_sorted.txt" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid_sig_sorted.txt"
	output:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid_genes_overlap.txt" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid_genes_overlap.txt"
	params:
		anno1 = config["ANNOTATION_PARENT_1"],
		anno2 = config["ANNOTATION_PARENT_2"]
	conda:
		"envs/environment_downstream.yaml"
	shell:
		"bedtools intersect -a {params.anno1} -b {input} -wo > {output}" if (sum(samples.origin == "parent1") > 0) else ("bedtools intersect -a {params.anno2} -b {input} -wo > {output}" if (sum(samples.origin == "parent2") > 0) else "bedtools intersect -a {params.anno1} -b {input} -wo > {output}; bedtools intersect -a {params.anno2} -b {input} -wo >> {output}")

# The third and final rule for downstream analyses. With overlap information and DMR information, we generate a summary file including gene ID of the genes overlapping with DMRs, all ranges for gene regions and DMRs and methylation status. Methylation status is based on the statistics given by the dmrseq output, taking condition 'A' as reference (i.e. increase means increase compared to condition 'A')

rule dmr_genes:
	input:
		i1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo_genes_overlap.txt",
		i2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo_genes_overlap.txt",
		dm1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo.txt",
		dm2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo.txt"
	output:
		o1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent1_v_allo_{context}.txt",
		o2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent2_v_allo_{context}.txt"
	params:
		geneID1 = config["GENE_ID_PARENT_1"],
		geneID2 = config["GENE_ID_PARENT_2"],
		o1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent1_v_allo_{context}",
		o2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent2_v_allo_{context}"
	conda:
		"envs/environment_downstream.yaml"
	shell:
		"Rscript scripts/DMGeneSummary.R {input.i1} {input.dm1} {params.geneID1} {params.o1};"
		"Rscript scripts/DMGeneSummary.R {input.i2} {input.dm2} {params.geneID2} {params.o2}"

# The third downstream rules for special modes: diploid vs diploid (first below) or polyploid vs polyploid (second below). Read above for a short explanation of the rule.

rule dmr_genes_special_diploid:
	input:
		i1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid_genes_overlap.txt",
		dm1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid.txt"
	output:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_diploid_{context}.txt"
	params:
		geneID1 = config["GENE_ID_PARENT_1"],
		geneID2 = config["GENE_ID_PARENT_2"],
		o1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_diploid_{context}"
	conda:
		"envs/environment_downstream.yaml"
	shell:
		"Rscript scripts/DMGeneSummary.R {input.i1} {input.dm1} {params.geneID1} {params.o1}" if (sum(samples.origin == "parent1") > 0) else "Rscript scripts/DMGeneSummary.R {input.i1} {input.dm1} {params.geneID2} {params.o1}"

rule dmr_genes_special_polyploid:
	input:
		i1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid_genes_overlap.txt",
		dm1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid.txt"
	output:
		first = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_1.txt",
		second = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_2.txt"
	params:
		geneID1 = config["GENE_ID_PARENT_1"],
		geneID2 = config["GENE_ID_PARENT_2"],
		o1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_diploid_{context}_1",
		o2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_diploid_{context}_2"
	conda:
		"envs/environment_downstream.yaml"
	shell:
		"Rscript scripts/DMGeneSummary.R {input.i1} {input.dm1} {params.geneID1} {params.o1};"
		"Rscript scripts/DMGeneSummary.R {input.i1} {input.dm1} {params.geneID2} {params.o2}"
