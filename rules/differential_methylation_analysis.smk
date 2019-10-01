# This file includes all rules related to differential methylation analysis
# DMR = Differentially Methylated Region
# Tools used in the rules: R scripts, dmrseq

# R scripts used: CoverageFileGeneratorComplete.R, dmrseq.R

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
