# This file includes all functions to define inputs to other rules
# SE = Single-end
# PE = Paired-end

## Define a function to create an input for MultiQC to include all the settings specified in the config file.

def multiqc_input(wildcards):
	input = []
	if config["CONVERSION_CHECK"]:
		if config["IS_PAIRED"]:
			if config["RUN_TRIMMING"]:
				## paired and trimmed
				input.extend(expand(OUTPUT_DIR + "Conversion_efficiency/{sample}/cc.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.bam", sample = samples.name[samples.type == 'PE'].values.tolist()))
			else:
				## paired and not trimmed
				input.extend(expand(OUTPUT_DIR + "Conversion_efficiency/{sample}/cc.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.bam", sample = samples.name[samples.type == 'PE'].values.tolist()))
		else:
			if config["RUN_TRIMMING"]:
				## not paired and trimmed
				input.extend(expand(OUTPUT_DIR + "Conversion_efficiency/{sample}/cc.{sample}_trimmed_bismark_bt2.bam", sample = samples.name[samples.type == 'SE'].values.tolist()))
			else:
				## not paired and not trimmed
				input.extend(expand(OUTPUT_DIR + "Conversion_efficiency/{sample}/cc.{sample}_bismark_bt2.bam", sample = samples.name[samples.type == 'SE'].values.tolist()))
	if config["IS_PAIRED"]:
		input.extend(expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_1"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
		input.extend(expand(OUTPUT_DIR + "FastQC/{sample}_" + str(config["PAIR_2"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
		if config["CONVERSION_CHECK"]:
			if config["RUN_TRIMMING"]:
				input.extend(expand(OUTPUT_DIR + "Conversion_efficiency/{sample}/cc.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.bam", sample = samples.name[samples.type == 'PE'].values.tolist()))
			else:
				input.extend(expand(OUTPUT_DIR + "Conversion_efficiency/{sample}/cc.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.bam", sample = samples.name[samples.type == 'PE'].values.tolist()))
	else:
		input.extend(expand(OUTPUT_DIR + "FastQC/{sample}_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist()))
		if config["CONVERSION_CHECK"]:
			if config["RUN_TRIMMING"]:
				input.extend(expand(OUTPUT_DIR + "Conversion_efficiency/{sample}/cc.{sample}_trimmed_bismark_bt2.bam", sample = samples.name[samples.type == 'PE'].values.tolist()))
			else:
				input.extend(expand(OUTPUT_DIR + "Conversion_efficiency/{sample}/cc.{sample}_bismark_bt2.bam", sample = samples.name[samples.type == 'PE'].values.tolist()))
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
			if config["RUN_TRIMMING"]:
				## paired and trimmed
				## alignment
				input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent2')].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent1')].values.tolist()))
				## deduplication
				input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent2')].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent1')].values.tolist()))
			else:
				## paired and not trimmed
				## alignment
				input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent2')].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "Bismark/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent1')].values.tolist()))
				## deduplication
				input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent2')].values.tolist()))
				input.extend(expand(OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam", sample = samples.name[(samples.type == 'PE') & (samples.origin != 'parent1')].values.tolist()))
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
	if config["CONVERSION_CHECK"]:
		param.append(OUTPUT_DIR + "Conversion_efficiency")
	if config["RUN_TRIMMING"]:
		param.append(OUTPUT_DIR + "FASTQtrimmed")
	if config["RUN_BISMARK"]:
		param.append(OUTPUT_DIR + "Bismark")
		param.append(OUTPUT_DIR + "qualimap")
	return param

# General rule to run all analyses depending on the config settings

def dmr_input(wildcards):
	input = []
	if config["RUN_DOWNSTREAM"]:
		if config["ONLY_CG_CONTEXT"]:
			if config["POLYPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_1.txt", context = ["CG_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_2.txt", context = ["CG_context"]))
			elif config["DIPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_diploid_{context}.txt", context = ["CG_context"]))
			else:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent1_v_allo_{context}.txt", context = ["CG_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent2_v_allo_{context}.txt", context = ["CG_context"]))
		else:
			if config["POLYPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_1.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_2.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
			elif config["DIPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_diploid_{context}.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
			else:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent1_v_allo_{context}.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent2_v_allo_{context}.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
	if config["RUN_DMR_ANALYSIS"]:
		if config["ONLY_CG_CONTEXT"]:
			if config["POLYPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid.txt", context = ["CG_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid.txt", context = ["CG_context"]))
			elif config["DIPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid.txt", context = ["CG_context"]))
			else:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo.txt", context = ["CG_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo.txt", context = ["CG_context"]))
		else:
			if config["POLYPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
			elif config["DIPLOID_ONLY"]:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
			else:
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
				input.extend(expand(OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo.txt", context = ["CG_context", "CHG_context", "CHH_context"]))
	return input

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
