# This file includes all rules related to downstream analyses of DMRs
# DMR = Differentially Methylated Region
# Tools used in the rules: R scripts, bedtools

# R scripts used: significantGenesToBed.R, DMGeneSummary.R

# The first downstream rule takes the dmrseq output and creates a bed file with all significant DMRs (specifically DMRs with q-values < 0.05)

rule dm_regions_bed:
	input:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent{one_or_two}_v_allo.txt"
	output:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent{one_or_two}_v_allo_sig_sorted.bed"
	params:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent{one_or_two}_v_allo_sig"
	conda:
		"../envs/environment_downstream.yaml"
	shell:
		"Rscript scripts/significantGenesToBed.R {input} {params}"

# The first downstream rule for special mode: diploid vs diploid or polyploid vs polyploid. Read above for a short explanation of the rule.

rule dm_regions_bed_special:
	input:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid.txt" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid.txt"
	output:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid_sig_sorted.bed" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid_sig_sorted.bed"
	params:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid_sig" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid_sig"
	conda:
		"../envs/environment_downstream.yaml"
	shell:
		"Rscript scripts/significantGenesToBed.R {input} {params}"

# The second downsteam rule checks the gene regions provided in the annotation file and finds any overlaps with DMRs. The output includes all genes showing an overlap with the overlapping DMR.

rule bedtools_intersect_1:
	input:
		i1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo_sig_sorted.bed"
	output:
		o1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo_genes_overlap.txt"
	params:
		anno1 = config["ANNOTATION_PARENT_1"]
	conda:
		"../envs/environment_downstream.yaml"
	shell:
		"bedtools intersect -a {params.anno1} -b {input.i1} -wo > {output.o1}"

rule bedtools_intersect_2:
	input:
		i2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo_sig_sorted.bed"
	output:
		o2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo_genes_overlap.txt"
	params:
		anno2 = config["ANNOTATION_PARENT_2"]
	conda:
		"../envs/environment_downstream.yaml"
	shell:
		"bedtools intersect -a {params.anno2} -b {input.i2} -wo > {output.o2}"

# The second downstream rule for special mode: diploid vs diploid or polyploid vs polyploid. Read above for a short explanation of the rule.

rule bedtools_intersect_special:
	input:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid_sig_sorted.bed" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid_sig_sorted.bed"
	output:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_diploid_genes_overlap.txt" if config["DIPLOID_ONLY"] else OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid_genes_overlap.txt"
	params:
		anno1 = config["ANNOTATION_PARENT_1"],
		anno2 = config["ANNOTATION_PARENT_2"]
	conda:
		"../envs/environment_downstream.yaml"
	shell:
		"bedtools intersect -a {params.anno1} -b {input} -wo > {output}" if (sum(samples.origin == "parent1") > 0) else ("bedtools intersect -a {params.anno2} -b {input} -wo > {output}" if (sum(samples.origin == "parent2") > 0) else "bedtools intersect -a {params.anno1} -b {input} -wo > {output}; bedtools intersect -a {params.anno2} -b {input} -wo >> {output}")

# The third and final rule for downstream analyses. With overlap information and DMR information, we generate a summary file including gene ID of the genes overlapping with DMRs, all ranges for gene regions and DMRs and methylation status. Methylation status is based on the statistics given by the dmrseq output, taking condition 'A' as reference (i.e. increase means increase compared to condition 'A')

rule dmr_genes_1:
	input:
		i1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo_genes_overlap.txt",
		dm1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent1_v_allo.txt"
	output:
		o1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent1_v_allo_{context}.txt"
	params:
		geneID1 = config["GENE_ID_PARENT_1"],
		o1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent1_v_allo_{context}"
	conda:
		"../envs/environment_downstream.yaml"
	shell:
		"Rscript scripts/DMGeneSummary.R {input.i1} {input.dm1} {params.geneID1} {params.o1}"

rule dmr_genes_2:
	input:
		i2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo_genes_overlap.txt",
		dm2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/parent2_v_allo.txt"
	output:
		o2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent2_v_allo_{context}.txt"
	params:
		geneID2 = config["GENE_ID_PARENT_2"],
		o2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_parent2_v_allo_{context}"
	conda:
		"../envs/environment_downstream.yaml"
	shell:
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
		"../envs/environment_downstream.yaml"
	shell:
		"Rscript scripts/DMGeneSummary.R {input.i1} {input.dm1} {params.geneID1} {params.o1}" if (sum(samples.origin == "parent1") > 0) else "Rscript scripts/DMGeneSummary.R {input.i1} {input.dm1} {params.geneID2} {params.o1}"

rule dmr_genes_special_polyploid_1:
	input:
		i1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid_genes_overlap.txt",
		dm1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid.txt"
	output:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_1.txt"
	params:
		geneID1 = config["GENE_ID_PARENT_1"],
		o1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_1"
	conda:
		"../envs/environment_downstream.yaml"
	shell:
		"Rscript scripts/DMGeneSummary.R {input.i1} {input.dm1} {params.geneID1} {params.o1}"

rule dmr_genes_special_polyploid_2:
	input:
		i1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid_genes_overlap.txt",
		dm1 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/A_v_B_polyploid.txt"
	output:
		OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_2.txt"
	params:
		geneID2 = config["GENE_ID_PARENT_2"],
		o2 = OUTPUT_DIR + "DMR_analysis/dmrseq/{context}/DM_genes_A_v_B_polyploid_{context}_2"
	conda:
		"../envs/environment_downstream.yaml"
	shell:
		"Rscript scripts/DMGeneSummary.R {input.i1} {input.dm1} {params.geneID2} {params.o2}"
