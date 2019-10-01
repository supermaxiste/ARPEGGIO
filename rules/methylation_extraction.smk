# This file includes all rules related to methylation extraction after alignment
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: Bismark


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
		"../envs/environment.yaml"
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
		"../envs/environment.yaml"
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
		"../envs/environment.yaml"
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
		"../envs/environment.yaml"
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
		"../envs/environment.yaml"
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
		"../envs/environment.yaml"
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
		"../envs/environment.yaml"
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
		"../envs/environment.yaml"
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
		"../envs/environment.yaml"
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
		"../envs/environment.yaml"
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
		"../envs/environment.yaml"
	shell:
		"""
		coverage2cytosine -CX --genome_folder {params.genome1} -o {params.filename1} {input.f1}
		coverage2cytosine -CX --genome_folder {params.genome2} -o {params.filename2} {input.f2}
		"""
