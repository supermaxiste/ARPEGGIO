# This file includes all rules related to methylation extraction after alignment
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: Bismark


## Run Bismark methylation extraction on SE bam files for parent species 1

rule methylation_extraction_SE_parent_1:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_bismark_bt2.deduplicated.bismark.cov.gz",
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CpG_context_1.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CpG_context_1.{sample}_bismark_bt2.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CHG_context_1.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CHG_context_1.{sample}_bismark_bt2.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CHH_context_1.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CHH_context_1.{sample}_bismark_bt2.deduplicated.txt")
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
		"bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config["UNFINISHED_GENOME"] else "bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"


## Run Bismark methylation extraction on SE bam files for parent species 2

rule methylation_extraction_SE_parent_2:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_bismark_bt2.deduplicated.bismark.cov.gz",
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CpG_context_2.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CHG_context_2.{sample}_bismark_bt2.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CHG_context_2.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CpG_context_2.{sample}_bismark_bt2.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CHH_context_2.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CHH_context_2.{sample}_bismark_bt2.deduplicated.txt")
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
		"bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config["UNFINISHED_GENOME"] else "bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on SE bam files for allopolyploid species (GENOME_PARENT_1)

rule methylation_extraction_SE_allo_1:
	input:
		OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified1.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_se/{sample}_classified1.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bismark.cov.gz",
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CpG_context_{sample}_classified1.ref.txt") if config["RUN_READ_SORTING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CpG_context_1.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CpG_context_1.{sample}_bismark_bt2.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHG_context_{sample}_classified1.ref.txt") if config["RUN_READ_SORTING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHG_context_1.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHG_context_1.{sample}_bismark_bt2.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHH_context_{sample}_classified1.ref.txt") if config["RUN_READ_SORTING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHH_context_1.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHH_context_1.{sample}_bismark_bt2.deduplicated.txt")
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
		"bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config["UNFINISHED_GENOME"] else "bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on SE bam files for allopolyploid species (GENOME_PARENT_2)

rule methylation_extraction_SE_allo_2:
	input:
		OUTPUT_DIR + "read_sorting/{sample}_se/{sample}_classified2.ref.bam" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_se/{sample}_classified2.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bismark.cov.gz",
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CpG_context_{sample}_classified2.ref.txt") if config["RUN_READ_SORTING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CpG_context_2.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CpG_context_2.{sample}_bismark_bt2.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHG_context_{sample}_classified2.ref.txt") if config["RUN_READ_SORTING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHG_context_2.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHG_context_2.{sample}_bismark_bt2.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHH_context_{sample}_classified2.ref.txt") if config["RUN_READ_SORTING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHH_context_2.{sample}_trimmed_bismark_bt2.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_se/CHH_context_2.{sample}_bismark_bt2.deduplicated.txt")
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
		"bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config["UNFINISHED_GENOME"] else "bismark_methylation_extractor -s -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for parent species 1

rule methylation_extraction_PE_parent_1:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bedGraph.gz",
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CpG_context_1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CpG_context_1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CHG_context_1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CHG_context_1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CHH_context_1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p1/CHH_context_1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt")
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
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config["UNFINISHED_GENOME"] else "bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for parent species 2

rule methylation_extraction_PE_parent_2:
	input:
		OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bedGraph.gz",
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CpG_context_2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CpG_context_2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CHG_context_2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CHG_context_2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CHH_context_2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_p2/CHH_context_2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt")
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
		"bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config["UNFINISHED_GENOME"] else "bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for allopolyploid species (GENOME_PARENT_1)

rule methylation_extraction_PE_allo_1:
	input:
		OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified1.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_classified1.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_classified1.ref.bedGraph.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bedGraph.gz",
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_1/CpG_context_{sample}_classified1.ref.txt") if config["IS_PAIRED"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_1/CpG_context_1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_1/CpG_context_1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_1/CHG_context_{sample}_classified1.ref.txt") if config["IS_PAIRED"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_1/CHG_context_1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_1/CHG_context_1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_1/CHH_context_{sample}_classified1.ref.txt") if config["IS_PAIRED"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_1/CHH_context_1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_1/CHH_context_1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt")
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
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config["UNFINISHED_GENOME"] else "bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark methylation extraction on PE bam files for allopolyploid species (GENOME_PARENT_2)

rule methylation_extraction_PE_allo_2:
	input:
		 OUTPUT_DIR + "read_sorting/{sample}/{sample}_classified2.ref.bam" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bam"
	output:
		OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_classified2.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bismark.cov.gz",
		OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_classified2.ref.bedGraph.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bedGraph.gz",
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_2/CpG_context_{sample}_classified2.ref.txt") if config["IS_PAIRED"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_2/CpG_context_2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_2/CpG_context_2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_2/CHG_context_{sample}_classified2.ref.txt") if config["IS_PAIRED"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_2/CHG_context_2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_2/CHG_context_2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt"),
		temp(OUTPUT_DIR + "Bismark/extraction/{sample}_2/CHH_context_{sample}_classified2.ref.txt") if config["IS_PAIRED"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_2/CHH_context_2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.txt") if config["RUN_TRIMMING"] else temp(OUTPUT_DIR + "Bismark/extraction/{sample}_2/CHH_context_2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.txt")
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
		"bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config["UNFINISHED_GENOME"] else "bismark_methylation_extractor -p -o {params.output} --genome_folder {params.genome} --multicore {BISMARK_CORES} --no_overlap --comprehensive --bedGraph --CX {input}"

## Run Bismark coverage2cytosine on extraction output to obtain a single file with information about all cytosines (parent 1)

rule coverage2cytosine_1:
	input:
		f1 = OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p1/1.{sample}_bismark_bt2.deduplicated.bismark.cov.gz")
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
		f2 = OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_p2/2.{sample}_bismark_bt2.deduplicated.bismark.cov.gz")
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

rule coverage2cytosine_allo_1:
	input:
		f1 = OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}_classified1.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/extraction/{sample}_se/{sample}_classified1.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bismark.cov.gz"
	output:
		o1 = OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}.CX_report.txt"
	benchmark:
		OUTPUT_DIR + "benchmark/c2c_allo_{sample}.txt"
	params:
		genome1 = config["GENOME_PARENT_1"],
		filename1 = OUTPUT_DIR + "Bismark/extraction/{sample}_1/{sample}"
	conda:
		"../envs/environment.yaml"
	shell:
		"coverage2cytosine -CX --genome_folder {params.genome1} -o {params.filename1} {input.f1}"

## Run Bismark coverage2cytosine on extraction output to obtain a single file with information about all cytosines (allopolyploid)

rule coverage2cytosine_allo_2:
	input:
		f2 = OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}_classified2.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] and config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["IS_PAIRED"] and config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_bismark_bt2_pe.deduplicated.bismark.cov.gz" if config["IS_PAIRED"] else OUTPUT_DIR + "Bismark/extraction/{sample}_se/{sample}_classified2.ref.bismark.cov.gz" if config["RUN_READ_SORTING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bismark.cov.gz" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/extraction/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bismark.cov.gz"
	output:
		o2 = OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}.CX_report.txt"
	benchmark:
		OUTPUT_DIR + "benchmark/c2c_allo_{sample}.txt"
	params:
		genome2 = config["GENOME_PARENT_2"],
		filename2 = OUTPUT_DIR + "Bismark/extraction/{sample}_2/{sample}"
	conda:
		"../envs/environment.yaml"
	shell:
		"coverage2cytosine -CX --genome_folder {params.genome2} -o {params.filename2} {input.f2}"
