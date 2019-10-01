# This file includes all rules related to read deduplication after alignment
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: Bismark

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
