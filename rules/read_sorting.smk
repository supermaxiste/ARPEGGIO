# This file includes all rules related to read sorting after alignment
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: EAGLE-RC

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
		reads1 = OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam",
		reads2 = OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_" + str(config["PAIR_1"]) + "_val_1_bismark_bt2_pe.deduplicated.bam"
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
