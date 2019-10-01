# This file includes all rules related to checking the quality of the data
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: FastQC,

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

# sort bam files for multiqc (parent1)

rule bam_sorting_p1:
	input:
		p1 = OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated.bam")
	output:
		o1 = OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_R1_val_1_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_trimmed_bismark_bt2.deduplicated_sorted.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_1/1.{sample}_bismark_bt2.deduplicated_sorted.bam")
	conda:
		"envs/environment.yaml"
	shell:
		"samtools sort {input.p1} > {output.o1}"

# sort bam files for multiqc (parent2)

rule bam_sorting_p2:
	input:
		p2 = OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated.bam"),
	output:
		o2 = OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_R1_val_1_bismark_bt2_pe.deduplicated_sorted.bam" if config["IS_PAIRED"] else (OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_trimmed_bismark_bt2.deduplicated_sorted.bam" if config["RUN_TRIMMING"] else OUTPUT_DIR + "Bismark/deduplication/{sample}_2/2.{sample}_bismark_bt2.deduplicated_sorted.bam")
	conda:
		"envs/environment.yaml"
	shell:
		"samtools sort {input.p2} > {output.o2}"

# sort bam files for multiqc (allopolyploid)

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
