# This file includes all rules related to read trimming
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: Trim Galore!

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
		"../envs/environment.yaml"
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
		"../envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log};"
		"trim_galore -q 20 --clip_R1 {params.trim_5_r1} --clip_R2 {params.trim_5_r2} --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2}" if config["RUN_TRIMMING"] and config["TRIM_5_ONLY"] else ("trim_galore -q 20 --three_prime_clip_R1 {params.trim_3_r1} --three_prime_clip_R2 {params.trim_3_r2} --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2}" if config["RUN_TRIMMING"] and config["TRIM_3_ONLY"] else ("trim_galore -q 20 --clip_R1 {params.trim_5_r1} --clip_R2 {params.trim_5_r2} --three_prime_clip_R1 {params.trim_3_r1} --three_prime_clip_R2 {params.trim_3_r2} --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2}" if config["RUN_TRIMMING"] else "trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2}"))
