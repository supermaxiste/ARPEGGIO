import os.path

# This file includes all rules related to read trimming
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: Trim Galore!

## Run TrimGalore to trim reads, there are two rules depending on paired or single-end status
## For single-end reads


rule trim_galore_se:
    input:
        fastq1=f"{RAW_DATA_DIR}{{sample}}.{str(config['RAW_DATA_EXTENSION'])}.gz",
    output:
        o1=f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_trimmed.fq.gz",
    log:
        f"{OUTPUT_DIR}logs/trim_galore_{{sample}}_se.log",
    params:
        FASTQtrimmeddir=lambda w, output: os.path.split(output.o1)[0],
        trim_5_r1=config["CLIP_5_R1"],
        trim_3_r1=config["CLIP_3_R1"],
        trim_cores=TRIM_CORES,
    threads: (TRIM_CORES * 3 + 3) * 2 if config["IS_PAIRED"] else (TRIM_CORES * 3 + 3)
    benchmark:
        f"{OUTPUT_DIR}benchmark/trim_se_{{sample}}.txt"
    conda:
        "../envs/environment.yaml"
    shell:
        "trim_galore -q 20 --clip_R1 {params.trim_5_r1} --phred33 --length 20 -j {params.trim_cores} -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq1} 2> {log}" if config[
        "RUN_TRIMMING"
        ] and config[
            "TRIM_5_ONLY"
        ] else (
            "trim_galore -q 20 --clip_R1 {params.trim_5_r1}  --three_prime_clip_R1 {params.trim_3_r1} --phred33 --length 20 -j {params.trim_cores} -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq1} 2> {log}"
            if config["RUN_TRIMMING"] and config["TRIM_3_ONLY"]
            else (
                "trim_galore -q 20 --clip_R1 {params.trim_5_r1}  --three_prime_clip_R1 {params.trim_3_r1} --phred33 --length 20 -j {params.trim_cores} -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq1} 2> {log}"
                if config["RUN_TRIMMING"]
                else "trim_galore -q 20 --phred33 --length 20 -j {params.trim_cores} -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq1} 2> {log}"
            )
        )


## For paired-end reads


rule trim_galore_pe:
    input:
        fastq1=f"{RAW_DATA_DIR}{{sample}}_{str(config['PAIR_1'])}.{str(config['RAW_DATA_EXTENSION'])}.gz",
        fastq2=f"{RAW_DATA_DIR}{{sample}}_{str(config['PAIR_2'])}.{str(config['RAW_DATA_EXTENSION'])}.gz",
    output:
        o1=f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_1'])}_val_1.fq.gz",
        o2=f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_2'])}_val_2.fq.gz",
    log:
        f"{OUTPUT_DIR}logs/trim_galore_{{sample}}_pe.log",
    params:
        FASTQtrimmeddir=lambda w, output: os.path.split(output.o1)[0],
        trim_5_r1=config["CLIP_5_R1"],
        trim_5_r2=config["CLIP_5_R2"],
        trim_3_r1=config["CLIP_3_R1"],
        trim_3_r2=config["CLIP_3_R2"],
        trim_cores=TRIM_CORES,
    threads: (TRIM_CORES * 3 + 3) * 2 if config["IS_PAIRED"] else (TRIM_CORES * 3 + 3)
    benchmark:
        f"{OUTPUT_DIR}benchmark/trim_pe_{{sample}}.txt"
    conda:
        "../envs/environment.yaml"
    shell:
        "trim_galore -q 20 --clip_R1 {params.trim_5_r1} --clip_R2 {params.trim_5_r2} --phred33 --length 20 -j {params.trim_cores} -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2} 2> {log}" if config[
        "RUN_TRIMMING"
        ] and config[
            "TRIM_5_ONLY"
        ] else (
            "trim_galore -q 20 --three_prime_clip_R1 {params.trim_3_r1} --three_prime_clip_R2 {params.trim_3_r2} --phred33 --length 20 -j {params.trim_cores} -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2} 2> {log}"
            if config["RUN_TRIMMING"] and config["TRIM_3_ONLY"]
            else (
                "trim_galore -q 20 --clip_R1 {params.trim_5_r1} --clip_R2 {params.trim_5_r2} --three_prime_clip_R1 {params.trim_3_r1} --three_prime_clip_R2 {params.trim_3_r2} --phred33 --length 20 -j {params.trim_cores} -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2} 2> {log}"
                if config["RUN_TRIMMING"]
                else "trim_galore -q 20 --phred33 --length 20 -j {params.trim_cores} -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2} 2> {log}"
            )
        )
