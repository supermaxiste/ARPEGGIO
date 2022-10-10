import os.path

# This file includes all rules related to read deduplication after alignment
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: Bismark

## Run deduplication of the alignments to remove duplicated reads for SE reads (GENOME_PARENT_1)


rule deduplication_SE:
    input:
        f"{OUTPUT_DIR}Bismark/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_trimmed_bismark_bt2.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_bismark_bt2.bam",
    output:
        f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_trimmed_bismark_bt2.deduplicated.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_bismark_bt2.deduplicated.bam",
    log:
        f"logs/deduplication_{{sample}}_{{one_or_two}}_SE.log",
    params:
        f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/",
    benchmark:
        f"{OUTPUT_DIR}benchmark/dedup_se{{one_or_two}}_{{sample}}.txt"
    conda:
        "../envs/environment.yaml"
    shell:
        "deduplicate_bismark -s --output_dir {params} --bam {input}"


## Run deduplication of the alignments to remove duplicated reads for PE reads (GENOME_PARENT_1)


rule deduplication_PE:
    input:
        f"{OUTPUT_DIR}Bismark/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.bam",
    output:
        f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bam",
    log:
        f"logs/deduplication_{{sample}}_{{one_or_two}}_PE.log",
    params:
        f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/",
    benchmark:
        f"{OUTPUT_DIR}benchmark/dedup_pe{{one_or_two}}_{{sample}}.txt"
    conda:
        "../envs/environment.yaml"
    shell:
        "deduplicate_bismark -p --output_dir {params} --bam {input}"
