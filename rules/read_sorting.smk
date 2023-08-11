import os.path

# This file includes all rules related to read sorting after alignment
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: EAGLE-RC

## Run EAGLE-RC to classify reads to the most probable genome for SE reads

EAGLE = ".eagle/bin/eagle-rc"


rule read_sorting_SE:
    input:
        eagle_bin=EAGLE,
        reads1=f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_trimmed_bismark_bt2.deduplicated.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_bismark_bt2.deduplicated.bam",
        reads2=f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_trimmed_bismark_bt2.deduplicated.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_bismark_bt2.deduplicated.bam",
    output:
        o1=f"{OUTPUT_DIR}read_sorting/{{sample}}_se/{{sample}}_classified1.ref.bam",
        o2=f"{OUTPUT_DIR}read_sorting/{{sample}}_se/{{sample}}_classified2.ref.bam",
    log:
        f"{OUTPUT_DIR}logs/read_sorting_{{sample}}_SE.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/readsorting_se_{{sample}}.txt"
    params:
        list=f"{OUTPUT_DIR}read_sorting/{{sample}}_se/{{sample}}_classified_reads.list",
        phred="--phred64" if config["PHRED_SCORE_64"] else "",
        genome1=config["ASSEMBLY_PARENT_1"],
        genome2=config["ASSEMBLY_PARENT_2"],
        output=lambda w, output: os.path.splitext(os.path.splitext(output.o1)[0])[0][
            :-1
        ],
    shell:
        "{input.eagle_bin} --ngi {params.phred} --ref1={params.genome1} --bam1={input.reads1} --ref2={params.genome2} --bam2={input.reads2} -o {params.output} --bs=3 > {params.list} 2>&1 {log}"


## Run EAGLE-RC to classify reads to the most probable genome for PE reads


rule read_sorting_PE:
    input:
        eagle_bin=EAGLE,
        reads1=f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bam",
        reads2=f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bam",
    output:
        o1=f"{OUTPUT_DIR}read_sorting/{{sample}}/{{sample}}_classified1.ref.bam",
        o2=f"{OUTPUT_DIR}read_sorting/{{sample}}/{{sample}}_classified2.ref.bam",
    log:
        f"{OUTPUT_DIR}logs/read_sorting_{{sample}}_PE.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/readsorting_pe_{{sample}}.txt"
    params:
        list=f"{OUTPUT_DIR}read_sorting/{{sample}}/{{sample}}_classified_reads.list",
        phred="--phred64" if config["PHRED_SCORE_64"] else "",
        genome1=config["ASSEMBLY_PARENT_1"],
        genome2=config["ASSEMBLY_PARENT_2"],
        output=lambda w, output: os.path.splitext(os.path.splitext(output.o1)[0])[0][
            :-1
        ],
    shell:
        "{input.eagle_bin} --ngi --paired {params.phred} --ref1={params.genome1} --bam1={input.reads1} --ref2={params.genome2} --bam2={input.reads2} -o {params.output} --bs=3 > {params.list} 2>&1 {log}"
