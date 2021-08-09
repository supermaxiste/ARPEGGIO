import os.path

# This file includes all rules related to methylation extraction after alignment
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: Bismark


## Run Bismark methylation extraction on SE bam files for parent species 1


rule methylation_extraction_SE_parent_1:
    input:
        f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_trimmed_bismark_bt2.deduplicated.bam" if config[
            "RUN_TRIMMING"
        ] else (
            f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_bismark_bt2.deduplicated.bam"
        ),
    output:
        sample=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/1.{{sample}}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
            if config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/1.{{sample}}_bismark_bt2.deduplicated.bismark.cov.gz"
            )
        ),
        tmp1=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CpG_context_1.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CpG_context_1.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                )
            )
        ),
        tmp2=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CHG_context_1.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CHG_context_1.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                )
            )
        ),
        tmp3=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CHH_context_1.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CHH_context_1.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                )
            )
        ),
    log:
        f"logs/methXtract_{{sample}}_SE_p1.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/extraction_se_p1_{{sample}}.txt"
    params:
        output=lambda w, output: os.path.split(output.sample)[0],
        genome=config["GENOME_PARENT_1"],
        bismark_cores=BISMARK_CORES,
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "bismark_methylation_extractor -s -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config[
            "UNFINISHED_GENOME"
        ] else (
            "bismark_methylation_extractor -s -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --CX {input}"
        )


## Run Bismark methylation extraction on SE bam files for parent species 2


rule methylation_extraction_SE_parent_2:
    input:
        f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_trimmed_bismark_bt2.deduplicated.bam" if config[
            "RUN_TRIMMING"
        ] else (
            f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_bismark_bt2.deduplicated.bam"
        ),
    output:
        sample=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/2.{{sample}}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
            if config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/2.{{sample}}_bismark_bt2.deduplicated.bismark.cov.gz"
            )
        ),
        tmp1=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CpG_context_2.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CHG_context_2.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                )
            )
        ),
        tmp2=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CHG_context_2.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CpG_context_2.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                )
            )
        ),
        tmp3=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CHH_context_2.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CHH_context_2.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                )
            )
        ),
    log:
        f"logs/methXtract_{{sample}}_SE_p2.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/extraction_se_p2_{{sample}}.txt"
    params:
        output=lambda w, output: os.path.split(output.sample)[0],
        genome=config["GENOME_PARENT_2"],
        bismark_cores=BISMARK_CORES,
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "bismark_methylation_extractor -s -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config[
            "UNFINISHED_GENOME"
        ] else (
            "bismark_methylation_extractor -s -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --CX {input}"
        )


## Run Bismark methylation extraction on SE bam files for allopolyploid species (GENOME_PARENT_1)


rule methylation_extraction_SE_allo_1:
    input:
        f"{OUTPUT_DIR}read_sorting/{{sample}}_se/{{sample}}_classified1.ref.bam" if config[
            "RUN_READ_SORTING"
        ] else (
            f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_trimmed_bismark_bt2.deduplicated.bam"
            if config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_bismark_bt2.deduplicated.bam"
            )
        ),
    output:
        sample=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/{{sample}}_classified1.ref.bismark.cov.gz"
            if config["RUN_READ_SORTING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/1.{{sample}}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
                if config["RUN_TRIMMING"]
                else (
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/1.{{sample}}_bismark_bt2.deduplicated.bismark.cov.gz"
                )
            )
        ),
        tmp1=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CpG_context_{{sample}}_classified1.ref.txt.gz"
            )
            if config["RUN_READ_SORTING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CpG_context_1.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CpG_context_1.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                    )
                )
            )
        ),
        tmp2=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHG_context_{{sample}}_classified1.ref.txt.gz"
            )
            if config["RUN_READ_SORTING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHG_context_1.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHG_context_1.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                    )
                )
            )
        ),
        tmp3=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHH_context_{{sample}}_classified1.ref.txt.gz"
            )
            if config["RUN_READ_SORTING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHH_context_1.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHH_context_1.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                    )
                )
            )
        ),
    log:
        f"logs/methXtract_{{sample}}_SE_allo1.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/extraction_se_allo1_{{sample}}.txt"
    params:
        output=lambda w, output: os.path.split(output.sample)[0],
        genome=config["GENOME_PARENT_1"],
        bismark_cores=BISMARK_CORES,
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "bismark_methylation_extractor -s -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config[
            "UNFINISHED_GENOME"
        ] else (
            "bismark_methylation_extractor -s -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --CX {input}"
        )


## Run Bismark methylation extraction on SE bam files for allopolyploid species (GENOME_PARENT_2)


rule methylation_extraction_SE_allo_2:
    input:
        f"{OUTPUT_DIR}read_sorting/{{sample}}_se/{{sample}}_classified2.ref.bam" if config[
            "RUN_READ_SORTING"
        ] else (
            f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_trimmed_bismark_bt2.deduplicated.bam"
            if config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_bismark_bt2.deduplicated.bam"
            )
        ),
    output:
        sample=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/{{sample}}_classified2.ref.bismark.cov.gz"
            if config["RUN_READ_SORTING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/2.{{sample}}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
                if config["RUN_TRIMMING"]
                else (
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/2.{{sample}}_bismark_bt2.deduplicated.bismark.cov.gz"
                )
            )
        ),
        tmp1=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CpG_context_{{sample}}_classified2.ref.txt.gz"
            )
            if config["RUN_READ_SORTING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CpG_context_2.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CpG_context_2.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                    )
                )
            )
        ),
        tmp2=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHG_context_{{sample}}_classified2.ref.txt.gz"
            )
            if config["RUN_READ_SORTING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHG_context_2.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHG_context_2.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                    )
                )
            )
        ),
        tmp3=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHH_context_{{sample}}_classified2.ref.txt.gz"
            )
            if config["RUN_READ_SORTING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHH_context_2.{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/CHH_context_2.{{sample}}_bismark_bt2.deduplicated.txt.gz"
                    )
                )
            )
        ),
    log:
        f"logs/methXtract_{{sample}}_SE_allo2.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/extraction_se_allo2_{{sample}}.txt"
    params:
        output=lambda w, output: os.path.split(output.sample)[0],
        genome=config["GENOME_PARENT_2"],
        bismark_cores=BISMARK_CORES,
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "bismark_methylation_extractor -s -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config[
            "UNFINISHED_GENOME"
        ] else (
            "bismark_methylation_extractor -s -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --CX {input}"
        )


## Run Bismark methylation extraction on PE bam files for parent species 1


rule methylation_extraction_PE_parent_1:
    input:
        f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bam" if config[
            "RUN_TRIMMING"
        ] else (
            f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bam"
        ),
    output:
        sample1=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
            if config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bismark.cov.gz"
            )
        ),
        sample2=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
            if config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bedGraph.gz"
            )
        ),
        tmp1=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CpG_context_1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CpG_context_1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                )
            )
        ),
        tmp2=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CHG_context_1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CHG_context_1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                )
            )
        ),
        tmp3=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CHH_context_1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/CHH_context_1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                )
            )
        ),
    log:
        f"logs/methXtract_{{sample}}_PE_p1.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/extraction_pe_p1_{{sample}}.txt"
    params:
        output=lambda w, output: os.path.split(output.sample1)[0],
        genome=config["GENOME_PARENT_1"],
        bismark_cores=BISMARK_CORES,
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config[
            "UNFINISHED_GENOME"
        ] else (
            "bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --CX {input}"
        )


## Run Bismark methylation extraction on PE bam files for parent species 2


rule methylation_extraction_PE_parent_2:
    input:
        f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bam" if config[
            "RUN_TRIMMING"
        ] else (
            f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bam"
        ),
    output:
        sample1=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
            if config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bismark.cov.gz"
            )
        ),
        sample2=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
            if config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bedGraph.gz"
            )
        ),
        tmp1=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CpG_context_2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CpG_context_2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                )
            )
        ),
        tmp2=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CHG_context_2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CHG_context_2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                )
            )
        ),
        tmp3=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CHH_context_2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
            )
            if config["RUN_TRIMMING"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/CHH_context_2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                )
            )
        ),
    log:
        f"logs/methXtract_{{sample}}_PE_p2.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/extraction_pe_p2_{{sample}}.txt"
    params:
        output=lambda w, output: os.path.split(output.sample1)[0],
        genome=config["GENOME_PARENT_2"],
        bismark_cores=BISMARK_CORES,
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config[
            "UNFINISHED_GENOME"
        ] else (
            "bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --CX {input}"
        )


## Run Bismark methylation extraction on PE bam files for allopolyploid species (GENOME_PARENT_1)


rule methylation_extraction_PE_allo_1:
    input:
        f"{OUTPUT_DIR}read_sorting/{{sample}}/{{sample}}_classified1.ref.bam" if config[
            "RUN_READ_SORTING"
        ] and config["IS_PAIRED"] else (
            f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bam"
            if config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bam"
            )
        ),
    output:
        sample1=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/{{sample}}_classified1.ref.bismark.cov.gz"
            if config["RUN_READ_SORTING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
                if config["RUN_TRIMMING"]
                else (
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bismark.cov.gz"
                )
            )
        ),
        sample2=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/{{sample}}_classified1.ref.bedGraph.gz"
            if config["RUN_READ_SORTING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
                if config["RUN_TRIMMING"]
                else (
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bedGraph.gz"
                )
            )
        ),
        tmp1=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/CpG_context_{{sample}}_classified1.ref.txt.gz"
            )
            if config["IS_PAIRED"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/CpG_context_1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/CpG_context_1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                    )
                )
            )
        ),
        tmp2=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/CHG_context_{{sample}}_classified1.ref.txt.gz"
            )
            if config["IS_PAIRED"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/CHG_context_1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/CHG_context_1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                    )
                )
            )
        ),
        tmp3=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/CHH_context_{{sample}}_classified1.ref.txt.gz"
            )
            if config["IS_PAIRED"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/CHH_context_1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/CHH_context_1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                    )
                )
            )
        ),
    log:
        f"logs/methXtract_{{sample}}_PE_allo1.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/extraction_pe_allo1_{{sample}}.txt"
    params:
        output=lambda w, output: os.path.split(output.sample1)[0],
        genome=config["GENOME_PARENT_1"],
        bismark_cores=BISMARK_CORES,
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config[
            "UNFINISHED_GENOME"
        ] else (
            "bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --CX {input}"
        )


## Run Bismark methylation extraction on PE bam files for allopolyploid species (GENOME_PARENT_2)


rule methylation_extraction_PE_allo_2:
    input:
        f"{OUTPUT_DIR}read_sorting/{{sample}}/{{sample}}_classified2.ref.bam" if config[
            "RUN_READ_SORTING"
        ] and config["IS_PAIRED"] else (
            f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bam"
            if config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bam"
            )
        ),
    output:
        sample1=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/{{sample}}_classified2.ref.bismark.cov.gz"
            if config["RUN_READ_SORTING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
                if config["RUN_TRIMMING"]
                else (
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bismark.cov.gz"
                )
            )
        ),
        sample2=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/{{sample}}_classified2.ref.bedGraph.gz"
            if config["RUN_READ_SORTING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
                if config["RUN_TRIMMING"]
                else (
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bedGraph.gz"
                )
            )
        ),
        tmp1=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/CpG_context_{{sample}}_classified2.ref.txt.gz"
            )
            if config["IS_PAIRED"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/CpG_context_2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/CpG_context_2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                    )
                )
            )
        ),
        tmp2=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/CHG_context_{{sample}}_classified2.ref.txt.gz"
            )
            if config["IS_PAIRED"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/CHG_context_2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/CHG_context_2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                    )
                )
            )
        ),
        tmp3=(
            temp(
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/CHH_context_{{sample}}_classified2.ref.txt.gz"
            )
            if config["IS_PAIRED"]
            else (
                temp(
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/CHH_context_2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.txt.gz"
                )
                if config["RUN_TRIMMING"]
                else (
                    temp(
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/CHH_context_2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.txt.gz"
                    )
                )
            )
        ),
    log:
        f"logs/methXtract_{{sample}}_PE_allo2.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/extraction_pe_allo2_{{sample}}.txt"
    params:
        output=lambda w, output: os.path.split(output.sample1)[0],
        genome=config["GENOME_PARENT_2"],
        bismark_cores=BISMARK_CORES,
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --scaffolds --CX {input}" if config[
            "UNFINISHED_GENOME"
        ] else (
            "bismark_methylation_extractor -p -o {params.output} --gzip --genome_folder {params.genome} --multicore {params.bismark_cores} --no_overlap --comprehensive --bedGraph --CX {input}"
        )


## Run Bismark coverage2cytosine on extraction output to obtain a single file with information about all cytosines (parent 1)


rule coverage2cytosine_1:
    input:
        f1=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
            if config["IS_PAIRED"] and config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bismark.cov.gz"
                if config["IS_PAIRED"]
                else (
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/1.{{sample}}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
                    if config["RUN_TRIMMING"]
                    else (
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/1.{{sample}}_bismark_bt2.deduplicated.bismark.cov.gz"
                    )
                )
            )
        ),
    output:
        o1=f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/{{sample}}.CX_report.txt",
    log:
        f"logs/c2c_{{sample}}_p1.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/c2c_1_{{sample}}.txt"
    params:
        genome1=config["GENOME_PARENT_1"],
        filename1=lambda w, output: os.path.splitext(os.path.splitext(output.o1)[0])[
            0
        ],
    conda:
        "../envs/environment.yaml"
    shell:
        "coverage2cytosine -CX --genome_folder {params.genome1} -o {params.filename1} {input.f1}"


## Run Bismark coverage2cytosine on extraction output to obtain a single file with information about all cytosines (parent 2)


rule coverage2cytosine_2:
    input:
        f2=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
            if config["IS_PAIRED"] and config["RUN_TRIMMING"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bismark.cov.gz"
                if config["IS_PAIRED"]
                else (
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/2.{{sample}}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
                    if config["RUN_TRIMMING"]
                    else (
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/2.{{sample}}_bismark_bt2.deduplicated.bismark.cov.gz"
                    )
                )
            )
        ),
    output:
        o2=f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/{{sample}}.CX_report.txt",
    log:
        f"logs/c2c_{{sample}}_p2.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/c2c_2_{{sample}}.txt"
    params:
        genome2=config["GENOME_PARENT_2"],
        filename2=lambda w, output: os.path.splitext(os.path.splitext(output.o2)[0])[
            0
        ],
    conda:
        "../envs/environment.yaml"
    shell:
        "coverage2cytosine -CX --genome_folder {params.genome2} -o {params.filename2} {input.f2}"


## Run Bismark coverage2cytosine on extraction output to obtain a single file with information about all cytosines (allopolyploid)


rule coverage2cytosine_allo_1:
    input:
        f1=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/{{sample}}_classified1.ref.bismark.cov.gz"
            if config["RUN_READ_SORTING"] and config["IS_PAIRED"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
                if config["IS_PAIRED"] and config["RUN_TRIMMING"]
                else (
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bismark.cov.gz"
                    if config["IS_PAIRED"]
                    else (
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/{{sample}}_classified1.ref.bismark.cov.gz"
                        if config["RUN_READ_SORTING"]
                        else (
                            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/1.{{sample}}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
                            if config["RUN_TRIMMING"]
                            else (
                                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/1.{{sample}}_bismark_bt2.deduplicated.bismark.cov.gz"
                            )
                        )
                    )
                )
            )
        ),
    output:
        o1=f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/{{sample}}.CX_report.txt",
    log:
        f"logs/c2c_{{sample}}_allo1.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/c2c_allo_{{sample}}.txt"
    params:
        genome1=config["GENOME_PARENT_1"],
        filename1=lambda w, output: os.path.splitext(os.path.splitext(output.o1)[0])[
            0
        ],
    conda:
        "../envs/environment.yaml"
    shell:
        "coverage2cytosine -CX --genome_folder {params.genome1} -o {params.filename1} {input.f1}"


## Run Bismark coverage2cytosine on extraction output to obtain a single file with information about all cytosines (allopolyploid)


rule coverage2cytosine_allo_2:
    input:
        f2=(
            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/{{sample}}_classified2.ref.bismark.cov.gz"
            if config["RUN_READ_SORTING"] and config["IS_PAIRED"]
            else (
                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
                if config["IS_PAIRED"] and config["RUN_TRIMMING"]
                else (
                    f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bismark.cov.gz"
                    if config["IS_PAIRED"]
                    else (
                        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_se/{{sample}}_classified2.ref.bismark.cov.gz"
                        if config["RUN_READ_SORTING"]
                        else (
                            f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/2.{{sample}}_bismark_bt2.deduplicated.bismark.cov.gz"
                            if config["RUN_TRIMMING"]
                            else (
                                f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/2.{{sample}}_bismark_bt2.deduplicated.bismark.cov.gz"
                            )
                        )
                    )
                )
            )
        ),
    output:
        o2=f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/{{sample}}.CX_report.txt",
    log:
        f"logs/c2c_{{sample}}_allo2.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/c2c_allo_{{sample}}.txt"
    params:
        genome2=config["GENOME_PARENT_2"],
        filename2=lambda w, output: os.path.splitext(os.path.splitext(output.o2)[0])[
            0
        ],
    conda:
        "../envs/environment.yaml"
    shell:
        "coverage2cytosine -CX --genome_folder {params.genome2} -o {params.filename2} {input.f2}"
