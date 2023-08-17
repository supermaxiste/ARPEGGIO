import os.path

# This file includes all rules related to checking the quality of the data
# SE = Single-end
# PE = Paired-end
# Tools used in the rules: FastQC,

## Run FastQC on raw untrimmed reads for quality control


rule quality_control:
    input:
        fastq=f"{RAW_DATA_DIR}{{sample}}.{str(config['RAW_DATA_EXTENSION'])}.gz",
    output:
        o1=f"{OUTPUT_DIR}FastQC/{{sample}}_fastqc.zip",
    params:
        FastQC=lambda w, output: os.path.split(output.o1)[0],
    log:
        f"{OUTPUT_DIR}logs/qc_{{sample}}.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/qc_{{sample}}.txt"
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "fastqc -o {params.FastQC} -t {threads} {input.fastq} 2> {log}"


## Run FastQC on SE trimmed reads resulting from trim_galore


rule quality_control_trimmed_SE:
    input:
        fastq=f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_trimmed.fq.gz",
    output:
        o1=f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_trimmed_fastqc.zip",
    params:
        FastQC=lambda w, output: os.path.split(output.o1)[0],
    log:
        f"{OUTPUT_DIR}logs/qc_{{sample}}_trimmedSE.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/qc_trim_se_{{sample}}.txt"
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "fastqc -o {params.FastQC} -t {threads} {input.fastq} 2> {log}"


## Run FastQC on PE trimmed reads resulting from trim_galore


rule quality_control_trimmed_PE:
    input:
        f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_1'])}_val_1.fq.gz",
        f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_2'])}_val_2.fq.gz",
    output:
        o1=f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_1'])}_val_1_fastqc.zip",
        o2=f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_2'])}_val_2_fastqc.zip",
    params:
        FastQC=lambda w, output: os.path.split(output.o1)[0],
    log:
        f"{OUTPUT_DIR}logs/qc_{{sample}}_trimmedPE.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/qc_trim_pe_{{sample}}.txt"
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "fastqc -o {params.FastQC} -t {threads} {input} 2> {log}"


## Run Bismark to prepare synthetic bisulfite converted control genome to check conversion efficiency


rule bismark_prepare_control:
    input:
        control=config["CONTROL_GENOME"],
    output:
        f"{CONTROL_GENOME}Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
    log:
        f"{OUTPUT_DIR}logs/bismark_prepare_control.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/prepare_control_genome.txt"
    conda:
        "../envs/environment.yaml"
    shell:
        "bismark_genome_preparation {input.control} 2> {log}"


## Run Bismark to perform alignment to the control genome to check conversion efficiency if reads are single-end.


rule bismark_alignment_SE_control:
    input:
        control=f"{CONTROL_GENOME}Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        fastq=f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_trimmed.fq.gz"
        if config["RUN_TRIMMING"]
        else f"{RAW_DATA_DIR}{{sample}}.{str(config['RAW_DATA_EXTENSION'])}.gz",
    output:
        sample=f"{OUTPUT_DIR}Conversion_efficiency/{{sample}}/cc.{{sample}}_trimmed_bismark_bt2.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Conversion_efficiency/{{sample}}/cc.{{sample}}_bismark_bt2.bam",
        report=f"{OUTPUT_DIR}Conversion_efficiency/{{sample}}/cc.{{sample}}_trimmed_bismark_bt2_SE_report.txt"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Conversion_efficiency/{{sample}}/cc.{{sample}}_bismark_bt2_SE_report.txt",
    log:
        f"{OUTPUT_DIR}logs/bismark_alignment_SE_{{sample}}_control.log",
    params:
        output=lambda w, output: os.path.split(output.sample)[0],
        control=lambda w, input: os.path.split(
            os.path.split(os.path.split(input.control)[0])[0]
        )[0],
        prefix="cc",
        bismark_cores=BISMARK_CORES,
    benchmark:
        f"{OUTPUT_DIR}benchmark/Conversion_efficiency_{{sample}}.txt"
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "bismark --prefix {params.prefix} --multicore {params.bismark_cores}  -o {params.output} --temp_dir {params.output} --genome {params.control} {input.fastq} 2> {log}"


rule bismark_alignment_PE_control:
    input:
        control=f"{CONTROL_GENOME}Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        fastq1=f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_1'])}_val_1.fq.gz"
        if config["RUN_TRIMMING"]
        else f"{RAW_DATA_DIR}{{sample}}_{str(config['PAIR_1'])}.{str(config['RAW_DATA_EXTENSION'])}.gz",
        fastq2=f"{OUTPUT_DIR}FASTQtrimmed/{{sample}}_{str(config['PAIR_2'])}_val_2.fq.gz"
        if config["RUN_TRIMMING"]
        else f"{RAW_DATA_DIR}{{sample}}_{str(config['PAIR_2'])}.{str(config['RAW_DATA_EXTENSION'])}.gz",
    output:
        sample=f"{OUTPUT_DIR}Conversion_efficiency/{{sample}}/cc.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Conversion_efficiency/{{sample}}/cc.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.bam",
        report=f"{OUTPUT_DIR}Conversion_efficiency/{{sample}}/cc.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_PE_report.txt"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Conversion_efficiency/{{sample}}/cc.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_PE_report.txt",
    log:
        f"{OUTPUT_DIR}logs/bismark_alignment_PE_{{sample}}_control.log",
    params:
        output=lambda w, output: os.path.split(output.sample)[0],
        control=lambda w, input: os.path.split(
            os.path.split(os.path.split(input.control)[0])[0]
        )[0],
        prefix="cc",
        bismark_cores=BISMARK_CORES,
    benchmark:
        f"{OUTPUT_DIR}benchmark/Conversion_efficiency_{{sample}}.txt"
    conda:
        "../envs/environment.yaml"
    threads: CORES
    shell:
        "bismark --prefix {params.prefix} --multicore {params.bismark_cores} --genome {params.control} -1 {input.fastq1} -2 {input.fastq2} -o {params.output} --temp_dir {params.output} 2> {log}"


# sort bam files for multiqc (parent1)


rule bam_sorting_p1:
    input:
        p1=f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bam"
        if config["IS_PAIRED"] and config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bam"
        if config["IS_PAIRED"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_trimmed_bismark_bt2.deduplicated.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_bismark_bt2.deduplicated.bam",
    output:
        o1=f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated_sorted.bam"
        if config["IS_PAIRED"] and config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated_sorted.bam"
        if config["IS_PAIRED"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_trimmed_bismark_bt2.deduplicated_sorted.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_bismark_bt2.deduplicated_sorted.bam",
    log:
        f"{OUTPUT_DIR}logs/bam_sorting_{{sample}}_p1.log",
    conda:
        "../envs/environment.yaml"
    shell:
        "samtools sort {input.p1} > {output.o1}"


# sort bam files for multiqc (parent2)


rule bam_sorting_p2:
    input:
        p2=f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bam"
        if config["IS_PAIRED"] and config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bam"
        if config["IS_PAIRED"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_trimmed_bismark_bt2.deduplicated.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_bismark_bt2.deduplicated.bam",
    output:
        o2=f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated_sorted.bam"
        if config["IS_PAIRED"] and config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated_sorted.bam"
        if config["IS_PAIRED"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_trimmed_bismark_bt2.deduplicated_sorted.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_bismark_bt2.deduplicated_sorted.bam",
    log:
        f"{OUTPUT_DIR}logs/bam_sorting_{{sample}}_p2.log",
    conda:
        "../envs/environment.yaml"
    shell:
        "samtools sort {input.p2} > {output.o2} 2> {log}"


# sort bam files for multiqc (allopolyploid)


rule bam_sorting_allo:
    input:
        f"{OUTPUT_DIR}read_sorting/{{sample}}/{{sample}}_classified{{one_or_two}}.ref.bam"
        if config["RUN_READ_SORTING"] and config["IS_PAIRED"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bam"
        if config["IS_PAIRED"] and config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bam"
        if config["IS_PAIRED"]
        else f"{OUTPUT_DIR}read_sorting/{{sample}}_se/{{sample}}_classified{{one_or_two}}.ref.bam"
        if config["RUN_READ_SORTING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_trimmed_bismark_bt2.deduplicated.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_bismark_bt2.deduplicated.bam",
    output:
        f"{OUTPUT_DIR}read_sorting/{{sample}}/{{sample}}_classified{{one_or_two}}_sorted.ref.bam"
        if config["RUN_READ_SORTING"] and config["IS_PAIRED"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated_sorted_allo.bam"
        if config["IS_PAIRED"] and config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated_sorted_allo.bam"
        if config["IS_PAIRED"]
        else f"{OUTPUT_DIR}read_sorting/{{sample}}_se/{{sample}}_classified{{one_or_two}}_sorted.ref.bam"
        if config["RUN_READ_SORTING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_trimmed_bismark_bt2.deduplicated_sorted_allo.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_bismark_bt2.deduplicated_sorted_allo.bam",
    log:
        f"{OUTPUT_DIR}logs/bam_sorting_{{sample}}_{{one_or_two}}allo.log",
    conda:
        "../envs/environment.yaml"
    shell:
        "samtools sort {input} > {output}"


## Qualimap rule to get statistics about bam files for parent1


rule qualimap_p1:
    input:
        f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated_sorted.bam"
        if config["IS_PAIRED"] and config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated_sorted.bam"
        if config["IS_PAIRED"]
        else (
            f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_trimmed_bismark_bt2.deduplicated_sorted.bam"
            if config["RUN_TRIMMING"]
            else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_1/1.{{sample}}_bismark_bt2.deduplicated_sorted.bam"
        ),
    output:
        o1=f"{OUTPUT_DIR}qualimap/{{sample}}_p1/qualimapReport.html",
    log:
        f"{OUTPUT_DIR}logs/qualimap_{{sample}}_p1.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/qualimap_p1_{{sample}}.txt"
    params:
        output=lambda w, output: os.path.split(output.o1)[0],
    conda:
        "../envs/environment2.yaml"
    threads: CORES
    shell:
        "qualimap bamqc -bam {input} -outdir {params.output} -nt {threads} --java-mem-size=4G 2> {log}"


## Qualimap rule to get statistics about bam files for parent2


rule qualimap_p2:
    input:
        f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated_sorted.bam"
        if config["IS_PAIRED"] and config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated_sorted.bam"
        if config["IS_PAIRED"]
        else (
            f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_trimmed_bismark_bt2.deduplicated_sorted.bam"
            if config["RUN_TRIMMING"]
            else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_2/2.{{sample}}_bismark_bt2.deduplicated_sorted.bam"
        ),
    output:
        o1=f"{OUTPUT_DIR}qualimap/{{sample}}_p2/qualimapReport.html",
    log:
        f"{OUTPUT_DIR}logs/qualimap_{{sample}}_p2.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/qualimap_p2_{{sample}}.txt"
    params:
        output=lambda w, output: os.path.split(output.o1)[0],
    conda:
        "../envs/environment2.yaml"
    threads: CORES
    shell:
        "qualimap bamqc -bam {input} -outdir {params.output} -nt {threads} --java-mem-size=4G 2> {log}"


## Qualimap rule to get statistics about bam files for allopolyploid SE reads


rule qualimap_allo_se:
    input:
        genome=f"{OUTPUT_DIR}read_sorting/{{sample}}_se/{{sample}}_classified{{one_or_two}}_sorted.ref.bam"
        if config["RUN_READ_SORTING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_trimmed_bismark_bt2.deduplicated_sorted_allo.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_bismark_bt2.deduplicated_sorted_allo.bam",
    output:
        o1=f"{OUTPUT_DIR}qualimap/{{sample}}_allo_se_{{one_or_two}}/qualimapReport.html",
    log:
        f"{OUTPUT_DIR}logs/qualimap_{{sample}}_{{one_or_two}}_alloSE.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/qualimap_allo_se_{{sample}}_{{one_or_two}}.txt"
    params:
        output=lambda w, output: os.path.split(output.o1)[0],
    conda:
        "../envs/environment2.yaml"
    threads: CORES
    shell:
        "qualimap bamqc -bam {input.genome} -outdir {params.output} -nt {threads} --java-mem-size=4G 2> {log}"


## Qualimap rule to get statistics about bam files for allopolyploid PE reads


rule qualimap_allo_pe:
    input:
        genome=f"{OUTPUT_DIR}read_sorting/{{sample}}/{{sample}}_classified{{one_or_two}}_sorted.ref.bam"
        if config["RUN_READ_SORTING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_{str(config['PAIR_1'])}_val_1_bismark_bt2_pe.deduplicated.bam"
        if config["RUN_TRIMMING"]
        else f"{OUTPUT_DIR}Bismark/deduplication/{{sample}}_{{one_or_two}}/{{one_or_two}}.{{sample}}_{str(config['PAIR_1'])}_bismark_bt2_pe.deduplicated.bam",
    output:
        o1=f"{OUTPUT_DIR}qualimap/{{sample}}_allo_pe_{{one_or_two}}/qualimapReport.html",
    log:
        f"{OUTPUT_DIR}logs/qualimap_{{sample}}_{{one_or_two}}_alloPE.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/qualimap_allo_pe_{{sample}}_{{one_or_two}}.txt"
    params:
        output=lambda w, output: os.path.split(output.o1)[0],
    conda:
        "../envs/environment2.yaml"
    threads: CORES
    shell:
        "qualimap bamqc -bam {input.genome} -outdir {params.output} -nt {threads} --java-mem-size=4G 2> {log}"


## Run MultiQC to combine all the outputs from QC, trimming and alignment in a single nice report


rule multiqc:
    input:
        multiqc_input,
    output:
        o1=f"{OUTPUT_DIR}MultiQC/multiqc_report.html",
    log:
        f"{OUTPUT_DIR}logs/multiqc.log",
    params:
        inputdir=multiqc_params,
        multiqcdir=lambda w, output: os.path.split(output.o1)[0],
    benchmark:
        f"{OUTPUT_DIR}benchmark/multiqc.txt"
    conda:
        "../envs/environment.yaml"
    shell:
        "multiqc {params.inputdir} -f -o {params.multiqcdir} 2> {log}"
