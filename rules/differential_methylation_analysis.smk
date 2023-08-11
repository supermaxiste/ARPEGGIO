import os.path

# This file includes all rules related to differential methylation analysis
# DMR = Differentially Methylated Region
# Tools used in the rules: R scripts, dmrseq

# R scripts used: CoverageFileGeneratorComplete.R, dmrseq.R

## Run R script to remove all rows with no cytosines and generate three files for each cytosine context for parent 1


rule context_separation_parent_1:
    input:
        p1=f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p1/{{sample}}.CX_report.txt",
    output:
        o1=f"{OUTPUT_DIR}DMR_analysis/context_separation/parent1/{{sample}}_CG.cov",
        o2=f"{OUTPUT_DIR}DMR_analysis/context_separation/parent1/{{sample}}_CHG.cov",
        o3=f"{OUTPUT_DIR}DMR_analysis/context_separation/parent1/{{sample}}_CHH.cov",
    log:
        f"{OUTPUT_DIR}/logs/context_separation_{{sample}}_p1.log",
    params:
        sample_name=f"{{sample}}",
        output=lambda w, output: os.path.split(output.o1)[0],
    conda:
        "../envs/environment_downstream.yaml"
    shell:
        "Rscript scripts/CoverageFileGeneratorComplete.R {input.p1} {params.output} {params.sample_name} 2>&1 {log}"


## Run R script to remove all rows with no cytosines and generate three files for each cytosine context for parent 2


rule context_separation_parent_2:
    input:
        p1=f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_p2/{{sample}}.CX_report.txt",
    output:
        o1=f"{OUTPUT_DIR}DMR_analysis/context_separation/parent2/{{sample}}_CG.cov",
        o2=f"{OUTPUT_DIR}DMR_analysis/context_separation/parent2/{{sample}}_CHG.cov",
        o3=f"{OUTPUT_DIR}DMR_analysis/context_separation/parent2/{{sample}}_CHH.cov",
    log:
        f"{OUTPUT_DIR}/logs/context_separation_{{sample}}_p2.log",
    params:
        sample_name=f"{{sample}}",
        output=lambda w, output: os.path.split(output.o1)[0],
    conda:
        "../envs/environment_downstream.yaml"
    shell:
        "Rscript scripts/CoverageFileGeneratorComplete.R {input.p1} {params.output} {params.sample_name} 2>&1 {log}"


## Run R script to remove all rows with no cytosines and generate three files for each cytosine context for parent 1


rule combine_context_allo:
    input:
        p1=f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_1/{{sample}}.CX_report.txt",
        p2=f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_2/{{sample}}.CX_report.txt",
    output:
        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_12/{{sample}}.CX_report.txt",
    log:
        f"{OUTPUT_DIR}/logs/combine_context_{{sample}}_12.log",
    shell:
        "cat {input.p1} {input.p2} > {output}"


rule context_separation_allo:
    input:
        f"{OUTPUT_DIR}Bismark/extraction/{{sample}}_12/{{sample}}.CX_report.txt",
    output:
        o1=f"{OUTPUT_DIR}DMR_analysis/context_separation/allopolyploid/{{sample}}_CG.cov",
        o2=f"{OUTPUT_DIR}DMR_analysis/context_separation/allopolyploid/{{sample}}_CHG.cov",
        o3=f"{OUTPUT_DIR}DMR_analysis/context_separation/allopolyploid/{{sample}}_CHH.cov",
    log:
        f"{OUTPUT_DIR}/logs/context_separation_{{sample}}_allo.log",
    params:
        sample_name=f"{{sample}}",
        output=lambda w, output: os.path.split(output.o1)[0],
    conda:
        "../envs/environment_downstream.yaml"
    shell:
        "Rscript scripts/CoverageFileGeneratorComplete.R {input} {params.output} {params.sample_name} 2>&1 {log}"


# Run dmrseq for CG context for parent1 subgenome


rule dmrseq_CG_1:
    input:
        p1=dmrseq_CG_input_p1,
        allo=dmrseq_CG_input_allo,
    output:
        comparison1=f"{OUTPUT_DIR}DMR_analysis/dmrseq/CG_context/parent1_v_allo.txt",
    log:
        f"{OUTPUT_DIR}/logs/dmrseq_CG_1.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/dmrseq_CG.txt"
    params:
        n_samples_p1=n_samples_p1,
        n_samples_allo=n_samples_allo,
        script="scripts/dmrseq.R",
    threads: CORES
    conda:
        "../envs/environment_R.yaml"
    shell:
        "Rscript scripts/dmrseq.R {params.n_samples_p1} {params.n_samples_allo} {output.comparison1} {threads} {input.p1} {input.allo} 2>&1 {log}"


# Run dmrseq for CG context for parent2 subgenome


rule dmrseq_CG_2:
    input:
        p2=dmrseq_CG_input_p2,
        allo=dmrseq_CG_input_allo,
    output:
        comparison2=f"{OUTPUT_DIR}DMR_analysis/dmrseq/CG_context/parent2_v_allo.txt",
    log:
        f"{OUTPUT_DIR}/logs/dmrseq_CG_2.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/dmrseq_CG.txt"
    params:
        n_samples_p2=n_samples_p2,
        n_samples_allo=n_samples_allo,
        script="scripts/dmrseq.R",
    threads: CORES
    conda:
        "../envs/environment_R.yaml"
    shell:
        "Rscript scripts/dmrseq.R {params.n_samples_p2} {params.n_samples_allo} {output.comparison2} {threads} {input.p2} {input.allo} 2>&1 {log}"


# Run dmrseq for CHG context for parent1 subgenome


rule dmrseq_CHG_1:
    input:
        p1=dmrseq_CHG_input_p1,
        allo=dmrseq_CHG_input_allo,
    output:
        comparison1=f"{OUTPUT_DIR}DMR_analysis/dmrseq/CHG_context/parent1_v_allo.txt",
    log:
        f"{OUTPUT_DIR}/logs/dmrseq_CHG_1.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/dmrseq_CHG.txt"
    params:
        n_samples_p1=n_samples_p1,
        n_samples_allo=n_samples_allo,
        script="scripts/dmrseq.R",
    threads: CORES
    conda:
        "../envs/environment_R.yaml"
    shell:
        "Rscript scripts/dmrseq.R {params.n_samples_p1} {params.n_samples_allo} {output.comparison1} {threads} {input.p1} {input.allo} 2>&1 {log}"


# Run dmrseq for CHG context for parent2 subgenome


rule dmrseq_CHG_2:
    input:
        p2=dmrseq_CHG_input_p2,
        allo=dmrseq_CHG_input_allo,
    output:
        comparison2=f"{OUTPUT_DIR}DMR_analysis/dmrseq/CHG_context/parent2_v_allo.txt",
    log:
        f"{OUTPUT_DIR}/logs/dmrseq_CHG_2.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/dmrseq_CHG.txt"
    params:
        n_samples_p2=n_samples_p2,
        n_samples_allo=n_samples_allo,
        script="scripts/dmrseq.R",
    threads: CORES
    conda:
        "../envs/environment_R.yaml"
    shell:
        "Rscript scripts/dmrseq.R {params.n_samples_p2} {params.n_samples_allo} {output.comparison2} {threads} {input.p2} {input.allo} 2>&1 {log}"


# Run dmrseq for CHH context for parent1 subgenome


rule dmrseq_CHH_1:
    input:
        p1=dmrseq_CHH_input_p1,
        allo=dmrseq_CHH_input_allo,
    output:
        comparison1=f"{OUTPUT_DIR}DMR_analysis/dmrseq/CHH_context/parent1_v_allo.txt",
    log:
        f"{OUTPUT_DIR}/logs/dmrseq_CHH_1.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/dmrseq_CHH.txt"
    params:
        n_samples_p1=n_samples_p1,
        n_samples_allo=n_samples_allo,
        script="scripts/dmrseq.R",
    threads: CORES
    conda:
        "../envs/environment_R.yaml"
    shell:
        "Rscript scripts/dmrseq.R {params.n_samples_p1} {params.n_samples_allo} {output.comparison1} {threads} {input.p1} {input.allo} 2>&1 {log}"


# Run dmrseq for CHH context for parent2 subgenome


rule dmrseq_CHH_2:
    input:
        p2=dmrseq_CHH_input_p2,
        allo=dmrseq_CHH_input_allo,
    output:
        comparison2=f"{OUTPUT_DIR}DMR_analysis/dmrseq/CHH_context/parent2_v_allo.txt",
    log:
        f"{OUTPUT_DIR}/logs/dmrseq_CHH_2.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/dmrseq_CHH.txt"
    params:
        n_samples_p2=n_samples_p2,
        n_samples_allo=n_samples_allo,
        script="scripts/dmrseq.R",
    threads: CORES
    conda:
        "../envs/environment_R.yaml"
    shell:
        "Rscript scripts/dmrseq.R {params.n_samples_p2} {params.n_samples_allo} {output.comparison2} {threads} {input.p2} {input.allo} 2>&1 {log}"


## Rules for dmrseq comparing polyploids to polyploids or diploids to diploids based on conditions given in the metadata file


rule dmrseq_CG_special:
    input:
        condA=dmrseq_CG_special_input_A,
        condB=dmrseq_CG_special_input_B,
    output:
        comparison=f"{OUTPUT_DIR}DMR_analysis/dmrseq/CG_context/A_v_B_diploid.txt"
        if config["DIPLOID_ONLY"]
        else f"{OUTPUT_DIR}DMR_analysis/dmrseq/CG_context/A_v_B_polyploid.txt",
    log:
        f"{OUTPUT_DIR}/logs/dmrseq_CG_special.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/dmrseq_CG_special.txt"
    params:
        samples_A=n_samples_A,
        samples_B=n_samples_B,
    threads: CORES
    conda:
        "../envs/environment_R.yaml"
    shell:
        "Rscript scripts/dmrseq.R {params.samples_B} {params.samples_A} {output.comparison} {threads} {input.condB} {input.condA} 2>&1 {log}"


rule dmrseq_CHG_special:
    input:
        condA=dmrseq_CHG_special_input_A,
        condB=dmrseq_CHG_special_input_B,
    output:
        comparison=f"{OUTPUT_DIR}DMR_analysis/dmrseq/CHG_context/A_v_B_diploid.txt"
        if config["DIPLOID_ONLY"]
        else f"{OUTPUT_DIR}DMR_analysis/dmrseq/CHG_context/A_v_B_polyploid.txt",
    log:
        f"{OUTPUT_DIR}/logs/dmrseq_CHG_special.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/dmrseq_CHG_special.txt"
    params:
        samples_A=n_samples_A,
        samples_B=n_samples_B,
    threads: CORES
    conda:
        "../envs/environment_R.yaml"
    shell:
        "Rscript scripts/dmrseq.R {params.samples_B} {params.samples_A} {output.comparison} {threads} {input.condB} {input.condA} 2>&1 {log}"


rule dmrseq_CHH_special:
    input:
        condA=dmrseq_CHH_special_input_A,
        condB=dmrseq_CHH_special_input_B,
    output:
        comparison=f"{OUTPUT_DIR}DMR_analysis/dmrseq/CHH_context/A_v_B_diploid.txt"
        if config["DIPLOID_ONLY"]
        else f"{OUTPUT_DIR}DMR_analysis/dmrseq/CHH_context/A_v_B_polyploid.txt",
    log:
        f"{OUTPUT_DIR}/logs/dmrseq_CHH_special.log",
    benchmark:
        f"{OUTPUT_DIR}benchmark/dmrseq_CHH_special.txt"
    params:
        samples_A=n_samples_A,
        samples_B=n_samples_B,
    threads: CORES
    conda:
        "../envs/environment_R.yaml"
    shell:
        "Rscript scripts/dmrseq.R {params.samples_B} {params.samples_A} {output.comparison} {threads} {input.condB} {input.condA} 2>&1 {log}"
