################## BASE IMAGE ######################
FROM ubuntu:20.04

################## METADATA ######################
LABEL base.image="ubuntu:20.04"
LABEL version="1"
LABEL software="ARPEGGIO"
LABEL software.version="3.0.0"
LABEL description="A SnakeMake workflow to analyse whole genome bisulfite sequencing data from allopolyploids."
LABEL website="https://github.com/supermaxiste/ARPEGGIO"
LABEL license="https://github.com/supermaxiste/ARPEGGIO/blob/master/LICENSE"
LABEL maintainer="Stefan Milosavljevic"
LABEL maintainer.email="stefan.milos.srb.ch@gmail.com"

################## INSTALLATION ######################

FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="ae3e1ce0bbee983c5c5cd02ab1241276924b60c4061813705b75ce8a3ca58ace"

# INSTALL DEPENDENCIES
# Step 1: Retrieve conda environments

# Conda environment:
#   source: envs/build_eagle.yaml
#   prefix: /conda-envs/764031fc805b9f21beb1f94250215b4d
#   name: build_eagle
#   channels:
#     - conda-forge
#   dependencies:
#     - wget=1.20.1
#     - tar=1.32
#     - bzip2=1.0.8
#     - make=4.2.1
#     - autoconf=2.69
#     - zlib=1.2.11
#     - libcurl=7.65.3
#     - openssl=1.1.1d
#     - backports.lzma=0.0.14
#     - gcc_linux-64=7.3.0
RUN mkdir -p /conda-envs/764031fc805b9f21beb1f94250215b4d
COPY envs/build_eagle.yaml /conda-envs/764031fc805b9f21beb1f94250215b4d/environment.yaml

# Conda environment:
#   source: envs/environment.yaml
#   prefix: /conda-envs/ffe4901fcdc7a9ea44cb41f11f379297
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - trim-galore=0.6.6
#     - fastqc=0.11.9
#     - bismark=0.23.0
#     - samtools=1.11
#     - multiqc=1.9
#     - pigz=2.3.4
RUN mkdir -p /conda-envs/ffe4901fcdc7a9ea44cb41f11f379297
COPY envs/environment.yaml /conda-envs/ffe4901fcdc7a9ea44cb41f11f379297/environment.yaml

# Conda environment:
#   source: envs/environment2.yaml
#   prefix: /conda-envs/89a15cf85889e37b88cb866b9ea166ab
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - qualimap=2.2.2d
RUN mkdir -p /conda-envs/89a15cf85889e37b88cb866b9ea166ab
COPY envs/environment2.yaml /conda-envs/89a15cf85889e37b88cb866b9ea166ab/environment.yaml

# Conda environment:
#   source: envs/environment_R.yaml
#   prefix: /conda-envs/521188b4452f95cba75575c660c372a2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bioconductor-dmrseq=1.10.0
#     - r-base=4.0.3
RUN mkdir -p /conda-envs/521188b4452f95cba75575c660c372a2
COPY envs/environment_R.yaml /conda-envs/521188b4452f95cba75575c660c372a2/environment.yaml

# Conda environment:
#   source: envs/environment_downstream.yaml
#   prefix: /conda-envs/11e671c19c824bbf0d6b8afbe55df0c7
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-base=4.0.3
#     - r-data.table=1.13.6
#     - bedtools=2.29.2
#     - r-stringr=1.4.0
#     - r-stringi=1.5.3
#     - libgfortran=3.0.0
RUN mkdir -p /conda-envs/11e671c19c824bbf0d6b8afbe55df0c7
COPY envs/environment_downstream.yaml /conda-envs/11e671c19c824bbf0d6b8afbe55df0c7/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/764031fc805b9f21beb1f94250215b4d --file /conda-envs/764031fc805b9f21beb1f94250215b4d/environment.yaml && \
    mamba env create --prefix /conda-envs/ffe4901fcdc7a9ea44cb41f11f379297 --file /conda-envs/ffe4901fcdc7a9ea44cb41f11f379297/environment.yaml && \
    mamba env create --prefix /conda-envs/89a15cf85889e37b88cb866b9ea166ab --file /conda-envs/89a15cf85889e37b88cb866b9ea166ab/environment.yaml && \
    mamba env create --prefix /conda-envs/521188b4452f95cba75575c660c372a2 --file /conda-envs/521188b4452f95cba75575c660c372a2/environment.yaml && \
    mamba env create --prefix /conda-envs/11e671c19c824bbf0d6b8afbe55df0c7 --file /conda-envs/11e671c19c824bbf0d6b8afbe55df0c7/environment.yaml && \
    mamba clean --all -y
