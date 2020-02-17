import os.path

BUILD_PREFIX = ".eagle-build"
INSTALL_PREFIX = ".eagle"

ENV_PATH = os.path.join("..", "envs", "build_eagle.yaml")

HTSLIB_VERSION = "1.9"
HTSLIB_DIR_NAME = f"htslib-{HTSLIB_VERSION}"
HTSLIB_TAR_NAME = f"{HTSLIB_DIR_NAME}.tar.bz2"
HTSLIB_DIR_PATH = os.path.join(BUILD_PREFIX, HTSLIB_DIR_NAME)
HTSLIB_TAR_PATH = os.path.join(BUILD_PREFIX, HTSLIB_TAR_NAME)
HTSLIB_MK_PATH = os.path.join(HTSLIB_DIR_PATH, "Makefile")

EAGLE_VERSION = "1.1.1"
EAGLE_DIR_NAME = f"eagle-{EAGLE_VERSION}"
EAGLE_TAR_NAME = f"{EAGLE_DIR_NAME}.tar.gz"
EAGLE_DIR_PATH = os.path.join(BUILD_PREFIX, EAGLE_DIR_NAME)
EAGLE_TAR_PATH = os.path.join(BUILD_PREFIX, EAGLE_TAR_NAME)
EAGLE_MK_PATH = os.path.join(EAGLE_DIR_PATH, "Makefile")
EAGLE_MK_HTSDIR = os.path.join("..", HTSLIB_DIR_NAME)
EAGLE_MK_PREFIX = os.path.join("..", "..", INSTALL_PREFIX)
EAGLE_BIN_PATH = os.path.join(EAGLE_DIR_PATH, "eagle-rc")
EAGLE_BIN_INSTALL_PATH = os.path.join(INSTALL_PREFIX, "bin", "eagle-rc")


# For whatever reason conda's gcc ignores LIBRARY_PATH
# so we must try and override make's LDFLAGS
# eagle decided to use the non standard name LFLAGS
INCLUDES = os.path.join("$CONDA_PREFIX", "include")
LIBS = os.path.join("$CONDA_PREFIX", "lib")
GCC_PREFIX = "x86_64-conda_cos6-linux-gnu-"
CC = f"{GCC_PREFIX}cc"
AR = f"{GCC_PREFIX}ar"
RANLIB= f"{GCC_PREFIX}ranlib"
HTSLIB_MK_FLAGS = f"LDFLAGS=-L{LIBS}"
EAGLE_MK_FLAGS = f'CC={CC} AR={AR} RANLIB={RANLIB} HTSDIR="{EAGLE_MK_HTSDIR}" PREFIX="{EAGLE_MK_PREFIX}" LFLAGS="-L{EAGLE_MK_HTSDIR} -L{LIBS}" LDLIBS="-lm -lz -llzma -lbz2 -lpthread -lcurl -lssl -lcrypto" MAKE="make {HTSLIB_MK_FLAGS}"'
EAGLE_MK_ENV = f'C_INCLUDE_PATH="{INCLUDES}" LIBRARY_PATH="{LIBS}"'

rule download_eagle:
    output:
        eagle_tar = EAGLE_TAR_PATH
    conda:
        ENV_PATH
    shell:
        "wget https://github.com/tony-kuo/eagle/archive/v{EAGLE_VERSION}.tar.gz -O {output.eagle_tar}"

rule extract_eagle:
    input:
        eagle_tar = EAGLE_TAR_PATH
    output:
       eagle_mk = EAGLE_MK_PATH
    conda:
        ENV_PATH
    shell:
        "tar -C {BUILD_PREFIX} -vxf {input.eagle_tar}"

rule download_htslib:
    output:
       htslib_tar = HTSLIB_TAR_PATH
    conda:
        ENV_PATH
    shell:
        "wget https://github.com/samtools/htslib/releases/download/{HTSLIB_VERSION}/{HTSLIB_TAR_NAME} -O {output.htslib_tar}"

rule extract_htslib:
    input:
        htslib_tar = HTSLIB_TAR_PATH
    output:
        htslib_mk = HTSLIB_MK_PATH
    conda:
        ENV_PATH
    shell:
        "tar -C {BUILD_PREFIX} -vxf {input.htslib_tar}"

rule build_eagle:
    input:
        eagle_mk = EAGLE_MK_PATH,
        htslib_mk = HTSLIB_MK_PATH
    output:
        eagle_bin = EAGLE_BIN_PATH
    conda:
        ENV_PATH
    threads:
        CORES
    shell:
        "{EAGLE_MK_ENV} make -j {threads} -C {EAGLE_DIR_PATH} {EAGLE_MK_FLAGS}"

rule install_eagle:
    input:
        eagle_bin = EAGLE_BIN_PATH
    output:
        eagle_bin = EAGLE_BIN_INSTALL_PATH
    conda:
        ENV_PATH
    shell:
        "{EAGLE_MK_ENV} make -C {EAGLE_DIR_PATH} install {EAGLE_MK_FLAGS}"
