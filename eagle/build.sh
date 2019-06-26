#!/bin/bash

git clone https://github.com/samtools/htslib.git
sed -i -e "s|PREFIX = /usr/local|PREFIX = $PREFIX|g" Makefile
make
make install
