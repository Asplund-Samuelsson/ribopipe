#!/bin/bash

# Ribosome profiling pipeline
#
# Authors:
# Jan Karlsen
# Johannes Asplund-Samuelsson

# Input variables
EXPERIMENT_NAME=""
THREADS="10"

# Input files

# Intended directory structure
#
# PLACEHOLDER!
#
# main_folder/
# ├── raw/
# │   ├── one.fastq.gz
# │   ├── two.fastq.gz
# │   ├── three.fastq.gz
# │   └── four.fastq.gz
# ├── fastqc/
# │   ├── report.xml
# │   └── standard_out.txt
# └── mapped/
#     └── bowtie1_output.map
