#!/bin/bash
################################################################################
# FastQC Analysis
################################################################################
# create output directory
mkdir -p fastqc_output

# run FastQC on all FASTQ files
fastqc *.fastq.gz \
-o fastqc_output

echo "FastQC completed"
