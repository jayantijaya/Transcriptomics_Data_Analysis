#!/bin/bash

################################################################################
# SAW Pipeline
################################################################################

SAW="/path/to/saw"

SAMPLE="your_sample_id"

CHIP="your_chip_file"

MASK="your_chip_mask_file"

FASTQS="/path/to/fastq_directory"

REFERENCE="/path/to/reference"

OMICS="RNA"

KIT_VERSION=""

SEQUENCING_TYPE=""

OUTPUT_DIR="results/saw/${SAMPLE}"

LOG="logs/saw_${SAMPLE}.log"


# -----------------------------
# Create output directories
# -----------------------------

mkdir -p "${OUTPUT_DIR}"
mkdir -p "$(dirname ${LOG})"


# -----------------------------
# Run SAW count
# -----------------------------

${SAW} count \
    --sn "${CHIP}" \
    --id "${SAMPLE}" \
    --chip-mask "${MASK}" \
    --fastqs "${FASTQS}" \
    --reference "${REFERENCE}" \
    --omics "${OMICS}" \
    --kit-version "${KIT_VERSION}" \
    --sequencing-type "${SEQUENCING_TYPE}" \
    --no-bam \
    --output "${OUTPUT_DIR}" \
    > "${LOG}" 2>&1


echo "SAW count completed for ${SAMPLE}"
