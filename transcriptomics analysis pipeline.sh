#!/bin/bash
set -euo pipefail

# ============================================================
# User-configurable paths
# ============================================================

# Project directory
PROJECT_DIR="$(pwd)"

# Input/output directories
IN_DIR="${PROJECT_DIR}/fastq"
OUT_DIR="${PROJECT_DIR}/results"

# Reference files
GENOME_FA="${PROJECT_DIR}/reference/GRCh38.primary_assembly.genome.fa"
GTF_FILE="${PROJECT_DIR}/reference/gencode.v49.basic.annotation.gtf"

# Software
FASTQC="fastqc"
FASTP="fastp"
HISAT2="hisat2"
SAMTOOLS="samtools"
FEATURECOUNTS="featureCounts"

# HISAT2 genome index
INDEX="${PROJECT_DIR}/reference/hisat2_index"


# ============================================================
# Quality control
# ============================================================

echo "Quality control of FASTQ files"

mkdir -p "${OUT_DIR}/QC"

for file in "${IN_DIR}"/*.fastq
do
    "${FASTQC}" \
        -o "${OUT_DIR}/QC" \
        "$file"
done

echo "Quality control completed"


# ============================================================
# Read trimming
# ============================================================

read -p "Do you want to perform trimming (y/n)? " choice

while true
do
    if [[ "$choice" == "Y" || "$choice" == "y" ]]; then

        mkdir -p "${OUT_DIR}/trimmed"
        mkdir -p "${OUT_DIR}/QC/trimmed"

        for read1 in "${IN_DIR}"/*_R1.fastq
        do
            sample=$(basename "$read1" _R1.fastq)
            read2="${IN_DIR}/${sample}_R2.fastq"

            "${FASTP}" \
                --thread 16 \
                -i "$read1" \
                -I "$read2" \
                -o "${OUT_DIR}/trimmed/${sample}_R1_trimmed.fastq" \
                -O "${OUT_DIR}/trimmed/${sample}_R2_trimmed.fastq" \
                --cut_front \
                --cut_tail \
                --cut_window_size 4 \
                --cut_mean_quality 20 \
                --qualified_quality_phred 20 \
                --length_required 36 \
                --html "${OUT_DIR}/QC/trimmed/${sample}_fastp.html" \
                --json "${OUT_DIR}/QC/trimmed/${sample}_fastp.json"
        done

        echo "Trimming completed"

        for file in "${OUT_DIR}"/trimmed/*_trimmed.fastq
        do
            "${FASTQC}" \
                -o "${OUT_DIR}/QC" \
                "$file"
        done

        echo "Trimmed FastQC completed"

        break

    elif [[ "$choice" == "N" || "$choice" == "n" ]]; then

        echo "Proceeding to alignment..."
        break

    else

        read -p "Please enter y or n: " choice

    fi
done


# ============================================================
# Alignment
# ============================================================

echo "Alignment process"

mkdir -p "${OUT_DIR}/alignment"

if [[ "$choice" == "Y" || "$choice" == "y" ]]; then
    READS="${OUT_DIR}/trimmed"
    SUFFIX="_trimmed.fastq"
else
    READS="${IN_DIR}"
    SUFFIX=".fastq"
fi


for read1 in "${READS}"/*_R1"${SUFFIX}"
do
    sample=$(basename "$read1" "_R1${SUFFIX}")
    read2="${READS}/${sample}_R2${SUFFIX}"

    echo "Aligning sample: ${sample}"

    "${HISAT2}" \
        -p 8 \
        -x "${INDEX}" \
        -1 "$read1" \
        -2 "$read2" \
        -S "${OUT_DIR}/alignment/${sample}.sam" \
        -un-conc "${OUT_DIR}/alignment/${sample}_unaligned.fa" \
        --summary-file "${OUT_DIR}/alignment/${sample}.hisat2.summary"

    "${SAMTOOLS}" view \
        -bS \
        -o "${OUT_DIR}/alignment/${sample}.bam" \
        "${OUT_DIR}/alignment/${sample}.sam"

    echo "${sample} alignment completed"
done


# ============================================================
# BAM sorting
# ============================================================

echo "Sorting BAM files"

mkdir -p "${OUT_DIR}/sort_files"

for file in "${OUT_DIR}"/alignment/*.bam
do
    base=$(basename "$file" .bam)

    "${SAMTOOLS}" sort \
        -@ 8 \
        -o "${OUT_DIR}/sort_files/${base}_sorted.bam" \
        "$file"
done

echo "BAM sorting completed"


# ============================================================
# BAM indexing
# ============================================================

echo "Indexing BAM files"

mkdir -p "${OUT_DIR}/index_bam"

for file in "${OUT_DIR}"/sort_files/*_sorted.bam
do
    "${SAMTOOLS}" index \
        -@ 8 \
        "$file"
done

echo "BAM indexing completed"


# ============================================================
# Feature counting
# ============================================================

echo "Creating count matrix"

mkdir -p "${OUT_DIR}/counts_matrix"

"${FEATURECOUNTS}" \
    -T 8 \
    -p \
    -a "${GTF_FILE}" \
    -o "${OUT_DIR}/counts_matrix/all_samples_counts.csv" \
    "${OUT_DIR}"/sort_files/*_sorted.bam

echo "Count matrix completed"

echo "=========================================="
echo "Transcriptomics pipeline completed"
echo "=========================================="