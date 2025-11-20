#!/bin/bash

#
# This script automates the assembly of short reads for all samples in the current directory.
# It expects paired-end reads in a directory named 'reads' ending in _R1.fastq.gz and _R2.fastq.gz.
#
# Workflow:
# 1. FastQC: Raw read quality control.
# 2. fastp: Adapter and quality trimming.
# 3. SPAdes: Genome assembly.
# 4. seqtk: Filter contigs by minimum length.
# 5. CheckM2: Quality assessment of all assemblies.
# 6. assembly-stats: Generate assembly statistics for all assemblies.
#

# --- Stop script on any error ---
# The following line is changed from 'set -e' to 'set +e' to allow the script to continue after an error.
set +e
set -o pipefail

# --- Activate illumina conda environment ---
eval "$(conda shell.bash hook)"
conda activate ShortReadAssembler

# ==============================================================================
#  CONFIGURATION (EDIT THESE VARIABLES)
# ==============================================================================
THREADS=120                                         # Number of CPU threads to use for each step
SPADES_THREADS=6                                   # Number of CPU threads for SPADES                                         
MEMORY_GB=200                                       # Max memory for SPAdes in Gigabytes
LENGTH_FILTER=500                                   # Minimum contig length to keep
checkm2_database=/home/mcknid01/softwareDependencies/CheckM2_database/uniref100.KO.1.dmnd  # CheckM2 database path
# ==============================================================================

echo "Starting the assembly pipeline..."

# --- Create output directories ---
mkdir -p 01_fastqc 02_trimmed_reads 03_assemblies 04_final_assemblies
touch failed.txt

# --- Main processing loop for each sample ---
for R1 in 00_reads/*_R1.fastq.gz; do
    # --- Define sample names and files ---
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    ID=$(basename "$R1" _R1.fastq.gz)

    echo "================================================="
    echo "Processing Sample: $ID"
    echo "================================================="

    # create per-sample subdirs for each step
    mkdir -p "01_fastqc/$ID" "02_trimmed_reads/$ID" "03_assemblies/$ID"

    # --- Step 1: FastQC ---
    echo "[1/4] Running FastQC on raw reads..."
    if [ -f "01_fastqc/$ID/${ID}_R1_fastqc.html" ]; then
        echo "     -> FastQC report for $ID already exists. Skipping."
    else
        fastqc -t $THREADS "$R1" "$R2" -o "01_fastqc/$ID"
    fi

    # --- Step 2: fastp Trimming ---
    echo "[2/4] Trimming reads with fastp..."
    TRIMMED_R1="02_trimmed_reads/$ID/${ID}_R1.trimmed.fastq.gz"
    TRIMMED_R2="02_trimmed_reads/$ID/${ID}_R2.trimmed.fastq.gz"
    if [ -f "$TRIMMED_R1" ]; then
        echo "     -> Trimmed reads for $ID already exist. Skipping."
    else
        fastp \
            -i "$R1" -I "$R2" \
            -o "$TRIMMED_R1" -O "$TRIMMED_R2" \
            -h "02_trimmed_reads/$ID/${ID}.fastp.html" \
            -j "02_trimmed_reads/$ID/${ID}.fastp.json" \
            -w $THREADS
    fi

    # --- Step 3: SPAdes Assembly ---
    echo "[3/4] Assembling with SPAdes..."
    ASSEMBLY_DIR="03_assemblies/$ID/${ID}_spades"
    ASSEMBLY_FASTA="03_assemblies/$ID/${ID}.fasta"
    if [ -f "$ASSEMBLY_FASTA" ]; then
        echo "     -> Assembly for $ID already exists. Skipping."
    else
        spades.py \
            --isolate \
            -1 "$TRIMMED_R1" -2 "$TRIMMED_R2" \
            -o "$ASSEMBLY_DIR" \
            -t $SPADES_THREADS -m $MEMORY_GB
        
        # Check the exit status of the SPAdes command
        if [ $? -ne 0 ]; then
            echo "     -> SPAdes assembly for $ID failed. Skipping to the next sample."
            echo "$ID" >> failed.txt
            continue # Move to the next iteration of the loop
        fi

        # Copy final contigs to a cleaner filename
        cp "${ASSEMBLY_DIR}/contigs.fasta" "$ASSEMBLY_FASTA"
    fi

    # --- Step 4: seqtk Length Filtering ---
    echo "[4/4] Filtering contigs by minimum length (${LENGTH_FILTER}bp)..."
    FINAL_ASSEMBLY="04_final_assemblies/${ID}.fasta"
    
    if [ -f "$FINAL_ASSEMBLY" ]; then
        echo "     -> Final assembly for $ID already exists. Skipping."
    else
        if [ -f "$ASSEMBLY_FASTA" ] && [ -s "$ASSEMBLY_FASTA" ]; then
            seqtk seq -L $LENGTH_FILTER "$ASSEMBLY_FASTA" > "$FINAL_ASSEMBLY"
            echo "Assembly of $ID is complete."
        else
            echo "     -> Error: No valid assembly found for $ID"
            echo "$ID" >> failed.txt
            continue
        fi
    fi
done

# --- Step 5: Run CheckM2 on all assemblies ---
echo "================================================="
echo "[5/6] Running CheckM2 on all final assemblies..."
echo "================================================="

# Switch to checkm2 conda environment
conda activate checkm2

# Run CheckM2 on the final assemblies directory
checkm2 predict \
    --database_path $checkm2_database \
    --input 04_final_assemblies \
    --output-directory 04_final_assemblies/checkm2_results \
    --threads "${THREADS}" \
    -x fasta \
    --force

conda deactivate

# --- Step 6: Run assembly-stats on all assemblies ---
echo "================================================="
echo "[6/6] Running assembly-stats on all final assemblies..."
echo "================================================="

# Switch back to illumina conda environment
conda activate illumina

# Run assembly-stats on all final assemblies and save to file
assembly-stats 04_final_assemblies/*.fasta > 04_final_assemblies/assembly_stats_report.txt

conda deactivate

echo "================================================="
echo "Pipeline completed!"
echo "================================================="
echo "Final assemblies are in: 04_final_assemblies/"
echo "CheckM2 results are in: 04_final_assemblies/checkm2_results/"
echo "Assembly statistics are in: 04_final_assemblies/assembly_stats_report.txt"
echo "Failed samples, if any, are listed in: failed.txt"
