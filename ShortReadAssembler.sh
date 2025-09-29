#!/bin/bash

#
# This script automates the assembly of short reads for all samples in the current directory.
# It expects paired-end reads in a directory named 'reads' ending in _R1.fastq.gz and _R2.fastq.gz.
#
# Workflow:
# 1. FastQC: Raw read quality control.
# 2. fastp: Adapter and quality trimming.
# 3. SPAdes: Genome assembly.
# 4. Racon: Polish assemblies with trimmed reads.
# 5. seqtk: Filter contigs by minimum length.
# 6. CheckM2: Quality assessment of all assemblies.
#

# --- Stop script on any error ---
# The following line is changed from 'set -e' to 'set +e' to allow the script to continue after an error.
set +e
set -o pipefail

# --- Activate illumina conda environment ---
eval "$(conda shell.bash hook)"
conda activate illumina

# ==============================================================================
#  CONFIGURATION (EDIT THESE VARIABLES)
# ==============================================================================
THREADS=120                                         # Number of CPU threads to use for each step
MEMORY_GB=800                                       # Max memory for SPAdes in Gigabytes
LENGTH_FILTER=500                                   # Minimum contig length to keep
RACON_ROUNDS=2                                      # Number of Racon polishing rounds
checkm2_database=/home/mcknid01/softwareDependencies/CheckM2_database/uniref100.KO.1.dmnd  # CheckM2 database path
# ==============================================================================

echo "Starting the assembly pipeline..."

# --- Create output directories ---
mkdir -p 01_fastqc 02_trimmed_reads 03_assemblies 04_racon_polish 05_final_assemblies
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
    mkdir -p "01_fastqc/$ID" "02_trimmed_reads/$ID" "03_assemblies/$ID" "04_racon_polish/$ID"

    # --- Step 1: FastQC ---
    echo "[1/6] Running FastQC on raw reads..."
    if [ -f "01_fastqc/$ID/${ID}_R1_fastqc.html" ]; then
        echo "     -> FastQC report for $ID already exists. Skipping."
    else
        fastqc -t $THREADS "$R1" "$R2" -o "01_fastqc/$ID"
    fi

    # --- Step 2: fastp Trimming ---
    echo "[2/6] Trimming reads with fastp..."
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
    echo "[3/6] Assembling with SPAdes..."
    ASSEMBLY_DIR="03_assemblies/$ID/${ID}_spades"
    ASSEMBLY_FASTA="03_assemblies/$ID/${ID}.fasta"
    if [ -f "$ASSEMBLY_FASTA" ]; then
        echo "     -> Assembly for $ID already exists. Skipping."
    else
        spades.py \
            --isolate \
            -1 "$TRIMMED_R1" -2 "$TRIMMED_R2" \
            -o "$ASSEMBLY_DIR" \
            -t $THREADS -m $MEMORY_GB
        
        # Check the exit status of the SPAdes command
        if [ $? -ne 0 ]; then
            echo "     -> SPAdes assembly for $ID failed. Skipping to the next sample."
            echo "$ID" >> failed.txt
            continue # Move to the next iteration of the loop
        fi

        # Copy final contigs to a cleaner filename
        cp "${ASSEMBLY_DIR}/contigs.fasta" "$ASSEMBLY_FASTA"
    fi

    # --- Step 4: Racon Polish ---
    echo "[4/6] Polishing assembly with Racon (${RACON_ROUNDS} rounds)..."
    FINAL_ASSEMBLY="05_final_assemblies/${ID}.fasta"
    
    if [ -f "$FINAL_ASSEMBLY" ]; then
        echo "     -> Final assembly for $ID already exists. Skipping."
    else
        # Combine trimmed reads for consistent polishing
        echo "     -> Combining paired-end reads..."
        COMBINED_READS="04_racon_polish/$ID/${ID}_combined.fastq.gz"
        if [ ! -f "$COMBINED_READS" ]; then
            cat "$TRIMMED_R1" "$TRIMMED_R2" > "$COMBINED_READS"
        fi
        
        # Initialize polishing with original assembly
        CURRENT_ASSEMBLY="$ASSEMBLY_FASTA"
        
        # Perform multiple rounds of Racon polishing
        for ((round=1; round<=RACON_ROUNDS; round++)); do
            echo "     -> Racon polishing round $round/$RACON_ROUNDS..."
            
            POLISHED_ASSEMBLY="04_racon_polish/$ID/${ID}_polished_round${round}.fasta"
            PAF_FILE="04_racon_polish/$ID/${ID}_aligned_round${round}.paf"
            
            # Skip if this round already completed
            if [ -f "$POLISHED_ASSEMBLY" ] && [ -s "$POLISHED_ASSEMBLY" ]; then
                echo "        -> Round $round already completed. Skipping."
                CURRENT_ASSEMBLY="$POLISHED_ASSEMBLY"
                continue
            fi
            
            # Align combined reads to current assembly using minimap2
            echo "        -> Aligning reads to assembly with minimap2..."
            minimap2 -x sr -t $THREADS "$CURRENT_ASSEMBLY" "$COMBINED_READS" > "$PAF_FILE"
            
            # Check if alignment was successful
            if [ ! -s "$PAF_FILE" ]; then
                echo "        -> Warning: No alignments found in round $round for $ID"
                break
            fi
            
            # Run Racon for polishing
            echo "        -> Running Racon polishing..."
            racon -t $THREADS "$COMBINED_READS" "$PAF_FILE" "$CURRENT_ASSEMBLY" > "$POLISHED_ASSEMBLY"
            
            # Check if polishing was successful
            if [ ! -s "$POLISHED_ASSEMBLY" ]; then
                echo "        -> Warning: Racon polishing failed in round $round for $ID"
                break
            fi
            
            # Update current assembly for next round
            CURRENT_ASSEMBLY="$POLISHED_ASSEMBLY"
            
            # Clean up intermediate alignment file to save space
            rm -f "$PAF_FILE"
        done
        
        # --- Step 5: seqtk Length Filtering ---
        echo "[5/6] Filtering contigs by minimum length (${LENGTH_FILTER}bp)..."
        FILTERED_ASSEMBLY="04_racon_polish/$ID/${ID}_filtered.fasta"
        
        # Use the final polished assembly for filtering
        if [ -f "$CURRENT_ASSEMBLY" ] && [ -s "$CURRENT_ASSEMBLY" ]; then
            seqtk seq -L $LENGTH_FILTER "$CURRENT_ASSEMBLY" > "$FILTERED_ASSEMBLY"
            
            # Copy filtered assembly to final directory
            cp "$FILTERED_ASSEMBLY" "$FINAL_ASSEMBLY"
            
            # Clean up combined reads file to save space
            rm -f "$COMBINED_READS"
            echo "Assembly of $ID is complete."
        else
            echo "     -> Error: No valid polished assembly found for $ID"
            echo "$ID" >> failed.txt
            continue
        fi
    fi
done

# --- Step 6: Run CheckM2 on all assemblies ---
echo "================================================="
echo "[6/6] Running CheckM2 on all final assemblies..."
echo "================================================="

# Switch to checkm2 conda environment
conda activate checkm2

# Run CheckM2 on the final assemblies directory
checkm2 predict \
    --database_path $checkm2_database \
    --input 05_final_assemblies \
    --output-directory 05_final_assemblies/checkm2_results \
    --threads "${THREADS}" \
    -x fasta\
    --force

conda deactivate
echo "================================================="
echo "Pipeline completed!"
echo "================================================="
echo "Final assemblies are in: 05_final_assemblies/"
echo "CheckM2 results are in: 05_final_assemblies/checkm2_results/"
echo "Failed samples, if any, are listed in: failed.txt"
