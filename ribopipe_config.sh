# Configuration file for the ribosome profiling pipeline

# Universal options
EXPERIMENT_NAME="ExperimentName"
OUTDIR="YYYY.MM.DD.${EXPERIMENT_NAME}"
THREADS="16"

# Step-specific options

# Step 1: Import data to working directory
# None

# Step 2: Make one fastq.gz file per sample
# None

# Step 3: Assess read quality
# None

# Step 4: Trim away adapter sequences (PARALLEL)

# Step 5: Assess footprint-length distribution

# Step 6: Trim ends with low-quality base calls (PARALLEL)

# Step 7: Remove rRNA and tRNA sequences

# Step 8: Map reads to the genome (PARALLEL)

# Step 9: Count the number of reads on read-occupied positions in genome (PARALLEL)

# Step 10: Calculate total number of mapped reads

# Step 11: Calculate RPM on read-occupied positions in genome

# Step 12: Complete RPM list by assigning “0” to all unoccupied positions

# Step 13: Count the number of reads on every gene
