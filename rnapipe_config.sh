# Configuration file for the RNAseq pipeline

# Universal options
EXPERIMENT_NAME="ExperimentName"
OUTDIR="YYYY.MM.DD.${EXPERIMENT_NAME}"
THREADS="16"
START_STEP=1
# Step-specific options

# Step 1: Import data to working directory
# None

# Step 2: Make one fastq.gz file per sample
# None

# Step 3: Trim away adapter sequences (PARALLEL)
cutadapt_a="AGATCGGAAGAGCACACGTCT"
cutadapt_O="6"
cutadapt_m="6"
cutadapt_n="3"
cutadapt_e="0.15"

# Step 4: Trim ends with low-quality base calls (PARALLEL)
sickle_t="sanger"
sickle_q="20"
sickle_l="6"

# Step 5: Remove rRNA and tRNA sequences
bowtie_5_n="2"
bowtie_5_l="28"
bowtie_5_p="10"
bowtie_5_ref="/ssd/common/tools/bowtie1-1.1.2/indexes/rRNAtRNA.syn6803.NC000911"

# Step 6: Map reads to the genome
reverse_complement=true # Make input reads reverse-complement for e.g. NEBNext Directional
bowtie_6_m="1"
bowtie_6_n="2"
bowtie_6_l="28"
bowtie_6_p="10"
bowtie_6_ref="/ssd/common/tools/bowtie1-1.1.2/indexes/syn6803.NC_000911"

# Step 7: Count the number of reads on read-occupied positions in genome (PARALLEL)
readCountScript="/ssd/common/tools/ribopipe/python/readCountScript.max48.sup2.RNAseq_mod.py"

# Step 8: Calculate total number of mapped reads
totalNbrMappedReadsScript="/ssd/common/tools/ribopipe/python/totalNbrMappedReadsScript.sup3.py"

# Step 9: Calculate RPM on read-occupied positions in genome
RPMscript="/ssd/common/tools/ribopipe/python/RPMscript.sup6.py"

# Step 10: Complete RPM list by assigning “0” to all unoccupied positions
RPMcompleteScript="/ssd/common/tools/ribopipe/python/RPMcompleteScript.NC000911.sup7.py"

# Step 11: Count the number of reads on every gene
readsPerGeneScript="/ssd/common/tools/ribopipe/python/readsPerGeneScript.sup4.edited.py"
genelistP="/ssd/jan/auxillaryData/genelists/syn6803.NC000911.genelist_p"
genelistM="/ssd/jan/auxillaryData/genelists/syn6803.NC000911.genelist_m"
