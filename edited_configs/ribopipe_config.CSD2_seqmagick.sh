# Configuration file for the ribosome profiling pipeline

# Universal options
EXPERIMENT_NAME="CSD2_seqmagick"
OUTDIR="/hdd/common/proj/RibosomeProfiling/results/2018-04-16/${EXPERIMENT_NAME}"
THREADS="16"
START_STEP=1 # Auto-stop after step 5 unless skipped. Continue from step 6.
HALT_S5=0
GENOME_FASTA="/ssd/common/tools/ribopipe/reference_sequences/syn_PCC6803/NC_000911.1_chr_7plasmids.fasta"

# Step-specific options

# Step 1: Import data to working directory
# None

# Step 2: Make one fastq.gz file per sample
# None

# Step 3: Assess read quality
# None

# Step 4: Trim away adapter sequences (PARALLEL)
cutadapt_a="AGATCGGAAGAGCACACGTCT"
cutadapt_O="6"
cutadapt_m="6"
cutadapt_n="3"
cutadapt_e="0.15"

# Step 5: Assess footprint-length distribution
# None

# Step 6: Filter reads with low-quality base calls (PARALLEL)
filter_q="25"
filter_l="6"

# Step 7: Remove rRNA and tRNA sequences
bowtie_7_n="2"
bowtie_7_l="28"
bowtie_7_p="10"
bowtie_7_ref="/ssd/common/tools/ribopipe/reference_sequences/syn_PCC6803/rRNAtRNA.syn6803.NC000911"

# Step 8: Map reads to the genome
bowtie_8_m="1"
bowtie_8_n="2"
bowtie_8_l="28"
bowtie_8_p="10"
bowtie_8_ref="/ssd/common/tools/ribopipe/reference_sequences/syn_PCC6803/NC_000911.1_chr_7plasmids"

# Step 9: Count the number of reads on read-occupied positions in genome (PARALLEL)
readCountScript="/ssd/common/tools/ribopipe/python/readCountScript.max48.sup2.3prime_mod.py"
min_length="25"
max_length="100"

# Step 10: Calculate total number of mapped reads
totalNbrMappedReadsScript="/ssd/common/tools/ribopipe/python/totalNbrMappedReadsScript.sup3.py"

# Step 11: Calculate RPM on read-occupied positions in genome
RPMscript="/ssd/common/tools/ribopipe/python/RPMscript.sup6.py"

# Step 12: Complete RPM list by assigning “0” to all unoccupied positions
RPMcompleteScript="/ssd/common/tools/ribopipe/python/RPMcompleteScript.sup7.edited.py"

# Step 13: Count the number of reads on every gene
readsPerGeneScript="/ssd/common/tools/ribopipe/python/readsPerGeneScript.sup4.edited.py"
genelistP="/ssd/common/tools/ribopipe/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist.p"
genelistM="/ssd/common/tools/ribopipe/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist.m"
