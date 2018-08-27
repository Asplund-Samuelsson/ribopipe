# Configuration file for the ribosome profiling pipeline

# Universal options
EXPERIMENT_NAME="ExperimentName"
OUTDIR="YYYY.MM.DD.${EXPERIMENT_NAME}"
THREADS="16"
START_STEP=1 # Auto-stop after step 5 unless skipped. Continue from step 6.
HALT_S5=1
RIBODIR="/ssd/common/tools/ribopipe/"
GENOME_FASTA="${RIBODIR}/reference_sequences/syn_PCC6803/NC_000911.1_chr_7plasmids.fasta"
GENE_LIST="${RIBODIR}/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist_full.tab"
SHIFT="-12"

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
bowtie_7_ref="${RIBODIR}/reference_sequences/syn_PCC6803/rRNAtRNA.syn6803.NC000911"

# Step 8: Map reads to the genome
bowtie_8_m="1"
bowtie_8_n="2"
bowtie_8_l="28"
bowtie_8_p="10"
bowtie_8_ref="${RIBODIR}/reference_sequences/syn_PCC6803/NC_000911.1_chr_7plasmids"

# Step 9: Count the number of reads on read-occupied positions in genome (PARALLEL)
readCountScript="${RIBODIR}/python/readCountScript.max48.sup2.3prime_mod.py"
min_length="25"
max_length="100"

# Step 10: Calculate total number of mapped reads
totalNbrMappedReadsScript="${RIBODIR}/python/totalNbrMappedReadsScript.sup3.py"

# Step 11: Calculate RPM on read-occupied positions in genome
RPMscript="${RIBODIR}/python/RPMscript.sup6.py"

# Step 12: Complete RPM list by assigning “0” to all unoccupied positions
RPMcompleteScript="${RIBODIR}/python/RPMcompleteScript.sup7.edited.py"

# Step 13: Count the number of reads on every gene
readsPerGeneScript="${RIBODIR}/python/readsPerGeneScript.sup4.edited.py"
genelistP="${RIBODIR}/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist.p"
genelistM="${RIBODIR}/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist.m"

# Step 14: Create CDS RPKM table
table_script="${RIBODIR}/R/create_gene_rpkm_table.R"

# Step 15: Plot average gene profiles
profile_script="${RIBODIR}/R/average_gene_profile.R"
