# Configuration file for the RNAseq pipeline

# Universal options
EXPERIMENT_NAME="ExperimentName"
OUTDIR="YYYY.MM.DD.${EXPERIMENT_NAME}"
THREADS="16"
START_STEP=1
RIBODIR="/ssd/common/tools/ribopipe/"
GENOME_FASTA="${RIBODIR}/reference_sequences/syn_PCC6803/NC_000911.1_chr_7plasmids.fasta"
GENE_LIST="${RIBODIR}/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist_full.tab"

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
bowtie_5_ref="${RIBODIR}/reference_sequences/syn_PCC6803/rRNAtRNA.syn6803.NC000911"

# Step 6: Map reads to the genome
reverse_complement=true # Make input reads reverse-complement for e.g. NEBNext Directional
bowtie_6_m="1"
bowtie_6_n="2"
bowtie_6_l="28"
bowtie_6_p="10"
bowtie_6_ref="${RIBODIR}/reference_sequences/syn_PCC6803/NC_000911.1_chr_7plasmids"

# Step 7: Count the number of reads on read-occupied positions in genome (PARALLEL)
readCountScript="${RIBODIR}/python/readCountScript.max48.sup2.RNAseq_mod.py"

# Step 8: Calculate total number of mapped reads
totalNbrMappedReadsScript="${RIBODIR}/python/totalNbrMappedReadsScript.sup3.py"

# Step 9: Calculate RPM on read-occupied positions in genome
RPMscript="${RIBODIR}/python/RPMscript.sup6.py"

# Step 10: Complete RPM list by assigning “0” to all unoccupied positions
RPMcompleteScript="${RIBODIR}/python/RPMcompleteScript.sup7.edited.py"

# Step 11: Count the number of reads on every gene
readsPerGeneScript="${RIBODIR}/python/readsPerGeneScript.sup4.edited.py"
genelistP="${RIBODIR}/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist.p"
genelistM="${RIBODIR}/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist.m"

# Step 12: Create CDS RPKM table
table_script="${RIBODIR}/R/create_gene_rpkm_table.rnapipe.R"
pca_script="${RIBODIR}/R/simple_3D-PCA_clustering_of_samples.R"
