# Configuration file for the ribosome profiling data analysis pipeline

# Universal options
EXPERIMENT_NAME="ExperimentName"
OUTDIR="YYYY.MM.DD.${EXPERIMENT_NAME}"
START_STEP=1
RIBODIR="/home/johannes/proj/ribo/tools/ribopipe"
GENOME_FASTA="${RIBODIR}/reference_sequences/syn_PCC6803/NC_000911.1_chr_7plasmids.fasta"
GENE_LIST="${RIBODIR}/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist_full.tab"
ANNOTATION_FILE="${RIBODIR}/genelists/syn_PCC6803/20180315_all_annotations_Micha.csv"
rnadir="/hdd/common/proj/RNAseq/results/2018-03-29/CO2_lim_with_plasmids"

# Step-specific options

# Step 1: Find codon positions
script_01="${RIBODIR}/data_analysis/step1.codon_positions.py"

# Step 2: Codon PS distribution I
script_02="${RIBODIR}/data_analysis/step2.codon_pause_score_distribution_1.R"

# Step 3: Gene Grubbs test
script_03="${RIBODIR}/data_analysis/step3.gene_Grubbs.R"

# Step 4: Stop codon Grubbs test
script_04="${RIBODIR}/data_analysis/step4.stop_Grubbs.R"

# Step 5: RNAseq vs RP
script_05="${RIBODIR}/data_analysis/step5.rnaseq_vs_ribprof.R"

# Step 6: Average gene profile
script_06="${RIBODIR}/data_analysis/step6.average_gene_profile.R"

# Step 7: Codon PS by sample
script_07="${RIBODIR}/data_analysis/step7.pause_score_A_vs_E.R"

# Step 8: Codon PS distributions II
script_08="${RIBODIR}/data_analysis/step8.codon_pause_score_distribution_2.R"

# Step 9: Codon PS distributions III
script_09="${RIBODIR}/data_analysis/step9.pause_score_ordered_codon_abundances.R"

# Step 10: Codon PS boxplots
script_10="${RIBODIR}/data_analysis/step10.codon_pause_score_boxplot.R"

# Step 11: High 5' PS
script_11a="${RIBODIR}/data_analysis/step11a.gene_profile_high_5prime_PS.R"
script_11b="${RIBODIR}/data_analysis/step11b.5prime_UTR_sampling_test.R"

# Step 12: 5' boxplots
script_12="${RIBODIR}/data_analysis/step12.average_gene_profile_5prime_boxplots.R"

# Step 13: All shifts codon medians
script_13="${RIBODIR}/data_analysis/step13.codon_pause_score_distribution_3.R"

# Step 14: Extract and align 5' UTRs
script_14a="${RIBODIR}/data_analysis/step14a.5prime_UTR_positions.R"
script_14b="${RIBODIR}/data_analysis/step14b.5prime_UTR_extraction.py"
script_14c="${RIBODIR}/data_analysis/step14c.5prime_UTRs_nucleotide_count.R"

# Step 15: Plot RPKM per process
script_15="${RIBODIR}/data_analysis/step15.RPKM_per_process.R"

# Step 16: Create translational efficiency table
script_16="${RIBODIR}/data_analysis/step16.create_TE_table.R"

# Step 17: Save CDS RPKM table without filtering
script_17="${RIBODIR}/data_analysis/step17.all_genes_RPKM_no_filter.R"
