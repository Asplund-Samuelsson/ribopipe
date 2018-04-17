# Configuration file for the ribosome profiling data analysis pipeline

# Universal options
EXPERIMENT_NAME="ExperimentName"
OUTDIR="YYYY.MM.DD.${EXPERIMENT_NAME}"
START_STEP=1 # Auto-stop after step 5 unless skipped. Continue from step 6.
RIBODIR="/home/johannes/proj/ribo/tools/ribopipe"
GENOME_FASTA="${RIBODIR}/reference_sequences/syn_PCC6803/NC_000911.1_chr_7plasmids.fasta"
GENE_LIST="${RIBODIR}/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist_full.tab"

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
rnadir="/hdd/common/proj/RNAseq/results/2018-03-29/CO2_lim_with_plasmids"

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
script_11="${RIBODIR}/data_analysis/step11.gene_profile_high_5prime_PS.R"

# Step 12: 5' boxplots
script_12="${RIBODIR}/data_analysis/step12.average_gene_profile_5prime_boxplots.R"

# Step 13: All shifts codon medians
script_13="${RIBODIR}/data_analysis/step13.codon_pause_score_distribution_3.R"
