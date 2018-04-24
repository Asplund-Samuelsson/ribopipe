#!/usr/bin/env bash

# Ribosome profiling data analysis pipeline
#
# Authors:
# Johannes Asplund-Samuelsson, KTH

# Input variables
# Read input variables from config file specified via command line
# e.g. "./analysis.sh analysis_config.sh"
source $1

# Work in output directory
cd $OUTDIR

# Create analysis results directory
mkdir analysis

################################################################################
# Step 1: Find codon positions
S=1

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Identifying codon positions...\e[0m\n"

    # Create CDS file
    grep -P "\tCDS\t" $GENE_LIST | cut -f 1,3,4,6,8 > analysis/CDS.tab

    # Identify codon positions
    $script_01 --inGen $GENOME_FASTA --inCDS analysis/CDS.tab > \
    analysis/codon_seqpos.tab

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 2: Codon PS distribution I
S=2

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Codon PS distribution (1)...\e[0m\n"

    # Run script
    $script_02 analysis/codon_seqpos.tab

    # Gzip big file
    pigz analysis/codon_PauseScore.tab

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 3: Gene Grubbs test
S=3

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Gene Grubbs test...\e[0m\n"

    # Run script
    $script_03

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 4: Stop codon Grubbs test
S=4

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Stop codon Grubbs test...\e[0m\n"

    # Run script
    $script_04

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 5: RNAseq vs RP
S=5

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: RNAseq vs RP...\e[0m\n"

    # Run script
    $script_05 $rnadir

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"

  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 6: Average gene profile
S=6

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Calculating average gene profile...\e[0m\n"

    # Run scripts
    $script_06

    # Gzip big file
    pigz analysis/gene_PauseScore_profiles.tab

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 7: Codon PS by sample
S=7

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Plotting codon PS by sample...\e[0m\n"

    # Run script
    $script_07 ${RIBODIR}/data_analysis/translation_table_11.tab \
    ${RIBODIR}/data_analysis/colours.txt

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 8: Codon PS distributions II
S=8

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Codon PS distributions (2)...\e[0m\n"

    # Run script
    $script_08 analysis/codon_seqpos.tab

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 9: Codon PS distributions III
S=9

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Codon PS distributions (3)...\e[0m\n"

    # Count codons
    cut -f 3 analysis/codon_seqpos.tab | sort | uniq -c | sed -e 's/^ \+//' | \
    sed -e 's/ \+/\t/' > analysis/count_codon.tab

    # Run script
    $script_09 ${RIBODIR}/data_analysis/translation_table_11.tab \
    ${RIBODIR}/data_analysis/colours.txt

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 10: Codon PS boxplots
S=10

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Creating codon PS boxplots...\e[0m\n"

    # Run script
    $script_10 ${RIBODIR}/data_analysis/translation_table_11.tab \
    ${RIBODIR}/data_analysis/colours.txt

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 11: High 5' PS genes
S=11

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Identifying high 5 prime PS genes...\e[0m\n"

    # Run scripts
    $script_11a
    $script_11b

    # Plot the genes
    mkdir analysis/high_5prime_PS_genes
    mkdir analysis/high_5prime_PS_genes/spike
    mkdir analysis/high_5prime_PS_genes/close
    mkdir analysis/high_5prime_PS_genes/closer
    mkdir analysis/high_5prime_PS_genes/significant

    cat analysis/high_5prime_PS_genes.txt | parallel --no-notice --jobs 7 \
    '/ssd/common/tools/ribopipe/plot_genes.R . {} operon -12 analysis/high_5prime_PS_genes/spike/{}.high_5prime_PS.spike.pdf'

    cat analysis/similar_5prime_PS_genes.txt | parallel --no-notice --jobs 7 \
    '/ssd/common/tools/ribopipe/plot_genes.R . {} operon -12 analysis/high_5prime_PS_genes/close/{}.similar_5prime_PS.pdf'

    cat analysis/25nt_5prime_PS_ge10_genes.txt | parallel --no-notice --jobs 7 \
    '/ssd/common/tools/ribopipe/plot_genes.R . {} operon -12 analysis/high_5prime_PS_genes/closer/{}.25nt_5prime_PS_ge10.pdf'

    cat analysis/significant_5prime_PS_genes.txt | parallel --no-notice --jobs 7 \
    '/ssd/common/tools/ribopipe/plot_genes.R . {} operon -12 analysis/high_5prime_PS_genes/closer/{}.significant_5prime_PS.pdf'

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 12: 5' boxplots
S=12

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Creating 5 prime boxplots...\e[0m\n"

    # Run script
    $script_12

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 13: All shifts codon medians
S=13

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Calculating codon medians for all shifts...\e[0m\n"

    # Run script
    $script_13 ${RIBODIR}/data_analysis/translation_table_11.tab \
    ${RIBODIR}/data_analysis/colours.txt

    # Gzip big file
    pigz analysis/codon_PauseScore.all_shifts.tab

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi
