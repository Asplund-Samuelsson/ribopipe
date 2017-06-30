#!/usr/bin/env bash

# Ribosome profiling pipeline
#
# Authors:
# Jan Karlsen, KTH
# Johannes Asplund-Samuelsson, KTH

# Input variables
# Read input variables from config file specified via command line
# e.g. "./ribopipe.sh ribopipe_config.sh"
source $1

# Work in output directory
cd $OUTDIR

# Manually place input fastq.gz files in rawfastqgz subdirectory
# Store input filenames in array
INPUT_FILES=(rawfastqgz/*fastq.gz)

################################################################################
# Step 1: Import data to working directory
S=1

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Setting up output directory...\e[0m\n"

    # Set up remaining output subdirectories automatically
    mkdir qualityCheck cutadapt lengthDistr highQuality
    mkdir tANDrRNAremoval mapped readcount RPM readsPerGene

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 2: Make one fastq.gz file per sample
S=2

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Concatenating fastq files per sample\e[0m\n"

    # Concatenate input files
    # Iterate over all infiles
    for infile in "${INPUT_FILES[@]}"
      do
        # Extract the sample name from each infile
        sample_name=$(echo "$infile" | rev | cut -f 1 -d \/ | rev | cut -f 1 -d _)
        # Construct an outfile name
        outfile="rawfastqgz/${EXPERIMENT_NAME}.${sample_name}.fastq"
        # Unzip the infile, keeping it intact, and add data to end of outfile
        gunzip -c $infile >> $outfile
      done

    # Gzip concatenated files
    pigz rawfastqgz/*fastq

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 3: Assess read quality
S=3

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: FastQC read quality assessment\e[0m\n"

    # Calculate number of threads to use
    # Count concatenated input files
    number_files=$(ls rawfastqgz/${EXPERIMENT_NAME}.*.fastq.gz | wc -l)
    # Pick lowest of number_files and THREADS
    threads_fastqc=$(dc -e "[${THREADS}]sM ${number_files}d ${THREADS}<Mp")

    # Run fastqc on all concatenated input files
    fastqc -t $threads_fastqc -o qualityCheck \
    rawfastqgz/${EXPERIMENT_NAME}.*.fastq.gz

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 4: Trim away adapter sequences (PARALLEL)
S=4

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Trim adapter sequences\e[0m\n"

    # Store input filenames in array (specific for step 4; files from step 2)
    input_files_4=(rawfastqgz/${EXPERIMENT_NAME}.*.fastq.gz)

    # Define run_cutadapt function to be used with GNU parallel application
    run_cutadapt() {
      # The infile name is stored in positional argument 1
      infile=$1
      # Extract the sample name from the infile
      sample_name=$(echo "$infile" | rev | cut -f 1 -d \/ | rev | cut -f 2 -d \.)
      # Construct outfile names
      prefix="cutadapt/${EXPERIMENT_NAME}.${sample_name}"
      out_cutadapt="${prefix}.cutadapt.fastq.gz"
      out_tooShort="${prefix}.tooShort.fastq.gz"
      out_untrimmed="${prefix}.untrimmed.fastq.gz"
      out_report="${prefix}.cutadapt.report.txt"
      # Run cutadapt with the supplied options
      cutadapt -a $cutadapt_a -O $cutadapt_O -m $cutadapt_m -n $cutadapt_n \
      -e $cutadapt_e --too-short-output=$out_tooShort \
      --untrimmed-output=$out_untrimmed -o $out_cutadapt $infile > $out_report
    }

    # Export function and variables so that each subprocess can access them
    export -f run_cutadapt
    export cutadapt_a cutadapt_O cutadapt_m cutadapt_n cutadapt_e
    export EXPERIMENT_NAME

    # Run cutadapt in parallel for the input files
    parallel --no-notice --jobs $THREADS run_cutadapt ::: ${input_files_4[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 5: Assess footprint-length distribution
S=5

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Assess footprint-length distribution\e[0m\n"

    # Calculate number of threads to use
    # Count concatenated input files
    number_files=$(ls cutadapt/${EXPERIMENT_NAME}.*.cutadapt.fastq.gz | wc -l)
    # Pick lowest of number_files and THREADS
    threads_fastqc=$(dc -e "[${THREADS}]sM ${number_files}d ${THREADS}<Mp")

    # Run fastqc on all concatenated input files
    fastqc -t $threads_fastqc -o lengthDistr \
    cutadapt/${EXPERIMENT_NAME}.*.cutadapt.fastq.gz

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"

    # Halt script
    echo -e "\n\e[33mStep $S: Halting. Check results and modify inputs.\e[0m\n"
    exit 1
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 6: Trim ends with low-quality base calls (PARALLEL)
S=6

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Trim ends with low-quality base calls\e[0m\n"

    # Store input filenames in array (specific for step 6; files from step 4)
    input_files_6=(cutadapt/${EXPERIMENT_NAME}.*.cutadapt.fastq.gz)

    # Define run_sickle function to be used with GNU parallel application
    run_sickle() {
      # The infile name is stored in positional argument 1
      infile=$1
      # Extract the sample name from the infile
      sample_name=$(echo "$infile" | rev | cut -f 1 -d \/ | rev | cut -f 2 -d \.)
      # Construct outfile names
      prefix="highQuality/${EXPERIMENT_NAME}.${sample_name}"
      out_quality="${prefix}.quality.fastq"
      out_report="${prefix}.quality.report.txt"
      # Run cutadapt with the supplied options
      sickle se -f $infile -t $sickle_t -q $sickle_q -l $sickle_l \
      -o $out_quality > $out_report
    }

    # Export function and variables so that each subprocess can access them
    export -f run_sickle
    export sickle_t sickle_q sickle_l
    export EXPERIMENT_NAME

    # Run sickle in parallel for the input files
    parallel --no-notice --jobs $THREADS run_sickle ::: ${input_files_6[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 7: Remove rRNA and tRNA sequences
S=7

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Remove rRNA and tRNA sequences\e[0m\n"

    # Calculate number of threads to use
    # Use lowest number out of THREADS and bowtie_p
    bowtie_7_p=$(dc -e "[${THREADS}]sM ${bowtie_7_p}d ${THREADS}<Mp")

    # Run bowtie1 for each input file
    ls highQuality/*quality.fastq | while read infile
      do
        # Create output filenames for each input file
        out_fastq=$(echo $infile | sed -e 's/highQuality/tANDrRNAremoval/' | \
        sed -e 's/quality/tANDrRNAdeplete/')
        out_sam=$(echo $out_fastq | sed -e 's/deplete\.fastq/.sam/')
        out_report=$(echo $out_fastq | sed -e 's/\.fastq/.report.txt/')
        # Run bowtie1
        bowtie -a --best --strata -t -n $bowtie_7_n -l $bowtie_7_l -S \
        -p $bowtie_7_p --un $out_fastq $bowtie_7_ref $infile $out_sam \
        2> $out_report
      done

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 8: Map reads to the genome
S=8

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Map reads to the genome\e[0m\n"

    # Calculate number of threads to use
    # Use lowest number out of THREADS and bowtie_p
    bowtie_8_p=$(dc -e "[${THREADS}]sM ${bowtie_8_p}d ${THREADS}<Mp")

    # Run bowtie1 for each input file
    ls tANDrRNAremoval/*tANDrRNAdeplete.fastq | while read infile
      do
        # Create output filenames for each input file
        out_not=$(echo $infile | sed -e 's/tANDrRNAremoval/mapped/' | \
        sed -e 's/tANDrRNAdeplete/notmapped/')
        out_more=$(echo $out_not | sed -e 's/\.notmapped\./.moremapped./')
        out_mapped=$(echo $out_not | sed -e 's/\.fastq/.bwt1/' | \
        sed -e 's/\.notmapped\./.mapped./')
        out_report=$(echo $out_not | sed -e 's/\.notmapped\.fastq/.report.txt/')
        # Run bowtie1
        bowtie -a --best --strata -m $bowtie_8_m -n $bowtie_8_n -l $bowtie_8_l \
        -p $bowtie_8_p --un $out_not --max $out_more $bowtie_8_ref $infile \
        $out_mapped 2> $out_report
      done

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 9: Count the number of reads on read-occupied positions in genome (PARALLEL)
S=9

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Count reads on read-occupied positions in genome\e[0m\n"



    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 10: Calculate total number of mapped reads
S=10

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Count total number of mapped reads\e[0m\n"



    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 11: Calculate RPM on read-occupied positions in genome
S=11

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Calculate RPM on read-occupied positions in genome\e[0m\n"



    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 12: Complete RPM list by assigning “0” to all unoccupied positions
S=12

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Complete RPM list with zeros\e[0m\n"



    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 13: Count the number of reads on every gene
S=13

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Count the number of reads on every gene\e[0m\n"



    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi
