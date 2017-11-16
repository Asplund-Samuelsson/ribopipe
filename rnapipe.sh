#!/usr/bin/env bash

# RNAseq pipeline (based on "Ribosome profiling pipeline")
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
    mkdir \
    cutadapt highQuality tANDrRNAremoval mapped readcount RPM readsPerGene

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
    echo -e "\n\e[94mStep $S: Concatenating fastq files per sample...\e[0m\n"

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
# Step 3: Trim away adapter sequences (PARALLEL)
S=3

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Trimming adapter sequences...\e[0m\n"

    # Store input filenames in array (specific for step 3; files from step 2)
    input_files_3=(rawfastqgz/${EXPERIMENT_NAME}.*.fastq.gz)

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
      out_report="${prefix}.cutadapt.report.txt"
      # Run cutadapt with the supplied options
      cutadapt -a $cutadapt_a -O $cutadapt_O -m $cutadapt_m -n $cutadapt_n \
      -e $cutadapt_e --too-short-output=$out_tooShort \
      -o $out_cutadapt $infile > $out_report
    }

    # Export function and variables so that each subprocess can access them
    export -f run_cutadapt
    export cutadapt_a cutadapt_O cutadapt_m cutadapt_n cutadapt_e
    export EXPERIMENT_NAME

    # Run cutadapt in parallel for the input files
    parallel --no-notice --jobs $THREADS run_cutadapt ::: ${input_files_3[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 4: Trim ends with low-quality base calls (PARALLEL)
S=4

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Trimming ends with low-quality base calls...\e[0m\n"

    # Store input filenames in array (specific for step 4; files from step 3)
    input_files_4=(cutadapt/${EXPERIMENT_NAME}.*.cutadapt.fastq.gz)

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
      # Run sickle with the supplied options
      sickle se -f $infile -t $sickle_t -q $sickle_q -l $sickle_l \
      -o $out_quality > $out_report
    }

    # Export function and variables so that each subprocess can access them
    export -f run_sickle
    export sickle_t sickle_q sickle_l
    export EXPERIMENT_NAME

    # Run sickle in parallel for the input files
    parallel --no-notice --jobs $THREADS run_sickle ::: ${input_files_4[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 5: Remove rRNA and tRNA sequences
S=5

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Removing rRNA and tRNA sequences...\e[0m\n"

    # Calculate number of threads to use
    # Use lowest number out of THREADS and bowtie_p
    bowtie_5_p=$(dc -e "[${THREADS}]sM ${bowtie_5_p}d ${THREADS}<Mp")

    # Run bowtie1 for each input file
    ls highQuality/*quality.fastq | while read infile
      do
        # Create output filenames for each input file
        out_fastq=$(echo $infile | sed -e 's/highQuality/tANDrRNAremoval/' | \
        sed -e 's/quality/tANDrRNAdeplete/')
        out_sam=$(echo $out_fastq | sed -e 's/deplete\.fastq/.sam/')
        out_report=$(echo $out_fastq | sed -e 's/\.fastq/.report.txt/')
        # Run bowtie1
        bowtie -a --best --strata -t -n $bowtie_5_n -l $bowtie_5_l -S \
        -p $bowtie_5_p --un $out_fastq $bowtie_5_ref $infile $out_sam \
        2> $out_report
      done

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 6: Map reads to the genome
S=6

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Mapping reads to the genome...\e[0m\n"

    # Calculate number of threads to use
    # Use lowest number out of THREADS and bowtie_p
    bowtie_6_p=$(dc -e "[${THREADS}]sM ${bowtie_6_p}d ${THREADS}<Mp")

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
        bowtie -a --best --strata -m $bowtie_6_m -n $bowtie_6_n -l $bowtie_6_l \
        -p $bowtie_6_p --un $out_not --max $out_more $bowtie_6_ref $infile \
        $out_mapped 2> $out_report
      done

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 7: Count the number of reads on read-occupied positions in genome (PARALLEL)
S=7

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Counting reads on read-occupied positions in genome...\e[0m\n"

    # Store input filenames in array (specific for step 7; files from step 6)
    input_files_7=(mapped/${EXPERIMENT_NAME}.*.mapped.bwt1)

    # Define run_count function to be used with GNU parallel application
    run_count() {
      # The infile name is stored in positional argument 1
      infile=$1
      # Extract the sample name from the infile
      sample_name=$(echo "$infile" | rev | cut -f 1 -d \/ | rev | cut -f 2 -d \.)
      # Construct outfile names
      prefix="readcount/${EXPERIMENT_NAME}.${sample_name}"
      out_p="${prefix}.readCount.p"
      out_m="${prefix}.readCount.m"
      # Run the read count script with the supplied options
      $readCountScript -i $infile --outP $out_p --outM $out_m
    }

    # Export function and variables so that each subprocess can access them
    export -f run_count
    export readCountScript min_length max_length
    export EXPERIMENT_NAME

    # Run read count in parallel for the input files
    parallel --no-notice --jobs $THREADS run_count ::: ${input_files_7[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 8: Calculate total number of mapped reads (PARALLEL)
S=8

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Counting total number of mapped reads...\e[0m\n"

    # Store input filenames in array (specific for step 8; files from step 7)
    input_files_8=(readcount/${EXPERIMENT_NAME}.*.readCount.p)

    # Define run_total function to be used with GNU parallel application
    run_total() {
      # The "p" infile name is stored in positional argument 1
      infile_p=$1
      # Reconstruct the corresponding "m" infile name
      infile_m=$(echo $infile_p | sed -e 's/\.p$/.m/')
      # Extract the sample name from the infile
      sample_name=$(echo "$infile_p" | rev | cut -f 1 -d \/ | rev | cut -f 2 -d \.)
      # Construct outfile name
      prefix="readcount/${EXPERIMENT_NAME}.${sample_name}"
      out="${prefix}.totalNbrMappedReads"
      # Run the total read count script with the supplied options
      $totalNbrMappedReadsScript --inP $infile_p --inM $infile_m --out $out
    }

    # Export function and variables so that each subprocess can access them
    export -f run_total
    export totalNbrMappedReadsScript
    export EXPERIMENT_NAME

    # Run total count in parallel for the input files
    parallel --no-notice --jobs $THREADS run_total ::: ${input_files_8[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 9: Calculate RPM on read-occupied positions in genome (PARALLEL)
S=9

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Calculating RPM on read-occupied positions in genome...\e[0m\n"

    # Store input filenames in array (specific for step 9; files from step 7)
    input_files_9=(readcount/${EXPERIMENT_NAME}.*.readCount.p)

    # Define run_rpm function to be used with GNU parallel application
    run_rpm() {
      # The "p" infile name is stored in positional argument 1
      infile_p=$1
      # Reconstruct the corresponding "m" infile name
      infile_m=$(echo $infile_p | sed -e 's/\.p$/.m/')
      # Reconstruct the corresponding "total" infile name (from step 10)
      infile_tot=$(echo $infile_p | sed -e 's/\.readCount\.p$/.totalNbrMappedReads/')
      # Extract the sample name from the infile
      sample_name=$(echo "$infile_p" | rev | cut -f 1 -d \/ | rev | cut -f 2 -d \.)
      # Construct outfile names
      prefix="RPM/${EXPERIMENT_NAME}.${sample_name}"
      out_p="${prefix}.RPM.p"
      out_m="${prefix}.RPM.m"
      # Run the read count script with the supplied options
      $RPMscript --inP $infile_p --inM $infile_m --number $infile_tot \
      --outP $out_p --outM $out_m
    }

    # Export function and variables so that each subprocess can access them
    export -f run_rpm
    export RPMscript
    export EXPERIMENT_NAME

    # Run total count in parallel for the input files
    parallel --no-notice --jobs $THREADS run_rpm ::: ${input_files_9[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 10: Complete RPM list by assigning “0” to all unoccupied positions (PARALLEL)
S=10

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Completing RPM list with zeros...\e[0m\n"

    # Store input filenames in array (specific for step 10; files from step 9)
    input_files_10=(RPM/${EXPERIMENT_NAME}.*.RPM.p)

    # Define run_complete function to be used with GNU parallel application
    run_complete() {
      # The "p" infile name is stored in positional argument 1
      infile_p=$1
      # Reconstruct the corresponding "m" infile name
      infile_m=$(echo $infile_p | sed -e 's/\.p$/.m/')
      # Extract the sample name from the infile
      sample_name=$(echo "$infile_p" | rev | cut -f 1 -d \/ | rev | cut -f 2 -d \.)
      # Construct outfile names
      prefix="RPM/${EXPERIMENT_NAME}.${sample_name}"
      out_p="${prefix}.RPM0.p"
      out_m="${prefix}.RPM0.m"
      # Run the zero count completion script with the supplied options
      $RPMcompleteScript --inP $infile_p --inM $infile_m \
      --outP $out_p --outM $out_m
    }

    # Export function and variables so that each subprocess can access them
    export -f run_complete
    export RPMcompleteScript
    export EXPERIMENT_NAME

    # Run zero count completion in parallel for the input files
    parallel --no-notice --jobs $THREADS run_complete ::: ${input_files_10[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 11: Count the number of reads on every gene
S=11

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Counting the number of reads on every gene...\e[0m\n"

    # Store input filenames in array (specific for step 11; files from step 7)
    input_files_11=(readcount/${EXPERIMENT_NAME}.*.readCount.p)

    # Define run_genes function to be used with GNU parallel application
    run_genes() {
      # The "p" infile name is stored in positional argument 1
      infile_p=$1
      # Reconstruct the corresponding "m" infile name
      infile_m=$(echo $infile_p | sed -e 's/\.p$/.m/')
      # Extract the sample name from the infile
      sample_name=$(echo "$infile_p" | rev | cut -f 1 -d \/ | rev | cut -f 2 -d \.)
      # Construct outfile names
      prefix="readsPerGene/${EXPERIMENT_NAME}.${sample_name}"
      out_p="${prefix}.readsPerGene.p"
      out_m="${prefix}.readsPerGene.m"
      # Run the reads per gene script with the supplied options
      $readsPerGeneScript --inP $infile_p --inM $infile_m \
      --listP $genelistP --listM $genelistM --outP $out_p --outM $out_m
    }

    # Export function and variables so that each subprocess can access them
    export -f run_genes
    export readsPerGeneScript genelistP genelistM
    export EXPERIMENT_NAME

    # Run reads per gene counting in parallel for the input files
    parallel --no-notice --jobs $THREADS run_genes ::: ${input_files_11[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi
