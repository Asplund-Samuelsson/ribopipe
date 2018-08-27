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
    mkdir analysis

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
# Step 3: Assess read quality
S=3

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Assessing read quality with FastQC...\e[0m\n"

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
    echo -e "\n\e[94mStep $S: Trimming adapter sequences...\e[0m\n"

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
    echo -e "\n\e[94mStep $S: Assessing footprint-length distribution...\e[0m\n"

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


    if [[ HALT_S5 -eq 1 ]]
	  then
	    # Halt script
   	    echo -e "\n\e[33mStep $S: Halting. Check results and modify inputs.\e[0m\n"
	    exit 1
    fi
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 6: Filter reads with low-quality base calls (PARALLEL)
S=6

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Removing reads with low-quality base calls...\e[0m\n"

    # Store input filenames in array (specific for step 6; files from step 4)
    input_files_6=(cutadapt/${EXPERIMENT_NAME}.*.cutadapt.fastq.gz)

    # Define run_sickle function to be used with GNU parallel application
    run_filter() {
      # The infile name is stored in positional argument 1
      infile=$1
      # Extract the sample name from the infile
      sample_name=$(echo "$infile" | rev | cut -f 1 -d \/ | rev | cut -f 2 -d \.)
      # Construct outfile names
      prefix="highQuality/${EXPERIMENT_NAME}.${sample_name}"
      out_quality="${prefix}.quality.fastq"
      out_report="${prefix}.quality.report.txt"
      # Run seqmagick quality-filter with the supplied options
      seqmagick quality-filter \
      --min-mean-quality $filter_q --min-length $filter_l \
      $infile $out_quality > $out_report
    }

    # Export function and variables so that each subprocess can access them
    export -f run_filter
    export filter_q filter_l
    export EXPERIMENT_NAME

    # Run sickle in parallel for the input files
    parallel --no-notice --jobs $THREADS run_filter ::: ${input_files_6[@]}

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
    echo -e "\n\e[94mStep $S: Removing rRNA and tRNA sequences...\e[0m\n"

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
        bowtie -a --best --strata -t -n $bowtie_7_n -l $bowtie_7_l \
        -p $bowtie_7_p --un $out_fastq $bowtie_7_ref $infile /dev/null \
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
    echo -e "\n\e[94mStep $S: Mapping reads to the genome...\e[0m\n"

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
    echo -e "\n\e[94mStep $S: Counting reads on read-occupied positions in genome...\e[0m\n"

    # Store input filenames in array (specific for step 9; files from step 8)
    input_files_9=(mapped/${EXPERIMENT_NAME}.*.mapped.bwt1)

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
      $readCountScript -i $infile --min $min_length --max $max_length \
      --outP $out_p --outM $out_m
    }

    # Export function and variables so that each subprocess can access them
    export -f run_count
    export readCountScript min_length max_length
    export EXPERIMENT_NAME

    # Run read count in parallel for the input files
    parallel --no-notice --jobs $THREADS run_count ::: ${input_files_9[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 10: Calculate total number of mapped reads (PARALLEL)
S=10

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Counting total number of mapped reads...\e[0m\n"

    # Store input filenames in array (specific for step 10; files from step 9)
    input_files_10=(readcount/${EXPERIMENT_NAME}.*.readCount.p)

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
    parallel --no-notice --jobs $THREADS run_total ::: ${input_files_10[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 11: Calculate RPM on read-occupied positions in genome (PARALLEL)
S=11

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Calculating RPM on read-occupied positions in genome...\e[0m\n"

    # Store input filenames in array (specific for step 11; files from step 9)
    input_files_11=(readcount/${EXPERIMENT_NAME}.*.readCount.p)

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
    parallel --no-notice --jobs $THREADS run_rpm ::: ${input_files_11[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 12: Complete RPM list by assigning “0” to all unoccupied positions (PARALLEL)
S=12

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Completing RPM list with zeros...\e[0m\n"

    # Store input filenames in array (specific for step 12; files from step 11)
    input_files_12=(RPM/${EXPERIMENT_NAME}.*.RPM.p)

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
      $RPMcompleteScript --inP $infile_p --inM $infile_m --inG $GENOME_FASTA \
      --outP $out_p --outM $out_m
    }

    # Export function and variables so that each subprocess can access them
    export -f run_complete
    export RPMcompleteScript
    export EXPERIMENT_NAME
    export GENOME_FASTA

    # Run zero count completion in parallel for the input files
    parallel --no-notice --jobs $THREADS run_complete ::: ${input_files_12[@]}

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
    echo -e "\n\e[94mStep $S: Counting the number of reads on every gene...\e[0m\n"

    # Store input filenames in array (specific for step 13; files from step 9)
    input_files_13=(readcount/${EXPERIMENT_NAME}.*.readCount.p)

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
      --listP $genelistP --listM $genelistM --outP $out_p --outM $out_m \
      --inG $GENOME_FASTA
    }

    # Export function and variables so that each subprocess can access them
    export -f run_genes
    export readsPerGeneScript genelistP genelistM GENOME_FASTA
    export EXPERIMENT_NAME

    # Run reads per gene counting in parallel for the input files
    parallel --no-notice --jobs $THREADS run_genes ::: ${input_files_13[@]}

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 14: Create CDS RPKM table
S=14

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Creating CDS RPKM table...\e[0m\n"

    # Run script
    $table_script $GENE_LIST $SHIFT

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi


################################################################################
# Step 15: Plot average gene profiles
S=15

# Check starting step
if [[ START_STEP -le S ]]
  then
    # Report progress
    echo -e "\n\e[94mStep $S: Plotting average gene profiles...\e[0m\n"

    # Run script
    $profile_script $GENE_LIST $SHIFT

    # Report step done
    echo -e "\n\e[92mStep $S: Done.\e[0m\n"
  else
    # Report skip
    echo -e "\n\e[31mStep $S: Skipping...\e[0m\n"
fi
