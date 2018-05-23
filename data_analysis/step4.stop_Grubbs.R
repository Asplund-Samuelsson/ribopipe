#!/usr/bin/env Rscript

### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Define input file
codon_file="analysis/codon_PauseScore.tab.gz"

### LOAD DATA ##################################################################

# Load the files
library(data.table)
codon_rpm = as.data.frame(fread(
  paste(c("gzip -dc ", codon_file), collapse=""),
  header=T, stringsAsFactors=F, sep="\t"
  ))

# Subset codon PauseScores to stop codons
codon_rpm = subset(codon_rpm, Codon %in% c("TAA","TAG","TGA"))

# Reshape data
library(reshape2)

codon_rpm$Name = gsub("'", '', codon_rpm$Name)
gen_wide = dcast(codon_rpm[,c("Sample","Name","PauseScore")], Name ~ Sample, value.var="PauseScore")

nrow(gen_wide) == length(unique(gen_wide$Name)) # No genes with >1 stop codon

### PERFORM GRUBBS ONE OUTLIER TESTS ###########################################

# Find sets that contain one outlier data point
library(foreach)
library(doMC)
registerDoMC(16) # Use 16 cores for "foreach"

library(outliers)

# Perform

t_data_raw = foreach(i=1:nrow(gen_wide)) %dopar% {
  gen_wide_sub = asinh(gen_wide[i,][,2:dim(gen_wide)[2]])
  grubb_data = grubbs.test(as.numeric(gen_wide_sub))
  alt = unlist(strsplit(grubb_data$alternative, " "))
  outlier_value = alt[!is.na(suppressWarnings(as.numeric(alt)))]
  outlier = colnames(gen_wide_sub)[as.character(gen_wide_sub[1,]) == outlier_value]
  c(outlier, grubb_data$p.value)
}

# Format Grubbs test data
t_data = data.frame(
  Name=gen_wide$Name,
  Outlier=unlist(lapply(t_data_raw, '[[', 1)),
  p_val=as.numeric(unlist(lapply(t_data_raw, '[[', 2)))
  )

# Adjust p-values
t_data$p_adj = p.adjust(t_data$p_val, method="bonferroni")

# Add significance stars
t_data$sign = ifelse(
  t_data$p_adj < 0.05, ifelse(
    t_data$p_adj < 0.01, ifelse(
      t_data$p_adj < 0.001, "***", "**"
      ),
    "*"
    ),
  ""
  )

t_data$sign[is.na(t_data$sign)] = "" # Add no significance to NA p-values

# Add mean of four first time points and value of last time point
early_PS = aggregate(PauseScore ~ Name, subset(codon_rpm, Sample != "E"), mean)
late_PS = gen_wide[,c("Name","E")]
colnames(late_PS)[2] = "PS_late"
colnames(early_PS)[2] = "PS_early"

t_data = merge(t_data, merge(early_PS, late_PS))

# Calculate change
t_data$abs_change = t_data$PS_late - t_data$PS_early
t_data$log2_change = log2(t_data$PS_late / t_data$PS_early)

# Save Bonferroni adjusted test results table
write.table(
  t_data,
  "analysis/stop_Grubbs.Bonferroni_adjusted.tab",
  quote=F, row.names=F, col.names=T, sep="\t"
  )

# Adjust p-values
t_data$p_adj = p.adjust(t_data$p_val, method="BH")

# Add significance stars
t_data$sign = ifelse(
  t_data$p_adj < 0.05, ifelse(
    t_data$p_adj < 0.01, ifelse(
      t_data$p_adj < 0.001, "***", "**"
      ),
    "*"
    ),
  ""
  )

t_data$sign[is.na(t_data$sign)] = "" # Add no significance to NA p-values

# Use merge to order columns
t_data = merge(t_data, merge(early_PS, late_PS))

# Save Bonferroni adjusted test results table
write.table(
  t_data,
  "analysis/stop_Grubbs.BH_adjusted.tab",
  quote=F, row.names=F, col.names=T, sep="\t"
  )
