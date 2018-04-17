#!/usr/bin/env Rscript

### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Define input file
infile="analysis/gene_RPKM.tab"

### LOAD DATA ##################################################################

# Load the files
gen = read.table(infile, header=T, stringsAsFactors=F, sep="\t")

# Reshape data
library(reshape2)

gen$Name = gsub("'", '', gen$Name)

# Transform RPKM values
gen$RPKM_asinh = asinh(gen$RPKM)

# Cast into wide format
gen_wide = dcast(gen[,c("Sample","Name","RPKM_asinh")], Name ~ Sample, value.var="RPKM_asinh")

### PERFORM GRUBBS ONE OUTLIER TESTS ###########################################

# Find sets that contain one outlier data point
library(foreach)
library(doMC)
registerDoMC(16) # Use 16 cores for "foreach"

library(outliers)

# Perform

t_data_raw = foreach(i=1:nrow(gen_wide)) %dopar% {
  gen_wide_sub = log(gen_wide[i,][,2:dim(gen_wide)[2]])
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
early_RPKM = aggregate(RPKM ~ Name, subset(gen, Sample != "E"), mean)
late_RPKM = gen_wide[,c("Name","E")]
colnames(late_RPKM)[2] = "RPKM_late"
colnames(early_RPKM)[2] = "RPKM_early"

t_data = merge(t_data, merge(early_RPKM, late_RPKM))

# Calculate change
t_data$abs_change = t_data$RPKM_late - t_data$RPKM_early
t_data$log2_change = log2(t_data$RPKM_late / t_data$RPKM_early)

# Save Bonferroni adjusted test results table
write.table(
  t_data,
  "analysis/gene_Grubbs.Bonferroni_adjusted.asinh.tab",
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
t_data = merge(t_data, merge(early_RPKM, late_RPKM))

# Save Bonferroni adjusted test results table
write.table(
  t_data,
  "analysis/ribo_Grubbs.BH_adjusted.asinh.tab",
  quote=F, row.names=F, col.names=T, sep="\t"
  )
