
### LOAD DATA ##################################################################

# Define input file
infile = "results/2018-04-04/gene_PauseScore_profiles.with_plasmids.1.tab.gz"

# Load the file
library(data.table)
gene_profiles = as.data.frame(fread(
  paste("gzip -dc ", infile), sep="\t", stringsAsFactors=F
  ))

### ANALYZE ####################################################################

# Clean up data
gene_profiles = gene_profiles[is.finite(gene_profiles$PauseScore),]

# Calculate distance from end of gene, and subset for 5 and 3 prime ends
gene_profiles$RelPosFromEnd = gene_profiles$RelPos - gene_profiles$Length + 1
gene_prof_5prime = subset(gene_profiles, RelPos %in% -50:0)

# Calculate number of zero and non-zero
gene_prof_5prime_zeros = aggregate(
  ifelse(PauseScore == 0, "0", ">0") ~ RelPos + Sample, gene_prof_5prime, table
  )
gene_prof_5prime_zeros = data.frame(
  RelPos = gene_prof_5prime_zeros$RelPos,
  Sample = gene_prof_5prime_zeros$Sample,
  Zero = gene_prof_5prime_zeros[,3][,1],
  AboveZero = gene_prof_5prime_zeros[,3][,2]
  )

# Should melt that
library(reshape2)
gene_prof_5prime_zeros = melt(gene_prof_5prime_zeros, id.vars=c("Sample", "RelPos"))
gene_prof_5prime_zeros$variable = factor(
  gene_prof_5prime_zeros$variable, levels = c("AboveZero", "Zero")
  )

### PLOT #######################################################################

library(ggplot2)
library(egg)

make_box = function(sample){
  df = subset(gene_prof_5prime, Sample == sample & PauseScore > 0)
  gp = ggplot(df, aes(RelPos, PauseScore, group=RelPos))
  gp = gp + geom_boxplot(outlier.size=0.4)
  gp = gp + theme_bw()
  gp = gp + scale_y_log10()
  gp = gp + ggtitle(paste("Sample ", sample))
  return(gp)
}

make_bar = function(sample){
  df = subset(gene_prof_5prime_zeros, Sample == sample)
  gp = ggplot(df, aes(RelPos, value, group=variable, fill=variable))
  gp = gp + geom_bar(stat="identity")
  gp = gp + theme_bw()
  gp = gp + scale_fill_manual(values=c("#6f6f6f", "#eeeeee"))
  gp = gp + theme(axis.title.x=element_blank())
  gp = gp + scale_x_continuous(position = "top")
  return(gp)
}

S = c("A","B","C","D","E")

outfile = "art/2018-04-04/5prime_PauseScore_boxplots.with_plasmids.1.pdf"
pdf(outfile, width=210/25.4, height=750/25.4, onefile=FALSE)
ggarrange(
  plots=c(rbind(lapply(S, make_box), lapply(S, make_bar))),
  ncol=1, heights=rep(c(3,1), 5)
  )
dev.off()
