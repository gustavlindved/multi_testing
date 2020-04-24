# Load data 
setwd('~/Google Drev/DTU/10. Semester/Applied Statistics/Project 3/')
array <- read.delim('rawdata.tab', sep = '\t')
rownames(array) <- array$geneID

# Set parameters
m <- dim(array)[1]
alpha <- 0.05

# Install DESeq2
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")
library('DESeq2')

# Load into DESeq dataframe
# Choose two conditions (ABS and BS)
array <- array[,c('ABS1','ABS2','ABS3','ABS4','ABS5','ABS6', 'BS1','BS2','BS3','BS4','BS5', 'BS6')]

# Make coldata matrix for dds object
coldata <- data.frame(c(rep('ABS',6), rep('BS', 6)))
rownames(coldata) <- colnames(array)
colnames(coldata) <- 'condition'

# Make dds object
dds <- DESeqDataSetFromMatrix(countData = array, colData = coldata, design = ~condition)

# Run DESeq test and get results
dds <- DESeq(dds)
res <- results(dds)

# Create matrix for testing statistics
tests <- data.frame(rep(NA, m))
row.names(tests) <- rownames(array)
colnames(tests) <- 'pvalue'

# Add p-values from DESeq
tests$pvalue_de <- res$pvalue

# Add p-values from Mann-whitney
tests$pvalue_mw <- rep(NA, m)
for (i in 1:m){
  tests$pvalue_mw[i] <- wilcox.test(as.numeric(array[i,1:6]), as.numeric(array[i,7:12]))$p.value
}

# Significant without bonferroni
tests$no_correction_de <- tests$pvalue_de < alpha
tests$no_correction_mw <- tests$pvalue_mw < alpha

# Significant tests with bonferroni correction
tests$p_bonferroni_de <- p.adjust(tests$pvalue_de, method = 'bonferroni')
tests$sig_bonferroni_de <- tests$p_bonferroni_de < alpha

tests$p_bonferroni_mw <- p.adjust(tests$pvalue_mw, method = 'bonferroni')
tests$sig_bonferroni_mw <- tests$p_bonferroni_mw < alpha

# Significant tests with holm correction
tests$p_holm_de <- p.adjust(tests$pvalue_de, method = 'holm')
tests$sig_holm_de <- tests$p_holm_de < alpha

tests$p_holm_mw <- p.adjust(tests$pvalue_mw, method = 'holm')
tests$sig_holm_mw <- tests$p_holm_mw < alpha

# Significant tests with BH correction
tests$p_BH_de <- p.adjust(tests$pvalue_de, method = 'BH')
tests$sig_BH_de <- tests$p_BH_de < alpha

tests$p_BH_mw <- p.adjust(tests$pvalue_mw, method = 'BH')
tests$sig_BH_mw <- tests$p_BH_mw < alpha

# Significant tests with Storeys q-value
#BiocManager::install("qvalue")
library('qvalue')
qobj_de <- qvalue(p = tests$pvalue_de)
tests$q_value_de <- qobj_de$qvalues
tests$sig_qvalue_de <- tests$q_value_de < alpha

qobj_mw <- qvalue(p = tests$pvalue_mw)
tests$q_value_mw <- qobj_mw$qvalues
tests$sig_qvalue_mw <- tests$q_value_mw < alpha

# How many significant?
cat('DESeq')
cat('Number of significant genes without correction is', sum(tests$no_correction_de, na.rm = T))
cat('Number of significant genes with bonferroni correction is', sum(tests$sig_bonferroni_de, na.rm = T))
cat('Number of significant genes with Holm correction is', sum(tests$sig_holm_de, na.rm = T))
cat('Number of significant genes with Benjamini-Hochberg correction is', sum(tests$sig_BH_de, na.rm = T))
cat('Number of significant genes with q-value correction is', sum(tests$sig_qvalue_de, na.rm = T))

cat('Mann Whitney')
cat('Number of significant genes without correction is', sum(tests$no_correction_mw, na.rm = T))
cat('Number of significant genes with bonferroni correction is', sum(tests$sig_bonferroni_mw, na.rm = T))
cat('Number of significant genes with Holm correction is', sum(tests$sig_holm_mw, na.rm = T))
cat('Number of significant genes with Benjamini-Hochberg correction is', sum(tests$sig_BH_mw, na.rm = T))
cat('Number of significant genes with q-value correction is', sum(tests$sig_qvalue_mw, na.rm = T))

