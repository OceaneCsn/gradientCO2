library(DIANE)
library(tidyverse)
library(patchwork)
library(splines)
library(edgeR)

# expression dataset
load('rdata/expression_data_no_mismatch.rdata')
# list of regulators, same as in the integration manuscript, 
# touttaken from ATFDB et plnTFDB
load('rdata/regulators.rdata')

#normalization
conditions <- str_split_fixed(colnames(data), '_', 2)[,1]
tcc <- normalize(data, norm_method = 'tmm', iteration = FALSE, 
                 conditions = conditions)
tcc <- filter_low_counts(tcc, 10*length(conditions))
normalized_counts <- TCC::getNormalizedData(tcc)
draw_PCA(normalized_counts)

# Differential expression via splines
DF <- 2
FDR <- 0.001

############### CO2*N design
CO2 <- as.numeric(str_split_fixed(conditions, '\\.', 2)[,1])
CO2_spline <- ns(CO2, df = DF)

N <- str_split_fixed(conditions, '\\.', 2)[,2]
design <- model.matrix(~CO2_spline*N)

dge <- edgeR::DGEList(
  counts = normalized_counts,
  norm.factors = rep(1,ncol(normalized_counts)))

normalized_counts_kno3 <- normalized_counts[,str_detect(colnames(normalized_counts), "KNO3")]

############## or CO2 design on KNO3 only
# dge <- edgeR::DGEList(
#   counts = normalized_counts_kno3,
#   norm.factors = rep(1,ncol(normalized_counts_kno3)))
# conditions <- str_split_fixed(colnames(normalized_counts_kno3), '_', 2)[,1]
# CO2 <- as.numeric(str_split_fixed(conditions, '\\.', 2)[,1])
# CO2_spline <- ns(CO2, df = DF)
# design <- model.matrix(~CO2_spline)


############## tests
y <- edgeR::estimateDisp(dge, design)
fit <- edgeR::glmQLFit(y, design)


CO2_coeffs <- 2:(1+DF)
CO2_dea <- edgeR::glmQLFTest(fit, coef = CO2_coeffs)$table %>%
  filter(PValue < 0.005)

tfs <- intersect(regulators, rownames(CO2_dea))

annot <- gene_annotations$`Arabidopsis thaliana`
CO2_dea[,c("label", "description")] <- 
  annot[rownames(CO2_dea), c("label", "description")]


CO2_responsive_genes <- list("counts" = normalized_counts_kno3[rownames(CO2_dea),], 
                             "genes" = rownames(CO2_dea),
                             "tfs" = tfs,
                             "counts_all_N" = normalized_counts)

save(CO2_responsive_genes, file = "rdata/CO2_degs_expression.rdata")
