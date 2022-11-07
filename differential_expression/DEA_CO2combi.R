library(stringr)
library(DIANE)
library(tictoc)
library(ggplot2)

data <- read.csv("data/raw_expression.csv", h = T, row.names = "Gene")
annotate_condition <- function(condition) {
  co2 <- substr(condition, 1, 1)
  nitrate <- substr(condition, 2, 2)
  iron <- substr(condition, 3, 3)
  if (co2 == 'c')
    res <- "ACO2"
  else
    res <- "ECO2"
  if (nitrate == 'n')
    res <- paste(res, "Low-Nitrate", sep = '.')
  else
    res <- paste(res, "High-Nitrate", sep = '.')
  if (iron == 'f')
    res <- paste(res, "Iron-starvation", sep = '.')
  else
    res <- paste(res, "Iron-supply", sep = '.')
  return(res)
}

colnames(data) <-
  paste0(
    sapply(str_split_fixed(colnames(data), '_', 2)[, 1],
           annotate_condition),
    '_',
    str_split_fixed(colnames(data), '_', 2)[, 2]
  )

########## normalization with tmm and filtering < 10 avg counts

tcc_object <-
  normalize(data, str_split_fixed(colnames(data), '_', 2)[, 1], iteration = FALSE)
threshold = 10 * ncol(data)
tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
normalized_counts <- TCC::getNormalizedData(tcc_object)
normalized_counts <-
  normalized_counts[, order(colnames(normalized_counts))]


###################################### DEA

fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions =
                                   str_split_fixed(colnames(data), '_', 2)[, 1])

topTags <-
  DIANE::estimateDEGs(
    fit,
    reference = "ACO2.Low-Nitrate.Iron-supply",
    perturbation = "ECO2.Low-Nitrate.Iron-supply",
    p.value = 0.05,
    lfc = 0
  )
genes <- topTags$table$genes


##################################### Clustering of CO2*N conditions


#removing iron starvation
normalized_counts <-
  normalized_counts[, !str_detect(colnames(normalized_counts), 'starvation')]


###################### Gene grouping
# list of regulators, same as in the integration manuscript, 
# touttaken from ATFDB et plnTFDB
load('rdata/regulators.rdata')

tfs <- intersect(genes, regulators)

DIANE::draw_PCA(data = normalized_counts[genes,])


CO2_responsive_genes <- list("counts" = normalized_counts, 
                             "genes" = genes,
                             "tfs" = tfs)

# r <- group_regressors(genes = genes, normalized.count = 
#                         normalized_counts[genes,], regressors = regressors, 
#                       corr_thr = 0.95)
# 
# visNetwork::visNetwork(r$correlated_regressors_graph$nodes,
#                        r$correlated_regressors_graph$edges)
# 
# data <- r$counts
# grouped_genes <- r$grouped_genes
# grouped_regressors <- r$grouped_regressors
save(CO2_responsive_genes, file = "rdata/CO2combi_degs_expression.rdata")
