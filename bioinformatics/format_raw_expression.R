# format htseq count output for usage in DIANE

###### read htseq count directly
# data <- read.table(paste0('quantif_files/', 241, '.txt'), row.names ="V1")
# colnames(data) = "241"
# for(sample in 242:280){
#   df <- read.table(paste0('quantif_files/', sample, '.txt'), row.names ="V1")
#   colnames(df) <- as.character(sample)
#   data <- cbind.data.frame(data, df)
# }

library(tidyverse)
library(corrplot)

data <- read.csv(header = TRUE, row.names = "Gene", file = 'data/OC_08_raw_expression.txt', 
                check.names = F) 

data <- data[!str_detect(row.names(data), '__'),]

corrplot(cor(data), method = 'color', order = 'alphabet', col.lim = c(0,1), addrect = T)


annot <- read.csv('data/OC_08_annotation.csv', sep=';') %>%
  filter(Tube > 240) %>%
  mutate(label = paste(CO2,N, sep = '.')) %>%
  mutate(label = paste0(label, '_', rep(1:4, 10))) %>%
  select(Tube, label)

colnames(data) <- annot[match(colnames(data), annot$Tube), "label"]
data <- data[order(colnames(data))]

save(data, file = "rdata/expression_data.rdata")

# DIANE input, has to add manually "Gene" in first column
# write.table(data, file = "data/OC_08_CO2_gradient_KNO3_Mix_Ecotron.csv", sep = ',', quote = F)


dc2 <- data
colnames(dc2) <- colnames(data) %>%
  str_replace("400.KNO3_3", "tmp") %>%
  str_replace("900.Mix_3", "400.KNO3_3") %>%
  str_replace("tmp", "900.Mix_3") %>%
  str_replace("400.KNO3_4", "tmp") %>%
  str_replace("900.Mix_4", "400.KNO3_4") %>%
  str_replace("tmp", "900.Mix_4") %>%
  str_replace("400.Mix_1", "tmp") %>%
  str_replace("775.Mix_3", "400.Mix_1") %>%
  str_replace("tmp", "775.Mix_3") %>%
  str_replace("900.KNO3_4", "tmp") %>%
  str_replace("775.Mix_4", "900.KNO3_4") %>%
  str_replace("tmp", "775.Mix_4")

colnames(dc2) <- colnames(dc2) %>%
  str_replace("775.KNO3_2", "tmp") %>%
  str_replace("400.Mix_1", "775.KNO3_2") %>%
  str_replace("tmp", "400.Mix_1")%>%
  str_replace("775.KNO3_4", "tmp") %>%
  str_replace("400.Mix_2", "775.KNO3_4") %>%
  str_replace("tmp", "400.Mix_2")%>%
  str_replace("900.KNO3_4", "tmp") %>%
  str_replace("400.Mix_4", "900.KNO3_4") %>%
  str_replace("tmp", "400.Mix_4")%>%
  str_replace("775.KNO3_4", "tmp") %>%
  str_replace("400.Mix_3", "775.KNO3_4") %>%
  str_replace("tmp", "400.Mix_3")


corrplot(cor(dc2), method = 'color', order = 'alphabet', col.lim = c(0,1), addrect = T)


corrplot(cor(dc2), method = 'color', order = 'hclust', col.lim = c(0,1), addrect = T)

# corrplot(cor(data), method = 'color', order = 'alphabet', col.lim = c(0,1), addrect = T)
# 
# corrplot(cor(data), method = 'color', order = 'hclust', col.lim = c(0,1), addrect = T)
# 
# corrplot(cor(data), method = 'color', order = "AOE", col.lim = c(0,1), addrect = T)
# 
# corrplot(cor(data), method = 'color', order = 'FPC')




write.table(dc2, file = "data/OC_08_CO2_gradient_KNO3_Mix_Ecotron_corrected_for_mismatch.csv", sep = ',', quote = F)

data <- dc2
save(data, file = "rdata/expression_data_no_mismatch.rdata")


corplot(data)

library(DIANE)
data("abiotic_stresses")

corrplot(cor(abiotic_stresses$normalized_counts), method = 'color', order = 'alphabet')

pca <- ade4::dudi.pca(df = t(log(data+1)), scannf = FALSE, nf = 3, center = T, scale = T)
pca$li$cond <- str_split_fixed(rownames(pca$li), '_', 2)[,1]


DIANE::draw_PCA(data)

ggplot(pca$li, aes(x=Axis1, y=Axis2, color=cond,label = rownames(pca$li)))+ geom_point() + geom_text()


pca$li$kmeans <- as.factor(kmeans$cluster[rownames(pca$li)])
ggplot(pca$li, aes(x=Axis1, y=Axis2, color=kmeans,label = rownames(pca$li)))+ geom_point() + geom_text()

kmeans <- kmeans(t(log(data+1)), 20)
kmeans$cluster
