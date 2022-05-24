# format htseq count output for usage in DIANE

###### read htseq count directly
data <- read.table(paste0('quantif_files/', 241, '.txt'), row.names ="V1")
colnames(data) = "241"
for(sample in 242:280){
  df <- read.table(paste0('quantif_files/', sample, '.txt'), row.names ="V1")
  colnames(df) <- as.character(sample)
  data <- cbind.data.frame(data, df)
}

library(tidyverse)

#data <- read.csv(header = TRUE, row.names = "Gene", file = 'data/OC_08_raw_expression.txt', 
#                check.names = F) 

data <- data[!str_detect(row.names(data), '__'),]

annot <- read.csv('data/OC_08_annotation.csv', sep=';') %>%
  filter(Tube > 240) %>%
  mutate(label = paste(CO2,N, sep = '.')) %>%
  mutate(label = paste0(label, '_', rep(1:4, 10))) %>%
  select(Tube, label)

colnames(data) <- annot[match(colnames(data), annot$Tube), "label"]
data <- data[order(colnames(data))]

save(data, file = "rdata/expression_data.rdata")

# DIANE input, has to add manually "Gene" in first column
write.table(data, file = "data/OC_08_CO2_gradient_KNO3_Mix_Ecotron.csv", sep = ',', quote = F)


data <- data[rowSums(data) >40,]
pca <- ade4::dudi.pca(df = t(log(data+1)), scannf = FALSE, nf = 3, center = T, scale = T)
pca$li$cond <- str_split_fixed(rownames(pca$li), '_', 2)[,1]


DIANE::draw_PCA(data)

ggplot(pca$li, aes(x=Axis1, y=Axis2, color=cond,label = rownames(pca$li)))+ geom_point() + geom_text()


pca$li$kmeans <- as.factor(kmeans$cluster[rownames(pca$li)])
ggplot(pca$li, aes(x=Axis1, y=Axis2, color=kmeans,label = rownames(pca$li)))+ geom_point() + geom_text()

kmeans <- kmeans(t(log(data+1)), 20)
kmeans$cluster
