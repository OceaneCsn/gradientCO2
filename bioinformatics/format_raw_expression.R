# format htseq count output for usage in DIANE and downstream scripts

# tests if there is a problem with my merging of quantif files, but no (same as python script)

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

# making the link between tube IDs and conditions
annot <- read.csv('data/OC_08_annotation.csv', sep=';') %>%
  filter(Tube > 240) %>%
  mutate(label = paste(CO2,N, sep = '.')) %>%
  mutate(label = paste0(label, '_', rep(1:4, 10))) %>%
  select(Tube, label)
colnames(data) <- annot[match(colnames(data), annot$Tube), "label"]
data <- data[order(colnames(data))]

# drawing corplot
corrplot(cor(data), method = 'color', order = 'alphabet', col.lim = c(0,1), addrect = T)

save(data, file = "rdata/expression_data.rdata")

# DIANE input, has to add manually "Gene" in first column
# write.table(data, file = "data/OC_08_CO2_gradient_KNO3_Mix_Ecotron.csv", sep = ',', quote = F)

# obviously, there are sample mix-ups, so I manually rearranged the annotation :(

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
  str_replace("tmp", "400.Mix_3")%>%
  str_replace("775.KNO3_2", "tmp") %>%
  str_replace("900.KNO3_4", "775.KNO3_2") %>%
  str_replace("tmp", "900.KNO3_4")%>%
  str_replace("775.Mix_2", "tmp") %>%
  str_replace("900.Mix_2", "775.Mix_2") %>%
  str_replace("tmp", "900.Mix_2")


corrplot(cor(dc2), method = 'color', order = 'alphabet', col.lim = c(0,1), addrect = T) 

dc2 <- dc2[order(colnames(dc2))]

# reshape2::melt(cor(dc2)) %>%
#   ggplot(aes(x=Var1, y=Var2, fill = -value)) + geom_tile()+ scale_fill_distiller(palette = "Spectral") +
#   theme(axis.text.x = element_text(angle=90))

corrplot(cor(dc2), method = 'color', order = 'hclust', col.lim = c(0,1), addrect = T)

# export for DIANE (need to add Gene manually)
write.table(dc2, file = "data/OC_08_CO2_gradient_KNO3_Mix_Ecotron_corrected_for_mismatch.csv", sep = ',', quote = F)

# save rdata file
data <- dc2
save(data, file = "rdata/expression_data_no_mismatch.rdata")


# to see how a clean corplot should look like :
library(DIANE)
data("abiotic_stresses")
corrplot(cor(abiotic_stresses$normalized_counts), method = 'color', order = 'alphabet')