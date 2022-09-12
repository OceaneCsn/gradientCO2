########## PWM occurrences in N response dataset
load(file = "rdata/CO2_degs_expression.rdata")
load('rdata/pwm_prom_jaspar_dap.rdata')

genes <- CO2_responsive_genes$genes
tfs <- CO2_responsive_genes$tfs
counts <- CO2_responsive_genes$counts

pwm_prom_ <- pwm_prom[pwm_prom$target %in% genes | 
                        pwm_prom$TF %in% tfs, ]
pwm_prom_ <- pwm_prom_[pwm_prom_$pval < 1e-4,]
known_tfs <- unique(pwm_prom_$TF)

intersect(known_tfs, tfs)
# construction de la matrice de présence du score
pwm_occurrence <- matrix(NA, nrow = length(genes),
                         ncol = length(tfs),
                         dimnames = list(genes, tfs))

for(tf in tfs){
  print(tf)
  if(length(intersect(tf, known_tfs)) >=1){
    pwm_prom_tf <- pwm_prom_[pwm_prom_$TF %in% tf,]
    for(gene in genes){
      pwm_occurrence[gene, tf] <-
        nrow(pwm_prom_tf[pwm_prom_tf$target %in% gene,])/(length(gene)*length(tf))
    }
  }
}

save(pwm_occurrence, file = "rdata/pwm_occurrences_CO2_response.rdata")