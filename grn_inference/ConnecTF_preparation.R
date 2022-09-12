########### how many tfs were studied in ConnceTF?

load("rdata/pwm_occurrences_CO2_response.rdata")
load(file = "rdata/CO2_degs_expression.rdata")

genes <- CO2_responsive_genes$genes
tfs <- CO2_responsive_genes$tfs
counts <- CO2_responsive_genes$counts


library(AraNetBench)
validated_edges <- validated_edges[validated_edges$from %in% tfs & validated_edges$to %in% genes,]
save(validated_edges, file = "rdata/connectf_CO2_genes.rdata")


intersect(unique(validated_edges$from), tfs)
target <- validated_edges[validated_edges$type == "TARGET",]
chip <- validated_edges[validated_edges$type == "CHIPSeq",]
dap <- validated_edges[validated_edges$type == "DAPSeq",]


DIANE::draw_venn(list("target" = intersect(unique(target$from), tfs), 
                      "chip" = intersect(unique(chip$from), tfs),
                      "dap" = intersect(unique(dap$from), tfs)))
