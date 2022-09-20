library(igraph)
library(visNetwork)

rank_nodes <- function(graph, regulators, gene_info = NULL){
  data <- visNetwork::toVisNetworkData(graph)
  
  degree <- igraph::degree(graph)
  degree_out <- igraph::degree(graph, mode = "out")
  degree_in <- igraph::degree(graph, mode = "in")
  
  # degree computaton
  data$nodes$degree <- degree[match(data$nodes$id, names(degree))]
  data$nodes$degree_out <- degree_out[match(data$nodes$id, names(degree_out))]
  data$nodes$degree_in <- degree_in[match(data$nodes$id, names(degree_in))]
  
  # modules computation
  
  data$nodes$group <- ifelse(data$nodes$id %in% regulators, "Regulator", 
                             ifelse(grepl("mean_", data$nodes$id), 
                                    "Grouped Regulators", "Target Gene"))
  data$nodes$gene_type <- data$nodes$group

  # adding additional infos
  if(!is.null(gene_info)){
    data$nodes[,colnames(gene_info)] <- 
      gene_info[match(data$nodes$id, rownames(gene_info)), ]
    
  }
  return(data$nodes)
}
