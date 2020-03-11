library(ggplot2)
library(gridExtra)
suppressMessages(library(coseq, warn.conflicts = F, quietly = T))
library(plotly)

load("./Data/normalized.count_At.RData")
load("./Data/OntologyAllGenes.RData")
load("./Data/filteredData.RData")
load("./Data/NitrateGenes.RData")

# DEGs <- DEGs[["cnF CnF"]]
# genes <- unique(unlist(DEGs))
# cluster <- clustering(genes, data, nb_clusters = 7:14, norm = "TMM")

load("~/Documents/ShinyApp/ClusteringAnalysis/Clusterings/AmbientCO2_LowNitrateFe-ElevatedCO2_LowNitrateFe.RData")


plotProfile <- function(cluster, k="none"){
  # plot all the profiles or the profile of cluster k
  results <- cluster[[2]]
  clusters <- cluster[[1]]
  
  #plots <- plot(results, conds = groups, collapse_reps="average", graphs = c("ICL", "boxplots", "profiles", "probapost_barplots"))
  profiles <- data.frame(results@y_profiles)
  profiles$gene <- rownames(profiles)
  d <- melt(profiles)
  d$group <- str_split_fixed(d$variable, '_', 2)[,1]
  d$cluster <- clusters[match(d$gene, names(clusters))]
  if(k=="none"){
    ggplot(data = d, aes(x=group, y=value)) + geom_boxplot(alpha=0.7, aes( color = group), fill = "grey")  + geom_jitter(width = 0.1, alpha=0.01) +facet_wrap(~cluster) +
      ylim(0, 0.3) + scale_colour_discrete("Condition")
  }
  else{
    ggplot(data = d[d$cluster==k,], aes(x=group, y=value)) + geom_boxplot(alpha=0.7, aes(color = group), fill = "grey")  + geom_jitter(width = 0.1, alpha=0.02) +
      scale_colour_discrete("Condition")
  }
}


findNitrateGenes <- function(k, cluster){
  genesK <- names(cluster[[1]][cluster[[1]]==k])
  res <- data.frame(Gene = genesK)
  res$Wang_2004 = ifelse(res$Gene %in% nGenes$Wang_2004, 1, 0)
  res$Marchive_2013 = ifelse(res$Gene %in% nGenes$Marchive_2013_20min, 1, 0)
  res$Widiez_2011 = ifelse(res$Gene %in% nGenes$HN_induced_Widiez_2011, 1, 0)
  res[,c("Name", "Description")] <- ontologies[match(res$Gene, ontologies$ensembl_gene_id),c("external_gene_name", "description")]
  res$NitrateScore <- rowSums(res[,grepl("_", colnames(res))])
  res <- res[order(-res$NitrateScore),]
  return(res)
}


#findNitrateGenes(5,cluster)

####### generation des ontologies de genes nitrates
# nitrateGenes <- unique(c(as.vector(nGenes$Wang_2004), as.vector(nGenes$Marchive_2013_20min), as.vector(nGenes$HN_induced_Widiez_2011)))
# nitrateGenes <- nitrateGenes[grepl("AT", nitrateGenes)]
# ontologies <- OntologyProfile(ids = AGIToEntrez$ensembl_gene_id)
# save(ontologies, file = "./Data/OntologyAllGenes.RData")
