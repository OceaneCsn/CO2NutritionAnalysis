# Those generic functions are made to be used for Differential Expression Analysis
# So scripts for the analysis on a particular dataset are lighter and less redundant

suppressMessages(library(biomaRt))
suppressMessages(library(ggplot2))
suppressMessages(library(org.At.tair.db, warn.conflicts = F, quietly = T))
suppressMessages(library(enrichplot, warn.conflicts = F, quietly = T))
suppressMessages(library(coseq, warn.conflicts = F, quietly = T))
suppressMessages(library(clusterProfiler, warn.conflicts = F, quietly = T))

mart = useMart(biomart="plants_mart",host="plants.ensembl.org", dataset = "athaliana_eg_gene")



########################################################## Ontology

OntologyProfile <- function(ids, specie="At", plot = T){
  #Plot ontology enrichment stats of a given a set of entrezgene IDs
  # only for Arabidopsis
  if(specie == "At"){
    results <- getBM( filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "description", "external_gene_name", "entrezgene_id"),
                      values = ids, mart = mart)
    results <- results[!rownames(results) %in% which(duplicated(results$ensembl_gene_id)), ]
    kable(results)
    return(results)
  }
}

OntologyEnrich <- function(ids, universe, plot = T){
  # ids and universe must be entrez gene ids
  ego <- enrichGO(gene = ids,
                  OrgDb = org.At.tair.db,
                  ont = "BP",
                  universe = universe,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable = TRUE)
  
  # Elimine les redondances, en fusionnant les GO terms dont la similarite excede le cutoff
  simpOnt <- clusterProfiler::simplify(ego, cutoff=0.65, by="p.adjust", select_fun=min)
  result <- simpOnt@result
  if(plot){
    print(barplot(simpOnt, showCategory = 40, font.size = 10))
    print(emapplot(simpOnt, layout = "kk"))
  }
  return(simpOnt)
}


compareOnt <- function(idsList, universe){
  ckreg <- compareCluster(geneCluster = idsList, fun = "enrichGO", OrgDb = org.At.tair.db, ont = "BP", pAdjustMethod = "BH", 
                          pvalueCutoff = 0.01, qvalueCutoff = 0.05, universe = universe)
  ckreg@compareClusterResult
  simCk <- clusterProfiler::simplify(ckreg, cutoff=0.65, by="p.adjust", select_fun=min)
  resCk <- simCk@compareClusterResult
  print(dotplot(simCk, x = ~Cluster, showCategory = 15, font.size = 10))
  return(resCk)
}

########################################################## Format

writeGenes <- function(export, labels, ont, DEresult){
  #write.table(x = ont, file = export, quote = F, sep = '\t', row.names = F)
  ont$lfc <- DEresult[match(ont$ensembl_gene_id,DEresult$gene_id),]$m.value
  ont$adj.p.value <- DEresult[match(ont$ensembl_gene_id,DEresult$gene_id),]$q.value
  ont <- ont[order(ont$adj.p.value),]
  write.table(x = ont, file = paste0(export, translateToOSX(labels), '.txt'), quote = F, sep = '\t', row.names = F)
  return(ont)
}

translateToOSX <- function(labels){
  # Parce que ce concon de system de fichier n'est pas case sensitive
  RES = ""
  for (text in labels){
    res = ""
    if(grepl("c", text)){res = paste0(res, "AmbientCO2_")}
    else{res = paste0(res, "ElevatedCO2_")}
    if(grepl("N", text)){res = paste0(res, "HightNitrate_")}
    else{res = paste0(res, "LowNitrate")}
    if(grepl("f", text)){res = paste0(res, "FeStarvation")}
    else{res = paste0(res, "Fe")}
    RES = paste0(RES, res, '-')
  }
  return(substr(RES, 1, nchar(RES)-1))
}

######################## Poisson mixture model for gene clustering on expression

clustering <- function(DEgenes, data, nb_clusters = 2:12){
  conds = colnames(data)
  groups <- str_split_fixed(conds, '_', 2)[,1]
  dataC <- data[DEgenes,conds]
  run_pois <- coseq(dataC, conds=groups, K=nb_clusters, model="Poisson", iter = 5, transformation = "none")
  print(coseq::plot(run_pois, conds = groups, collapse_reps="average", graphs = c("ICL", "boxplots", "profiles", "probapost_barplots")))
  print(summary(run_pois))
  clusters_per_genes <- coseq::clusters(run_pois)
  return(clusters_per_genes)
}