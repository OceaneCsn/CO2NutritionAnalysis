# Those generic functions are made to be used for Differential Expression Analysis
# So scripts for the analysis on a particular dataset are lighter and less redundant

suppressMessages(library(biomaRt))
suppressMessages(library(ggplot2))
suppressMessages(library(org.At.tair.db, warn.conflicts = F, quietly = T))
suppressMessages(library(enrichplot, warn.conflicts = F, quietly = T))
suppressMessages(library(coseq, warn.conflicts = F, quietly = T))
suppressMessages(library(clusterProfiler, warn.conflicts = F, quietly = T))
suppressMessages(library(TCC, warn.conflicts = F, quietly = T))


mart = useMart(biomart="plants_mart",host="plants.ensembl.org", dataset = "athaliana_eg_gene")

########################################################## Ontology

OntologyProfile <- function(ids, specie="At", plot = T){
  #Plot ontology enrichment stats of a given a set of ensembl_gene_ids
  # only for Arabidopsis
  if(specie == "At"){
    results <- getBM( filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "description", "external_gene_name", "entrezgene_id"),
                      values = ids, mart = mart)
    results <- results[!rownames(results) %in% which(duplicated(results$ensembl_gene_id)), ]
    return(results)
  }
}

OntologyEnrich <- function(ids, universe, plot = T, simCutoff = 0.8){
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
  simpOnt <- clusterProfiler::simplify(ego, cutoff=simCutoff, by="p.adjust", select_fun=min)
  result <- simpOnt@result
  if(plot){
    print(barplot(simpOnt, showCategory = 40, font.size = 10))
    print(emapplot(simpOnt, layout = "kk"))
  }
  return(simpOnt)
}


compareOnt <- function(idsList, universe, simCutoff = 0.8){
  ckreg <- compareCluster(geneCluster = idsList, fun = "enrichGO", OrgDb = org.At.tair.db, ont = "BP", pAdjustMethod = "BH", 
                          pvalueCutoff = 0.01, qvalueCutoff = 0.05, universe = universe)
  ckreg@compareClusterResult
  simCk <- clusterProfiler::simplify(ckreg, cutoff=simCutoff, by="p.adjust", select_fun=min)
  resCk <- simCk@compareClusterResult
  print(dotplot(simCk, x = ~Cluster, showCategory = 25, font.size = 7))
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

clustering <- function(DEgenes, data, nb_clusters = 2:12, norm = "none"){
  conds = colnames(data)
  groups <- str_split_fixed(conds, '_', 2)[,1]
  dataC <- data[DEgenes,conds]
  run_pois <- coseq(dataC, conds=groups, K=nb_clusters, model="Poisson", iter = 5, transformation = "none", normFactors =norm, parallel = TRUE)
  print(coseq::plot(run_pois, conds = groups, collapse_reps="average", graphs = c("ICL", "boxplots", "profiles", "probapost_barplots")))
  print(summary(run_pois))
  clusters_per_genes <- coseq::clusters(run_pois)
  return(list(clusters_per_genes, run_pois))
}


writeExpression <- function(comp){
  filename = paste0(path,translateToOSX(comp), ".txt")
  At <- read.csv(filename, h=T, sep = "\t")
  # mean expression for the three replicates
  for(c in comp){
    At[,c] <- rowMeans(normalized.count[match(At$ensembl_gene_id, rownames(normalized.count)),grepl(c, colnames(normalized.count))])
  }
  At$MeanNormalizedExpression <- rowMeans(At[,comp])
  write.table(x = At, file = filename, quote = F, sep = '\t', row.names = F)
}


dualDE <- function(data, labels, pval=0.01, method="edger", lfc_filter = 0, plot=T){
  # selecting the right labels for pairwise comparison
  headers <- c(colnames(data)[(grepl(labels[1], colnames(data)))] , colnames(data)[grepl(labels[2], colnames(data))])
  data <- data[,headers]
  group <- str_split_fixed(colnames(data), "_", 2)[,1]
  group <- factor(group)
  group <- relevel(group, labels[1])
  # tcc object
  tcc <- new("TCC", data, group)
  print(model.matrix(~group))
  print(colnames(data))
  #2 steps normalisation
  tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 1, FDR = pval, floorPDEG = 0.05)
  print(tcc$norm.factors)
  tcc$DEGES$execution.time
  s <- sample(rownames(tcc$count), size = 200)
  normalized.count <- getNormalizedData(tcc)
  if(plot){
    heatmap(as.matrix(tcc$count[s,]), main = "Before normalisation")
    heatmap(as.matrix(normalized.count[s,]), main = "After normalisation")
  }
  #DEtest
  tcc <- estimateDE(tcc, test.method = method, FDR = pval, design = model.matrix(~group))
  result <- getResult(tcc, sort = TRUE)
  
  DEgenes <- subset(result,estimatedDEG==1 & abs(m.value) > lfc_filter)
  DEgenes$upreg = ifelse(DEgenes$m.value > 0, 1, 0)
  print(paste(dim(DEgenes)[1], " genes DE"))
  head(result)
  if(plot){
    print(ggplot(data = result, aes(a.value, m.value, color=factor(estimatedDEG))) +
            scale_fill_discrete("Set2") + geom_point(alpha=0.7) + ggtitle(paste0("M.A Plot : ", labels[2], " vs ", labels[1])) + xlab("Average expression") + ylab("Log Fold Change")+ theme(
              plot.title = element_text(size = 20, face="bold")) + labs(color = "Is DE"))
    
    print(ggplot(data = result, aes(m.value, -log10(q.value), color=factor(estimatedDEG))) +
            scale_fill_discrete("Set2") + geom_point(alpha=0.7) + ggtitle(paste0("Vulcano Plot : ", labels[2], " vs ", labels[1])) +
            xlab("Log Fold Change") + ylab("-Log(adj.pvalue)") + theme(
              plot.title = element_text(size = 20, face="bold")) + labs(color = "Is DE"))
    plotMDS(normalized.count, main="Multidimensional scaling plot of distances between gene expression profiles")
    heatmap(normalized.count[DEgenes$gene_id,])
  }
  return(DEgenes)
}
