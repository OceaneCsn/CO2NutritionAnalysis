library(coseq)
library(stringr)

load("./Data/filteredData.RData")
load("./Data/DEGsListsFilteredNoIronStarv.RData")
genes <- DEGs[["cnF CnF"]]


conds = colnames(data)
groups <- str_split_fixed(conds, '_', 2)[,1]
dataC <- data[genes,conds]
run_pois <- coseq(dataC, conds=groups, K=6:12, model="Poisson", iter = 5, transformation = "none", normFactors ="TMM", parallel = TRUE)
print(coseq::plot(run_pois, conds = groups, collapse_reps="average", graphs = c("ICL", "boxplots", "profiles", "probapost_barplots")))
print(summary(run_pois))
clusters_per_genes <- coseq::clusters(run_pois)
run_pois@metadata

Co2 <- str_split_fixed(groups, "", 3)[,1]
nitrate <- str_split_fixed(groups, "", 3)[,2]
fer <- str_split_fixed(groups, "", 3)[,3]
model.matrix(~ Co2 * nitrate * fer)
