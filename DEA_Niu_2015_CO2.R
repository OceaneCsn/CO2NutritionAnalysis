source("./Imports/ToolFunctions.R")
library(HTSFilter)
library(stringr)
library(ggVennDiagram)
data <- read.csv("./Data/Niu_2015_CO2.txt", h=T, sep='\t')
genes <- data$Tracking_id
rownames(data) <- genes

cols <- grepl("CR", colnames(data))
data <- data[,cols]
conds = colnames(data)
groups <- str_split_fixed(conds, '_', 2)[,1]


data <- data[data>10,]
data <- na.omit(data)

head(data)


heatmap(as.matrix(data[sample(rownames(data), 400),]) )


res <- dualDE(data, c("ACR", "ECR"), plot = T)
genesNiu <- res
save(genesNiu, file = "./Data/DEG_Nui_2015_CO2.RData")


load("Data/DEGsListsFiltered.RData")
DEG_CO2_controle <- DEGs[["cNF CNF"]]


ggVennDiagram(list("Niu_2015" = genesNiu$gene_id, "Nous" = DEG_CO2_controle ))
common <- OntologyProfile(ids = intersect(genesNiu$gene_id, DEG_CO2_controle))

ggVennDiagram(list("Niu_2015" = genesNiu$gene_id, "Nous" = DEGs[["cnF CnF"]] ))
common <- OntologyProfile(ids = intersect(genesNiu$gene_id, DEGs[["cnF CnF"]]))


# on ne retrouve que 6 gènes en commun avec mes 68 gènes détectés en CO2 elevé et controles
# 54 en commun avec ceux en faible nitrate sur les 1300 c'est tout nul