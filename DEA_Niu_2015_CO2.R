source("./Imports/ToolFunctions.R")
library(HTSFilter)
library(stringr)
library(ggVennDiagram)
library(plotly)
data <- read.csv("./Data/Niu_2015_CO2.txt", h=T, sep='\t')
genes <- data$Tracking_id
rownames(data) <- genes

cols <- grepl("CR", colnames(data))
data <- data[,cols]
conds = colnames(data)
groups <- str_split_fixed(conds, '_', 2)[,1]

data <- data[rowSums(data) > 40,]
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

niu <- OntologyProfile(genesNiu$gene_id)

# on ne retrouve que 6 gènes en commun avec mes 68 gènes détectés en CO2 elevé et controles
# 54 en commun avec ceux en faible nitrate sur les 1300 c'est tout nul

data <- read.csv("./Data/Jauregui_2015.txt", h=T, sep='\t')
data <- na.omit(data)
Jauregui_2015 <- toupper(data$gen.TAIR.10)
ggVennDiagram(list("Niu_2016" = genesNiu$gene_id, "Jauregui_2015" = Jauregui_2015, "cNF-CNF" = DEG_CO2_controle, "cnF-CnF" = DEGs[["cnF CnF"]])) + scale_fill_distiller(palette = "Spectral") + ggtitle("Genes in common for eCO2 response")
OntologyProfile(intersect( DEGs[["cnF CnF"]],DEGs[["cNF CNF"]] ))
ggVennDiagram(list("Niu_2016" = genesNiu$gene_id, "Jauregui_2015" = Jauregui_2015, "cNF-CNF" = DEG_CO2_controle, "cnF-CnF" = DEGs[["cnF CnF"]])) + ggtitle("Genes in common for eCO2 response")
ggVennDiagram(list("Niu_2016" = genesNiu$gene_id, "Jauregui_2015" = Jauregui_2015, "cnF-CnF" = DEGs[["cnF CnF"]])) + ggtitle("Genes in common for eCO2 response")
