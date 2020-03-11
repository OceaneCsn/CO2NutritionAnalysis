library(ggplot2)
library(gridExtra)

source("./Imports/ToolFunctions.R")
source("./Imports/VisuFunctions.R")

load("./Data/normalized.count_At.RData")
data <- read.table("./Data/AmbientCO2_LowNitrateFe-ElevatedCO2_LowNitrateFe.txt", h = T, stringsAsFactors = F)

genes <- data$ensembl_gene_id
heatmapPerso(normalized.count = normalized.count, conds = c("cNF", "cnF"), genes = genes)

activators <- data[data$Status == "Activator",]$ensembl_gene_id
repressors <- data[data$Status == "Repressor",]$ensembl_gene_id
par(mfrow = c(6,6))

act <- list()
for (gene in as.vector(activators)){
  print(getExpression(normalized.count = normalized.count, conds = c("cnF", "CnF"), gene))
}

conds = c("cnF", "CnF")
conds = grepl(conds[1], colnames(normalized.count)) | grepl(conds[2], colnames(normalized.count))
df <- normalized.count[c(activators, repressors), conds]
d<- melt(df, silent=T)
d$Role <- ifelse(d$Var1 %in% activators, "Activator", "Repressor")
d$group = str_split_fixed(d$Var2, "_", 2)[,1]
d$Name = data[match(d$Var1, data$ensembl_gene_id),"external_gene_name"]


d <- transform(d, Role=factor(Role,levels=c("Activator","Repressor")))

ggplot(data = d[d$Role=="Repressor",], aes(x=group, y=value, fill=group)) + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 2.5, alpha = 0.8) +
  scale_fill_manual(name = "Conditions", labels = c("Ambient CO2", "Elevated CO2"), values = c('#006666','#D35400')) + facet_wrap(~Name, nrow=3, scales="free") + ggtitle("Repressors of nitrate pathways") +
  theme(strip.text.x = element_text(size = 20), plot.title = element_text(size=25, face="bold"),
        legend.title = element_text(size = 25, face="bold"), legend.text = element_text(size=24),
        axis.text.y = element_text(size = 8, angle = 30), axis.text.x = element_text(size = 0),
        axis.title=element_text("",size=17)) + ylab("Normalized counts") + xlab("") 

#### See how many genes in common the DEGs have with the kown nitrate related genes

load("Data/DEGsListsFiltered.RData")
nGenes <- read.table("./Data/N_regulated_genes.txt", h=TRUE, sep = '\t')
nGenes$Wang_2004 <- toupper(nGenes$Wang_2004)
nGenes$HN_induced_Widiez_2011 <- toupper(nGenes$HN_induced_Widiez_2011)

library(ggVennDiagram)
ggVennDiagram(nGenes) + scale_fill_distiller(palette = "Spectral") + ggtitle("Common genes of nitrate pathways")
lengths(DEGs)
df <- data.frame(Comparison = names(DEGs), DEGs = lengths(DEGs))



getCommonGenes <- function(comp, littG){
  return(length(intersect(DEGs[[as.character(comp)]], nGenes[,littG]))/length(DEGs[[as.character(comp)]])*100)
}
getCommonGenes(comp = "cnF CnF", littG = "Wang_2004")

getCommonGenes(comp = "cnF CnF", littG = "Wang_2004")

df$Wang_2004 <- sapply(X = df$Comparison,getCommonGenes, littG ="Wang_2004")
df$Marchive_2013_20min <- sapply(X = df$Comparison,getCommonGenes, littG ="Marchive_2013_20min")
df$HN_induced_Widiez_2011 <- sapply(X = df$Comparison,getCommonGenes, littG ="HN_induced_Widiez_2011")

nGenes$HN_induced_Widiez_2011

library(reshape2)
df$Factor <- c(rep("CO2", 4),rep("Nitrate", 4),rep("Iron", 4))
d<- melt(df[,!grepl("DEGs", colnames(df))], silent=T)
ggplot(data <- d, aes(x=Comparison, y = value)) + geom_bar(aes(fill = variable), position = "dodge", stat="identity", alpha = 0.7) + facet_wrap(~Factor, nrow = 1) +theme(strip.text.x = element_text(size = 26), plot.title = element_text(size=22, face="bold"),
        legend.title = element_text(size = 25, face="bold"), legend.text = element_text(size=20),
        axis.text.y = element_text(size = 18, angle = 30), axis.text.x = element_text(size = 26, angle = 320, hjust = 0, colour = "grey50"),
        axis.title=element_text(size=17)) + ylab("Percentage of common genes") +
  ggtitle(paste("% of genes in each list for each DE analysis"))
