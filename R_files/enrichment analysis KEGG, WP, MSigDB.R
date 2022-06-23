###Script for pathway enrichment analysis using the KEGG, WP and MSig database 


setwd("N:/")

library(readxl)
library(clusterProfiler)
library(enrichplot)
library(devEMF)
library(extrafont)
library(ggplot2)
library(org.Hs.eg.db)

##Load DEGs; Exemplary shown for the DEGs of cluster NSC1a
d <- read.csv("degs_NSC1a.csv") #Supplementary Table 5
gene <- d$gene

##Gene identifier changed from HGNC symbol to EntrezID
hs <- org.Hs.eg.db
genesEntrez <- select(hs, 
               keys = gene,
               columns = c("ENTREZID"),
               keytype = "SYMBOL")
geneEntrez <- genesEntrez$ENTREZID


####KEGG terms####
ekegg <- enrichKEGG(gene         = geneEntrez,
                    pvalueCutoff = 0.05)

#Export Table
write.table(ekegg,"ekegg.csv", row.names = FALSE, sep = ",")

#Generate Treeplot
treeplot(pairwise_termsim(ekegg), showCategory = 10,  nCluster = 3)

#Generate Emapplot
emapplot(pairwise_termsim(ekegg))


####Wikipathway terms####
ewp <- enrichWP(geneEntrez, organism = "Homo sapiens") 

#Export Table
write.table(ewp,"ewp.csv", row.names = FALSE, sep = ",")

#Generate Treaplot
treeplot(pairwise_termsim(ewp),showCategory = 10, nCluster = 4)

#Generate Emapplot
emapplot(pairwise_termsim(ewp))


####MsigDB terms####
#curated gene sets downloaded from https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C2 and stored in the folder extdata of the package clusterProfiler
#Import curated gene sets
gmtfile_c2_cp <- system.file("extdata", "c2.cp.v7.4.symbols.gmt", package="clusterProfiler")
c2_cp <- read.gmt(gmtfile_c2_cp) 

#Enrichment analysis
try(egmt_c2_cp  <- enricher(gene, TERM2GENE=c2_cp))

#Export Table
write.table(egmt_c2_cp,"egmt_c2_cp.csv", row.names = FALSE, sep = ",")

#Generate Emapplot
try(emapplot(pairwise_termsim(egmt_c2_cp)))
