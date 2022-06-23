###Script to generate the Volcano plots depicted in Supplementary Fig. 4 and 11


setwd("N:/")

library(devEMF)
library(EnhancedVolcano)

##Vulcano plot for genes/DEGs of cluster NSC1a
###load genes
NSC1a <- read.csv("full-genelist_NSC1a.csv")

###list of genes that should be labeled
label <- c("SLC1A2", "LINGO2", "BBS5", "GET4",  "SRCAP", "SFRP4", "GLI1", "GLI2", "GLI3", "FOXA2", "GBP1", "SHOX2", "PTCH1", "PTCH2", "CFHR1", "GRB14", "PITX3")


###Vulcano plot
EnhancedVolcano(toptable = NSC1a,
                title = "",
                subtitle = "",
                lab = NSC1a$gene, 
                x = "log2fc", y = "qval", 
                ylim = c(0,15), xlim = c(-4,4),
                pCutoff = 0.05,
                FCcutoff = 0.26,
                pointSize = 2.0,
                gridlines.major = FALSE,
                colConnectors = "black",
                selectLab = label,
                labCol = "black",
                drawConnectors = TRUE
)
