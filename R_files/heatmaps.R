###Script to generate the heaptmaps depicted in Fig. 4 and Supplementary Fig. 5


setwd("N:/")

library(readxl)
library(scales)
library(RColorBrewer)
library(gplots)
library(devEMF)

##label top degs
l_NSC1a <- c("TTR", "PLEKHD1", "SCML4", "ADAMTS8", "SLC22A13",
             "CCDC134", "PSD4", "PRIMA1", "ZNF285", "SRCAP")
l_NSC2a <- c("CATSPERD", "PCDHB1", "STOML3", "LHX4", "CAVIN2",
             "SPATA4", "SMTNL2", "RYR1", "SYNGR3", "CYB561D1")
l_NSC1b <- c("COL19A1", "ARX", "GRIN2A", "MOV10L1","C11orf88",
             "RNASE6", "PANX2", "ELOVL3", "ZNF737", "ATP10A")
l_NSC2b <- c("FAM20A", "POLM", "PHLDB3", "TMPPE", "SERHL2",
             "HIST1H2BE", "DDTL", "RNF212", "KCNC1", "CLIP4")
l_NCSC <- c("ART5","MUSK", "ZNF283", "GABRB3", "CYTIP",
            "ARMC5", "GM2A", "ERCC5", "A1BG", "MRPS24")
l_Glial <- c("PIF1", "PCDH8", "VPREB3", "CENPA", "UBE2C",
             "ATG2A", "MYL4", "TUBGCP6", "ADRA2A", "PNPLA3")
l_Neurons <- c("TACR3", "CDH19", "AKAP3", "C21orf62", "VGLL2",
               "HERC3", "GRAMD4", "TBC1D10A", "SNAPC4", "FAM227B")
label <- c(l_NSC1a, l_NSC2a, l_NSC1b, l_NSC2b, l_NCSC, l_Glial, l_Neurons)

####extract values and create matrix
NSC1a <- read.csv("full-genelist_NSC1a.csv")
NSC2a <- read.csv("full-genelist_NSC2a.csv")
NSC1b <- read.csv("full-genelist_NSC1b.csv")
NSC2b <- read.csv("full-genelist_NSC2b.csv")
NCSC <- read.csv("full-genelist_NCSC.csv")
Glial <- read.csv("full-genelist_Glial_precursors.csv")
Neurons <- read.csv("full-genelist_Immature_neurons.csv")

a <- data.frame(gene = label)
NSC1a_sub <- NSC1a[NSC1a$gene %in% label, c(2,4)] 
NSC2a_sub <- NSC2a[NSC2a$gene %in% label,c(2,4)]
NSC1b_sub <- NSC1b[NSC1b$gene %in% label,c(2,4)] 
NSC2b_sub <- NSC2b[NSC2b$gene %in% label,c(2,4)] 
NCSC_sub <- NCSC[NCSC$gene %in% label,c(2,4)]
Glial_sub <- Glial[Glial$gene %in% label,c(2,4)]
Neurons_sub <- Neurons[Neurons$gene %in% label,c(2,4)]

colnames(NSC1a_sub) <- c("gene", "NSC1a")
colnames(NSC2a_sub) <- c("gene", "NSC2a")
colnames(NSC1b_sub) <- c("gene", "NSC1b")
colnames(NSC2b_sub) <- c("gene", "NSC2b")
colnames(NCSC_sub) <- c("gene", "NCSC")
colnames(Glial_sub) <- c("gene", "Glial")
colnames(Neurons_sub) <- c("gene", "Neurons")


c <- merge(merge(merge(merge(merge(merge(merge(
  a,
  NSC1a_sub , all = TRUE), 
NSC2a_sub, all = TRUE), 
NSC1b_sub, all = TRUE), 
NSC2b_sub, all = TRUE), 
NCSC_sub, all = TRUE), 
Glial_sub, all = TRUE), 
Neurons_sub, all = TRUE) 
c <- c[match(label, c$gene),]

row.names(c) <- c$gene
c <- c[,-1]

col <- brewer_pal(palette = "RdYlBu", direction = -1)(10)
#or 
#col <- bluered(100)

heatmap.2(as.matrix(c), scale = "row", col = col,
          trace = "none", density.info = "none", dendrogram = "col", Rowv = FALSE,
          hclustfun=function(x) hclust(x,method = "complete"))



##load PD and mt genes
pd <- read_xlsx("PD-genes_for_heatmap.xlsx")
mt <- read_xlsx("mt-genes_for_heatmap.xlsx")

label <- mt$gene
#or
#label <- pd$gene


heatmap.2(as.matrix(c), scale = "row", col = col,
          trace = "none", density.info = "none", dendrogram = "col", Rowv = FALSE,
          hclustfun=function(x) hclust(x,method = "complete"))
