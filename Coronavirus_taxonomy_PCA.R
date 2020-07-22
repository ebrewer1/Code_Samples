#Phylogram and Principle Components Analysis Code
#Code Template provided by Dr. T.F. Khang and Drew Kerkhoff. See references for citations

#library and setup
library(seqinr)
library(Biostrings)
library(ape)
library(textmineR)
library(phangorn)
library(ggplot2)
library(ggrepel)
setwd("/Users/Ned/BMMB_852/Genomics_Final/viral_genomes")

#read into phyDat, create taxonomy, create color scheme
corona.aln.phydat=read.phyDat("/Users/Ned/BMMB_852/Genomics_Final/clustalo_msa.maf", format="fasta", type = "DNA")
summary(corona.aln.phydat)[1:5,]
corona.names <- c("SARS-Cov-2", "SARS-Cov", "MERS-Cov", "HCov-HKU1", "HCov-229E", "HCov-OC43", "HCov-NL63")
names(corona.aln.phydat)=corona.names
names(corona.aln.phydat)
taxonomy <- data.frame(names(corona.aln.phydat), 
                       c("Beta Sarbecovirus", "Beta Sarbecovirus", "Beta Merbecovirus", "Beta Embecovirus", 
                         "Alpha Duvinacovirus","Beta Embecovirus", "Alpha Setracovirus"))
colnames(taxonomy) <- c("Virus","Family")
tipcolor <- c("red","blue","darkviolet","green4","black")[unclass(taxonomy$Family)]

#perform distance analysis, neighbor joining, and plot the phylogram
corona.phydat.dist=dist.ml(corona.aln.phydat)
corona.NJ=NJ(corona.phydat.dist)
plot.phylo(corona.NJ, x.lim = 15, y.lim = NULL, use.edge.length=FALSE, cex=1, tip.color = tipcolor)

#normalize phyDat data, perform principle components analysis, and plot pairs and final graph
M <- do.call(rbind, corona.aln.phydat)
M_norm <- t(apply(M, 1, function(k) k/sum(k)))
pca <- prcomp(M_norm)
summary(pca) 
pairs(pca$x[,1:5], pch=16, cex=2, col=tipcolor)
pc1_v_2 <- data.frame(pca$x[,1], pca$x[,2])
colnames(pc1_v_2) <- c("PC1", "PC2")
ggplot(data = pc1_v_2, aes(x = PC1, y = PC2, color = taxonomy$Family)) +
  geom_point(color = tipcolor, size = 3, aes(color = taxonomy$Family)) +
  geom_label_repel(aes(label = rownames(pc1_v_2))) +
  labs(x = "PC1 (39.3%)", y = "PC2 (27.7%)", title = "Figure 2: Normalized PCA of Human-infecting coronaviruses") +
  scale_color_manual(name = "Virus Family",
                     breaks = c("Alpha Duvinacovirus", "Alpha Setracovirus", "Beta Embecovirus", 
                                "Beta Merbecovirus", "Beta Sarbecovirus"),
                     values = c("Alpha Duvinacovirus" = "red", "Alpha Setracovirus" = "blue",
                                "Beta Embecovirus" = "darkviolet", "Beta Merbecovirus" = "green4", 
                                "Beta Sarbecovirus" = "black")) 
