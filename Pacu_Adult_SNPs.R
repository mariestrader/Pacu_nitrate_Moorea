library(vcfR)
library(adegenet)
library(pcadapt)
library(ggplot2)
library(pegas)
library (ape)
library(hierfstat)
library(pheatmap)
library(wesanderson)
library(dartR)

#loading all samples together, and all species separate
setwd("/Users/mes0192/Dropbox/Moorea/RNA/writing/scriptsforsubmission")
setwd("/Users/mariestrader/Dropbox/Moorea/RNA/writing/scriptsforsubmission")

pacu<-read.vcfR('Pacu_snps_filtered_adults.recode.vcf') #Just the adult samples with their replicate

##Adegenet DAPC and Fstat
pacu<-read.vcfR('Pacu_snps_filtered_adults.recode.vcf')
pacu_gl=vcfR2genlight(pacu) #converts the matrix to Genlight format for adegenet
pacu_gl@ind.names<- sub("_[^_]+$", "", pacu_gl@ind.names)
pacu_gl@ind.names<- sub("_[^_]+$", "", pacu_gl@ind.names)
pacu_gl@ind.names<- sub("_[^_]+$", "", pacu_gl@ind.names)

pop(pacu_gl) <- c(rep("A",2),rep("C",1),rep("E",1),rep("B",2),rep("A",3),rep("D",1),rep("A",4),rep("F",2),rep("D",2),rep("A",2))

pc <- gl.pcoa(pacu_gl, nfactors=5)
#Performing a PCoA, individuals as entities, SNP loci as attributes
#Ordination yielded 5 informative dimensions from 19 original dimensions
#  PCoA Axis 1 explains 24.3 % of the total variance
#  PCoA Axis 1 and 2 combined explain 34.3 % of the total variance
#  PCoA Axis 1-3 combined explain 43.2 % of the total variance

gl.pcoa.plot(pc, pacu_gl, labels="pop", xaxis=1, yaxis=2) +theme_bw() 
gl.pcoa.scree(pc)

######################################## identity by state
library(WGCNA) # for coloring clonal groups
library(sparcl) # for ColorDendrogram

# reading list of bam files = order of samples in IBS matrix

bams=read.table("bams",header=F)[,1]
bams=sub("\\.fastq.*","",bams,perl=T)

bams <- sub("_[^_]+$", "", bams)
bams <- sub("_[^_]+$", "", bams)
bams <- sub("_[^_]+$", "", bams)

# reading IBS matrix based on SNPs with allele frequency >= 0.05:

ma = as.matrix(read.table("ibs05.ibsMat"))

dimnames(ma)=list(bams,bams)


# plotting hierarchical clustering tree

hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.6)




