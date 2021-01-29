

setwd("/Users/mes0192/Dropbox/Moorea/RNA/writing/scriptsforsubmission")
setwd("/Users/mariestrader/Dropbox/Moorea/RNA/writing/scriptsforsubmission")

library('DESeq2')
library('arrayQualityMetrics')
library('vegan')
library('rgl')
library('ape')
library('pheatmap')
library('dplyr')
library('adegenet')
library('VennDiagram')
library('stringr')
library(tidyr)
library(readr)
library(Rmisc)
library(ggplot2)
library(wesanderson)

############################ upload gene counts file
gcountsALL=read.table("enrichmentCounts_pdamtranscriptome_host.txt", header=T)

#clean sample names
colnames(gcountsALL)<-sub("X","", colnames(gcountsALL))
colnames(gcountsALL) <- sub("_[^_]+$", "", colnames(gcountsALL))
colnames(gcountsALL) <- sub("_[^_]+$", "", colnames(gcountsALL))

length(gcountsALL[,1]) #23681 genes

summary(colSums(gcountsALL))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 491600  623900  684700  696300  758900  899300 

############################ remove genes with low mean counts

mns = apply(gcountsALL, 1, mean)
gcounts=gcountsALL[mns>10,] 
table(mns > 10)
#FALSE  TRUE 
#15361  8322 
dim(gcounts) 
#8322   34

############################ set up metadata from sample IDs

colData <- data.frame(row.names= colnames(gcounts), sample= colnames(gcounts))
colData$stage=str_sub(colData$sample, -1,-1)
colData$treat = str_sub(colData$sample, -2, -2)
colData$colony = as.numeric(gsub("([0-9]+).*$", "\\1", colData$sample))
colData$colony=as.factor(colData$colony)
colData$clone=as.factor(c(rep("A",4),rep("C",1),rep("E",2),rep("B",2),rep("A",6),rep("D",1),rep("A",7),rep("F",3),rep("D",4),rep("A",4)))
colData$group <- factor(paste0(colData$treat,colData$clone))

############################ search for outliers. none removed. 

dds<-DESeqDataSetFromMatrix(countData=gcounts, colData=colData, design= ~ treat+stage+geno)
vsd=varianceStabilizingTransformation(dds, blind=T)
e=ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("treat"), force=T, outdir= "report_for_genes_sub")

############################ subset the datasets between stages

# l is for larvae, a is for adults
gcountsL=gcounts[,c(2,4,7,11,13,15,18,21,23,25,28,30,32,34)]
gcountsA=gcounts[,-c(2,4,7,11,13,15,18,21,23,25,28,30,32,34)]

colDataL=colData[colData$stage=="L",]
colDataA=colData[colData$stage=="A",]

############################ model interaction using group for adult samples

dds<-DESeqDataSetFromMatrix(gcountsA,
	colData = colDataA, 
	design = formula(~ group))

# create rlog file
rld=rlog(dds)
rldAG.df = assay(rld)
ddsA<-DESeq(dds, minReplicatesForReplace=Inf) 

#test for role of nitrate on clonal group A (colonies 10,9,6,5,3,4C)
resA1<-results(ddsA, contrast=c('group', 'NA', 'CA')) 
mcols(resA1,use.names=TRUE)
summary(resA1)
#out of 8322 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 9, 0.11% 
#LFC < 0 (down)   : 6, 0.072% 
#outliers [1]     : 0, 0% 
#low counts [2]   : 0, 0% 
table(resA1$padj < 0.05)
#FALSE  TRUE 
# 8321     1 
rownames(resA1)<-sub("isogroup","pdam_", rownames(resA1))
resA1=resA1[!is.na(resA1$pvalue),]
resA1=data.frame(cbind("gene"=row.names(resA1),"stat"=resA1$stat))
head(resA1)
write.csv(resA1,file="GO/resA1_stat.csv",quote=F, row.names=F)
### this file used for treatment specific GO enrichment using GOMWU, scripts found at https://github.com/z0on/GO_MWU

#test for role of nitrate on clonal group B (colony 2)
resA4<-results(ddsA, contrast=c('group', 'NB', 'CB'))
mcols(resA4,use.names=TRUE)
summary(resA4)
#out of 8322 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 20, 0.24%
#LFC < 0 (down)     : 14, 0.17%
#outliers [1]       : 0, 0%
#low counts [2]     : 1291, 16%
table(resA4$padj < 0.05)
#FALSE  TRUE 
#  7014    17  
rownames(resA4)<-sub("isogroup","pdam_", rownames(resA4))
resA4=resA4[!is.na(resA4$pvalue),]
resA4=data.frame(cbind("gene"=row.names(resA4),"stat"=resA4$stat))
head(resA4)
write.csv(resA4,file="GO/resA4_stat.csv",quote=F, row.names=F)
### this file used for treatment specific GO enrichment using GOMWU, scripts found at https://github.com/z0on/GO_MWU

#test for role of nitrate on clonal group D (Colonies 8 + 4N)
resA5<-results(ddsA, contrast=c('group', 'ND', 'CD')) 
mcols(resA5,use.names=TRUE)
summary(resA5)
#out of 8322 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 0, 0% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 0, 0% 
#low counts [2]   : 0, 0% 
rownames(resA5)<-sub("isogroup","pdam_", rownames(resA5))
resA5=resA5[!is.na(resA5$pvalue),]
resA5=data.frame(cbind("gene"=row.names(resA5),"stat"=resA5$stat))
head(resA5)
write.csv(resA5,file="GO/resA5_stat.csv",quote=F, row.names=F)
### this file used for treatment specific GO enrichment using GOMWU, scripts found at https://github.com/z0on/GO_MWU

#test for role of nitrate on clonal group F (Colony 7)
resA6<-results(ddsA, contrast=c('group', 'NF', 'CF')) 
mcols(resA6,use.names=TRUE)
summary(resA6)
#out of 8322 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 0, 0% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 0, 0% 
#low counts [2]   : 0, 0% 
rownames(resA6)<-sub("isogroup","pdam_", rownames(resA6))
resA6 = resA6[!is.na(resA6 $pvalue),]
resA6 =data.frame(cbind("gene"=row.names(resA6),"stat"= resA6 $stat))
head(resA6)
write.csv(resA6,file="GO/resA6_stat.csv",quote=F, row.names=F)
### this file used for treatment specific GO enrichment using GOMWU, scripts found at https://github.com/z0on/GO_MWU

################################## model for treatment + clone in adults

dds<-DESeqDataSetFromMatrix(gcountsA,
	colData = colDataA, 
	design = formula(~ treat+clone))

rld=rlog(dds)
rldA.df = assay(rld)

ddsA<-DESeq(dds, minReplicatesForReplace=Inf) 
resA<-results(ddsA, contrast=c('treat', 'N', 'C')) 
mcols(resA,use.names=TRUE)
summary(resA)
#out of 8322 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 2, 0.024% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 0, 0% 
#low counts [2]   : 0, 0% 
resA = resA[order(resA$padj),]

rownames(resA)<-sub("isogroup","pdam_", rownames(resA))
resA=resA[!is.na(resA$pvalue),]
resA=data.frame(cbind("gene"=row.names(resA),"stat"=resA$stat))
head(resA)
write.csv(resA,file="GO/GO_allAdults_nitrate_stat.csv",quote=F, row.names=F)
### this file used for treatment specific GO enrichment using GOMWU, scripts found at https://github.com/z0on/GO_MWU


################################## model for group in larvae

dds<-DESeqDataSetFromMatrix(gcountsL,
	colData = colDataL, 
	design = formula(~ group))

rldL=rlog(dds)
rldLG.df = assay(rldL)

ddsL<-DESeq(dds, minReplicatesForReplace=Inf) 

#test for role of nitrate on clonal group A (colonies 10,9,6,5,3,4C)
resL1<-results(ddsL, contrast=c('group', 'NA', 'CA')) 
mcols(resL1,use.names=TRUE)
summary(resL1)
#out of 8320 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 0, 0% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 0, 0% 
#low counts [2]   : 0, 0% 
rownames(resL1)<-sub("isogroup","pdam_", rownames(resL1))
resL1 = resL1[!is.na(resL1 $pvalue),]
resL1 =data.frame(cbind("gene"=row.names(resL1),"stat"= resL1 $stat))
head(resL1)
write.csv(resL1,file="GO/resL1_stat.csv",quote=F, row.names=F)

#test for role of nitrate on clonal group D (Colonies 8 + 4N)
resL5<-results(ddsL, contrast=c('group', 'ND', 'CD')) 
mcols(resL5,use.names=TRUE)
summary(resL5)
#out of 8320 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2, 0.024%
#LFC < 0 (down)     : 0, 0%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
rownames(resL5)<-sub("isogroup","pdam_", rownames(resL5))
resL5 = resL5[!is.na(resL5 $pvalue),]
resL5 =data.frame(cbind("gene"=row.names(resL5),"stat"= resL5 $stat))
head(resL5)
write.csv(resL5,file="GO/resL5_stat.csv",quote=F, row.names=F)


################################# model for treatment in larvae

dds<-DESeqDataSetFromMatrix(gcountsL,
	colData = colDataL, 
	design = formula(~ treat + clone))

rld=rlog(dds)
rldL.df = assay(rld)

ddsL<-DESeq(dds, minReplicatesForReplace=Inf) 
resL<-results(ddsL, contrast=c('treat', 'N', 'C')) 
mcols(resL,use.names=TRUE)
summary(resL)
#out of 8320 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 0, 0% 
#LFC < 0 (down)   : 0, 0% 
#outliers [1]     : 0, 0% 
#low counts [2]   : 0, 0% 


rownames(resL)<-sub("isogroup","pdam_", rownames(resL))
resL=resL[!is.na(resL$pvalue),]
resL=data.frame(cbind("gene"=row.names(resL),"stat"=resL$stat))
head(resL)
write.csv(resL,file="GO_allLarvae_nitrate_stat.csv",quote=F, row.names=F)

save(rldAG.df, rldLG.df, colDataA, colDataL, resA1, resA4, resA5, resA6, resL1, resL5, file="Pacu_enrichment_BM10_groupALL_061020.Rdata")
save(rldA.df,rldL.df,colDataL,colDataA,resA,resL, file="Pacu_enrichment_BM10_treatALL_061020.Rdata")

################################# make datafile for all comparisons with annotations - adults

ll=load("Pacu_enrichment_BM10_groupALL_061020.Rdata") 

valsA1 <- cbind(resA1$stat, resA1$pvalue, resA1$padj)
colnames(valsA1)=c("stat_resA1","pval_resA1", "padj_resA1")

valsA4 <- cbind(resA4$stat, resA4$pvalue, resA4$padj)
colnames(valsA4)=c("stat_resA4","pval_resA4", "padj_resA4")

valsA5 <- cbind(resA5$stat, resA5$pvalue, resA5$padj)
colnames(valsA5)=c("stat_resA5","pval_resA5", "padj_resA5")

valsA6 <- cbind(resA6$stat, resA6$pvalue, resA6$padj)
colnames(valsA6)=c("stat_resA6","pval_resA6", "padj_resA6")

vsdpvalsA <- cbind(rldAG.df,valsA1, valsA4, valsA5, valsA6)
dim(vsdpvalsA)#8322   32

################################# for larvae
valsL1 <- cbind(resL1$stat, resL1$pvalue, resL1$padj)
colnames(valsL1)=c("stat_resL1","pval_resL1", "padj_resL1")

valsL5 <- cbind(resL5$stat, resL5$pvalue, resL5$padj)
colnames(valsL5)=c("stat_resL5","pval_resL5", "padj_resL5")

vsdpvalsL <- cbind(rldLG.df,valsL1, valsL5)
dim(vsdpvalsL) #8322   20

################################# add annotations
annot=read.table("pdam_iso2gene.tab",header=T,sep="\t",quote=NULL,fill=T)
names(annot)=c("iso","gene")

genesA = row.names(vsdpvalsA)
genesL = row.names(vsdpvalsL)

genes2annotA = match(genesA,annot$iso)
genes2annotL = match(genesL,annot$iso)

vsdpvalsA=data.frame(GeneID=row.names(vsdpvalsA), annot[genes2annotA,], vsdpvalsA)
vsdpvalsL=data.frame(GeneID=row.names(vsdpvalsL), annot[genes2annotL,], vsdpvalsL)

#save rldpvals with annotations#
save(vsdpvalsA,vsdpvalsL,file="vsdpvalsannots_host.Rdata")

#------------------ heatmap of sig adult genes
ll=load("vsdpvalsannots_host.Rdata")
colnames(vsdpvalsA)<-sub("X","", colnames(vsdpvalsA))
colnames(vsdpvalsL)<-sub("X","", colnames(vsdpvalsL))

row.names(vsdpvalsA)=make.names(vsdpvalsA$gene, unique=T)
row.names(vsdpvalsA)<-sub("Similar.to.","", row.names(vsdpvalsA))

row.names(vsdpvalsL)=make.names(vsdpvalsL$gene, unique=T)
row.names(vsdpvalsL)<-sub("Similar.to.","", row.names(vsdpvalsL))

# subset rld and pvals table to only those significant for each comparison #
vsdpvalsA1 = vsdpvalsA[vsdpvalsA$padj_resA1<= 0.1,]
dim(vsdpvalsA1) #15 35
selA1=as.data.frame(vsdpvalsA1[-c(1:3, 6:9,13,18:21,24:35) ]) #remove the other genotypes

vsdpvalsA4 = vsdpvalsA[vsdpvalsA$padj_resA4<= 0.1 & !is.na(vsdpvalsA$padj_resA4),]
dim(vsdpvalsA4) #45 35
selA4=as.data.frame(vsdpvalsA4[-c(1:7,10:35) ]) #remove the other genotypes

vsdpvalsL5 = vsdpvalsL[vsdpvalsL$padj_resL5<= 0.1 & !is.na(vsdpvalsL$padj_resL5),]
dim(vsdpvalsL5) #14 23
selL5=as.data.frame(vsdpvalsL5[-c(1:13, 16:23) ])

### heatmap plotting 
means=apply(selL5,1,mean) # means of rows
explc=selL5-means


heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1)(100)
heat.colors = colorRampPalette(wes_palette("Zissou1", 30, type = "continuous"), bias=1)(100)
pheatmap(explc,color=heat.colors,cluster_cols=T,border_color=NA,clustering_distance_rows="correlation")
pheatmap(explc,color=heat.colors,cluster_cols=T,border_color=NA)



################################# principal coordinate calculation 
ll=load("Pacu_enrichment_BM10_groupALL_061020.Rdata") 

conditions=colDataA

grp=rep("#228B22",ncol(rldAG.df))
grp[grep("C",conditions$treat)]="#A9A9A9"


fitpc=capscale(dist(t(rldAG.df),method="manhattan")~1)
summary(eigenvals(fitpc)) #37.8,  10.55 for adults, 21.26 and 10.56 for larvae

plot(fitpc$CA$u,pch=16, col="white",main="Pocillopora Adults", xlab="MDS1  21.8%", ylab="MDS2 13.5%")
abline(h=0,v=0, col="grey")
points(fitpc$CA$u,pch=16, col=grp)
ordispider(fitpc$CA$u,conditions$treat,col="black", label=F)
ordihull(fitpc$CA$u,conditions$treat,draw="polygon",label=T, col=grp)

plot(fitpc$CA$u,pch=16, col="white",main="Pocillopora Adults", xlab="MDS1  21.8%", ylab="MDS2 13.5%")
abline(h=0,v=0, col="grey")
points(fitpc$CA$u,pch=16, col=grp)
ordispider(fitpc$CA$u,conditions$clone,col="black", label=F)
ordiellipse(fitpc$CA$u,conditions$clone,draw="polygon",label=T)

ad=adonis(t(rldAG.df)~treat+clone+colony,data=conditions,method="manhattan")
#adonis(formula = t(rldAG.df) ~ treat + clone + colony, data = conditions,      method = "manhattan") 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treat      1   2050296 2050296  1.4149 0.04785  0.102    
#clone      5  20838175 4167635  2.8760 0.48632  0.001 ***
#colony     6   9815989 1635998  1.1290 0.22909  0.241    
#Residuals  7  10143832 1449119         0.23674           
#Total     19  42848292                 1.00000               

labs=c("Treatment","Clone","Residuals")
cols=c("#F21A00","#3B9AB2","grey80")
#cols = append(gg_color_hue(3), 'grey')
labs2 = paste(labs, round(ad$aov.tab$R2[1:3]*100, digits=1))
pie(ad$aov.tab$R2[1:3],labels=labs2,col=cols,main="Pocillopora Adults")

ad=adonis(t(rldLG.df)~treat+clone,data=conditions,method="manhattan")
#Call:
#adonis(formula = t(rldLG.df) ~ treat + clone, data = conditions,      method = "manhattan") 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#treat      1    797779  797779 0.81268 0.05203  0.700   
#clone      3   5699716 1899905 1.93540 0.37174  0.003 **
#Residuals  9   8834965  981663         0.57623          
#Total     13  15332459                 1.00000          

labs=c("Treatment","Clone","Residuals")
cols=c("#F21A00","#3B9AB2","grey80")
#cols = append(gg_color_hue(3), 'grey')
labs2 = paste(labs, round(ad$aov.tab$R2[1:3]*100, digits=1))
pie(ad$aov.tab$R2[1:3],labels=labs2,col=cols,main="Pocillopora Larvae")




