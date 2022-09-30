#CA 
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
library(affy)  
library(GEOquery) 
if (!dir.exists("GEO")) {
  dir.create("GEO")
}
#get "GSE71348"
gse71348 <-getGEO('GSE71348')
gse71348
class(gse71348)
gse71348 <- gse71348[[1]]
class(gse71348)
library(tidyverse)
library(oligo)
varLabels(gse71348)
gse71348$supplementary_file
pd <- pData(gse71348)
pd
class(pd)
class(pd['cel_file'])
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
pd$cel_file
class(pd$cel_file)
library(pd.mouse430.2)
gse71348_celdata <- read.celfiles(paste0('cel/',pd$cel_file),phenoData=phenoData(gse71348))
image(gse71348_celdata[,1],main="Image of GSE71348")
pData(gse71348_celdata)[,c("cell line:ch1","time point:ch1","treatment:ch1")]

#RMA
hist(gse71348_celdata,main="CEL file densities before RMA")
boxplot(exprs(gse71348_celdata),outline=T,main="boxplot before RMA")
gse71348_rma <- rma(gse71348_celdata)
gse71348_rma
head(exprs(gse71348_celdata))
head(exprs(gse71348_rma))
hist(gse71348_rma,main="CEL file densities after RMA")
boxplot(exprs(gse71348_rma), outline=T,main="boxplot after RMA")

#OR Normalization + Summarisation
oligo_normalised <- normalize(gse71348_celdata,method='quantile',which='pm')
oligo_normalised
hist(oligo_normalised,lwd=2,xlab='log intensity', which='pm',
     main="CEL file densities after quantile normalisation")
oligo_summarised <- rma(oligo_normalised,background=FALSE,normalize=FALSE)
oligo_summarised
head(exprs(oligo_summarised))
hist(oligo_summarised)
boxplot(exprs(oligo_summarised), outline=F)

#Identifying differentially expresses genes using Linear models
library(limma)
gse71348_rma[['treatment:ch1']]
design <-  model.matrix(~ gse71348_rma[['treatment:ch1']])
colnames(design)[2:5] <- c("CHIR Day4","CHIR Day6","RSPO Day4","RSPO Day6")
design
fit <- lmFit(gse71348_rma,design)
fitted.ebayes <- eBayes(fit)
topTable(fitted.ebayes)
names(fitted.ebayes)
summary(decideTests(fitted.ebayes[,c("CHIR Day4","CHIR Day6","RSPO Day4","RSPO Day6")],lfc = 1))

library(dplyr)
pd <- pData(gse71348_rma)
pd
pd$cell_type <- pd[["cell line:ch1"]]
pd$treatment <- pd[['treatment:ch1']]
pd$treatment <- as.factor(pd$treatment)
factor(pd$treatment)
levels(pd$treatment) <- c("Control Day0","CHIRO Day4","CHIRO Day6","Rspo3 Day4","Rspo3 Day6")
pd$group <- as.factor(paste(pd$cell_type,pd$treatment))
pd$group
levels(pd$group) <- c("CK35D6.Chir","CK35D6.Rspo3","E14D4.Chir","E14D0.Control","E14D4.Rspo3")

design <- model.matrix(~0 + pd$group)
design
colnames(design) <- levels(pd$group)
design
contrasts_matrix <- makeContrasts(Msgn1_in_Chir=E14D4.Chir - E14D0.Control,
                                  Msgn1_in_Rspo3=E14D4.Rspo3 - E14D0.Control,
                                  Pax3_in_Chir=CK35D6.Rspo3 - E14D0.Control,
                                  Pax3_in_Chir=CK35D6.Rspo3 - E14D0.Control,
                                  interaction=(E14D4.Rspo3 - E14D4.Chir) - (CK35D6.Rspo3 - CK35D6.Chir),
                                  levels=design)
 contrasts_matrix
gse71348_fit <- lmFit(gse71348_rma,design)
gse71348_fit2 <- contrasts.fit(gse71348_fit,contrasts=contrasts_matrix)
gse71348_fit2 <- eBayes(gse71348_fit2)
summary(decideTests(gse71348_fit2,lfc=1))

#Annotation of genomic data in AnnotationDbi interface
ps <- rownames(topTable(fitted.ebayes))
ps
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('mouse4302.db')
library(mouse4302.db)
list(mouse4302.db)
columns(mouse4302.db)
keytypes(mouse4302.db)
head(keys(mouse4302.db,keytype="PROBEID"))
AnnotationDbi::select(mouse4302.db,ps,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

#Volcano plots
volcanoplot(fitted.ebayes,coef=2)
interesting_genes <- topTable(fitted.ebayes,coef=2,number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fitted.ebayes, coef=2, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red' )

#Heatmaps
eset_of_interest <- gse71348_rma[rownames(interesting_genes),]
eset_of_interest
heatmap(exprs(eset_of_interest))
library(RColorBrewer)
heatmap(exprs(eset_of_interest),
        labCol=gse71348_rma[['treatment']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))
