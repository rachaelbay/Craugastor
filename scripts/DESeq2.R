library(tximport)
library(dplyr)
library(ggplot2)
library(DESeq2)

##Read in metadata
meta <- read.csv("data/frog_metadata_RAB.csv")
meta$Treatment <- as.factor(meta$Treatment)

##Read in expression data from Salmon separately for CRCR and CRFI
cr.samples <- list.files("data/salmon/crcr")
cr.names <- sapply(cr.samples,function(x) unlist(strsplit(x,split="[.]"))[1])
cr.files <- file.path("data/salmon","crcr",cr.samples)
names(cr.files) <- cr.names
txi.crcr <- tximport(cr.files,type='salmon',txOut=TRUE)

fi.samples <- list.files("data/salmon/crfi")
fi.names <- sapply(fi.samples,function(x) unlist(strsplit(x,split="[.]"))[1])
fi.files <- file.path("data/salmon","crfi",fi.samples)
names(fi.files) <- fi.names
txi.crfi <- tximport(fi.files,type='salmon',txOut=TRUE)

##Read in orthologroup table
ortho.counts <- read.delim("data/Orthogroups/Orthogroups.GeneCount.tsv")
singles <- which(ortho.counts$CRCR_fixednames==1 & ortho.counts$CRFI_fixednames==1)
length(singles) # These are 1:1 orthologs
nozeros <- which(ortho.counts$CRCR_fixednames>0 & ortho.counts$CRFI_fixednames>0)
length(nozeros) # These are all groups with at least one in each species

##Read in orthogroups
orthos <- read.delim("data/Orthogroups/Orthogroups_fixednames.tsv") #Had to find/replace "__" with "::"
single.orthos <- orthos[singles,]

##combine one-to-one orthologues into one frame
singles.crcr <- txi.crcr$counts[match(single.orthos$CRCR_fixednames,rownames(txi.crcr$counts)),]
singles.crfi <- txi.crfi$counts[match(single.orthos$CRFI_fixednames,rownames(txi.crfi$counts)),]
all.singles <- cbind(singles.crcr,singles.crfi)
rownames(all.singles) <- single.orthos[,1]
all.singles.rounded <- apply(all.singles,2,round) #This is janky but DESeq needs integers
#saveRDS(all.singles.rounded,"Singles_counts.rds")

##DESeq all together
##DON'T USE THIS, THIS IS WRONG!
# ordermeta <- meta[match(colnames(all.singles),meta$Sample_ID),c("Sample_ID","Species","Treatment")]
# dds <- DESeqDataSetFromMatrix(countData=all.singles.rounded,
#                               colData=ordermeta,
#                               design=~Species+Treatment)
# dds <- DESeq(dds)
# saveRDS(ordermeta,"Meta_ordered.rds")
# vsd <- vst(dds,blind=FALSE)
# plotPCA(vsd,intgroup=c("Species","Treatment"))
# 
# resultsNames(dds)
# SpecLFC <- lfcShrink(dds, coef="Species_CRFI_vs_CRCR")
# HeatLFC <- lfcShrink(dds, coef="Treatment_30_vs_20")
#saveRDS(SpecLFC,"Species_DDSeq_LFC.rds")
#saveRDS(HeatLFC,"Heat_DDSeq_LFC.rds")




###############################
##I also analyzed each species separately (because DESeq can do this 'right')
##Do the two species respond to heat stress the same? Analyze each species separately
##Note: This is the only 'correct' way to use DESeq for this dataset!
cr.meta <- meta[match(colnames(txi.crcr$counts),meta$Sample_ID),c("Sample_ID","Species","Treatment")]
fi.meta <- meta[match(colnames(txi.crfi$counts),meta$Sample_ID),c("Sample_ID","Species","Treatment")]

cr.dds <- DESeqDataSetFromTximport(txi.crcr,cr.meta,design=~Treatment)
cr.dds <- DESeq(cr.dds)
cr.res <- results(cr.dds)
cr.single.res <- cr.res[match(single.orthos$CRCR_fixednames,rownames(cr.res)),]

fi.dds <- DESeqDataSetFromTximport(txi.crfi,fi.meta,design=~Treatment)
fi.dds <- DESeq(fi.dds)
fi.res <- results(fi.dds)
fi.single.res <- fi.res[match(single.orthos$CRFI_fixednames,rownames(fi.res)),]

length(which(cr.single.res$padj<0.01))
length(which(fi.single.res$padj<0.01))
length(which(cr.single.res$padj<0.01 & fi.single.res$padj<0.01)) # Lots of overlap here

###Save dataframe with both
singles.frame <- data.frame(CR.mean=cr.single.res$baseMean,
                            CR.lfc=cr.single.res$log2FoldChange,
                            CR.padj=cr.single.res$padj,
                            FI.mean=fi.single.res$baseMean,
                            FI.lfc=fi.single.res$log2FoldChange,
                            FI.padj=fi.single.res$padj,
                            row.names=single.orthos$Orthogroup)
saveRDS(singles.frame,"processed/DEres_separate_species.rds")



plot(cr.single.res$log2FoldChange[cr.single.res$padj<0.01 | fi.single.res$padj<0.01],
     fi.single.res$log2FoldChange[cr.single.res$padj<0.01 | fi.single.res$padj<0.01],xlab="CRCR heat stress response",ylab="CRFI heat stress response")
points(cr.single.res$log2FoldChange[cr.single.res$padj<0.01 & fi.single.res$padj<0.01],
       fi.single.res$log2FoldChange[cr.single.res$padj<0.01 & fi.single.res$padj<0.01], col="red")
abline(0,1,col="blue") # No clear pattern here, probably just noise


