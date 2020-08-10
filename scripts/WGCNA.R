library(WGCNA)

##Read in gene expression and metadata
ge <- readRDS("data/Singles_TPM.rds")
datExpr <- t(ge)
meta <- read.csv("data/frog_metadata_RAB.csv")
meta$Treatment <- as.factor(meta$Treatment)
meta$Species <- as.factor(as.numeric(as.factor(meta$Species)))
ordermeta <- meta[match(rownames(datExpr),meta$Sample_ID),]
datTraits <- ordermeta[,c("Species","Days","Treatment")]
rownames(datTraits) <- ordermeta$Sample_ID

##Filter for genes with very low expression
low.thresh=1
up.thresh=10000
means <- colMeans(datExpr)
keep <- which(means>low.thresh & means<up.thresh)
datExpr <- datExpr[,keep]

plotClusterTreeSamples(datExpr=datExpr) #No obvious outliers

## Soft Thresholding (I ended up chosing 10)
options(stringsAsFactors = FALSE)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


## Blockwise network construction
bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
                         power = 10, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         verbose = 3)

# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], labels2colors(bwnet$colors)[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

## Relating modules to traits
nGenes=ncol(datExpr)
nSamples=nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, bwnet$colors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");

dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1),mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#Species modules=0*,1*,4,8,11,12
#Heat modules=5,8,6
MEs <- readRDS("WGCNA_MEs.rds")
par(mfrow=c(3,5),mar=c(5,2,2,2))
library(beeswarm)
for (i in 1:ncol(MEs)) {
beeswarm(MEs[,i]~datTraits$Treatment+datTraits$Species,main=names(MEs)[i],las=2,xlab="")  
boxplot(MEs[,i]~datTraits$Treatment+datTraits$Species,add=T,col=NA,names=NA,axes=F)
}

###Save important info
#saveRDS(MEs,"WGCNA_MEs.rds")
#saveRDS(cbind(bwnet$colors),"WGCNA_modules.rds")
