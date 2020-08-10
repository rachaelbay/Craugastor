library(topGO)

##Format for topGO
geneID2GO <- readMappings(file="data/annotation/GO_Orthogroups_07.21.20.txt")

###list of interest - this can be whatever you want it to be!
mods <- data.frame(readRDS("processed/WGCNA_modules.rds"))
allgenes <- rownames(mods) # This is your 'background' list of all possible genes. In this case, all genes that have been assigned to a module
myIG <- factor(as.integer(mods[,1]==8) # These is your list of genes of interest
names(myIG) <- allgenes


##Read into topGO object. The "BP" means "biological proces". Can also do "MF" and "CC"
GOdata <- new("topGOdata",ontology="BP",allGenes=myIG, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO) #

##Fisher test. There are a couple other possible tests too
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)

res <- GenTable(GOdata,classic=resultFisher,topNodes=500,numChar=100)
res$padj <- p.adjust(res$classic,method="fdr")
filt <- res[res$padj<0.01,]
filt
