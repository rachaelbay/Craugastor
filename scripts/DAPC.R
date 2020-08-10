library(adegenet)


norm <- readRDS("processed/Singles_TPM.rds")
meta <- readRDS("processed/Meta_ordered.rds")
groups <- paste(meta$Species,meta$Treatment,sep=".")

d <- dapc(t(norm),groups) # The choice of #PCs seems to matter quantitatively, but not qualitatively. I chose 10,3

##Plot DAPC
par(mfrow=c(2,1),mar=c(2,2,1,1))
scatter(d)
scatter(d,1,1)

