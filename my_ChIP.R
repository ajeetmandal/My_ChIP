##source("https://bioconductor.org/biocLite.R")
##biocLite("ChIPseeker")

## loading packages
require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(clusterProfiler)

files <- getSampleFiles()
print(files)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
## resample = 500 by default, here use 100 to speed up the compilation of this vignette.
plotAvgProf(tagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=100, facet="row")

##Peak heatmaps
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-3000, 3000), verbose=FALSE)
head(peakAnnoList)
summary(peakAnnoList)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

#Functional profiles comparison
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
summary (genes)
names(genes) = sub("_", "\n", names(genes), options(max.print=1000000) )

##Overlap of peaks and annotated genes
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)

##The following chunk of codes work
compGO <- compareCluster(geneCluster   = genes, 
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05, 
                         pAdjustMethod = "BH")
plot(compGO, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

##Unfortunately, the following codes doesn't work
compKEGG <- compareCluster(geneCluster   = genes, 
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
plot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")


