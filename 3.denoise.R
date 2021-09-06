#correct the batch effect using SoupX
sample_id <- c("B6","E4","E5")
soup_sample<-function()
{
	library(Seurat)
	library(SoupX)
		sampleList <- list()
	for(i in 1:length(sample_id)){
		cells10x<-load10X(paste0("data/",sample_id[i]))
		cells10x<-autoEstCont(cells10x)
		cells10x<-adjustCounts(cells10x)
		#then save 10x matrix data or run obj<-CreateSeuratObject(counts = cells10x)
	}
}
soup_sample()

