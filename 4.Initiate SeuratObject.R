library(Seurat)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(pheatmap)
library(scales)
library(ggplot2)

sample_id <- c("B6","E4","E5")
sample_name <- c("OVA", "Saline 1","Saline 2")
condition <- c("OVA", "Saline","Saline")
results_path <- "result"
#create Seurat object
initialization<-function()
{
	sampleList <- list()
	for(i in 1:length(sample_id)){
		#read SoupX treated data
		cells10x <- Read10X(paste0("data/",sample_id[i], "soup"))
		#read doublet scores calculated by scrublet.py
		score_doublets<-t(read.csv(paste0("data/",sample_id[i],'soup/doublet_scores.txt'),header=F))
		colnames(score_doublets)<-colnames(cells10x)
		sampleList[[i]] <- CreateSeuratObject(counts = cells10x, assay = "RNA", 
			project = sample_name[i], min.cells = 3, min.features = 200)
		sampleList[[i]] <- AddMetaData(object = sampleList[[i]], metadata = t(score_doublets), col.name = "score.doublets")	
		sampleList[[i]]$condition <- condition[i]
		sampleList[[i]][["percent.mt"]] <- PercentageFeatureSet(sampleList[[i]], pattern = "^mt-")	
		sampleList[[i]][["percent.rb"]] <- PercentageFeatureSet(sampleList[[i]], pattern = "^Rp[sl]")	
	}
	return(sampleList)
}
sampleList<-initialization()
