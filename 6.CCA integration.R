#integrate data
pre_process_int<-function(obj){
	sampleFiltered <- lapply(obj, FUN = function(x) {
		x<- NormalizeData(x, verbose = F)
		x<- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = F)
	})
	obj.anchors <- FindIntegrationAnchors(object.list = sampleFiltered, dims = 1:30,reduction = "cca")
	obj <- IntegrateData(anchorset = obj.anchors,features.to.integrate = rownames(tmp_obj) ,dims = 1:30)
	DefaultAssay(obj) <- "integrated"
	obj<-ScaleData(obj, verbose = T,features = rownames(tmp_obj),vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
}
#tmp_obj<-sampleFiltered[[1]]
#for(i in 2:length(sampleFiltered)){
#    tmp_obj<-merge(tmp_obj,sampleFiltered[[i]])
#}


asthma_obj<-pre_process_int(sampleFiltered)
DefaultAssay(asthma_obj)<-"integrated"
