

qcplot<-function(sampleObj){
	plot1<-VlnPlot(sampleObj, features = c("nFeature_RNA",  "percent.mt","score.doublets"), ncol = 3,pt.size = 0)
	return(plot1)
}
tmp_obj<-sampleList[[1]]
for(i in 2:length(sampleList)){
    tmp_obj<-merge(tmp_obj,sampleList[[i]])
}
# tmp_obj
# 21061 features across 48664 samples within 1 assay 

#FigS1
p1<-qcplot(tmp_obj)
ggsave(filename=paste0(results_path,"/FigS1_QC_before.pdf"),plot=p1,width=10,height=7)
#Quality control filtering
sampleFiltered <- lapply(sampleList, FUN = function(x) {
	x$limit_nFeature<-median(x$nFeature_RNA)+3*mad(x$nFeature_RNA)
	x$limit_doublet_score<-median(x$score.doublets)+3*mad(x$score.doublets)
	x$limit_percent_mt<-median(x$percent.mt)+3*mad(x$percent.mt)
	x <- subset(x, subset =nFeature_RNA < limit_nFeature
		&score.doublets<limit_doublet_score & percent.mt < limit_percent_mt
		&`Hbb-bt`<=0 & `Hba-a1`<=0 & `Hba-a2`<=0)
})

tmp_obj<-sampleFiltered[[1]]
for(i in 2:length(sampleFiltered)){
    tmp_obj<-merge(tmp_obj,sampleFiltered[[i]])
}
#FigS2
p1<-qcplot(tmp_obj)
ggsave(filename=paste0(results_path,"/FigS2_QC_after.pdf"),plot=p1,width=10,height=7)
