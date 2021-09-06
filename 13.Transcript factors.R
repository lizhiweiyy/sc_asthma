
#1.IRIS3
init_data_TF<-function(obj){
	subcls<-obj
	count_raw <- data.frame(subcls@assays$RNA@counts)
	write.csv(count_raw, 'neu_count.csv')
	meta_data <- cbind(rownames(subcls@meta.data), subcls@meta.data[,'cellType', drop=F]) 
	write.csv(meta_data, 'neu_cellType.csv',row.names = F)
}
init_data_TF(asthma_neu)
#
#Then uploaded to IRIS3(https://bmbl.bmi.osumc.edu/iris3/submit.php) for online analysis

#2.MetaCore
markers <- FindAllMarkers(asthma_neu,only.pos = T)
write.csv(markers,"markers.csv")
#Then uploaded to MetaCore software for online analysis

