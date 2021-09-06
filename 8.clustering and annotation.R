sc_process<-function(x,ndims=30){	
	x <- RunPCA(x, npcs = ndims, verbose = F)	
	x <- FindNeighbors(x, reduction = "pca", dims = 1:ndims)
	x <- RunUMAP(x, reduction = "pca", dims = 1:ndims)
	x <- RunTSNE(x, reduction = "pca", dims = 1:ndims)
}

#clustering
sc_clustree<-function(x){
	library(clustree)
	resuo<-c(0,0.01,0.05,0.1,0.2,0.5,1,2)
	for(i in 1:length(resuo)){
		x <- FindClusters(x, resolution =resuo[i]) 
	}
	x.experi<-as.SingleCellExperiment(x)
	p1<-clustree(x.experi, prefix = "integrated_snn_res.",exprs='logcounts')
}


asthma_obj<-sc_process(asthma_obj,60)
ElbowPlot(asthma_obj,ndims = 60)
asthma_obj <- JackStraw(asthma_obj, num.replicate = 20,dims=60)
asthma_obj<- ScoreJackStraw(asthma_obj, dims = 1:60)
JackStrawPlot(asthma_obj, dims = 1:60)
sc_clustree(asthma_obj)
asthma_obj <- FindClusters(object = asthma_obj, resolution = 1)
asthma_obj.markers <- FindAllMarkers(asthma_obj)

#  annotation
asthma_obj<-RenameIdents(asthma_obj, '0'='Endo','1'='TRAM','2'='Neu','3'='Endo','4'='B',
	'5'='B','6'='T','7'='Monocyte','8'='Endo','9'='Fibro','10'='T','11'='TRIM',
	'12'='Fibro','13'='MDM','14'='Endo','15'='Fibro','16'='NKT','17'='Endo','18'='T',
	'19'='Neu','20'='SMC','21'='T','22'='DC','23'='Endo','24'='AT2','25'='Endo',
	'26'='TRAM','27'='T','28'='DC','29'='Pericyte','30'='T','31'='lym-Endo','32'='DC',
	'33'='ILC2','34'='AT1','35'='35','36'='Ciliated','37'='Mesothelial','38'='Mast','39'='Endo')
#remove low quality cluster
asthma_obj<-subset(asthma_obj,idents=c("35"),invert=T)	
#SMC has two unique subsets
c20<-subset(asthma_obj,idents=c("SMC"))
c20<-sc_process(c20,60)
p1<-sc_clustree(c20)
c20<-FindClusters(c20,resolution = 0.05)
c20.markers <- FindAllMarkers(c20, only.pos =TRUE, min.pct =0.2,logfc.threshold = 0.3)
c20<-RenameIdents(c20, '0'='VSMC','1'='ASMC')
cells.use <- WhichCells(object = c20, idents = "ASMC")
asthma_obj <- SetIdent(object = asthma_obj, cells = cells.use, value = 'ASMC')
cells.use <- WhichCells(object = c20, idents = "VSMC")
asthma_obj <- SetIdent(object = asthma_obj, cells = cells.use, value = 'VSMC')
asthma_obj$cellType<-Idents(asthma_obj)
#Fig1G
DimPlot(asthma_obj,reduction = "tsne",label = T)
ggsave(filename=paste0(results_path,"/Fig1G.tiff"),width=8,height=7)
DimPlot(asthma_obj,label=F,reduction="tsne",group.by = "orig.ident")
ggsave(filename=paste0(results_path,"/Fig1H.tiff"),width=8,height=7)
#FigS1C
top <- asthma_obj.markers %>% group_by(cluster) %>% top_n(n =3, wt = avg_logFC)
p1<-DoHeatmap(subset(asthma_obj, downsample =300), features = unique(top$gene))
ggsave(filename=paste0(results_path,"/FigS1C_heatmap.tiff"),plot=p1,width=14,height=10)
#FigS2A
p1<-DimPlot(asthma_obj,reduction = "umap",label=T)
##FigS2B
clfeatures<-c("Ager","Sftpc","Cd3d","Ms4a1","Gzma","S100a8","Icos","Col1a2","Myl9","Cdh5","Cpa3","Cd68","Mrc1","Siglecf","Itgam","Cx3cr1","Scgb3a2","Prox1","F13a1","Msln")
FeaturePlot(asthma_obj,features = clfeatures,order = F, min.cutoff = "q1", max.cutoff ="q99",reduction="umap")+NoLegend()
