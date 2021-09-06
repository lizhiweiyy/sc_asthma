#subcluster analysis
sub_clustree<-function(x){
	library(clustree)
	resuo<-c(0,0.05,0.1,0.2)
	for(i in 1:length(resuo)){
		x <- FindClusters(x, resolution =resuo[i]) 
	}
	x.experi<-as.SingleCellExperiment(x) 
	p1<-clustree(x.experi, prefix = "RNA_snn_res.",exprs='logcounts')
}
sub_preprocess<-function(obj){
	obj<- NormalizeData(obj,verbose=F)
	obj<- FindVariableFeatures(obj,nfeatures = 2000,verbose=F)
	obj<-ScaleData(obj, verbose = T,vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
	obj <- RunPCA(obj, npcs = 50, verbose = T)	
}

sub_analysis<-function(subcls){
	DefaultAssay(subcls)<-"RNA"
	subcls<-sub_preprocess(subcls)
	subcls<-sc_process(subcls,22)
}


sub_heatmap<-function(obj,markers){
	top <- markers %>% group_by(cluster) %>% top_n(n =number, wt = avg_logFC)
	p1<-DoHeatmap(subset(obj, downsample =300), features = unique(top$gene))+
	scale_fill_gradientn(colors = colorRampPalette(c("#blue","white","red"))(30))
}

#trans： transpose of matrix 
sc_subpHeatmap<-function(obj,fileName,markgenes,nwidth=4,nheight=5,trans=FALSE){
		library(dplyr)
		library(pheatmap)
		library(scales)
		genes<-markgenes
		avg_obj<-AverageExpression(obj)
		mat<-avg_obj$RNA[markgenes,]
		coulrange = c('purple', 'black', 'yellow')
		palettebreaks <- seq(-2, 2, 0.05)
		mypalette <- colorRampPalette(coulrange, space = 'Lab')
		tiff(fileName,units = "in",width=nwidth,height=nheight)
		if(trans==F){
			pheatmap(mat,scale = "row",cluster_cols = F,cluster_rows = F,
			  color=mypalette(length(palettebreaks)-1),
			  breaks= palettebreaks, border_color = NA,
			  show_colnames = T,
			  )
		 }
		 else{
			pheatmap(t(mat),scale = "column",cluster_cols = F,cluster_rows = F,
			  color=mypalette(length(palettebreaks)-1),
			  breaks= palettebreaks, border_color  = NA,
			  show_colnames = T,
			  annotation_col=annotation_rows,
			  angle_col="90"
			  )
		 }
		dev.off()
}


###################
##Fig2
##################
asthma_neu<-subset(asthma_obj,idents=c("Neu"))
#The saline group have insufficient number of neutrophils for subclustering,
#to avoid batched effect, analysis the OVA group first, and use the Seruat's label transfer algorithm
ova_neu<-subset(asthma_neu,subset=condition=="OVA")
ova_neu<-sub_analysis(ova_neu,1)
#FigS5A clustree
p1<-sub_clustree(ova_neu)
ggsave(filename=paste0(results_path,"/FigS5A.tiff"),width=6.5,height=5.2)
ova_neu<-FindClusters(ova_neu,resolution = 0.2)
ova_neu$cellType<-Idents(ova_neu)

#Fig3 Neutrophil subset label transfer
ova_obj<-subset(asthma_obj,subset=condition=="OVA")
ova_obj<-SetIdent(object = ova_obj, cells = WhichCells(object = ova_neu), value = ova_neu$cellType)
ova_obj$cellType<-Idents(ova_obj)
ova_obj<- FindVariableFeatures(ova_obj, selection.method = "vst", nfeatures = 2000, verbose = F)
DefaultAssay(asthma_obj)<-"RNA"
saline_obj<-subset(asthma_obj,subset=condition=="Saline")
neu.anchors <- FindTransferAnchors(reference = ova_obj, query = saline_obj)
predictions <- TransferData(anchorset = neu.anchors, refdata = ova_obj$cellType)
saline_obj <- AddMetaData(saline_obj, metadata = predictions)
Idents(saline_obj)<-factor(saline_obj$predicted.id)
saline_neu<-subset(saline_obj,idents=c(0:3))
asthma_neu<-merge(saline_neu,ova_neu)
#Fig2A
p1<-DimPlot(ova_neu,reduction = "tsne")
ggsave(filename=paste0(results_path,"/Fig2A.tiff"),width=2.8,height=3)

asthma_neu.marker <- FindAllMarkers(asthma_neu,logfc.threshold = 0.2)

#Fig2B
p1<-sub_heatmap(asthma_neu,asthma_neu.marker)
ggsave(filename=paste0(results_path,"/Fig2B.tiff"),width=2.6,height=3.2)
#Fig2D
vlnfeatures<-c("Gda","Ifitm1","Ccl6","Hmgb2","Il1rn","Bhlhe40","Ccl3","Fcgr2b","Slfn5","Slfn4","Ifit1","Rsad2","Cd177","Ly6g","Ngp","Mmp8")
VlnPlot(asthma_neu,features = vlnfeatures,pt.size = 0)
ggsave(filename=paste0(results_path,"/Fig2D.tiff"),width=10,height=10)

#Fig2E
g_BM<-c("Camp","Ngp","Ltf","Lcn2","Ifitm6","Cd177","Serpinb1a","Ly6g","Chil3","Elane")
g_PBMC<-c("Il1b","Rps27rt","Gm9843","Ifitm1","Wfdc17","Csf3r","Fgl2")
g_Spleen<-c("Fos","Dusp1","Junb","Zfp36","Jund","Ccl5","Cebpb","Cxcl2","Cd74","Slc7a11")
genes<-c(g_BM,g_PBMC,g_Spleen)
sc_subpHeatmap(asthma_neu,"Fig2E.tiff",genes,nwidth=7.5,nheight=2.2,TRUE)
#Fig2F     Some genes have alias:Itgam(Cd11b),Ptprc(Cd45),Fcgr3(Cd16)
DotPlot(asthma_neu,features=c("Ptprc","Itgam","Ly6g","Cd33","Ly6c2","Cd14","Fcgr3"))

#FigS5C
DotPlot(asthma_neu,features = c("S100a8","S100a9","Csf3r","Il1b","Hdc","Il1rn","Cxcl2","Acod1","Marcksl1","Ptgs2","Ccrl2","Ccl3","Ccl4","Csf1","Slfn4","Rsad2","Ifit1","Cd177","Ly6g",
	"Ngp","Ifitm3","Ifitm6","Mcpt8","Cpa3","Fcer1a","F13a1","S100a4","Csf1r","Prg2","Epx","Prg3"))+coord_flip()
#FigS5D
c1_mast_genes<-intersect(asthma_neu.marker[asthma_neu.marker$cluster=="1"&asthma_neu.marker$avg_logFC>0.5,]$gene,asthma_obj.markers[asthma_obj.markers$cluster=="Mast"&asthma_obj.markers$avg_logFC>0.25,]$gene)
p1<-DotPlot(subset(asthma_obj,ident=c("C0","C1","C2","C3","Mast")),features = c1_mast_genes)+coord_flip()

#FigS6A
gsigs<-read.csv("G-sigs.csv")#gene list available in  TableS6
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=gsigs[gsigs$cluster %in% "G0",]$gene.symbol),  ctrl = 100,name = "G0")
p0<-VlnPlot(asthma_neu,features = c("G01"),pt.size = 0)+ggtitle("G0 score")
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=gsigs[gsigs$cluster %in% "G1",]$gene.symbol),  ctrl = 100,name = "G1")
p1<-VlnPlot(asthma_neu,features = c("G11"),pt.size = 0)+ggtitle("G1 score")
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=gsigs[gsigs$cluster %in% "G2",]$gene.symbol),  ctrl = 100,name = "G2")
p2<-VlnPlot(asthma_neu,features = c("G21"),pt.size = 0)+ggtitle("G2 score")
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=gsigs[gsigs$cluster %in% "G3",]$gene.symbol),  ctrl = 100,name = "G3")
p3<-VlnPlot(asthma_neu,features = c("G31"),pt.size = 0)+ggtitle("G3 score")
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=gsigs[gsigs$cluster %in% "G4",]$gene.symbol),  ctrl = 100,name = "G4")
p4<-VlnPlot(asthma_neu,features = c("G41"),pt.size = 0)+ggtitle("G4 score")
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=gsigs[gsigs$cluster %in% "G5a",]$gene.symbol),  ctrl = 100,name = "G5a")
p5<-VlnPlot(asthma_neu,features = c("G5a1"),pt.size = 0)+ggtitle("G5a score")
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=gsigs[gsigs$cluster %in% "G5b",]$gene.symbol),  ctrl = 100,name = "G5b")
p6<-VlnPlot(asthma_neu,features = c("G5b1"),pt.size = 0)+ggtitle("G5b1 score")
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=gsigs[gsigs$cluster %in% "G5c",]$gene.symbol),  ctrl = 100,name = "G5c")
p7<-VlnPlot(asthma_neu,features = c("G5c1"),pt.size = 0)+ggtitle("G5c score")

#FigS6B
p1<-DotPlot(asthma_neu,features =intersect(asthma_neu.marker[asthma_neu.marker$avg_logFC>0.5,]$gene,gsigs[gsigs$cluster %in% "G3"&gsigs$avg_logFC>0.5,]$gene.symbol))+coord_flip()
p2<-DotPlot(asthma_neu,features =intersect(asthma_neu.marker[asthma_neu.marker$avg_logFC>0.5,]$gene,gsigs[gsigs$cluster %in% "G5a",]$gene.symbol))+coord_flip()
p3<-DotPlot(asthma_neu,features =intersect(asthma_neu.marker[asthma_neu.marker$avg_logFC>0.5,]$gene,gsigs[gsigs$cluster %in% "G5b",]$gene.symbol))+coord_flip()
p4<-DotPlot(asthma_neu,features =intersect(asthma_neu.marker[asthma_neu.marker$avg_logFC>0.5,]$gene,gsigs[gsigs$cluster %in% "G5c",]$gene.symbol))+coord_flip()
