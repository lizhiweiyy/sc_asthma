DefaultAssay(asthma_obj)<-"RNA"
ova_obj<-subset(asthma_obj,subset=condition=="OVA")
#Fig6
library("CellChat")
#Assign neutrophil subsets names to celltype 
ova_obj$cellTypesub<-paste0("C",ova_obj$cellType)
#Set neutrophil subsets
ova_obj<-SetIdent(ova_obj,cells=colnames(asthma_neu),value=asthma_neu$cellTypesub)
ova_obj$cellTypesub<-Idents(ova_obj)

clchat_ova<-createCellChat(ova_obj,group.by = "cellTypesub")
clchat_ova@DB<-subsetDB(full_db,search = "Secreted Signaling")

cchat_process<-function(obj)
{
	obj <- subsetData(obj)
	obj <- identifyOverExpressedGenes(obj)
	obj <- identifyOverExpressedInteractions(obj)
	obj <- projectData(obj, PPI.mouse)
	obj <- computeCommunProb(obj)
	obj <- filterCommunication(obj, min.cells = 10)
	obj <- computeCommunProbPathway(obj)
	obj <- aggregateNet(obj)#
}
clchat_ova<-cchat_process(clchat_ova)

#Fig11A
groupSize <- as.numeric(table(clchat_ova@idents))
mat <- clchat_ova@net$weight
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[1:4, ] <- mat[1:4,]
tiff(file=paste0(results_path,"FigS11A.tiff"),width=4,height=4)
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat))
dev.off()
#FigS11B
netVisual_bubble(clchat_ova, sources.use = c(1:4), targets.use = c(15:23), remove.isolate = T)
#FigS11C
netVisual_bubble(clchat_ova, sources.use = c(1:4), targets.use = c(1:14), remove.isolate = T)
#FigS11D
netVisual_bubble(clchat_ova, sources.use = c(5:23), targets.use = c(1:4), remove.isolate = T)




#The expression of Th2/Th17 cytokines
#FigS12A
p1<-DotPlot(asthma_obj,assay = "RNA",features = c("Il4","Il13","Il5","Il17a","Il17f","Il6","Il1b","Il33"),scale.by = "size",dot.scale = 4)+coord_flip()
#FigS12B
p1<-DotPlot(asthma_obj,assay = "RNA",features = c("Il4","Il13","Il5","Il17a","Il17f","Il6","Il1b","Il33"),scale.by = "size",group.by = "condition")+coord_flip()



#compare saline and OVA neutrophils
saline_obj_neu<-subset(saline_obj,idents=c("C0","C1","C2","C3"))
saline_obj_neu$cellTypesub<-saline_obj_neu@active.ident
clchat_saline_neu<-createCellChat(saline_obj_neu,group.by = "cellTypesub")
clchat_saline_neu@DB<-secreted_db
clchat_saline_neu<-cchat_process(clchat_saline_neu)
ova_obj_sub<-subset(ova_obj,idents=c("C0","C1","C2","C3"))
ova_obj_sub$cellTypesub<-ova_obj_sub@active.ident
clchat_ova_neu<-createCellChat(ova_obj_sub,group.by = "cellTypesub"))
clchat_ova_neu@DB<-secreted_db
clchat_ova_neu<-cchat_process(clchat_ova_neu)

object.sublist <- list( Saline = clchat_saline_neu, OVA = clchat_ova_neu)
cellchat_asthma_neu <- mergeCellChat(object.sublist, add.names = names(object.sublist))
#Fig5C The expression of representative Th2/Th17 receptors
p1<-plotGeneExpression(cellchat_asthma_neu,split.by="condition",features =c("Il4ra","Il13ra1","Il17ra"))
##split group by condition and reordering
asthma_neu$cellType_condition<-factor(paste0(asthma_neu$cellType,"--",asthma_neu$condition))
cellTypes<-levels(asthma_neu$cellType)
asthma_neu$cellTypeID<-2*as.numeric(asthma_neu$cellType)
asthma_neu$cdtID<- 2-as.numeric(asthma_neu$condition)
asthma_neu$cellTypeID<-asthma_neu$cellTypeID+asthma_neu$cdtID
asthma_neu_splitByCondition<-asthma_neu
Idents(asthma_neu_splitByCondition)<-"cellType_condition"
asthma_neu_splitByCondition<-ReorderIdent(asthma_neu,var='cellTypeID')
#Fig5D
p2<-VlnPlot(asthma_neu_splitByCondition,features = c("Fcgr2b"),pt.size = 0,assay = "RNA")
#Fig5E
p1<-DotPlot(asthma_neu_splitByCondition,features=rev(c("Runx1","Irf1","Prdm1","Irf8","Nfe2l2","Bhlhe40","Egr2","Egr3","Rbpj","Xbp1","Cebpb")))+ggtitle("Transcription factors")+coord_flip()
ggsave(paste0(results_path,"/Fig5E.tiff"),plot = p1,pointsize=8,width=5.5,height=5)

#Fig5F
p1<-DotPlot(asthma_neu_splitByCondition,features=strsplit(go_neu1@result["GO:0016570","geneID"],split = '/')[[1]])+ggtitle("histone modification")+
	theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+ylab("")+coord_flip()
ggsave(paste0(results_path,"/Fig5F.tiff"),plot = p1,pointsize=8,width=6,height=5)

#Fig5G
p1<-DotPlot(asthma_neu_splitByCondition,features=rev(strsplit(go_neu1@result["GO:0002224","geneID"],split = '/')[[1]]))+ggtitle("toll-like receptor signaling pathway")+coord_flip()
ggsave(paste0(results_path,"/Fig5G.tiff"),plot = p1,pointsize=8,width=6,height=5)


#FigS13A
Idents(asthma_neu)<-"cellType"
p1<-DotPlot(asthma_neu,features=rev(strsplit(go_neu1@result["GO:0019933","geneID"],split = '/')[[1]]))+ggtitle("cAMP-mediated signaling")+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+ylab("")+coord_flip()
ggsave(paste0(results_path,"/FigS13A.tiff"),plot = p1,pointsize=8,width=6,height=6)
#FigS13B
p1<-DotPlot(asthma_neu,features=rev(strsplit(go_neu1@result["GO:0019882","geneID"],split = '/')[[1]]))+ggtitle("antigen processing and presentation")+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+ylab("")+coord_flip()
ggsave(paste0(results_path,"/FigS13B.tiff"),plot = p1,pointsize=8,width=6,height=6)

#FigS14A
p1<-VlnPlot(asthma_neu,features = c("Csf3r","Ly6g","Fcgr2b","Cd177"),pt.size = .1,ncol=4)
ggsave(paste0(results_path,"/FigS14A.tiff"),plot = p1,pointsize=8,width=8.4,height=3.6)

