############
#Fig2
#############
#GSEA analysis using Reactome database, DE:differentially expressed genes
GSEA_REACTOME<-function(DE,group,geneSetID)
{
	EntrezID<-bitr(DE$gene,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Mm.eg.db)
	DE<-merge(x=DE,y=EntrezID,by.x="gene",by.y="SYMBOL")
	geneList<-as.numeric(as.character(DE[DE$cluster==group,"avg_logFC"]))
	names(geneList) = as.character(DE[DE$cluster==group,"gene"])
	geneList= sort(geneList, decreasing = TRUE)
	wp2gene <- read.gmt("ReactomePathways.gmt")
	wpid2name <- wp2gene
	tmp<-HumanToMouse(wp2gene$gene)#HumanToMouse was defined in 7.OVA vs saline.R
	colnames(tmp)<-c("gene","ms_gene")
	wp2gene<-merge(wp2gene,tmp,by="gene",all=F)
	wp2gene<-wp2gene[,-1]
	wpid2name$gene<-wpid2name$ont
	gseaResult <- GSEA(geneList = geneList,nPerm= 1000,minGSSize= 10, TERM2GENE=wp2gene,TERM2NAME =wpid2name,pvalueCutoff = 0.5,verbose= T)
	write.csv(gseaResult@result,file=path)
	gseaplot2(gseaResult,geneSetID) 
}


#FigS4A 
eos_sig<-read.csv("Eos_asthma_signature.csv")#Gene list availabe in Table S5
eos_sig$pathway<-"Mouse EOS asthma"
neu_sig<-read.csv("Neu_asthma_signature.csv")#Gene list availabe in Table S5
neu_sig$pathway<-"Mouse Neu asthma"
#GSEA_SC was defined in 7.OVA vs saline.R
p2<-GSEA_SC(condmarker,"OVA",rbind(neu_sig[,c("pathway","Gene.symbol")],eos_sig[,c("pathway","Gene.symbol")]))
ggsave(paste0(results_path,"/FigS4A.tiff"),plot = p2,width=3.47,height=4.2)
#FigS4B
eos_pbmc_sig<-HumanToMouse(asthma_human_pbmc$gene_eos)
eos_pbmc_sig$Gene.symbol<-eos_pbmc_sig$MGI.symbol
eos_pbmc_sig$pathway<-"Human EOS asthma"
#NEU asthma
neu_pbmc_sig<-HumanToMouse(asthma_human_pbmc$gene_neu)
neu_pbmc_sig$Gene.symbol<-neu_pbmc_sig$MGI.symbol
neu_pbmc_sig$pathway<-"Human NEU asthma"
p3<-GSEA_SC(condmarker,"OVA",rbind(neu_pbmc_sig[,c("pathway","Gene.symbol")],eos_pbmc_sig[,c("pathway","Gene.symbol")]))
ggsave(paste0(results_path,"/FigS4B.tiff"),plot = p3,width=3.47,height=4.2)
#FigS4C
p1<-GSEA_REACTOME(condmarker,"OVA",c("Neutrophil degranulation","Interleukin-10 signaling"))
ggsave(paste0(results_path,"/FigS4C.tiff"),plot = p1,width=3.47,height=4.2)




sc_pHeatmap<-function(obj,fileName,markgenes){
		
		genes<-markgenes
		avg_obj<-AverageExpression(obj)
		mat<-avg_obj$integrated[markgenes,]
		coulrange = c('purple', 'black', 'yellow')
		palettebreaks <- seq(-4, 4, 0.1) #
		mypalette <- colorRampPalette(coulrange, space = 'Lab')
		tiff(fileName,units = "in",compression = "lzw",res=300,pointsize=8,width=17,height=12)
		pheatmap(mat,scale = "row",cluster_cols = F,cluster_rows = F,
		  color=mypalette(length(palettebreaks)-1),
		  breaks= palettebreaks, border_color      = NA)
		dev.off()
}
#FigS4D
sc_pHeatmap(asthma_obj,"FigS4D.tiff",c("Il1b","Il33","S100a8","S100a9","Csf3r","Retnla","Retnlg","Slc7a11","Cxcl2","Cd300lf","S100a11","Tnfaip2","Chil3","Ifitm1","AA467197","Lilr4b","Marcksl1","Clec4d","Spi1","Ftl1","Cyp4f18","Pglyrp1","Ccl9","Slfn4","Rsad2","Cd44","Tmsb4x","Fth1","Cd14","Abcg1","Pilra","Cldn5","Macf1","Gsn","Tyrobp","Fgl2","Mmp9","Mmp8","Slpi","Hp","Fcer1g","Cd33","Itgam","Ptgs2","Il1rn","Ccr1","Ptafr","Ccl4","Ccr2","Tnfrsf1b","Ccl3","Il1r2"))



