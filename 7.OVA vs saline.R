sc_condHeatmap<-function(obj,condmarker,filename){
	library(dplyr)
	library(ComplexHeatmap)
	genes<-condmarker[condmarker$p_val_adj<0.05&condmarker$avg_logFC>0,] %>% group_by(cluster) %>% top_n(n =500, wt = avg_logFC)
	genes<-genes[!(genes$gene %in% grep("^mt-",genes$gene,value=T)),]
	genes<-genes[!(genes$gene %in% grep("^Rp[sl]",genes$gene,value=T)),]
	obj@active.ident<-factor(obj$condition)
	subobj<-subset(obj , downsample = 500) 
	mat <- GetAssayData(subobj, slot = "counts",assay = "RNA")
	mat<-log2(mat+1)
	mat <- t(scale(t(mat)))
	cluster_info <- sort(obj$orig.ident)
	mat<-as.matrix(mat[unique(genes$gene),names(cluster_info)])
	markgenes<-condmarker[condmarker$p_val_adj<0.05&condmarker$avg_logFC>0,,]
	markgenes<-markgenes[!(markgenes$gene %in% grep("^mt-",markgenes$gene,value=T)),]
	markgenes<-markgenes[!(markgenes$gene %in% grep("^Rp[sl]",markgenes$gene,value=T)),]
	markgenes <- markgenes %>% group_by(cluster) %>% top_n(n=20,wt=avg_logFC)
	gene_pos <- which(rownames(mat) %in% markgenes$gene)
	row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, labels = markgenes$gene))
	tiff(filename,units = "in",res=900,pointsize=8,width=4.5,height=5.6)
	Heatmap(mat,cluster_rows = F,cluster_columns = T,show_column_names = F,show_row_names = F,cluster_column_slices = F,column_split = obj$orig.ident[colnames(mat)],name="Expression",
	right_annotation = row_anno,heatmap_legend_param = list(title="Expression",size = unit(0.1, "inches")))
	dev.off()
}
Idents(asthma_obj)<-"condition"
condmarker<-FindAllMarkers(asthma_obj)
#Fig1D
sc_condHeatmap(subset(obj , downsample = 5000),condmarker,"Fig1D_condition_heatmap.tiff")

reactome_pathway<-function(marker){
	library(ReactomePA)
	marker<-marker[marker$p_val_adj<0.05,]
	geneList<-unique(marker$gene)
	EntrezID<-bitr(geneList,fromType='SYMBOL',toType='ENTREZID',OrgDb='org.Mm.eg.db')
	ewp<-enrichPathway(gene=EntrezID$ENTREZID,organism='mouse',pvalueCutoff=0.05)
	return(ewp)
}

#Fig1E
ewp<-reactome_pathway(condmarker[condmarker$cluster %in% "OVA"&condmarker$avg_logFC>0.25,])
p1<-dotplot(ewp,showCategory=10)
tiff("Fig1E.tiff",units = "in",width=8,height=8.5)
p1
dev.off()

HumanToMouse<-function(hsGenes)
{
	library(biomaRt)
		mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
		human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
	geneMm <- getLDS(attributes = "hgnc_symbol",filters="hgnc_symbol",values=unique(hsGenes),mart=human,attributesL = "mgi_symbol",martL=mouse,uniqueRows = TRUE)
}
HumanToMouse(geneHs)

#Fig1F
CTD_asthma_human<-read.csv("CTD_asthma_genes_human.csv")#Gene list availabe in Table S5
CTD_asthma_mouse<-CTD_asthma_human
tmp<-HumanToMouse(CTD_asthma_mouse$Gene.Symbol)
colnames(tmp)<-c("Gene.Symbol","gene")
CTD_asthma_mouse<-merge(CTD_asthma_mouse,tmp,by="Gene.Symbol",all=T)
#Some gene names are NA after convertion,rename manually
CTD_asthma_mouse[CTD_asthma_mouse$Gene.Symbol=="Ccl5",]<-"Ccl5"
CTD_asthma_mouse[CTD_asthma_mouse$Gene.Symbol=="GSDMD",]<-"Gsdmb"
CTD_asthma_mouse[CTD_asthma_mouse$Gene.Symbol=="LTC4S",]<-"Ltc4s"
CTD_asthma_mouse[CTD_asthma_mouse$Gene.Symbol=="MMP10",]<-"Mmp10"


GSEA_SC<-function(DE,group,REF_Genes,geneSetID=1)
{
	EntrezID<-bitr(DE$gene,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Mm.eg.db)
	DE<-merge(x=DE,y=EntrezID,by.x="gene",by.y="SYMBOL")
	geneList<-as.numeric(as.character(DE[DE$cluster==group,"avg_logFC"]))
	names(geneList) = as.character(DE[DE$cluster==group,"gene"])
	geneList= sort(geneList, decreasing = TRUE)
	gseaResult <- GSEA(geneList = geneList,nPerm= 1000,minGSSize= 5, TERM2GENE=REF_Genes,TERM2NAME =REF_Genes[,c(1,1)],pvalueCutoff = 0.5,verbose= T)
	write.csv(gseaResult@result,file=path)
	gseaplot2(gseaResult,1:geneSetID)
}

CTD_asthma_mouse$pathway<-"CTD_asthma"
p1<-GSEA_SC(condmarker,"OVA",CTD_asthma_mouse[,c("pathway","gene")])
ggsave(filename="Fig1F.tiff",plot=p1,width=10,height=7)
