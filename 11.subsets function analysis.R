my_go<-function(markers,group){
	topgenes<-markers[markers$cluster %in% group,]
	EntrezID<-bitr(unique(topgenes$gene),fromType='SYMBOL',toType='ENTREZID',OrgDb='org.Mm.eg.db')
	gobp<-enrichGO(gene=EntrezID$ENTREZID,'org.Mm.eg.db',ont="BP",pvalueCutoff=0.05)
}
my_kegg<-function(markers,group){
	topgenes<-markers[markers$cluster %in% group,]
	EntrezID<-bitr(unique(topgenes$gene),fromType='SYMBOL',toType='ENTREZID',OrgDb='org.Mm.eg.db')
	kegg<-enrichKEGG(gene=EntrezID$ENTREZID,organism='mmu',pvalueCutoff=0.05)
}
my_react<-function(markers,group){
	library(ReactomePA)
	topgenes<-markers[markers$cluster %in% group,]
	EntrezID<-bitr(unique(topgenes$gene),fromType='SYMBOL',toType='ENTREZID',OrgDb='org.Mm.eg.db')
	react<-enrichPathway(gene=EntrezID$ENTREZID,organism='mouse',pvalueCutoff=0.05)
}

markers <- FindAllMarkers(asthma_neu,only.pos = T)
#Fig3A Gene Ontology Enrichment
mkclusters<-levels(markers$cluster)
markerList<-list()
for(i in 1:length(mkclusters))
{
	tmpmarkers<-markers[markers$cluster %in% mkclusters[i],]
	EntrezID<-bitr(unique(tmpmarkers$gene),fromType='SYMBOL',toType='ENTREZID',OrgDb='org.Mm.eg.db')
	markerList[[i]]<-EntrezID$ENTREZID
}
names(markerList)<-mkclusters
GO=compareCluster(markerList, fun='enrichGO',OrgDb=org.Mm.eg.db,ont="BP",readable=T)
dotplot(GO,includeAll=T,showCategory=5)
ggsave(filename=paste0(results_path,"/Fig3A.tiff"),width=2.4,height=4.6)


#Fig3B KEGG
markerList<-list()
for(i in 1:length(mkclusters))
{
	tmpmarkers<-markers[markers$cluster %in% mkclusters[i],]
	EntrezID<-bitr(unique(tmpmarkers$gene),fromType='SYMBOL',toType='ENTREZID',OrgDb='org.Mm.eg.db')
	markerList[[i]]<-EntrezID$ENTREZID
}
names(markerList)<-mkclusters
KEGG=compareCluster(markerList, fun='enrichKEGG',organism='mmu')
dotplot(KEGG,includeAll=T,showCategory=5)
ggsave(filename=paste0(results_path,"/Fig3B.tiff"),width=2.4,height=4.6)

my_gseaKEGG<-function(markers){

	EntrezID<-bitr(unique(markers$gene),fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Mm.eg.db)
	markers<-merge(x=markers,y=EntrezID,by.x="gene",by.y="SYMBOL")
	geneList<-as.numeric(as.character(markers[markers$cluster==group,"avg_logFC"]))
	names(geneList) = as.character(markers[markers$cluster==group,"ENTREZID"])
	geneList= sort(geneList, decreasing = TRUE)
	gkegg <- gseKEGG(geneList = geneList,organism= 'mmu',nPerm= 1000,minGSSize= 10, pvalueCutoff=0.5)
}

my_gseaGO<-function(markers,group){
	EntrezID<-bitr(unique(markers$gene),fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Mm.eg.db)
	markers<-merge(x=markers,y=EntrezID,by.x="gene",by.y="SYMBOL")
	geneList<-as.numeric(as.character(markers[markers$cluster==group,"avg_logFC"]))
	names(geneList) = as.character(markers[markers$cluster==group,"ENTREZID"])
	geneList= sort(geneList, decreasing = TRUE)
	gseaGO<- gseGO(geneList =geneList, 'org.Mm.eg.db', ont="BP", nPerm = 10000, minGSSize = 10, pvalueCutoff=1)
}

sc_GSEA_REACTOME<-function(DE)
{
	EntrezID<-bitr(DE$gene,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Mm.eg.db)
	DE<-merge(x=DE,y=EntrezID,by.x="gene",by.y="SYMBOL")
	geneList<-as.numeric(as.character(DE[DE$cluster==group,"avg_logFC"]))
	names(geneList) = as.character(DE[DE$cluster==group,"gene"])
	geneList= sort(geneList, decreasing = TRUE)
	wp2gene <- read.gmt("/path/ReactomePathways.gmt")
	wpid2name <- wp2gene
	tmp<-HumanToMouse(wp2gene$gene)
	colnames(tmp)<-c("gene","ms_gene")
	wp2gene<-merge(wp2gene,tmp,by="gene",all=F)
	wp2gene<-wp2gene[,-1]
	wpid2name$gene<-wpid2name$ont
	gseaResult <- GSEA(geneList = geneList,nPerm= 10000,minGSSize= 10, TERM2GENE=wp2gene,TERM2NAME =wpid2name,pvalueCutoff = 1,verbose= T)
}

gseGO_neu0<-my_gseaGO(asthma_neu.marker,group = "0")
gseGO_neu1<-my_gseaGO(asthma_neu.marker,group = "1")
gseGO_neu2<-my_gseaGO(asthma_neu.marker,group = "2")
gseGO_neu3<-my_gseaGO(asthma_neu.marker,group = "3")

#Fig3D
gseKEGG_neu0<-my_gseaKEGG(asthma_neu.marker,group = "0")
gseKEGG_neu1<-my_gseaKEGG(asthma_neu.marker,group = "1")
gseKEGG_neu2<-my_gseaKEGG(asthma_neu.marker,group = "2")
gseKEGG_neu3<-my_gseaKEGG(asthma_neu.marker,group = "3")
gseaplot2(gseKEGG_neu1,c("mmu04060","mmu04210"))
ggsave(paste0(results_path,"/Fig3D.tiff"),width=3.45,height=2.5)
#Fig3E
gseReact_neu0<-sc_GSEA_REACTOME(asthma_neu.marker,group = "0")
gseReact_neu1<-sc_GSEA_REACTOME(asthma_neu.marker,group = "1")
gseReact_neu2<-sc_GSEA_REACTOME(asthma_neu.marker,group = "2")
gseReact_neu3<-sc_GSEA_REACTOME(asthma_neu.marker,group = "3")
gseaplot2(gseReact_neu2,c("Interferon Signaling","Antiviral mechanism by IFN-stimulated genes")
ggsave(paste0(results_path,"/Fig3E.tiff"),width=3.45,height=2.5)
#Fig3F
gseaplot2(gseReact_neu3,c("Cytokine Signaling in Immune system","Neutrophil degranulation")
ggsave(paste0(results_path,"/Fig3F.tiff"),width=3.45,height=2.5)

#Fig3G
g_Apoptosis<-unlist(strsplit(gseKEGG_neu1@result["mmu04210","core_enrichment"],"/"))
g_Cytokine<-unlist(strsplit(gseKEGG_neu1@result["mmu04060","core_enrichment"],"/"))
g_Apoptosis<-setdiff(g_Apoptosis,g_Cytokine)
genes<-c(g_Apoptosis,g_Cytokine)
sc_subpHeatmap(asthma_neu,"Fig3G.tiff",genes,2.5,5.2)
#Fig3H
g_interferon<-unlist(strsplit(gseReact_neu2@result["Interferon Signaling","core_enrichment"],"/"))
g_antivirus<-unlist(strsplit(gseReact_neu2@result["Antiviral mechanism by IFN-stimulated genes","core_enrichment"],"/"))
genes<-c(g_interferon,g_antivirus)
sc_subpHeatmap(asthma_neu,"Fig3H.tiff",genes,2.5,5.2)
#Fig3I
g_degranulation<-unlist(strsplit(gseReact_neu3@result["Neutrophil degranulation","core_enrichment"],"/"))
g_Netosis<-unlist(strsplit(gseKEGG_neu3@result["mmu04613","core_enrichment"],"/"))
g_degranulation<-g_degranulation[-c(17:22)]
genes<-c(g_degranulation,g_Netosis)
sc_subpHeatmap(asthma_neu,"Fig3I.tiff",genes,2.5,4.4)





phago_genes<-bitr(gseKEGG_neu1@geneSets["mmu04145"][[1]],fromType='ENTREZID',toType='SYMBOL',OrgDb=org.Mm.eg.db)
tmp<-intersect(phago_genes$SYMBOL,rownames(asthma_neu))
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=tmp),  ctrl = 100,name = "Phagosome")
#Lysosome
Lysosome_genes<-bitr(gseKEGG_neu1@geneSets["mmu04142"][[1]],fromType='ENTREZID',toType='SYMBOL',OrgDb=org.Mm.eg.db)
tmp<-intersect(Lysosome_genes$SYMBOL,rownames(asthma_neu))
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=tmp),  ctrl = 100,name = "Lysosome")
#Apoptosis
Apoptosis_genes<-bitr(gseKEGG_neu1@geneSets["mmu04210"][[1]],fromType='ENTREZID',toType='SYMBOL',OrgDb=org.Mm.eg.db)
tmp<-intersect(Apoptosis_genes$SYMBOL,rownames(asthma_neu))
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=tmp),  ctrl = 100,name = "Apoptosis")
#Cytokine-cytokine receptor interaction
Cytokine_genes<-bitr(gseKEGG_neu1@geneSets["mmu04060"][[1]],fromType='ENTREZID',toType='SYMBOL',OrgDb=org.Mm.eg.db)
tmp<-intersect(Lysosome_genes$SYMBOL,rownames(asthma_neu))
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=tmp),  ctrl = 100,name = "Cytokine")
#neutrophil degranulation
neu_degranulation_genes<-bitr(gseGO_neu1@geneSets["GO:0043312"][[1]],fromType='ENTREZID',toType='SYMBOL',OrgDb=org.Mm.eg.db)
tmp<-intersect(neu_degranulation_genes$SYMBOL,rownames(asthma_neu))
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=tmp),  ctrl = 100,name = "neu_degranulation")
#response to molecule of bacterial origin
bacterialres_genes<-bitr(gseGO_neu1@geneSets["GO:0002237"][[1]],fromType='ENTREZID',toType='SYMBOL',OrgDb=org.Mm.eg.db)
tmp<-intersect(bacterialres_genes$SYMBOL,rownames(asthma_neu))
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=tmp),  ctrl = 100,name = "bacterialres")
#response to type I interferon
ifn1res_genes<-bitr(gseGO_neu1@geneSets["GO:0034340"][[1]],fromType='ENTREZID',toType='SYMBOL',OrgDb=org.Mm.eg.db)
tmp<-intersect(ifn1res_genes$SYMBOL,rownames(asthma_neu))
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=tmp),  ctrl = 100,name = "ifn1res")
#response to  interferon-gamma
ifngres_genes<-bitr(gseGO_neu1@geneSets["GO:0034341"][[1]],fromType='ENTREZID',toType='SYMBOL',OrgDb=org.Mm.eg.db)
tmp<-intersect(ifngres_genes$SYMBOL,rownames(asthma_neu))
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=tmp),  ctrl = 100,name = "ifngres")
#Neutrophil extracellular trap formation
NETosis_genes<-bitr(gseKEGG_neu1@geneSets["mmu04613"][[1]],fromType='ENTREZID',toType='SYMBOL',OrgDb=org.Mm.eg.db)
tmp<-intersect(NETosis_genes$SYMBOL,rownames(asthma_neu))
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=tmp),  ctrl = 100,name = "NETosis")
#defense response to virus	
defvirus_genes<-bitr(gseGO_neu1@geneSets["GO:0051607"][[1]],fromType='ENTREZID',toType='SYMBOL',OrgDb=org.Mm.eg.db)
tmp<-intersect(defvirus_genes$SYMBOL,rownames(asthma_neu))
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=tmp),  ctrl = 100,name = "defvirus")
#Interleukin-10 signaling
IL10_genes<-gseReact_neu1@geneSets["Interleukin-10 signaling"][[1]]
tmp<-intersect(IL10_genes,rownames(asthma_neu))
asthma_neu <- AddModuleScore(asthma_neu, features = list(cc=tmp),  ctrl = 100,name = "IL10")
#Fig3J-L
vlnfeatures<-c("Apoptosis1","defvirus1","neu_degranulation1","Cytokine1","ifngres1","NETosis1")
titles<-c("Apoptosis","defense response to virus","neutrophil degranulation","Cytokine-cytokine receptor interaction","response to type I interferon","Neutrophil extracellular ntrap formation")
pvfeatures<-list()
for(i in 1 :length(vlnfeatures)){
	pvfeatures[[i]]<-VlnPlot(asthma_neu,features = c(vlnfeatures[i]),pt.size = 0))+coord_flip()
}
cowplot::plot_grid(plotlist = pvfeatures,ncol = 3)


