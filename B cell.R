##B cell
Idents(ifnb1) <- "celltype"
pbmc_B <- subset(ifnb1, idents = "B cell")
rm(pbmc1)

pbmc_B<- NormalizeData(pbmc_B) 
pbmc_B <- FindVariableFeatures(pbmc_B)
pbmc_B <- ScaleData(pbmc_B)
pbmc_B <- RunPCA(pbmc_B, features = VariableFeatures(object = pbmc_B))
ElbowPlot(pbmc_B,ndims = 30)

pbmc1  <- FindNeighbors(pbmc_B, dims = 1:15)
pbmc1 <- FindClusters(pbmc1, resolution = 0.2) 
pbmc1 <- RunUMAP(pbmc1, dims = 1:15)
DimPlot(pbmc1, reduction = "umap", shuffle = T, label = T)
rm(pbmc_B2)


Idents(pbmc1) <- "seurat_clusters"
pbmc2 <- pbmc1

new.cluster.ids <- c("preB cell","mB cell", "immuB cell","preB cell", "proB cell", "immuB cell")
names(new.cluster.ids) <- levels(pbmc2)
pbmc2 <- RenameIdents(pbmc2, new.cluster.ids) 
pbmc2$sub_celltype <- pbmc2@active.ident


##UMAP
color1=c("#97A175","#A36F3D","#F2C847","#D3CBB8")
Idents(pbmc2) <- factor(Idents(pbmc2), levels = c("proB cell", "preB cell","immuB cell","mB cell"))
plot <- print(DimPlot(pbmc2, reduction = "umap", shuffle = T,cols = color1))
plot
ggsave(plot, filename = "umap_B.pdf", width=6, height = 5)

DimPlot(pbmc2, reduction = "umap", group.by ="group")
ggsave("umap_B_group.pdf", width=6, height = 5)

colourCount = length(unique(pbmc2@meta.data$celltype))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
celltype_colors <- getPalette(colourCount)
plotC <- table(pbmc2@meta.data$group, pbmc2@meta.data$sub_celltype) %>% melt()
colnames(plotC) <- c("Group", "CellType","Number")
write.csv(plotC, "Cellp_Neu.csv", row.names = T)
plotC <- data.table::fread("Cellp_B.csv", header = T)


## cell proportion
color2=c("#F2C847","#D3CBB8","#97A175","#A36F3D")
plotC$CellType <- factor(plotC$CellType,ordered = T,levels = c("proB cell", "preB cell","immuB cell","mB cell"))
pC2 <- ggplot(data = plotC, aes(x = Group, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=color1) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x=NULL,y=NULL)+
  scale_y_continuous(labels = percent)+  
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6))

pC2
ggsave(pC2,filename = "cellp_B.pdf",width = 3.5,height = 5)


#### marker dot
genes <- c("Vpreb1","Vpreb3",
           "Bcl7a","Dnajc7",
           "Iglc1","Ms4a1",
           "H2-Aa","H2-Eb1")
)
p1 <- DotPlot(pbmc2, features = genes)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x.top = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
  labs(x="Markers", y="Celltype")+guides(size=guide_legend(order=3))          
p1

##
Idents(pbmc2) <- factor(Idents(pbmc2), levels = c("proB cell","preB cell","immuB cell","mB cell"))
gene <- c("Isg15","Ifi213","Bst2","Ifi206","Trim30d","Eif2ak2","Trim30a","Ifi214","Stat1","Slfn8"
)
p1 <- DotPlot(pbmc2, features = gene)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x.top = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#F7DCD0", "#CF443B","#80181D"))+#颜色渐变设置  
  labs(x="Markers", y="Celltype")+guides(size=guide_legend(order=3))
p1  
ggsave(p1,filename="genedot1.pdf",width=7,height=3)


##top5 heatmap
pbmc1=pbmc2
rm(pbmc2)

pbmc1 <- JoinLayers(pbmc1)
markers <- FindAllMarkers(pbmc1, only.pos = T, logfc.threshold = 0.25)
write.csv(markers, "allmarkers_pos.csv", row.names = T)
markers <- read.csv("allmarkers_pos.csv",row.names = 1)
marker.sig <- markers %>% 
  mutate(pct.fc = pct.1-pct.2) %>%
  filter(p_val_adj < 0.05, 
         pct.1 > 0.5,
         pct.2 < 0.5,
         pct.fc > 0.2,
         avg_log2FC > 1)
write.csv(marker.sig, "allmarkers_pos_fil.csv", row.names = T)
setwd("H:/HBOC/data1/B")
marker.sig <- read.csv("allmarkers_pos_fil.csv",row.names = 1)

marker.sig1 <- markers %>% 
  filter(pct.1 >0.1 & p_val_adj < 0.05) %>%
  filter(abs(avg_log2FC) > 0.5)

top5 = marker.sig  %>% group_by(cluster) %>% top_n(n = 5, wt =avg_log2FC)
top20= marker.sig  %>% group_by(cluster) %>% top_n(n = 20, wt =avg_log2FC)
View (top20)

markerdata <- ScaleData(pbmc2,feature=as.character(unique(top20$gene)), 
                        assay = "RNA")

Idents(markerdata) <- "sub_celltype"
Idents(markerdata) <- factor(Idents(markerdata), 
                             levels = c("proB cell","preB cell","immuB cell","mB cell"))
DoHeatmap(markerdata,
          feature=as.character(unique(top5$gene)),
          assay= "RNA",
          group.colors=c("#97A175","#A36F3D","#F2C847","#D3CBB8"))+
  scale_fill_gradientn(colors=c("white","grey","firebrick3"))
ggsave("heatmap_celltype.pdf", width = 6,height =7.2)



##volcano plot
B_group <- FindMarkers(pbmc1,
                       threshold = 0.25,
                       only.pos = F,
                       ident.1 = "HBOC", ident.2 = "Control",
                       group.by = "group") %>%
  mutate(gene = rownames(.))

write.csv(B_group,file = "B_group.csv")


library(ggrepel)
B_group <- B_group %>%
  mutate(Difference = pct.1 - pct.2) 
B_group$change=0

for (i in 1:nrow(B_group)) {
  if(B_group$avg_log2FC[i]>=1){B_group$change[i]="up"}
  else if(B_group$avg_log2FC[i]<=-1){B_group$change[i]="down"}
  else{B_group$change[i]="ns"}
}


p1 <- ggplot(B_group, aes(x=Difference, y=avg_log2FC) )+ 
  geom_point(size=1.2,aes(color=change)) + 
  scale_color_manual("change",labels=c("down","ns","up"),
                     values=c("#0d1b46", "grey", "tomato"))+
  geom_label_repel(data=subset(B_group, 
                               avg_log2FC >= 2 & Difference >= 0.4 & p_val_adj <= 0.05), 
                   aes(label=gene), 
                   color="black", 
                   segment.colour = "black",
                   label.padding = 0.2, 
                   segment.size = 0.3,
                   size=4,
                   max.overlaps = 20) + 
  geom_label_repel(data=subset(B_group, 
                               avg_log2FC <= -1 & Difference <= -0.1 & p_val_adj <= 0.05), 
                   aes(label=gene), 
                   label.padding = 0.2, 
                   color="black",
                   segment.colour = "black",
                   segment.size = 0.3, size=4,
                   max.overlaps = 20) + 
  geom_vline(xintercept = 0,linetype = 2) +
  geom_hline(yintercept = 0,linetype = 2) +
  labs(x="Percentage difference",y="Log-Fold change") + 
  theme_bw()+
  xlim(-0.4,1.2)+
  ylim(-4,8)+
  ggtitle(label ="HBOC vs Control")
p1
ggsave(p1,filename="p1.pdf",width=6.5,height=5.5)


##GO/KEGG
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)

B_group_up <- B_group %>%
  filter(Difference >= 0.25 & p_val_adj < 0.05) %>%
  filter(avg_log2FC >=2)

B_group_up$gene <- rownames(B_group_up)
ids=bitr(B_group_up$gene, 'SYMBOL', 'ENTREZID', 'org.Mm.eg.db')
B_group_up=merge(B_group_up, ids,  by.x='gene', by.y='SYMBOL')
head(B_group_up)
##entichGo
B_group_up <- B_group_up[order(B_group_up$avg_log2FC,decreasing = T),]
group_degs_list <- as.numeric(B_group_up$avg_log2FC)
names(group_degs_list) <- B_group_up$ENTREZID
head(group_degs_list)

group_de <- names(group_degs_list)[group_degs_list > 2]
head(group_de)

group_ego <- enrichGO(group_de, OrgDb = "org.Mm.eg.db", ont = "BP", readable = T)
head(group_ego)
write.csv(group_ego, "B_group_ego_up.csv")
setwd("H:/HBOC/数据1/B")
group_ego1 <- read.csv("fil.csv",row.names = 1)

p1 <- ggplot(group_ego1, aes(y=Description,x=Count)) + 
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=p.adjust),stat="identity") + 
  scale_fill_gradient2(expression(p.adjust),low="#80181D", mid="#CF443B", high="#F7DCD0",midpoint = 0.00015)+
  ylab(NULL) + xlab("Gene count") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() )
p1
ggsave("go_bar.pdf",width=6,height=3)


##KEGG
group_ekg_up <- enrichKEGG(gene = group_de, organism = "mmu", pvalueCutoff = 0.05)
head(group_ekg_up)
write.csv(group_ekg_up, "B_group_kegg_up.csv")

group_kegg1 <- read.csv("fil1.csv",row.names = 1)

p1 <- ggplot(group_kegg1, aes(y=Description,x=Count)) + 
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=p.adjust),stat="identity") + 
  scale_fill_gradient2(expression(p.adjust),low="#80181D", mid="#CF443B", high="#F7DCD0",midpoint = 0.008)+
  ylab(NULL) + xlab("Gene count") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() )
p1
ggsave("kegg_bar.pdf",width=4.5,height=3)


##gene set score
# Differentiation
Differentiation=list(c("Itfg2","Cmtm7","Traf3ip2","Nckap1l","Foxp1",
                       "Card11","Igkj5","Abl1","Ada","Adam17",
                       "Ahr","Slc25a5","Aqp8","Atm","Bad",
                       "Bak1","Bax","Bcl2","Bcl3","Bcl6",
                       "Zfp36l1","Zfp36l2","Btk","Cd19","Ms4a1","Cd24a",
                       "Cd79a","Cdh17","Cebpg","Cr2","Dll1",
                       "Dpp4","Bcl11a",'Ezh2','Fas','Il4i1',
                       'Flt3','Fosl2','Fzd9','Lilrb4a','H2-Ab1',
                       "Slc39a7","Ptpn6","Hdac5","Hhex","Hmgb3",
                       "Hmga1","Onecut1","Irf8","Id2","Ifna1",
                       "Ifna11","Ifna2","Ifna4","Ifna5","Ifna6",
                       "Ifna7","Ifna9","Ifnab","Ifnb1","Cd79b",
                       "Ighg1","Ighm","Igkc","Il10","Il2",
                       "Il2rg","Il6","Il7","Il9","Il9r",
                       "Gimap1","Inpp5d","Itm2a","Jak3","Kit",
                       "Laptm5","Lfng","Lgals1","Zbtb7a","Lyl1",
                       "Mfng","Mmp14","Msh2","Myb","Nfatc1",
                       "Nkx2-3","Notch2","Ntrk1","Enpp1","Pik3r1",
                       "Pou1f1","Pnp","Pou2af1","Pou2f2","Prkdc",
                       "Ptk2b","Ptpn2","Ptprc","Ptprj","Rag1",
                       "Rag2","Rbpj","Spi1","Sfrp1","St3gal1",
                       "Sp3","Stat5a","Stat5b","Syk","Dock10",
                       "Tcf3","Fnip1","Tnfaip3","Cd27","Cd40lg",
                       "Top2b","Tpd52","Trp53","Tshr","Xbp1",
                       "Plcl2","Yy1","Dclre1c","Ikzf1","Ikzf3",
                       "Ifna13","Ifna16","Ifne","Pcid2","Ankle1",
                       "Plcg2","Muc19","Malt1","Tnfsf13b","Ifna15",
                       "Ifna12","Zbtb1","Irf2bp2","Tcirg1","Spib",
                       "Dnajb9","Ifnz","Atp11c","Dcaf1","Gpr183",
                       "Ep300","Ighe","Mir18","Mir20a","Mir150",
                       "Ifnk","Ifna14","Gm13271","Polm","Gm13283",
                       "Gm13272","Gm13276","Gm13277","Gm13275","Adgrg3",
                       "Lgals8","Hdac7","Gps2","Clcf1","Ppp2r3c",
                       "Il21","Rabl3","Mir19a","Mir17","Nfam1",
                       "Syvn1","Cyld","Slamf8","Mir19b-1","Mir92-1",
                       "Nhej1","Phf14","Dock11","Gon4l","6030468B19Rik",
                       "Hdac9","Nfkbiz","Tlr9"
                       
))

pbmc2 <- JoinLayers(pbmc2)
Idents(pbmc2) <- "sub_celltype"
m_sce <- AddModuleScore(pbmc2,features = Differentiation,
                        ctrl = 100,
                        name="Differentiation")
m_sce@meta.data <- m_sce@meta.data %>% rename("Differentiation1"="Differentiation")

library(ggpubr)
my_comparisons <- list(c("proB cell","preB cell"),c("preB cell","immuB cell"),c("immuB cell","mB cell"))

m_sce@meta.data$sub_celltype <- factor(m_sce@meta.data$sub_celltype, 
                                       levels = c("proB cell","preB cell","immuB cell","mB cell"))
p1 <- ggviolin(m_sce@meta.data,
               x="sub_celltype",
               y="Differentiation",
               width = 0.8,color = "black",
               fill="sub_celltype",
               xlab = F, 
               add = 'mean_sd', 
               bxp.errorbar=T,       
               bxp.errorbar.width=0.05,          
               size=0.5, 
               palette =c("#97A175","#A36F3D","#F2C847","#D3CBB8"), 
               legend = "right"
)+
  stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(angle=45))+
  ggtitle(label ="Differentiation")
p1
ggsave(p1,filename="Differentiation_score.pdf",width=5,height=4.5)

##Proliferation
Proliferation=list(c("Nckap1l","Ticam1","Card11","Abl1","Chrnb2",
                     "Pawr","Ada","Ahr","Tirap","Atm",
                     "Bax","Bcl2","Bcl6","Bmi1","Bst1",
                     "Btk","Casp3","Ctla4","Cd19","Cd22",
                     "Cd24a","Cd38","Cd79a","Cd81","Cdkn1a",
                     "Cdkn2a","Cr2","Ephb2","Fcgr2b","Fosl2",
                     "Cfb","Hspd1","Ifna1","Ifna11","Ifna2",
                     "Ifna4","Ifna5","Ifna6","Ifna7","Ifna9",
                     "Ifnab","Ifnb1","Ighm","Cd74","Il10",
                     "Il13","Il2","Il3","Il4","Il5",
                     "Il7","Il7r","Il9","Il9r","Inpp5d",
                     "Lef1","Cd180","Lyn","Mef2c","Mif",
                     "Nfatc2","Prkcd","Prlr","Pten",
                     "Ptprc","Rag2","Rasgrp1","Btla","Tcf3",
                     "Cd300a","Tlr4","Cd40","Cd27","Cd40lg",
                     "Cd70","Tfrc","Tsc2","Tnfrsf4","Tyrobp",
                     "Vpreb1a","Vpreb1b","Wnt3a","Plcl2",
                     "Slc39a10","Ikzf3","Shb","Ifna13","Ifna16",
                     "Ifne","Atad5","Gapt","Muc19","Tnfsf13b",
                     "Ifna15","Ifna12","Siglecg","Ifnz","Pkn1",
                     "Gpr183","Ighe","Ighd","Rc3h1","Irs2",
                     "Ifnk","Ifna14","Gm13271","Ctps1","Cd320",
                     "Gm13283","Gm13272","Gm13276","Gm13277","Gm13275",
                     "Clcf1","Vav3","Tnfrsf13b","Il21","Peli1",
                     "Mzb1","Tnfrsf13c","Sash3","Nfkbiz","Tlr9",
                     "Tnfrsf21"
))

pbmc2 <- JoinLayers(pbmc2)
Idents(pbmc2) <- "sub_celltype"
m_sce <- AddModuleScore(pbmc2,features = Proliferation,
                        ctrl = 100,
                        name="Proliferation")
m_sce@meta.data <- m_sce@meta.data %>% rename("Proliferation1"="Proliferation")

library(ggpubr)
my_comparisons <- list(c("proB cell","preB cell"),c("preB cell","immuB cell"),c("immuB cell","mB cell"))

m_sce@meta.data$sub_celltype <- factor(m_sce@meta.data$sub_celltype, 
                                       levels = c("proB cell","preB cell","immuB cell","mB cell"))
p1 <- ggviolin(m_sce@meta.data,
               x="sub_celltype",
               y="Proliferation",
               width = 0.8,color = "black",
               fill="sub_celltype", 
               xlab = F,  
               add = 'mean_sd', 
               bxp.errorbar=T,       
               bxp.errorbar.width=0.05,         
               size=0.5, 
               palette =c("#97A175","#A36F3D","#F2C847","#D3CBB8"), 
               legend = "right"
)+
  stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(angle=45))+
  ggtitle(label ="Proliferation")
p1
ggsave(p1,filename="Proliferation_score.pdf",width=5,height=4.5)


##Monocle
library(magrittr)
library(RColorBrewer)
library(reshape2)
library(Biobase)
library(ggsci)
library(ggpubr)
library(data.table)
library(monocle)

pbmc2 <- JoinLayers(pbmc2)
expr_matrix <- pbmc2@assays$RNA@layers$counts
sample_sheet <- pbmc2@meta.data
gene_annoation=data.frame(gene_short_name=rownames(pbmc1))
rownames(gene_annoation)=rownames(pbmc2)
pd <- new("AnnotatedDataFrame",data=sample_sheet)
fd <- new("AnnotatedDataFrame",data=gene_annoation)
cds <- newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,
                      expressionFamily = negbinomial.size())
save(cds,cds2,file = "B_mono.RData")
cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds,min_expr=0.1)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))

##方法一
diff_celltype <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~sub_celltype", cores=4)
head(diff_celltype)
write.csv(diff_celltype,"diff_celltype.csv")
diff_celltype <- diff_celltype[order(diff_celltype$qval),]
ordering_genes <- row.names(diff_celltype[1:1000,])
cds2 <- setOrderingFilter(cds,ordering_genes = ordering_genes)
plot_ordering_genes(cds2)
cds2 <- reduceDimension(cds2,max_components=2,method="DDRTree")
cds2 <- orderCells(cds2,reverse = T)
p1 <- plot_cell_trajectory(cds2,color_by ="State" ,cell_size =0.8 )+
  theme(text = element_text(size = 18))
p1
p2 <- plot_cell_trajectory(cds2,color_by = "Pseudotime" ,cell_size= 0.8)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#F7DCD0","#CF443B","#80181D"))
p2
ggsave(p2,filename="B_mon_pse.pdf",width = 3,height =3.5)

cds2$sub_celltype <- factor(cds2$sub_celltype,levels = c("proB cell", "preB cell","immuB cell","mB cell"))
p3 <- plot_cell_trajectory(cds2,color_by="sub_celltype",cell_size=0.8)+
  scale_color_manual(values=color1)+
  theme(legend.position = "top")+
  guides(color = guide_legend(nrow = 2))

p3
ggsave(p3,filename="B_mon_sub.pdf",width = 3.5,height =4)
p4 <- plot_cell_trajectory(cds2,color_by ="orig.sample" ,cell_size=0.8)+
  theme(text = element_text(size = 18))


library(viridis)
gene_pse <- differentialGeneTest(cds2[expressed_genes[1:1000],],
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                 cores=4)
gene_pse <- gene_pse[order(gene_pse$qval),]
gene_pse1 <- gene_pse[1:20,]

p <- plot_pseudotime_heatmap(cds2[rownames(gene_pse1),],
                             num_clusters = 4,
                             cores=4,
                             show_rownames = T,
                             return_heatmap = T,
                             hmcols = colorRampPalette(rev(brewer.pal(9,"PRGn")))(100)
                             
)

ggsave(p,filename="pse_heatmap.pdf",width = 4,height =4)
save(cds2,diff_celltype,gene_pse,top20g,file = "monocle.RData")


## gene correlation
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(ggstatsplot)
library(corrplot)
setwd("H:/HBOC/data1/B")
cluster.averages <- AverageExpression(pbmc2, group.by="orig.ident");

expr<-as.matrix(cluster.averages[["RNA"]])
head(expr)
data <- t(expr)
write.csv(data,"expr.csv")

dat1 <- read.csv("fil.csv",row.names = 1)

addcol <- colorRampPalette(c("CornflowerBlue","White","tomato"))
corr <- cor(dat1, method="pearson")
testRes <- cor.mtest(dat1, method="pearson",conf.level=0.95)
corrplot(corr, method = "circle",     
         type="lower",        
         col = addcol(100),        
         tl.col = "black",       
         tl.cex = 0.8,        
         tl.srt = 90,        
         tl.pos = "lt",        
         cl.pos = "r",         
         p.mat = testRes$p,         
         diag = F,         
         sig.level = c(0.0001,0.001, 0.01, 0.05),      
         pch.cex = 0.8,     
         insig = 'label_sig',         
         family="serif")
corrplot(corr, method = "circle", 
         col = addcol(100),          
         tl.col = "black", 
         tl.cex = 0.8, 
         tl.srt = 45,
         tl.pos = "lt",         
         p.mat = testRes$p, 
         diag = T, type = 'lower',         
         sig.level = c(0.0001,0.001, 0.01, 0.05), 
         pch.cex = 1.2,         
         insig = 'label_sig', pch.col = 'grey20', order = 'original')
corrplot(corr, method = "number", 
         type = "upper",col = addcol(100),          
         tl.col = "n", 
         tl.cex = 0.8,
         cl.pos = 'n', 
         tl.pos = "n",
         order = 'original',
         add=T,insig='blank')

save("corr.pdf", width=5,height=10)