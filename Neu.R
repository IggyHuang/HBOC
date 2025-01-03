##
Idents(ifnb1) <- "celltype"
pbmc_Neu <- subset(ifnb1, idents = "Neutrophil")
pbmc_Neu<- NormalizeData(pbmc_Neu) 
pbmc_Neu <- FindVariableFeatures(pbmc_Neu)
pbmc_Neu <- ScaleData(pbmc_Neu)
pbmc_Neu <- RunPCA(pbmc_Neu, features = VariableFeatures(object = pbmc_Neu))
ElbowPlot(pbmc_Neu,ndims = 30)

pbmc_N1  <- FindNeighbors(pbmc_Neu, dims = 1:15)
pbmc_N1 <- FindClusters(pbmc_N1, resolution = 0.4) 
pbmc_N1 <- RunUMAP(pbmc_N1, dims = 1:15)
DimPlot(pbmc_N1, reduction = "umap", shuffle = T,label = T)


##
pbmc_N2=pbmc_N2
Idents(pbmc_N2) <- "seurat_clusters"
new.cluster.ids <- c("mNeu_b","mNeu_a", "mNeu_b","immNeu", "immNeu", "mNeu_a","preNeu","immNeu","immNeu")
names(new.cluster.ids) <- levels(pbmc_N2)
pbmc_N2 <- RenameIdents(pbmc_N2, new.cluster.ids) 
DimPlot(pbmc_N2,label = T) +ggtitle("")
pbmc_N2$sub_celltype <- pbmc_N2@active.ident


##UMAP
color1=c("#636863","#8B4F8C","#BD7884","#AEBBDA","#D1A6C2")
Idents(pbmc_N1) <- "sub_celltype"
Idents(pbmc_N1) <- factor(Idents(pbmc_N1), levels = c("proNeu","preNeu","immNeu","mNeu_a","mNeu_b"))
plot <- print(DimPlot(pbmc_N1, reduction = "umap", shuffle = T,cols = color1))
plot
ggsave(plot, filename = "umap_Neu1.pdf", width=6, height = 5)
DimPlot(pbmc_N1, reduction = "umap", group.by ="group")
ggsave("umap_Neu_group1.pdf", width=6, height = 5)


##marker Dot
Idents(pbmc_N1) <- "sub_celltype"
Idents(pbmc_N1) <- factor(Idents(pbmc_N1), levels = c("proNeu","preNeu","immNeu","mNeu_a","mNeu_b"))
genes <- c("Elane","Mpo","Prtn3",
           "Chil3","Tuba1b","Fcnb",
           "Ltf","Ngp","Camp",
           "Ccl6","Gm16556", "Stfa2l1",
           "Isg15","Rsad2", "Ifit3")

p1 <- DotPlot(pbmc_N1, features = genes)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x.top = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
  labs(x="Markers", y="Celltype")+guides(size=guide_legend(order=3))
p1
ggsave(p1, filename = "Dot_Neu1.pdf", width=7, height = 4)

##
gene <- c("Isg15","Irf7","Zbp1","Ifih1","Oas3",
          "Oasl1","Rsad2","Ddx60","Trim30a")

p1 <- DotPlot(pbmc1, features = gene)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x.top = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#F7DCD0", "#CF443B","#80181D"))+#颜色渐变设置  
  labs(x="Markers", y="Celltype")+guides(size=guide_legend(order=3))
p1  
ggsave(p1,filename="genedot2.pdf",width=6.5,height=3)


##cell proportion
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)

pbmc1=pbmc_N1
colourCount = length(unique(pbmc1@meta.data$celltype))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
celltype_colors <- getPalette(colourCount)
plotC <- table(pbmc1@meta.data$group, pbmc1@meta.data$sub_celltype) %>% melt()
colnames(plotC) <- c("Group", "CellType","Number")
write.csv(plotC, "Cellp_Neu1.csv", row.names = T)
plotC <- data.table::fread("Cellp_Neu.csv", header = T)

plotC$CellType <- factor(plotC$CellType,ordered = T,levels = c("proNeu","preNeu","immNeu","mNeu_a","mNeu_b"))
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
ggsave(pC2,filename = "cellp_Neu1.pdf",width = 3.5,height =5)


##top5 heatmap
setwd("H:/HBOC/data1/Neu")
pbmc1=pbmc_N2
DefaultAssay(pbmc) <- "RNA"
pbmc1 <- JoinLayers(pbmc1)
markers <- FindAllMarkers(pbmc1, only.pos = T, logfc.threshold = 0.25)
write.csv(markers, "allmarkers_pos1.csv", row.names = T)
markers <- read.csv("allmarkers_pos.csv",row.names = 1)

marker.sig <- markers %>% 
  mutate(pct.fc = pct.1-pct.2) %>%
  filter(p_val_adj < 0.05, 
         pct.1 > 0.5,
         pct.2 < 0.5,
         pct.fc > 0.2,
         avg_log2FC > 1)
write.csv(marker.sig, "allmarkers_pos_fil1.csv", row.names = T)
marker.sig <- read.csv("allmarkers_pos_fil.csv",row.names = 1)

top5 = marker.sig  %>% group_by(cluster) %>% top_n(n = 5, wt =avg_log2FC)

write.csv(top5,file = "top5_1.csv")
top5 <- read.csv("top5_1.csv",row.names = 1)

markerdata <- ScaleData(pbmc1,feature=as.character(unique(top5$gene)), 
                        assay = "RNA")
Idents(markerdata) <- "sub_celltype"
Idents(markerdata) <- factor(Idents(markerdata), 
                             levels = c("proNeu","preNeu","immNeu","mNeu_a","mNeu_b"))
DoHeatmap(markerdata,
          feature=as.character(unique(top5$gene)),
          assay= "RNA",
          group.colors=color1)+
  scale_fill_gradientn(colors=c("white","grey","firebrick3"))
ggsave("heatmap_celltype1.pdf", width = 6,height =7.2)


##volcano plot
pbmc1=pbmc_N2
rm(pbmc_N2)
pbmc1 <- JoinLayers(pbmc1)
Neu_group <- FindMarkers(pbmc1,
                         threshold = 0.25,
                         only.pos = F,
                         ident.1 = "HBOC", ident.2 = "Control",
                         group.by = "group") %>%
  mutate(gene = rownames(.))

library(ggrepel)
Neu_group <- Neu_group %>%
  mutate(Difference = pct.1 - pct.2) 
Neu_group$change=0

for (i in 1:nrow(Neu_group)) {
  if(Neu_group$avg_log2FC[i]>=1){Neu_group$change[i]="up"}
  else if(Neu_group$avg_log2FC[i]<=-1){Neu_group$change[i]="down"}
  else{Neu_group$change[i]="ns"}
}
write.csv(Neu_group,file = "Neu_group1.csv")

p1 <- ggplot(Neu_group, aes(x=Difference, y=avg_log2FC) )+ 
  geom_point(size=1.2,aes(color=change)) + 
  scale_color_manual("change",labels=c("down","ns","up"),
                     values=c("#0d1b46", "grey", "tomato"))+
  geom_label_repel(data=subset(Neu_group, 
                               avg_log2FC >= 3 & Difference >= 0.5 & p_val_adj <= 0.005), 
                   aes(label=gene),  #添加label
                   color="black", #设置label中标签的颜色
                   segment.colour = "black",#设置label框的颜色
                   label.padding = 0.2, 
                   segment.size = 0.3,  #框的大小
                   size=4,
                   max.overlaps = 20) +  # 增加 max.overlaps 参数
  geom_label_repel(data=subset(Neu_group, 
                               avg_log2FC <= -3 & Difference <= -0.25 & p_val_adj <= 0.005), 
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
  xlim(-0.5,1.2)+
  ylim(-5,8)+
  ggtitle(label ="HBOC vs Control")
p1
ggsave(p1,filename="Neu_vocano1.pdf",width=7,height=6)


##GO/KEGG
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)

write.csv(Neu_group,"Neu_group.csv")
Neu_group_up <- Neu_group %>%
  filter(Difference >= 0.3 & p_val_adj < 0.05) %>%
  filter(avg_log2FC >=2)

Neu_group_up$gene <- rownames(Neu_group_up)
ids=bitr(Neu_group_up$gene, 'SYMBOL', 'ENTREZID', 'org.Mm.eg.db')
Neu_group_up=merge(Neu_group_up, ids,  by.x='gene', by.y='SYMBOL')
head(Neu_group_up)
##entichGo
Neu_group_up <- Neu_group_up[order(Neu_group_up$avg_log2FC,decreasing = T),]
group_degs_list <- as.numeric(Neu_group_up$avg_log2FC)
names(group_degs_list) <- Neu_group_up$ENTREZID
head(group_degs_list)

group_de <- names(group_degs_list)[group_degs_list > 2]
head(group_de)

group_ego <- enrichGO(group_de, OrgDb = "org.Mm.eg.db", ont = "BP", readable = T)
head(group_ego)
write.csv(group_ego, "neu_group_ego_up.csv")
setwd("H:/HBOC/data1/Neu/go1")

group_ego1 <- read.csv("fil.csv")

p1 <- ggplot(group_ego1, aes(y=Description,x=Count)) + 
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=p.adjust),stat="identity") + 
  scale_fill_gradient2(expression(p.adjust),low="#80181D", mid="#CF443B", high="#F7DCD0",midpoint = 0.0010)+
  ylab(NULL) + xlab("Gene count") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank() )
p1
ggsave("go_bar1.pdf",width=6,height=3)


##KEGG
group_ekg_down <- enrichKEGG(gene = group_de, organism = "mmu", pvalueCutoff = 0.05)
head(group_ekg_down)
write.csv(group_ekg_up, "neu_group_kegg_down.csv")

group_kegg1 <- read.csv("fil.csv",row.names = 1)

ggplot(group_kegg1, aes(y=Description,x=Count)) + 
  geom_bar(aes(y=Description,x=Count,fill=p.adjust),stat="identity") + 
  scale_fill_gradient2(expression(p.adjust),low="#80181D", mid="#CF443B", high="#F7DCD0",midpoint = 0.025)+
  ylab(NULL) + xlab("Gene count") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() )+
  scale_y_discrete(limits=group_kegg1$Description)

ggsave("kegg_bar.pdf",width=5.5,height=3)


##violin plot of gene score
##maturation
pbmc1=pbmc_N1
rm(pbmc_N1)
Maturation=list(c("Retnlg","Ccl6","S100a6","Clec4d","Prr13",
                  "Cebpb","Slpi","S100a11","Btg1","Cxcr2",
                  "Fth1","Grina","Mmp8","Fxyd5","Msrb1",
                  "H2-D1","Anxa2","Mmp9","Ftl1","Map1lc3b",
                  "Tmcc1","Sat1","Cyp4f18","Junb","Mxd1",
                  "Stk17b","Ypel3","Selplg","Il1f9","Dusp1",
                  "Slc16a3","Ccr1","Rdh12","Clec4e","Arg2",
                  "Cd300ld","Ctsd","Gda","Hacd4",
                  "Timp2","Fpr1","Ifi27l2a","Slc7a11","Stfa2l1",
                  "Il1b","Asprv1","Cxcl2","Ifitm1"
))
pbmc1 <- JoinLayers(pbmc1)
Idents(pbmc1) <- "sub_celltype"
m_sce <- AddModuleScore(pbmc1,features = Maturation,
                        ctrl = 100,
                        name="Maturation")
m_sce@meta.data <- m_sce@meta.data %>% rename("Maturation1"="Maturation")
m_sce@meta.data$sub_celltype <- factor(m_sce@meta.data$sub_celltype, 
                                       c("preNeu","proNeu","immNeu","mNeu_a","mNeu_b"))
library(ggpubr)
ggviolin(m_sce@meta.data,
         x="sub_celltype",
         y="Maturation",
         width = 0.8,color = "black",
         fill="sub_celltype", 
         xlab = F,  
         add = 'mean_sd', 
         bxp.errorbar=T,        
         bxp.errorbar.width=0.05,           
         size=0.5, 
         palette =color1, 
         legend = "right")+
  stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(angle=45))+
  ggtitle(label ="Maturation")
ggsave("maturation_score2.pdf",width=4.5,height=4)

##phagocytosis
Phagocytosis <- list(c("Abca1","Adgrb1","Aif1","Arhgap12","Arhgap25",
                       "Becn1","Bin2","Cdc42","Clcn3","Elmo1",
                       "Fcer1g","Fcgr1","Fcgr3","Gsn","Gulp1",
                       "Igll1","Itgam","Itgb2","Marco","Megf10",
                       "Mfge8","Msr1","Myh9","Rac1","Rac3",
                       "Rhobtb1","Rhobtb2","Sh3bp1","Sirpa","Thbs1",
                       "Trem2","Treml4","Vamp7","Xkr4","Xkr6",
                       "Xkr7","Xkr8","Xkr9"))

Idents(pbmc1) <- "sub_celltype"
n_sce <- AddModuleScore(pbmc1,features = Phagocytosis,
                        ctrl = 100,
                        name="Phagocytosis")
n_sce@meta.data <- n_sce@meta.data %>% rename("Phagocytosis1"="Phagocytosis")
n_sce@meta.data$sub_celltype <- factor(n_sce@meta.data$sub_celltype, 
                                       c("proNeu","preNeu","immNeu","mNeu_a","mNeu_b"))

ggviolin(n_sce@meta.data,
         x="sub_celltype",
         y="Phagocytosis",
         width = 0.8,color = "black",
         fill="sub_celltype", 
         xlab = F,  
         add = 'mean_sd', 
         bxp.errorbar=T,         
         bxp.errorbar.width=0.05,           
         size=0.5, 
         palette =color1, 
         legend = "right")+
  stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(angle=45))+
  ggtitle(label ="Phagocytosis")
ggsave("Phagocytosis_score1.pdf",width=5,height=4)

##chemotaxis
Chemotaxis <- list(c("C5ar1","Ccl1","Ccl12","Ccl17",
                     "Ccl19","Ccl2","Ccl22","Ccl24","Ccl25",
                     "Ccl3","Ccl4","Ccl5","Ccl6","Ccl7",
                     "Ccl8","Ccl9","Cklf","Csf3r","Cx3cl1",
                     "Cxadr","Cxcl10","Cxcl15",
                     "Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1",
                     "Cxcr2","Edn3","Fcer1g","Fcgr3","Gbf1",
                     "Ifng","Il1b",
                     "Il1rn","Itga1","Itga9","Itgam","Itgb2",
                     "Lgals3","Nckap1l","Pde4b","Pde4d","Pf4",
                     "Pla2g1b","Ppbp","Prex1","Prkca","S100a8",
                     "S100a9","Slc37a4","Spp1","Syk","Tgfb2",
                     "Trem1","Trem3","Vav1","Vav3","Xcl1"
))

Idents(pbmc1) <- "sub_celltype"
c_sce <- AddModuleScore(pbmc1,features = Chemotaxis,
                        ctrl = 100,
                        name="Chemotaxis")
c_sce@meta.data <- c_sce@meta.data %>% rename("Chemotaxis1"="Chemotaxis")
c_sce@meta.data$sub_celltype <- factor(c_sce@meta.data$sub_celltype, 
                                       c("proNeu","preNeu","immNeu","mNeu_a","mNeu_b"))

ggviolin(c_sce@meta.data,
         x="sub_celltype",
         y="Chemotaxis",
         width = 0.8,color = "black",
         fill="sub_celltype",
         xlab = F, 
         add = 'mean_sd', 
         bxp.errorbar=T,       
         bxp.errorbar.width=0.05,           
         size=0.5, 
         palette =color1, 
         legend = "right")+
  stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(angle=45))+
  ggtitle(label ="Chemotaxis")
ggsave("Chemotaxis_score.pdf",width=5,height=4)


#### monocle2
library(magrittr)
library(RColorBrewer)
library(reshape2)
library(Biobase)
library(ggsci)
library(ggpubr)
library(data.table)
library(monocle)

pbmc1 <- JoinLayers(pbmc1)
expr_matrix <- pbmc1@assays$RNA@layers$counts
sample_sheet <- pbmc1@meta.data
gene_annoation=data.frame(gene_short_name=rownames(pbmc1))
rownames(gene_annoation)=rownames(pbmc1)
pd <- new("AnnotatedDataFrame",data=sample_sheet)
fd <- new("AnnotatedDataFrame",data=gene_annoation)
cds <- newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,
                      expressionFamily = negbinomial.size())

cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds,min_expr=0.1)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))


diff_celltype <- differentialGeneTest(cds1[expressed_genes,],
                                      fullModelFormulaStr = "~sub_celltype", 
                                      cores=4)
head(diff_celltype)
write.csv(diff_celltype,"diff_celltype1.csv")
diff_celltype <- diff_celltype[order(diff_celltype$qval),]
ordering_genes <- row.names(diff_celltype[1:1000,])
cds1 <- setOrderingFilter(cds,ordering_genes = ordering_genes)
plot_ordering_genes(cds1)
cds1 <- reduceDimension(cds1,max_components=2,method="DDRTree")
cds1 <- orderCells(cds1,reverse=F)

cds1 <- orderCells(cds1,reverse=T)

p1 <- plot_cell_trajectory(cds1,color_by ="State" ,cell_size =0.8 )
p1
p2 <- plot_cell_trajectory(cds1,color_by = "Pseudotime" ,cell_size= 0.8)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#F7DCD0","#CF443B","#80181D"))
p2
ggsave(p2,filename="Neu_mono_pse2.pdf",width = 6,height =6)

##heatmap
library(viridis)
gene_pse <- differentialGeneTest(cds1[var.genes,],
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                 cores=4)
gene_pse <- gene_pse[order(gene_pse$qval),]
gene_pse1 <- gene_pse[1:20,]

p5 <- plot_pseudotime_heatmap(cds1[rownames(gene_pse1),],
                              num_clusters = 6,
                              cores=4,
                              show_rownames = T,
                              hmcols = colorRampPalette(rev(brewer.pal(9,"PRGn")))(100)
)
p5

ggsave(p5,filename="pse_heatmap.pdf",width = 8,height =8)


##gene correlation
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(ggstatsplot)
library(corrplot)
setwd("H:/HBOC/data1/Neu/gsva")
cluster.averages <- AverageExpression(pbmc1, group.by="orig.ident");

expr<-as.matrix(cluster.averages[["RNA"]])
head(expr)
data <- t(expr)
write.csv(data,"expr.csv")

dat1 <- read.csv("fil.csv",row.names = 1)

addcol <- colorRampPalette(c("CornflowerBlue","White","tomato"))
corr <- cor(dat1, method="pearson")
testRes <- cor.mtest(dat1, method="pearson",conf.level=0.95)

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
