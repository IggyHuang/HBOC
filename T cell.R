##T
Idents(ifnb1) <- "celltype"
pbmc_T <- subset(ifnb1, idents = "T cell")
rm(ifnb1)
pbmc_T=pbmc2
pbmc_T<- NormalizeData(pbmc_T) 
pbmc_T <- FindVariableFeatures(pbmc_T)
pbmc_T <- ScaleData(pbmc_T)
pbmc_T<- RunPCA(pbmc_T, features = VariableFeatures(object = pbmc_T))
ElbowPlot(pbmc_T,ndims = 30)

pbmc1 <- FindNeighbors(pbmc_T, dims = 1:15)
pbmc1 <- FindClusters(pbmc1, resolution = 0.2) 
pbmc1 <- RunUMAP(pbmc1, dims = 1:15)
plot <- print(DimPlot(pbmc1, reduction = "umap", shuffle = T, label = T))
DimPlot(pbmc1, reduction = "umap", group.by ="group")

pbmc2=pbmc1
Idents(pbmc2) <- "seurat_clusters"
new.cluster.ids <- c("Naive Cd8+ T cell","NKT cell", "NKT cell","NK cell", "Cd4+ Treg", "Proliferation T cell")
names(new.cluster.ids) <- levels(pbmc2)
pbmc2 <- RenameIdents(pbmc2, new.cluster.ids) 
pbmc2$sub_celltype <- pbmc2@active.ident


##umap
setwd("H:/HBOC/data1/T")
color1=c("#98DE84","#E7D70C","#45836B","#1D9DC3","#C7DEBD")
Idents(pbmc2) <- factor(Idents(pbmc2), levels = c("Naive Cd8+ T cell", "Cd4+ Treg","Proliferation T cell","NKT cell", "NK cell"))
plot <- print(DimPlot(pbmc2, reduction = "umap", shuffle = T,cols = color1, split.by = "group"))
plot
ggsave(plot, filename = "umap_T1.pdf", width=9, height = 5)

DimPlot(pbmc1, reduction = "umap", group.by ="group")
ggsave("umap_T_group.pdf", width=6.5, height = 5)


##marker dot
Idents(pbmc2) <- "sub_celltype"
Idents(pbmc2) <- factor(Idents(pbmc2), levels = c("Naive Cd8+ T cell", "Cd4+ Treg","Proliferation T cell","NKT cell", "NK cell" ))

genes <- c("Cd8a","Tcf7","Sell","Lef1",
           "Cd4","Foxp3","Il2ra",
           "Top2a","Mki67",
           "Ly6c2","Ccl5","Cd69",
           "Gzma","Klra8","Klra4")
p1 <- DotPlot(pbmc2, features = genes)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x.top = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
  labs(x="Markers", y="Celltype")+guides(size=guide_legend(order=3))          
p1
ggsave(p1, filename = "Dot_T.pdf", width=8, height = 3.5)

genes <- c("Ifit3","Isg15","Ifit1","Rtp4","Bst2","Irf7","Ifi213","Zbp1"
)

p1 <- DotPlot(pbmc2, features = genes)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x.top = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#F7DCD0", "#CF443B","#80181D"))+#颜色渐变设置  
  labs(x="Markers", y="Celltype")+guides(size=guide_legend(order=3))
p1  
ggsave(p1,filename="genedot1.pdf",width=7,height=3)


##cell proportion
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
setwd("H:/HBOC/data1/T/go")

colourCount = length(unique(pbmc2@meta.data$celltype))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
celltype_colors <- getPalette(colourCount)
plotC <- table(pbmc2@meta.data$group, pbmc2@meta.data$sub_celltype) %>% melt()
colnames(plotC) <- c("Group", "CellType","Number")
write.csv(plotC, "Cellp.csv", row.names = T)
plotC <- data.table::fread("Cellp_Neu.csv", header = T)

color2=c("#98DE84","#E7D70C","#45836B","#1D9DC3","#C7DEBD")
pC2 <- ggplot(data = plotC, aes(x = Group, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=color2) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x=NULL,y=NULL)+
  scale_y_continuous(labels = percent)+  
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6))
pC2
ggsave(pC2,filename = "cellp_T.pdf",width = 3.5,height = 4.5)


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
marker.sig <- read.csv("allmarkers_pos_fil.csv",row.names = 1)

marker.sig1 <- markers %>% 
  filter(pct.1 >0.1 & p_val_adj < 0.05) %>%
  filter(abs(avg_log2FC) > 0.5)

top5 = marker.sig  %>% group_by(cluster) %>% top_n(n = 5, wt =avg_log2FC)

markerdata <- ScaleData(pbmc1,feature=as.character(unique(top5$gene)), 
                        assay = "RNA")

Idents(markerdata) <- "sub_celltype"
Idents(markerdata) <- factor(Idents(markerdata), 
                             levels = c("Naive Cd8+ T cell", "Cd4+ Treg","Proliferation T cell","NKT cell", "NK cell"))
DoHeatmap(markerdata,
          feature=as.character(unique(top5$gene)),
          assay= "RNA",
          group.colors=color2)+
  scale_fill_gradientn(colors=c("white","grey","firebrick3"))
ggsave("heatmap_celltype.pdf", width = 6,height =7.2)


##volcano plot
pbmc2 <- JoinLayers(pbmc2)
group_degs_T = FindMarkers(pbmc2, logfc.threshold = 0.25,
                           only.pos = F,
                           ident.1 = "HBOC", ident.2 = "Control",
                           group.by = "group") %>%
  mutate(gene = rownames(.))
write.csv(group_degs_T, "group_degs_T.csv", row.names = T)
setwd("H:/HBOC/data1/T")
group_degs_T <- read.csv("group_degs_T.csv")
library(ggrepel)
group_degs_T <- group_degs_T %>%
  mutate(Difference = pct.1 - pct.2) 
group_degs_T$change=0

for (i in 1:nrow(group_degs_T)) {
  if(group_degs_T$avg_log2FC[i]>=1){group_degs_T$change[i]="up"}
  else if(group_degs_T$avg_log2FC[i]<=-1){group_degs_T$change[i]="down"}
  else{group_degs_T$change[i]="ns"}
}


p1 <- ggplot(group_degs_T, aes(x=Difference, y=avg_log2FC) )+ 
  geom_point(size=1.2,aes(color=change)) + 
  scale_color_manual("change",labels=c("down","ns","up"),
                     values=c("#0d1b46", "grey", "tomato"))+
  geom_label_repel(data=subset(group_degs_T, 
                               avg_log2FC >= 2 & Difference >= 0.4 & p_val_adj <= 0.05), 
                   aes(label=gene),  
                   color="black", 
                   segment.colour = "black",
                   label.padding = 0.2, 
                   segment.size = 0.3,  
                   size=4,
                   max.overlaps = 20) + 
  geom_label_repel(data=subset(group_degs_T, 
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

##GO
T_group_up <- group_degs_T %>%
  filter(Difference >= 0.3 & p_val_adj < 0.05) %>%
  filter(avg_log2FC >=2)

T_group_up$gene <- rownames(T_group_up)
ids=bitr(T_group_up$gene, 'SYMBOL', 'ENTREZID', 'org.Mm.eg.db')
T_group_up=merge(T_group_up, ids,  by.x='gene', by.y='SYMBOL')
head(T_group_up)

T_group_up <- T_group_up[order(T_group_up$avg_log2FC,decreasing = T),]
group_degs_list <- as.numeric(T_group_up$avg_log2FC)
names(group_degs_list) <- T_group_up$ENTREZID
head(group_degs_list)

group_de <- names(group_degs_list)[group_degs_list > 2]
head(group_de)

group_ego <- enrichGO(group_de, OrgDb = "org.Mm.eg.db", ont = "BP", readable = T)
head(group_ego)
write.csv(group_ego, "T_group_ego_up.csv")
setwd("H:/HBOC/data1/T")
group_ego_up <- read.csv("fil.csv",row.names = 1)

p1 <- ggplot(group_ego_up, aes(y=Description,x=Count)) + 
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=p.adjust),stat="identity") + 
  scale_fill_gradient2(expression(p.adjust),low="#80181D", mid="#CF443B", high="#F7DCD0",midpoint = 0.0010)+
  ylab(NULL) + xlab("Gene count") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank() )+
  scale_y_discrete(limits=group_ego_up$Description)

p1
ggsave("go_bar.pdf",width=5.5,height=3)

#KEGG
group_ekg_up <- enrichKEGG(gene = group_de, organism = "mmu", pvalueCutoff = 0.05)
head(group_ekg_up)
write.csv(group_ekg_up, "T_group_kegg_up.csv")

group_kegg_up <- read.csv("T_group_kegg_up.csv",row.names = 1)

p1 <- ggplot(group_kegg_up, aes(y=Description,x=Count)) + 
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=p.adjust),stat="identity") + 
  scale_fill_gradient2(expression(p.adjust),low="#80181D", mid="#CF443B", high="#F7DCD0",midpoint = 0.035)+
  ylab(NULL) + xlab("Gene count") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() )+
  scale_y_discrete(limits=group_kegg_up$Description)

p1
ggsave("kegg_bar.pdf",width=5.5,height=3)


##gene correaltion
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(ggstatsplot)
library(corrplot)
setwd("H:/HBOC/data1/T")
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