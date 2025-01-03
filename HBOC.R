##
setwd("H:/HBOC")
library(Seurat)
library(harmony)
library(tidyverse)
library(dplyr)
library(tidydr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(clustree)

counts1 <- Read10X("data/2404192_C1/filtered_feature_bc_matrix")
pbmc1 <- CreateSeuratObject(counts1, min.cells = 3, min.features = 300,project = 'C1')
counts2 <- Read10X("data/2404192_C2/filtered_feature_bc_matrix")
pbmc2 <- CreateSeuratObject(counts2, min.cells = 3, min.features = 300,project = 'C2')
counts3 <- Read10X("data/2404192_C3/filtered_feature_bc_matrix")
pbmc3 <- CreateSeuratObject(counts3, min.cells = 3, min.features = 300,project = 'C3')
counts4 <- Read10X("data/2404192_H1/filtered_feature_bc_matrix")
pbmc4 <- CreateSeuratObject(counts4, min.cells = 3, min.features = 300,project = 'H1')
counts5 <- Read10X("data/2404192_H2/filtered_feature_bc_matrix")
pbmc5 <- CreateSeuratObject(counts5, min.cells = 3, min.features = 300,project = 'H2')
counts6 <- Read10X("data/2404192_H3/filtered_feature_bc_matrix")
pbmc6 <- CreateSeuratObject(counts6, min.cells = 3, min.features = 300,project = 'H3')

pbmc <-  merge(x=pbmc1,y=c(pbmc2,pbmc3,pbmc4,pbmc5,pbmc6))
rm(pbmc1,pbmc2,pbmc3,pbmc4,pbmc5,pbmc6)
rm(counts1,counts2,counts3,counts4,counts5,counts6)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & nCount_RNA < 20000
               & percent.mt < 10)

pbmc<- NormalizeData(pbmc) 
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

ifnb1 <- RunHarmony(pbmc, group.by.vars = "group")
ElbowPlot(ifnb1, ndims = 30)
DimPlot(ifnb1, reduction = "harmony")

ifnb1 <- FindNeighbors(ifnb1, reduction = "harmony", dims = 1:15)
ifnb1 <- FindClusters(ifnb1, resolution = 1)
ifnb1 <- RunUMAP(ifnb1, dims = 1:15, reduction = "harmony")
DimPlot(ifnb1, reduction = "umap",label=T)
DimPlot(ifnb1, reduction = "umap", group.by = "group")
DimPlot(ifnb, reduction = "umap", split.by = "orig.ident")

Idents(ifnb1) <- "seurat_clusters"
new.cluster.ids <- c("Neutrophil", "NMP", "B cell", "Neutrophil", "B cell", "Neutrophil",
                     "Neutrophil", "Neutrophil", "Monocyte", "B cell", "Neutrophil", 
                     "T cell", "B cell","DC", "Erythrocyte", "B cell",
                     "Neutrophil", "Erythrocyte","Macrophage", "Monocyte", "Monocyte",
                     "B cell", "B cell","Macrophage", "Erythrocyte", "HSC",
                     "Neutrophil", "T cell","Mast cell", "MSC", "B cell",
                     "Erythrocyte", "Monocyte")
names(new.cluster.ids) <- levels(ifnb1)
ifnb1 <- RenameIdents(ifnb1, new.cluster.ids) 
DimPlot(ifnb1,label = T) +ggtitle("")
ifnb1$celltype <- ifnb1@active.ident

##UMAP
Idents(ifnb1) <- "celltype"
Idents(ifnb1) <- factor(Idents(ifnb1), levels = c("Neutrophil", "B cell","Monocyte", "T cell","DC","Erythrocyte","Macrophage","HSC","Mast cell","MSC"))
color=c("#E1CEE4","#E35A48","#7E2116","#8CCBB0","#3A81AA","#FCE770","#AE8867","#CDD86A","#FB9C42","#9966CC")
plot <- DimPlot(ifnb1, reduction = "umap", shuffle = T,cols = color) +coord_flip()
plot
ggsave(plot, filename = "umap.pdf", width=6, height = 5)
DimPlot(ifnb1, reduction = "umap", group.by ="group")+coord_flip()
DimPlot(ifnb1, reduction = "umap", split.by ="group",cols = color)+coord_flip()

ggsave("umap_group1.pdf", width=10, height = 5)

## top5 heatmap
setwd("H:/HBOC/data1")
DefaultAssay(pbmc) <- "RNA"
pbmc1 <- JoinLayers(pbmc1)
markers <- FindAllMarkers(pbmc, only.pos = T, logfc.threshold = 0.25)
write.csv(markers, "allmarkers_pos.csv", row.names = T)
marker.sig <- read.csv("allmarkers_pos_fil.csv")

all.markers =markers %>% 
  dplyr::select(gene, everything()) %>%  subset(p_val<0.05)

marker.sig <- all.markers %>% 
  mutate(Ratio = round(pct.1/pct.2,3),
         pct.fc = pct.1-pct.2) %>%
  filter(p_val_adj < 0.05, # 本条件为过滤统计学不显著的基因
         pct.1 > 0.5,
         pct.2 < 0.5,
         pct.fc > 0.2,
         avg_log2FC > 1)

top5 = marker.sig  %>% group_by(cluster) %>% top_n(n = 5, wt =avg_log2FC)
View (top10)
write.csv(top5, "top5.csv", row.names = T)
top5 <- read.csv("top5.csv", row.names = 1)
markers1 <- top5$gene

top5 <- as.data.frame(top5)

ifnb1 <- ifnb1@meta.data$celltype
markerdata <- ScaleData(ifnb1,feature=as.character(unique(top5$gene)),                       
                        assay = "RNA")
Idents(markerdata) <- "celltype"
Idents(markerdata) <- factor(Idents(markerdata), 
                             levels = c("Neutrophil", "B cell","Monocyte", "T cell","DC","Erythrocyte","Macrophage","HSC","Mast cell","MSC"))

color=c("#E1CEE4","#AE8867","#E35A48","#7E2116","#8CCBB0","#3A81AA","#FCE770","#FAEACE","#CDD86A","#FB9C42","#9966CC")

DoHeatmap(markerdata,
          feature=as.character(unique(top5$gene)),
          assay= "RNA",
          group.colors=color)+
  scale_fill_gradientn(colors=c("white","grey","firebrick3"))

ggsave("heatmap_celltype_.pdf", width = 12,height = 12)

## markers dot
setwd("H:/HBOC/数据1")
Idents(ifnb1) <- "celltype"
genes <- c("Cxcr2","Ly6g","Ltf",
           "Cd79a","Cd19","Ighm",
           "Fn1","Ccr2","F13a1",
           "Cd3d","Cd3e", "Nkg7",
           "Siglech","Irf8","Bst2",
           "Hba-a1", "Hbb-bs", "Hbb-bt",
           "C1qa", "Vcam1","Mrc1",
           "Cd34", "Adgrg1", "Cdk6",
           "Fcer1a","Gata2","Ms4a2",
           "Cxcl12","Col1a2","Lepr")

p2 <- DotPlot(ifnb1, features = genes)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x.top = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
  labs(x="Markers", y="Celltype")+guides(size=guide_legend(order=3))
p2
ggsave(p2,filename = "Dot.pdf",width = 9.5,height = 4)

##cell proportion
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
setwd("H:/HBOC/data1")
colourCount = length(unique(ifnb1@meta.data$celltype))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
celltype_colors <- getPalette(colourCount)
plotC <- table(ifnb1@meta.data$group, ifnb1@meta.data$celltype) %>% melt()
write.csv(plotC, "Cellp.csv", row.names = T)
colnames(plotC) <- c("Group", "CellType","Number")

plotC$CellType <- factor(plotC$CellType,ordered = T,levels = c("Neutrophil","B cell","Monocyte", "T cell","DC","Erythrocyte","Macrophage","HSC","Mast cell","MSC"))
p1 <- ggplot(data = plotC, aes(x = Group, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x=NULL,y=NULL)+
  scale_y_continuous(labels = percent)+  
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(linewidth = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6))
p1
ggsave(p1,filename = "cellproportion.pdf",width = 3.5,height = 5)


##
setwd("H:/HBOC/数据1")
ifnb1 <- JoinLayers(ifnb1)
group_deg <- FindMarkers(ifnb1,
                         threshold = 0.25,
                         only.pos = F,
                         ident.1 = "HBOC", ident.2 = "Control",
                         group.by = "group") %>%
  mutate(gene = rownames(.))

write.csv(group_deg,file = "group_deg.csv")

pbmc2 <- JoinLayers(pbmc2)
group_deg5 <- FindMarkers(pbmc2,
                          threshold = 0.25,
                          only.pos = F,
                          ident.1 = "HBOC", ident.2 = "Control",
                          group.by = "group") %>%
  mutate(gene = rownames(.))

write.csv(group_deg5,file = "group_deg5.csv")

##volcano plot
group_deg5 <- group_deg5 %>%
  mutate(Difference = pct.1 - pct.2) 

group_deg5_up <- group_deg5 %>%
  filter(Difference >= 0.3 & p_val_adj < 0.05) %>%
  filter(avg_log2FC >=2)

group_deg5 <- group_deg5 %>%
  mutate(Difference = pct.1 - pct.2) 
group_deg5$change=0

for (i in 1:nrow(group_deg5)) {
  if(group_deg5$avg_log2FC[i]>=1){group_deg5$change[i]="up"}
  else if(group_deg5$avg_log2FC[i]<=-1){group_deg5$change[i]="down"}
  else{group_deg5$change[i]="ns"}
}
write.csv(group_deg5,file = "group_deg51.csv")
p1 <- ggplot(group_deg5, aes(x=Difference, y=avg_log2FC) )+ 
  geom_point(size=1.2,aes(color=change)) + 
  scale_color_manual("change",labels=c("down","ns","up"),
                     values=c("#0d1b46", "grey", "tomato"))+
  geom_label_repel(data=subset(group_deg5, 
                               avg_log2FC >= 2.5 & Difference >= 0.5 & p_val_adj <= 0.005), 
                   aes(label=gene),  #添加label
                   color="black", #设置label中标签的颜色
                   segment.colour = "black",#设置label框的颜色
                   label.padding = 0.2, 
                   segment.size = 0.3,  #框的大小
                   size=4,
                   max.overlaps = 20) +  # 增加 max.overlaps 参数
  geom_label_repel(data=subset(group_deg5, 
                               avg_log2FC <= -1.5 & Difference <= -0.20 & p_val_adj <= 0.005), 
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
ggsave(p1,filename="five_vocano1.pdf",width=7,height=6)


##Go analysis
group_deg5_up$gene <- rownames(group_deg5_up)
ids=bitr(group_deg5_up$gene, 'SYMBOL', 'ENTREZID', 'org.Mm.eg.db')
group_deg5_up=merge(group_deg5_up, ids,  by.x='gene', by.y='SYMBOL')
head(group_deg5_up)
##entichGo
group_deg5_up <- group_deg5_up[order(group_deg5_up$avg_log2FC,decreasing = T),]
group_degs_list <- as.numeric(group_deg5_up$avg_log2FC)
names(group_degs_list) <- group_deg5_up$ENTREZID
head(group_degs_list)

group_de <- names(group_degs_list)[group_degs_list > 2]
head(group_de)

group_ego <- enrichGO(group_de, OrgDb = "org.Mm.eg.db", ont = "BP", readable = T)
head(group_ego)
write.csv(group_ego, "group5_ego_up.csv")

group_ego1 <- read.csv("fil.csv")

p1 <- ggplot(group_ego1, aes(y=Description,x=Count)) + 
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=p.adjust),stat="identity") + 
  scale_fill_gradient2(expression(p.adjust),low="#80181D", mid="#CF443B", high="#F7DCD0",midpoint = 7.5e-10)+
  ylab(NULL) + xlab("Gene count") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank() )
p1
ggsave("go_bar1.pdf",width=6,height=3)


