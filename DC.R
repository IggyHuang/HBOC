##DC
Idents(ifnb1) <- "celltype"
pbmc_Dc <- subset(ifnb1, idents = "DC")

pbmc_Dc<- NormalizeData(pbmc_Dc) 
pbmc_Dc <- FindVariableFeatures(pbmc_Dc)
pbmc_Dc <- ScaleData(pbmc_Dc)
pbmc_Dc <- RunPCA(pbmc_Dc, features = VariableFeatures(object = pbmc_Dc))
ElbowPlot(pbmc_Dc,ndims = 30)

pbmc1  <- FindNeighbors(pbmc_Dc, dims = 1:15)
pbmc1 <- FindClusters(pbmc1, resolution = 0.4) 
pbmc1 <- RunUMAP(pbmc1, dims = 1:15)
DimPlot(pbmc1, reduction = "umap", shuffle = T, label = T)

pbmc2=pbmc1
Idents(pbmc2) <- "seurat_clusters"
new.cluster.ids <- c("pDC","cDC2", "cDC1")
names(new.cluster.ids) <- levels(pbmc2)
pbmc2 <- RenameIdents(pbmc2, new.cluster.ids) 
pbmc2$sub_celltype <- pbmc2@active.ident


##UMAP
color1=c("#79A5E0","#615F9C","#C9D9EB")
Idents(pbmc2) <- factor(Idents(pbmc2), levels = c("pDC","cDC1", "cDC2"))
plot <- print(DimPlot(pbmc2, reduction = "umap", shuffle = T,cols = color1, split.by = "group"))
plot
setwd("H:/HBOC/data1/DC")
ggsave(plot, filename = "umap_DC1.pdf", width=7, height = 5)

DimPlot(pbmc2, reduction = "umap", group.by ="group")
ggsave("umap_DC_group.pdf", width=6, height = 5)

##cell proportion
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
colourCount = length(unique(pbmc2@meta.data$celltype))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
celltype_colors <- getPalette(colourCount)
plotC <- table(pbmc2@meta.data$group, pbmc2@meta.data$sub_celltype) %>% melt()
colnames(plotC) <- c("Group", "CellType","Number")
write.csv(plotC, "Cellp_DC.csv", row.names = T)
plotC <- data.table::fread("Cellp_B.csv", header = T)

plotC$CellType <- factor(plotC$CellType,ordered = T,levels = c("pDC","cDC1", "cDC2"))
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
ggsave(pC2,filename = "cellp_DC.pdf",width = 3.5,height = 5)

##marker dot
genes <- c("Tcf4","Runx2","Bcl11a","Siglech","Ctsl",
           "Cbfa2t3","Pa2g4","Arsb","Ifi205","Clec9a",
           "Cd209a","Cd300a","Itgb7","Tmem176b")

p2 <- DotPlot(pbmc2, features = genes)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x.top = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
  labs(x="Markers", y="Celltype")+guides(size=guide_legend(order=3))          
p2
ggsave(p2, filename = "Dot_Dc.pdf", width=6, height = 3.6)


p1 <- DotPlot(pbmc2, features = gene)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x.top = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#F7DCD0", "#CF443B","#80181D"))+#颜色渐变设置  
  labs(x="Markers", y="Celltype")+guides(size=guide_legend(order=3))+
  coord_flip()
p1  
ggsave(p1,filename="genedot1_Dc.pdf",width=4.5,height=5)


##top5 heatmap
pbmc_DC=pbmc2
pbmc2 <- JoinLayers(pbmc2)
markers <- FindAllMarkers(pbmc2, only.pos = T, logfc.threshold = 0.25)
setwd("H:/HBOC/data1/DC")

write.csv(markers, "allmarkers_pos.csv", row.names = T)
markers <- read.csv("allmarkers_pos.csv",row.names = 1)
marker.sig1 <- markers %>% 
  mutate(pct.fc = pct.1-pct.2) %>%
  filter(p_val_adj < 0.05, 
         pct.fc > 0.2,
         avg_log2FC > 1)
write.csv(marker.sig1, "allmarkers_pos_fil.csv", row.names = T)
marker.sig <- read.csv("allmarkers_pos_fil.csv",row.names = 1)

top5 = marker.sig  %>% group_by(cluster) %>% top_n(n = 5, wt =avg_log2FC)
View (top5)

markerdata <- ScaleData(pbmc2,feature=as.character(unique(top5$gene)), 
                        assay = "RNA")

Idents(markerdata) <- "sub_celltype"
Idents(markerdata) <- factor(Idents(markerdata), 
                             levels = c("pDC","cDC1", "cDC2"))
DoHeatmap(markerdata,
          feature=as.character(unique(top5$gene)),
          assay= "RNA",
          group.colors=color1)+
  scale_fill_gradientn(colors=c("white","grey","firebrick3"))
ggsave("heatmap_celltype.pdf", width = 6,height =7.2)


pbmc3 <- JoinLayers(pbmc3)
Dc_group <- FindMarkers(pbmc3,
                        threshold = 0.25,
                        only.pos = F,
                        ident.1 = "HBOC", ident.2 = "Control",
                        group.by = "group") %>%
  mutate(gene = rownames(.))


##volcano plot
setwd("H:/HBOC/data1/DC/pDC")
write.csv(Dc_group,file = "Dc_group.csv")

library(ggrepel)
Dc_group <- Dc_group %>%
  mutate(Difference = pct.1 - pct.2) 
Dc_group$change=0

for (i in 1:nrow(Dc_group)) {
  if(Dc_group$avg_log2FC[i]>=1){Dc_group$change[i]="up"}
  else if(Dc_group$avg_log2FC[i]<=-1){Dc_group$change[i]="down"}
  else{Dc_group$change[i]="ns"}
}
write.csv(Dc_group,file = )

p1 <- ggplot(Dc_group, aes(x=Difference, y=avg_log2FC) )+ 
  geom_point(size=1.2,aes(color=change)) + 
  scale_color_manual("change",labels=c("down","ns","up"),
                     values=c("#0d1b46", "grey", "tomato"))+
  geom_label_repel(data=subset(Dc_group, 
                               avg_log2FC >= 2 & Difference >= 0.4 & p_val_adj <= 0.05), 
                   aes(label=gene),  
                   color="black", 
                   segment.colour = "black",
                   label.padding = 0.2, 
                   segment.size = 0.3, 
                   size=4,
                   max.overlaps = 26) + 
  geom_label_repel(data=subset(Dc_group, 
                               avg_log2FC <= -1 & Difference <= -0.25 & p_val_adj <= 0.05), 
                   aes(label=gene), 
                   label.padding = 0.2, 
                   color="black",
                   segment.colour = "black",
                   segment.size = 0.3, size=4,
                   max.overlaps = 30) + 
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

Dc_group_up <- Dc_group %>%
  filter(Difference >= 0.3 & p_val_adj < 0.05) %>%
  filter(avg_log2FC >=1)

Dc_group_up$gene <- rownames(Dc_group_up)
ids=bitr(Dc_group_up$gene, 'SYMBOL', 'ENTREZID', 'org.Mm.eg.db')
Dc_group_up=merge(Dc_group_up, ids,  by.x='gene', by.y='SYMBOL')
head(Dc_group_up)
##entichGo
Dc_group_up <- Dc_group_up[order(Dc_group_up$avg_log2FC,decreasing = T),]
group_degs_list <- as.numeric(Dc_group_up$avg_log2FC)
names(group_degs_list) <- Dc_group_up$ENTREZID
head(group_degs_list)

group_de <- names(group_degs_list)[group_degs_list > 2]
head(group_de)

group_ego <- enrichGO(group_de, OrgDb = "org.Mm.eg.db", ont = "BP", readable = T)
head(group_ego)
write.csv(group_ego, "Dc_group_ego_up.csv")
group_ego1 <- read.csv("fil.csv",row.names = 1)

ggplot(group_ego1, aes(y=Description,x=Count)) + 
  geom_bar(aes(y=Description,x=Count,fill=p.adjust),stat="identity") + 
  scale_fill_gradient2(expression(p.adjust),low="#80181D", mid="#CF443B", high="#F7DCD0",midpoint =4e-04)+
  ylab(NULL) + xlab("Gene count") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank() )+
  scale_y_discrete(limits=group_ego1$Description)

ggsave("go_bar.pdf",width=5.5,height=3)


##KEGG
group_ekg_up <- enrichKEGG(gene = group_de, organism = "mmu", pvalueCutoff = 0.05)
head(group_ekg_up)
write.csv(group_ekg_up, "Dc_group_kegg_up.csv")

group_kegg1 <- read.csv("fil1.csv",row.names = 1)

ggplot(group_kegg1, aes(y=Description,x=Count)) + 
  geom_bar(aes(y=Description,x=Count,fill=p.adjust),stat="identity") + 
  scale_fill_gradient2(expression(p.adjust),low="#80181D", mid="#CF443B", high="#F7DCD0",midpoint =0.008)+
  ylab(NULL) + xlab("Gene count") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() )+
  scale_y_discrete(limits=group_kegg1$Description)

ggsave("kegg_bar.pdf",width=5.5,height=3)


##gene correlation
setwd("H:/HBOC/data1/DC")
cluster.averages <- AverageExpression(pbmc2, group.by="orig.ident");

expr<-as.matrix(cluster.averages[["RNA"]])
head(expr)
data <- t(expr)
write.csv(data,"expr.csv")

dat1 <- read.csv("fil2.csv",row.names = 1)

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