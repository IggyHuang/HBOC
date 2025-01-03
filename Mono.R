##
Idents(ifnb1) <- "celltype"
pbmc_Mo <- subset(ifnb1, idents = "Monocyte")

pbmc_Mo<- NormalizeData(pbmc_Mo) 
pbmc_Mo <- FindVariableFeatures(pbmc_Mo)
pbmc_Mo <- ScaleData(pbmc_Mo)
pbmc_Mo <- RunPCA(pbmc_Mo, features = VariableFeatures(object = pbmc_Mo))
ElbowPlot(pbmc_Mo,ndims = 30)

pbmc3  <- FindNeighbors(pbmc_Mo, dims = 1:15)
pbmc3 <- FindClusters(pbmc3, resolution = 0.4) 
pbmc3 <- RunUMAP(pbmc3, dims = 1:15)
DimPlot(pbmc1, reduction = "umap", shuffle = T, label = T)

pbmc4=pbmc3
new.cluster.ids <- c("Non-classical monocyte","Classical monocyte", "Monocyte progenitor","Non-classical monocyte","Classical monocyte", "Monocyte progenitor")
names(new.cluster.ids) <- levels(pbmc4)
pbmc4 <- RenameIdents(pbmc4, new.cluster.ids) 
pbmc4$sub_celltype <- pbmc4@active.ident

##marker dot
genes <- c( "Stmn1","Mki67",
            "Tmsb10","Crip1","Vcan",
            "Apoe","Cx3cr1","Hpgd")

Idents(pbmc4) <- factor(Idents(pbmc4), levels = c("Monocyte progenitor","Classical monocyte","Non-classical monocyte"))
p1 <- DotPlot(pbmc4, features = genes)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x.top = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.1),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
  labs(x="Markers", y="Celltype")+guides(size=guide_legend(order=3))
p1
ggsave(p1, filename = "Dot_Mono.pdf", width=7, height = 3)

##
genes <- c("Isg15","Ifi213","Ifi206","Ifit3","Ifit2",
           "Oasl2","Irf7","Oas3","Ifih1","Zbp1",
           "Apoe","Sulf2","Neurl3","Ltf","Cx3cr1")

p1 <- DotPlot(pbmc4, features = genes)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(angle = 45))+
  theme(axis.text.x.top = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#F7DCD0", "#CF443B","#80181D"))+#颜色渐变设置  
  labs(x="Markers", y="Celltype")+guides(size=guide_legend(order=3))
p1  
ggsave(p1, filename = "Dot_Mono2.pdf", width=8.5, height = 3.5)


##UMAP
color1=c("#8DDD7C","#D594E0","#A0C4EC")
plot <- print(DimPlot(pbmc4, reduction = "umap", split.by="group",shuffle = T,cols = color1))
plot
ggsave(plot, filename = "umap_Mo1.pdf", width=8, height = 3.5)

DimPlot(pbmc4, reduction = "umap", group.by ="group")
ggsave("umap_Mo_group.pdf", width=4.5, height = 3.5)


##cell proportion
colourCount = length(unique(pbmc4@meta.data$celltype))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
celltype_colors <- getPalette(colourCount)
plotC <- table(pbmc4@meta.data$group, pbmc4@meta.data$sub_celltype) %>% melt()
colnames(plotC) <- c("Group", "CellType","Number")
write.csv(plotC, "Cellp_Neu.csv", row.names = T)
plotC <- data.table::fread("Cellp_Neu.csv", header = T)

pC2 <- ggplot(data = plotC, aes(x = Group, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=color1) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x=NULL,y=NULL)+
  scale_y_continuous(labels = percent)+  ####用来将y轴移动位置
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6))
pC2
ggsave(pC2,filename = "cellp_Mo.pdf",width = 3.5,height = 3)


##top5 heatmap
pbmc1=pbmc4
rm(pbmc4)
pbmc1 <- JoinLayers(pbmc1)
markers <- FindAllMarkers(pbmc1, only.pos = T, logfc.threshold = 0.25)

write.csv(markers, "allmarkers_pos.csv", row.names = T)
markers <- read.csv("all_allmarkers_pos.csv",row.names = 1)
marker.sig <- markers %>% 
  mutate(pct.fc = pct.1-pct.2) %>%
  filter(p_val_adj < 0.05, # 本条件为过滤统计学不显著的基因
         pct.1 > 0.5,
         pct.2 < 0.5,
         pct.fc > 0.2,
         avg_log2FC > 1)
write.csv(marker.sig, "allmarkers_pos_fil.csv", row.names = T)

marker.sig <- read.csv("all_allmarkers_pos_fil.csv",row.names = 1)

top5 = marker.sig  %>% group_by(cluster) %>% top_n(n = 5, wt =avg_log2FC)
View (top5)


markerdata <- ScaleData(ifnb1,feature=as.character(unique(top5$gene)), 
                        assay = "RNA")

Idents(markerdata) <- "sub_celltype"
Idents(markerdata) <- factor(Idents(markerdata), 
                             levels = c("Monocyte progenitor","Classical monocyte","Non-classical monocyte"))
DoHeatmap(markerdata,
          feature=as.character(unique(top5$gene)),
          assay= "RNA",
          group.colors=color1)+
  scale_fill_gradientn(colors=c("white","grey","firebrick3"))
ggsave("heatmap_celltype.pdf", width = 6,height =7.2)


gene <- c("Esco2","Brca1","Fbxo5","Rrm2","Dtl",
          "Ifit2","Ifit3","Trim30c","Rsad2","Oasl1",
          "Apoe","Cx3cr1","Sirpb1c","Sirpb1b","Apoc2")

DoHeatmap(pbmc1,
          feature=gene,
          assay= "RNA",
          group.colors=color1)+
  scale_fill_gradientn(colors=c("white","grey","firebrick3"))
ggsave("heatmap1.pdf", width = 5.5,height =7)


##volcano plot
pbmc4 <- JoinLayers(pbmc4)
Mo_group <- FindMarkers(pbmc4,
                        threshold = 0.25,
                        only.pos = F,
                        ident.1 = "HBOC", ident.2 = "Control",
                        group.by = "group") %>%
  mutate(gene = rownames(.))

write.csv(Mo_group,file = "Mo_group.csv")

library(ggrepel)
Mo_group <- read.csv("Mo_group.csv")
Mo_group <- Mo_group %>%
  mutate(Difference = pct.1 - pct.2) 
Mo_group$change=0

for (i in 1:nrow(Mo_group)) {
  if(Mo_group$avg_log2FC[i]>=1){Mo_group$change[i]="up"}
  else if(Mo_group$avg_log2FC[i]<=-1){Mo_group$change[i]="down"}
  else{Mo_group$change[i]="ns"}
}

genes <- c("Isg15","Ifi213","Ifi206","Ifit3","Ifit2",
           "Oasl2","Irf7","Ifi205","Ifi47","Rtp4",
           "Apoe","Sulf2","Neurl3","Ltf","Cx3cr1"
           
)
p1 <- ggplot(Mo_group, aes(x=Difference, y=avg_log2FC) )+ 
  geom_point(size=1.2,aes(color=change)) + 
  scale_color_manual("change",labels=c("down","ns","up"),
                     values=c("#0d1b46", "grey", "tomato"))+
  geom_label_repel(data=subset(Mo_group, 
                               avg_log2FC >= 2 & Difference >= 0.6 & p_val_adj <= 0.05), 
                   aes(label=gene),  
                   color="black",
                   segment.colour = "black",
                   label.padding = 0.2, 
                   segment.size = 0.3, 
                   size=4,
                   max.overlaps = 30) + 
  geom_label_repel(data=subset(Mo_group, 
                               avg_log2FC <= -1 & Difference <= -0.2 & p_val_adj <= 0.05), 
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

##GO
Mo_group_up <- Mo_group %>%
  filter(Difference >= 0.5 & p_val_adj < 0.05) %>%
  filter(avg_log2FC >=2)

Mo_group_up$gene <- rownames(Mo_group_up)
ids=bitr(Mo_group_up$gene, 'SYMBOL', 'ENTREZID', 'org.Mm.eg.db')
Mo_group_up=merge(Mo_group_up, ids,  by.x='gene', by.y='SYMBOL')
head(Mo_group_up)

Mo_group_up <- Mo_group_up[order(Mo_group_up$avg_log2FC,decreasing = T),]
group_degs_list <- as.numeric(Mo_group_up$avg_log2FC)
names(group_degs_list) <- Mo_group_up$ENTREZID
head(group_degs_list)

group_de <- names(group_degs_list)[group_degs_list > 2]
head(group_de)

group_ego <- enrichGO(group_de, OrgDb = "org.Mm.eg.db", ont = "BP", readable = T)
head(group_ego)
write.csv(group_ego, "Mo_group_ego_up.csv")
setwd("H:/HBOC/data1")
group_ego1 <- read.csv("fil.csv",row.names = 1)

ggplot(group_ego1, aes(y=Description,x=Count)) + 
  geom_bar(aes(y=Description,x=Count,fill=p.adjust),stat="identity") + 
  scale_fill_gradient2(expression(p.adjust),low="#80181D", mid="#CF443B", high="#F7DCD0",midpoint =0.02 )+
  ylab(NULL) + xlab("Gene count") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank() )+
  scale_y_discrete(limits=group_ego1$Description)

ggsave("go_bar_up.pdf",width=5.5,height=3)

##KEGG
group_ekg_up <- enrichKEGG(gene = group_de, organism = "mmu", pvalueCutoff = 0.05)
head(group_ekg_up)
write.csv(group_ekg_up, "Mo_group_kegg_up.csv")
setwd()
group_kegg1 <- read.csv("fil.csv")

ggplot(group_kegg1, aes(y=Description,x=Count)) + 
  geom_bar(aes(y=Description,x=Count,fill=p.adjust),stat="identity") + 
  scale_fill_gradient2(expression(p.adjust),low="#80181D", mid="#CF443B", high="#F7DCD0",midpoint =0.02 )+
  ylab(NULL) + xlab("Gene count") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank() )+
  scale_y_discrete(limits=group_kegg1$Description)

ggsave("kegg_bar1.pdf",width=4.5,height=3)


##gene correlation
setwd("H:/HBOC/data1/Mono")
cluster.averages <- AverageExpression(pbmc4, group.by="orig.ident");

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
