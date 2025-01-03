##GSE264189
setwd("H:/HBOC/data1/RNA/GSE264189")
data1 <- read.csv("GSE264189_normalized_read_count.csv",header = T,sep = "\t")
str(data1)
data2 <- data1[,c(1,2,3,4,11,12,13,20,21,22,27,28,29,36,37,38)]
write.csv(data2,file="data2.csv")
row.names(data2) <- data2$X
data2 <- data2[,-1]
data2 <- log2(data2)
write.csv(data2,file="data.csv")

gene1 <- read.csv("fil1.csv")
write.csv(data1,file="data1.csv")

gene1$Gene <- factor(gene1$Gene,
                     "Isg15")
p <- ggplot(data=gene1, aes(x=Type,y=Expression)) +
  geom_boxplot(aes(fill=Type)) +  
  labs(x = NULL, y="Normalized expression") +  
  scale_fill_manual(values=c("#72b1b9","#f6ae00")) +  
  ylim(3.5,9.5)+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5), 
        text = element_text(size = 15),  
        axis.text.x = element_text(angle = 45, hjust = 1),  
        legend.position="top", legend.title = element_blank())
p 

ggsave(p,filename="gene1.pdf",width=5,height=6) 


##GSE267631
setwd("H:/HBOC/data1/RNA/GSE267631")
data1 <- read.csv("GSE267631_gene_count.csv")
data2 <- data1[,c(2,3,4,7,8,9,14)]
write.csv(data2,file="data2.csv")
data2 <- distinct(data2, gene_name,.keep_all = T)
row.names(data2) <- data2$gene_name
data2 <- data2[,-7]

list <- read.csv("fil1.csv")
list$group <- factor(list$group,levels = c("Control","SeV"))
list$group <- as.factor(list$group)
data1 <- data1[,gpl$sample]

data1$infe2 <- as.numeric(data1$infe2)
str(data1)

dds <- DESeqDataSetFromMatrix(countData = data2,
                              colData = list,
                              design=~group)

dds <- dds[rownames(counts(dds))>1,]
dds <- estimateSizeFactors(dds)
dds
dds <- DESeq(dds)
rld <- rlogTransformation(dds)
expr <- assay(rld)
write.csv(dds,file="data2.csv")

gene1 <- read.csv("fil.csv")
write.csv(gene1,file="gene.csv")

gene1$Gene <- factor(gene1$Gene,
                     c("Isg15"))
p <- ggplot(data=gene1, aes(x=Type,y=Expression)) +
  geom_boxplot(aes(fill=Type),linetype=1) +  
  labs(x = NULL, y="Normalized expression") +  
  scale_fill_manual(values=c("#72b1b9","#f6ae00")) + 
  ylim(7.5,15)+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5), 
        text = element_text(size = 15),  
        axis.text.x = element_text(angle = 45, hjust = 1),  
        legend.position="top", legend.title = element_blank())
p 

ggsave(p,filename="gene2.pdf",width=5,height=6)


##GSE200811
setwd("H:/HBOC/data1/RNA/GSE200811")
data1 <- read.csv("gene_counts_matrix.tsv",header = T,sep = "\t")
write.csv(data1,file="data1.csv")
data2 <- data1[,c(1,2,3,4,19,20,21)]
data2 <- distinct(data2, gene,.keep_all = T)
row.names(data2) <- data2$gene
data2 <- data2[,-1]

list <- read.csv("fil.csv")
list$group <- factor(list$group,levels = c("Control","HCV"))
list$group <- as.factor(list$group)
data1 <- data1[,gpl$sample]

dds <- DESeqDataSetFromMatrix(countData = data2,
                              colData = list,
                              design=~group)

dds <- dds[rownames(counts(dds))>1,]
dds <- estimateSizeFactors(dds)
dds
dds <- DESeq(dds)
rld <- rlogTransformation(dds)
expr <- assay(rld)
write.csv(expr,file="data2.csv")

gene1 <- read.csv("fil1.csv")
write.csv(gene1,file="gene.csv")

p <- ggplot(data=gene1, aes(x=Type,y=Expression)) +
  geom_boxplot(aes(fill=Type),linetype=1) +  
  labs(x = NULL, y="Normalized expression") +  
  scale_fill_manual(values=c("#72b1b9","#f6ae00")) +  
  ylim(11,14)+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5), 
        text = element_text(size = 15),  
        axis.text.x = element_text(angle = 45, hjust = 1),  
        legend.position="top", legend.title = element_blank())
p 

ggsave(p,filename="gene1.pdf",width=5,height=6)


##GSE163285
setwd("H:/HBOC/data1/RNA/GSE163285")
data1 <- read.csv("1.csv",header = T,sep = ",")
write.csv(data1,file="data1.csv")

data1 <- distinct(data1, gene_name,.keep_all = T)
row.names(data1) <- data1$gene_name
data1 <- data1[,-7]

list <- read.csv("fil.csv")
list$group <- factor(list$group,levels = c("Control","HBV"))
list$group <- as.factor(list$group)
data1 <- data1[,gpl$sample]

dds <- DESeqDataSetFromMatrix(countData = data1,
                              colData = list,
                              design=~group)

dds <- dds[rownames(counts(dds))>1,]
dds <- estimateSizeFactors(dds)
dds
dds <- DESeq(dds)
rld <- rlogTransformation(dds)
expr <- assay(rld)
write.csv(expr,file="data1.csv")

gene1 <- read.csv("fil.csv")
write.csv(gene1,file="gene.csv")
gene1$Type <- factor(gene1$Type,levels = c("Control","HBV"))

p <- ggplot(data=gene1, aes(x=Type,y=Expression)) +
  geom_boxplot(aes(fill=Type),linetype=1) +  
  labs(x = NULL, y="Normalized expression") +  
  scale_fill_manual(values=c("#72b1b9","#f6ae00")) + 
  ylim(9.8,11.8)+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5), 
        text = element_text(size = 15),  
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position="top", legend.title = element_blank())
p 

ggsave(p,filename="gene1.pdf",width=5,height=6)


##GSE211979
setwd("H:/HBOC/data1/RNA/GSE211979")
data1 <- read.csv("GSE211979_RNAseq_all_rpkm.txt",header = T,sep = "\t")

gene1 <- read.csv("fil1.csv")
write.csv(data1,file="data1.csv")

p <- ggplot(data=gene1, aes(x=Type,y=Expression)) +
  geom_boxplot(aes(fill=Type)) + 
  labs(x = NULL, y="Normalized expression") +  
  scale_fill_manual(values=c("#72b1b9","#f6ae00")) +  
  ylim(7,16)+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5), 
        text = element_text(size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position="top", legend.title = element_blank())
p 

ggsave(p,filename="gene2.pdf",width=5,height=6)  