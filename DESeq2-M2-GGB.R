library('DESeq2')
library('ggplot2')

########### mRNA
#data
mRNA=read.table('/projet/umr7139/ga/acormier/finalresult/lncRNA/HTSeq/mRNA.count', header=F, row.names=1)

#subclass
mRNA_GA_imm_f_vs_m <- mRNA[,c("V2","V3","V4","V5")]
mRNA_GA_mat_f_vs_m <- mRNA[,c("V6","V7","V8","V9")]

#design
mRNA_GA_imm_f_vs_m_design=data.frame(
  row.names = colnames(mRNA_GA_imm_f_vs_m),
  condition = c("GA_imm_f","GA_imm_f","GA_imm_m","GA_imm_m"),
  type = c(rep("single-end", 4)))

mRNA_GA_mat_f_vs_m_design=data.frame(
  row.names = colnames(mRNA_GA_mat_f_vs_m),
  condition = c("GA_mat_f","GA_mat_f","GA_mat_m","GA_mat_m"),
  type = c(rep("single-end", 4)))

#create matrix for each DE analysis
mRNA_GA_imm_f_vs_m_matrix <- DESeqDataSetFromMatrix(countData=mRNA_GA_imm_f_vs_m, colData=mRNA_GA_imm_f_vs_m_design, design = ~condition)
mRNA_GA_mat_f_vs_m_matrix <- DESeqDataSetFromMatrix(countData=mRNA_GA_mat_f_vs_m, colData=mRNA_GA_mat_f_vs_m_design, design = ~condition)

# run deseq2
mRNA_GA_imm_f_vs_m_deseq2 <- DESeq(mRNA_GA_imm_f_vs_m_matrix)
mRNA_GA_mat_f_vs_m_deseq2 <- DESeq(mRNA_GA_mat_f_vs_m_matrix)

#get the results
mRNA_GA_imm_f_vs_m_res <- results(mRNA_GA_imm_f_vs_m_deseq2)
mRNA_GA_mat_f_vs_m_res <- results(mRNA_GA_mat_f_vs_m_deseq2)

#filter
mRNA_GA_imm_f_vs_m_sig <- subset(mRNA_GA_imm_f_vs_m_res, padj < 0.1 & abs(log2FoldChange) >=1)
mRNA_GA_mat_f_vs_m_sig <- subset(mRNA_GA_mat_f_vs_m_res, padj < 0.1 & abs(log2FoldChange) >=1)

#reorder
mRNA_GA_imm_f_vs_m_resOrdered <- mRNA_GA_imm_f_vs_m_res[order(mRNA_GA_imm_f_vs_m_res$padj),]
mRNA_GA_mat_f_vs_m_resOrdered <- mRNA_GA_mat_f_vs_m_res[order(mRNA_GA_mat_f_vs_m_res$padj),]

#summary of the results
summary(mRNA_GA_imm_f_vs_m_res)
summary(mRNA_GA_mat_f_vs_m_res)

#number of DE using p-value
sum(mRNA_GA_imm_f_vs_m_res$padj < 0.1, na.rm=TRUE)
sum(mRNA_GA_mat_f_vs_m_res$padj < 0.1, na.rm=TRUE)

#make MA plot
plotMA(mRNA_GA_imm_f_vs_m_res, main="DESeq2", ylim=c(-2,2))
plotMA(mRNA_GA_mat_f_vs_m_res, main="DESeq2", ylim=c(-2,2))

#prepare for plotting
mRNA_GA_imm_f_vs_m_rld <- rlog(mRNA_GA_imm_f_vs_m_deseq2, blind=FALSE)
mRNA_GA_mat_f_vs_m_rld <- rlog(mRNA_GA_mat_f_vs_m_deseq2, blind=FALSE)

mRNA_GA_imm_f_vs_m_sampleDists <- dist(t(assay(mRNA_GA_imm_f_vs_m_rld)))
mRNA_GA_mat_f_vs_m_sampleDists <- dist(t(assay(mRNA_GA_mat_f_vs_m_rld)))


library("RColorBrewer")
mRNA_GA_imm_f_vs_m_sampleDistMatrix <- as.matrix(mRNA_GA_imm_f_vs_m_sampleDists)
mRNA_GA_mat_f_vs_m_sampleDistMatrix <- as.matrix(mRNA_GA_mat_f_vs_m_sampleDists)


rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# plot PCA samples
plotPCA(mRNA_GA_imm_f_vs_m_rld, intgroup=c("condition", "type"))
plotPCA(mRNA_GA_mat_f_vs_m_rld, intgroup=c("condition", "type"))

nrow(mRNA_GA_imm_f_vs_m_sig)
nrow(mRNA_GA_mat_f_vs_m_sig)

write.table(mRNA_GA_imm_f_vs_m_sig, sep="\t", quote=F, row.names = T, file="~/Bureau/mRNA_GA_imm_m-vs-f_sig.deseq2")
write.table(mRNA_GA_mat_f_vs_m_sig, sep="\t", quote=F, row.names = T, file="~/Bureau/mRNA_GA_mat_m-vs-f_sig.deseq2")

write.table(mRNA_GA_imm_f_vs_m_res, sep="\t", quote=F, row.names = T, file="~/Dropbox/Alex_thesis/Projects/lncRNA/DESeq2/mRNA/mRNA_GA_imm_m-vs-f_res.deseq2")
write.table(mRNA_GA_mat_f_vs_m_res, sep="\t", quote=F, row.names = T, file="~/Dropbox/Alex_thesis/Projects/lncRNA/DESeq2/mRNA/mRNA_GA_mat_m-vs-f_res.deseq2")










