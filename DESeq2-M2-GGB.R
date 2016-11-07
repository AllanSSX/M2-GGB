### chargement du package DESeq2
library('DESeq2')

# chargement des données
mRNA=read.table('C:/Users/Master/Desktop/mRNA-count_GA_imm-mat_m-f.count', header=F, row.names=1)

# sélection des colonnes pour les gamétophytes immatures
mRNA_GA_imm_f_vs_m <- mRNA[,c("V2","V3","V4","V5")]

# design expérimental du fichier chargé
mRNA_GA_imm_f_vs_m_design=data.frame(
  row.names = colnames(mRNA_GA_imm_f_vs_m),
  condition = c("GA_imm_f","GA_imm_f","GA_imm_m","GA_imm_m"),
  type = c(rep("single-end", 4)))

# création de la matrice pour DESeq2
mRNA_GA_imm_f_vs_m_matrix <- DESeqDataSetFromMatrix(countData=mRNA_GA_imm_f_vs_m, colData=mRNA_GA_imm_f_vs_m_design, design = ~condition)

# analyse DESeq2
mRNA_GA_imm_f_vs_m_deseq2 <- DESeq(mRNA_GA_imm_f_vs_m_matrix)

# obtention des résultats
mRNA_GA_imm_f_vs_m_res <- results(mRNA_GA_imm_f_vs_m_deseq2)

# filtre pour conserver les gènes différentiellement exprimés
mRNA_GA_imm_f_vs_m_sig <- subset(mRNA_GA_imm_f_vs_m_res, padj < 0.1 & abs(log2FoldChange) >=1)

# création du MA plot
plotMA(mRNA_GA_imm_f_vs_m_res, main="DESeq2", ylim=c(-2,2))

# enregistrement des résultats
write.table(mRNA_GA_imm_f_vs_m_sig, sep="\t", quote=F, row.names = T, file="C:/Users/Master/Desktop/mRNA_GA_imm_m-vs-f_sig.deseq2")
write.table(mRNA_GA_imm_f_vs_m_res, sep="\t", quote=F, row.names = T, file="C:/Users/Master/Desktop/mRNA_GA_imm_m-vs-f_res.deseq2")

