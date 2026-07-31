#oncoplot
library(maftools)

pdf("Oncoplot_Top20.pdf", width = 12, height = 8)

oncoplot(
  maf = maf,
  top = 20,
  removeNonMutated = TRUE,
  drawColBar = TRUE,
  drawRowBar = TRUE,
  sortByMutation = TRUE,
  showTumorSampleBarcodes = FALSE,
  fontSize = 0.8
)

#TOP 20 mutation
library(maftools)
library(ggplot2)
library(dplyr)

maf <- read.maf("your_file.maf")

top20 <- getGeneSummary(maf)$Hugo_Symbol[1:20]

res <- clinicalEnrichment(
  maf = maf,
  clinicalFeature = "Group"
)

df <- res$groupwise_comparision %>%
  filter(Hugo_Symbol %in% top20)

ggplot(df, aes(x = Hugo_Symbol, y = MutatedSamples, fill = Clinical_Group)) +
  geom_col(position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Mutated Samples", fill = "Group")

#Enrichment analysis
library(limma)
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggrepel)
library(dplyr)


# Differential expression analysis
expr <- log2(expr + 1)

design <- model.matrix(~0 + group$Group)

colnames(design) <- levels(group$Group)

fit <- lmFit(expr, design)

contrast.matrix <- makeContrasts(
  Treatment-Control,
  levels=design
)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

deg <- topTable(fit2,
                number = Inf,
                adjust.method = "BH")

deg$Gene <- rownames(deg)

deg$Regulation <- "NotSig"

deg$Regulation[
  deg$adj.P.Val < 0.05 &
    deg$logFC > 1
] <- "Up"

deg$Regulation[
  deg$adj.P.Val < 0.05 &
    deg$logFC < -1
] <- "Down"

write.csv(deg,
          "DEG_results.csv",
          row.names = FALSE)


# Volcano plot
p <- ggplot(deg,
            aes(logFC,
                -log10(adj.P.Val),
                color=Regulation))+
  
  geom_point(size=2)+
  
  scale_color_manual(values=c(
    Up="#E64B35",
    Down="#4DBBD5",
    NotSig="grey70"
  ))+
  
  geom_vline(
    xintercept=c(-1,1),
    linetype=2
  )+
  
  geom_hline(
    yintercept=-log10(0.05),
    linetype=2
  )+
  
  theme_bw()+
  
  labs(
    x="log2 Fold Change",
    y="-log10(FDR)"
  )

ggsave("Volcano_plot.pdf",
       p,
       width=7,
       height=6)


# Gene ID conversion

gene.df <- bitr(
  deg$Gene,
  fromType="SYMBOL",
  toType="ENTREZID",
  OrgDb=org.Hs.eg.db
)

deg2 <- merge(
  deg,
  gene.df,
  by.x="Gene",
  by.y="SYMBOL"
)

sig_gene <- deg2 %>%
  filter(adj.P.Val<0.05 & abs(logFC)>1)

geneList <- sig_gene$ENTREZID

# Prepare ranked gene list

gene_rank <- deg$logFC

names(gene_rank) <- deg$Gene

gene_rank <- sort(
  gene_rank,
  decreasing = TRUE
)

# Run GSEA
gsea_result <- GSEA(
  geneList = gene_rank,
  TERM2GENE = gmt,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# GSEA enrichment map

gsea_result2 <- pairwise_termsim(
  gsea_result
)


pdf(
  "GSEA_emapplot.pdf",
  width=10,
  height=8
)

emapplot(
  gsea_result2,
  showCategory=30
)

#ssGSEA score
library(GSVA)
library(GSEABase)
library(pheatmap)
library(ggplot2)
library(dplyr)

expr <- as.matrix(expr)

expr <- log2(
  expr + 1
)

# ssGSEA analysis

ssgsea_score <- gsva(
  expr,
  gene_sets,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = TRUE
)

