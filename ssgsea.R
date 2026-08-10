# tabelas de expressao para geneset obtidas de https://doi.org/10.1016/j.cell.2018.05.060
# (suplementar 2 e 3)

library(tidyverse)
library(openxlsx)

clusters1 <- read.xlsx('mmc2.xlsx', startRow = 3)
clusters2 <- read.xlsx('mmc3.xlsx', sheet = 2)

clusters1 <- clusters1[is.na(clusters1$Removed.from.Analysis),]

cluster <- clusters2 %>% left_join(clusters1[, c(1:3, 5,6)], join_by('cluster' == 'ClusterID'))
cluster <- na.omit(cluster)

library(biomaRt)
# realizar conversao de genes humanos para ortologos em ratos
hs_genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = unique(cluster$gene),
  mart = mart_hs
)

hs_orth <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "rnorvegicus_homolog_associated_gene_name",
    "rnorvegicus_homolog_orthology_type"
  ),
  filters = "ensembl_gene_id",
  values = hs_genes$ensembl_gene_id,
  mart = mart_hs
)

#filtrar df de genes humanos e df de ortologos
conv <- hs_genes %>%
  inner_join(hs_orth, by = "ensembl_gene_id") %>%
  filter(
    rnorvegicus_homolog_associated_gene_name != "",
    rnorvegicus_homolog_orthology_type == "ortholog_one2one"
  ) %>%
  dplyr::select(
    human_gene = hgnc_symbol,
    rat_gene   = rnorvegicus_homolog_associated_gene_name
  ) %>%
  distinct()

conv[!(conv$human_gene == toupper(conv$rat_gene)),]
#==== ssGSEA =====
#agrupar os clusters em um objeto com lista de gene para cada cluster
#esta sendo agrupado pelo nome do cluster, unindo diferentes cluster de uma mesma celula
genesets_human <- cluster %>%
  group_by(Final.Annotation.based.on.bulk.combined.with.differential.expressed.genes) %>%
  summarise(genes = list(unique(gene))) %>%
  deframe()

#gerar o geneset de ratos utilizando o geneset de humanos
#para cada cluster g, pegar o valor da coluna rat_gene corresponende à presença do gene em humanos naquele cluster
genesets_rat <- lapply(genesets_human, function(g) {
  conv$rat_gene[conv$human_gene %in% g]
})

# remover gene sets muito pequenos
genesets_rat <- genesets_rat[sapply(genesets_rat, length) >= 5]

## Para GSVA, se usa expressão normalizada
## nosso rlog(counts) contém o ENSEMBL como linha, mas precisamos do gene symbol
head(rld)
head(annot)

rld.gsva <- rld %>% as.data.frame() %>% rownames_to_column(., 'ENSEMBL') %>% left_join(., dplyr::select(annot, ENSEMBL, SYMBOL)) %>% na.omit
rld.gsva <- rld.gsva[!duplicated(rld.gsva$SYMBOL),]

rld.gsva$ENSEMBL <- NULL
rownames(rld.gsva) <- NULL
rld.gsva <- rld.gsva %>% column_to_rownames(., 'SYMBOL')

library(GSVA)
ssgsea_param <- ssgseaParam(
  exprData   = as.matrix(rld.gsva),
  geneSets   = genesets_rat
)


ssgsea_res <- gsva(ssgsea_param)
sapply(genesets_rat, length)

sapply(genesets_rat, function(g)
  sum(g %in% rownames(rld.gsva))
)


library(pheatmap)

pheatmap(
  ssgsea_res,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete"
)

