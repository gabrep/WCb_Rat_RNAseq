library(GSVA)
library(msigdbr)

#==========GSVA===========#
## Obter gene sets

# Hallmark (processos biológicos amplos)
rat_hallmark <- msigdbr(species = "Rattus norvegicus", category = "H")

# C2 (curated gene sets)
rat_reactome <- msigdbr(species = "Rattus norvegicus", category = "C2", subcategory = 'CP:REACTOME')
rat_bioc <- msigdbr(species = "Rattus norvegicus", category = "C2", subcategory = 'CP:BIOCARTA')
rat_wp <- msigdbr(species = "Rattus norvegicus", category = "C2", subcategory = 'CP:WIKIPATHWAYS')
rat_kegg <- msigdbr(species = "Rattus norvegicus", category = "C2", subcategory = 'CP:KEGG')

rat_c2 <- rbind(rat_reactome, rat_bioc, rat_wp, rat_kegg)

#============== C7 immune==============#
rat_imm<- msigdbr(species = "Rattus norvegicus", category = "C7", subcategory = 'IMMUNESIGDB')

# padrões que queremos manter
termos <- c("MACROPHAGE_vs", "NEUTROPHIL_VS", "NK_VS", "T_CELL_VS", "B_CELL_VS")

innatepat <- c("NKCELL_VS_.+_UP", "CD4_TCELL_VS_.+_UP", "CD8_TCELL_VS_.+_UP",
               "MACROPHAGE_VS_.+_UP", "NEUTROPHIL_VS_.+_UP"
              # "EOSINOPHIL_VS_.+_UP", "BASOPHIL_VS_.+_UP", "MAST_CELL_VS_.+_UP"
              # , "NEUTROPHIL_VS_.+_UP",
              # "CD4_TCELL_VS_.+_UP", "CD8_TCELL_VS_.+_UP", "BCELL_VS_.+_UP")
)

rat_c7_filt <- rat_imm %>%
  filter(grepl(paste(innatepat, collapse = "|"), gs_name, ignore.case = TRUE))

excludepat <- c("NAIVE", "LUPUS", "MYELOID")
rat_c7_filt <- rat_c7_filt %>%
  filter(!grepl(paste(excludepat, collapse = "|"), gs_name, ignore.case = TRUE))

## transformar para lista de genesets com os genes listados
genesets_h <- rat_hallmark %>%
  split(x = .$gene_symbol, f = .$gs_name)

genesets_c2 <- rat_c2 %>%
  split(x = .$gene_symbol, f = .$gs_name)

genesets_imm <- rat_c7_filt %>%
  split(x = .$gene_symbol, f = .$gs_name)

## Para GSVA, se usa expressão normalizada
## nosso rlog(counts) contém o ENSEMBL como linha, mas precisamos do gene symbol
head(rld)
head(annot)

rld.gsva <- rld %>% as.data.frame() %>% rownames_to_column(., 'ENSEMBL') %>% left_join(., dplyr::select(annot, ENSEMBL, SYMBOL)) %>% na.omit
rld.gsva <- rld.gsva[!duplicated(rld.gsva$SYMBOL),]

rld.gsva$ENSEMBL <- NULL
rownames(rld.gsva) <- NULL
rld.gsva <- rld.gsva %>% column_to_rownames(., 'SYMBOL')


## GSVA
gsva_res <- gsvaParam(
  exprData = as.matrix(rld.gsva),
  geneSets = genesets_imm,  minSize = 5 # ou genesets_c2
)
es <- gsva(gsva_res)

group.order[-1,]

model <- model.matrix(~ group, group.order[-1,])
model1 <- model.matrix(~1, group.order[-1,])

library(sva)
library(limma)
sv <- sva(es, model, model1) 

mod <- cbind(model, sv$sv)
fit <- lmFit(es, mod)
fit2 <- eBayes(fit, robust = T)
decideTests(fit2) %>% summary()
topTable(fit2, n=Inf) %>% View
gssizes <- geneSetSizes(es)
fit.eb.trend <- eBayes(fit, robust=TRUE, trend=gssizes)

topTable(fit.eb.trend, n=Inf) %>% View
decideTests(fit.eb.trend) %>% summary()
