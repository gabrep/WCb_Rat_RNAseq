library(pheatmap)
library(factoextra);library(FactoMineR)

## Functions used to plot dendrogram and to run enrichment analysis can be pulled from https://github.com/gabrep/funcoes_R

##dendrogram 
tpm <- read.csv('tpm.csv')
tpm <- filter(tpm, X %in% rownames(dds))
tpm <- column_to_rownames(tpm,'X')
colnames(tpm) <- group.order$sample_name
source('../../../Doutorado/Bioinformática/funcoes_R/color_dend.R')
color_dend(tpm[,-1], group.order$group[-1])
color_dend(counts[-1], group.order$group[-1])

##PCA
counts(rld)
pca <- PCA(t(tpm[,-1]), graph=F)

fviz_pca_ind(pca, mean.point=F,
             fill.ind = group.order$group[-1],
             addEllipses = T, ellipse.type='confidence')+
  theme_classic()+
  scale_fill_manual(values=cores, name='Grupos')+
  scale_color_manual(values=cores, name='Grupos')

##Heatmap

pheatmap(tpm[,-1],
         scale = 'row',
         show_rownames = F,
         treeheight_row = 0,
         cluster_cols = F)

## Differences between log2FC from res and shrinked log2FC
res %>% filter(abs(log2FoldChange) > 1, padj <0.05) %>% dim()
shr.res %>% filter(abs(log2FoldChange) > 1, padj <0.05) %>% dim()


BioVenn::draw.venn(list_x = filter(res, abs(log2FoldChange) > 1, padj < 0.05)$SYMBOL,
                   list_y = filter(shr.res, abs(log2FoldChange) > 1, padj < 0.05)$SYMBOL,list_z = NULL)

#Enrichment analysis
source('../../../Doutorado/Bioinformática/funcoes_R/run_GSEA.R')
source('../../../Doutorado/Bioinformática/funcoes_R/run_ORA.R')
library(org.Rn.eg.db)
ora <- run_ORA(res, lfc_threshold = 1, 
        OrgDb = org.Rn.eg.db, organismKEGG = 'rno', organismWP = 'Rattus norvegicus')

ora$go_up@result %>% View
ora$go_down@result %>% View

ora$kegg_up@result %>% View()
ora$wp_down@result %>% View

ora.shr <- run_ORA(shr.res, lfc_threshold = 1, 
                   OrgDb = org.Rn.eg.db, organismKEGG = 'rno', organismWP = 'Rattus norvegicus')

ora.shr$go_up@result %>% View
ora.shr$go_down@result %>% View

ora.shr$kegg_up@result %>% View
ora.shr$kegg_down@result %>% View
ora.shr$wp_up@result %>% View
ora.shr$wp_down@result %>% View

test <- enrichplot::pairwise_termsim(ora.shr$go_up)
test <- enrichplot::treeplot(test)

enrichplot::pairwise_termsim(go.mf$go_up) %>% enrichplot::treeplot(.)
plot(test)


gsea <- run_GSEA(res, org = 'Rattus norvegicus')



##color testing
'darkmagenta' 'aquamarine' 'springgreen'

x <- data.frame(valor=c(8, 12, 6, 5, 7, 8), grupo = c('a', 'b','c','d','e','f'))

x %>% ggplot(aes(grupo, valor, fill=grupo))+
  geom_col(color='black')+
  scale_fill_manual(values=viridis::mako(10, begin = .2))

"#6FD4ADFF" '#D46F96'


x %>% filter(grupo %in% c('a','b')) %>%  ggplot(aes(grupo, valor, fill=grupo))+
  geom_col(color='black')+
  scale_fill_manual(values=c("#6FD4ADFF",'#D46F96'))
