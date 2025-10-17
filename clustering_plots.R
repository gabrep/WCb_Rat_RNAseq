library(pheatmap)
library(factoextra);library(FactoMineR)

## Functions used to plot dendrogram and to run enrichment analysis can be pulled from https://github.com/gabrep/funcoes_R

cores <- c('#ED5A99', '#20DBA4')

##dendrogram 
tpm <- read.csv('tpm.csv')
tpm <- filter(tpm, X %in% rownames(dds))
tpm <- column_to_rownames(tpm,'X')
colnames(tpm) <- group.order$sample_name
source('../../../Doutorado/Bioinformática/funcoes_R/color_dend.R')

color_dend(tpm[,-1], group.order$group[-1], col_branches = F, cores_dend = cores)

##PCA
pca <- PCA(t(tpm[,-1]), graph=F)

fviz_pca_ind(pca, title='',
             geom = 'point', pointshape=21, pointsize=4,
             mean.point=F,
             fill.ind = group.order$group[-1],
             addEllipses = T, ellipse.type='confidence')+
  theme_classic()+
  scale_fill_manual(values=cores, name='Grupos')+
  scale_color_manual(values=cores, name='Grupos')
ggsave('figures/pca.png', width = 5, height = 4, dpi = 300)

##Heatmap
ann <- group.order %>% filter(!sample_name == 'W3') %>% dplyr::select(group, sample_name) %>% column_to_rownames('sample_name')
ann$group <- factor(ann$group)
levels(ann$group)
pheatmap(na.omit(tpm[1:500,-1]),
         scale = 'row',
         show_rownames = F,
         treeheight_row = 0,
         cluster_cols = F,
         annotation_col = ann,
         annotation_colors = list(group=c('W' = cores[1], 'WCb75' = cores[2])),

         color=colorRampPalette(c('deepskyblue3', 'lightblue', 'white','pink', 'firebrick1'))(100))

#Degs counts
shr.res %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% 
  mutate(reg = ifelse(log2FoldChange > 0, 'up', 'down')) %>% 
  ggplot(aes(reg, fill=reg))+
  geom_bar(color='black', width = .75)+
  stat_count(geom = "text", colour = c("black", 'white'), size = 4, aes(label = ..count..), vjust=1.5)+
  scale_fill_manual(values=c('white', 'navy'), name='Regulação')+
  theme_classic()+
  coord_cartesian(xlim=c(.5,2.5), expand = F, ylim = c(0,35))
ggsave('figures/degs.png', width = 3.5, height = 3.5, dpi = 300)
