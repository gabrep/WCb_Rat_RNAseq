library(pheatmap)
library(factoextra);library(FactoMineR)
library(EnhancedVolcano)
## Functions used to plot dendrogram and to run enrichment analysis can be pulled from https://github.com/gabrep/funcoes_R

cores <- c('#ED5A99', '#20DBA4')
rld <- rlog(counts)
colnames(rld) <- group.order$sample_name[-1] #-1 to remove the excluded sample

##dendrogram 
tpm <- read.csv('tpm.csv')
tpm <- filter(tpm, X %in% rownames(dds))
tpm <- column_to_rownames(tpm,'X')
colnames(tpm) <- group.order$sample_name
source('../../../Doutorado/Bioinformática/funcoes_R/color_dend.R')

color_dend(tpm[,-1], group.order$group[-1], col_branches = F, cores_dend = cores)

##PCA
pca <- PCA(t(rld), graph=F)

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
colnames(ann) <- 'Grupo'
pdf('figures/heatmap.pdf', width = 4, height = 5)
pheatmap(na.omit(rld),
         scale = 'row',
         show_rownames = F,
         treeheight_row = 0,
         cluster_cols = T, clustering_distance_cols = 'euclidean',
         annotation_col = ann,
         annotation_colors = list(Grupo=c('W' = cores[1], 'WCb75' = cores[2])),
         color=colorRampPalette(c('deepskyblue3', 'lightblue', 'white','pink', 'violetred'))(100))
dev.off()
#Degs counts
shr.res %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>% 
  mutate(reg = ifelse(log2FoldChange > 0, 'Up', 'Down')) %>% 
  ggplot(aes(reg, fill=reg))+
  geom_bar(color='black', width = .75)+
  labs(y=NULL, x=NULL)+
  stat_count(geom = "text", colour = c("black", 'white'), size = 4, aes(label = ..count..), vjust=1.5)+
  scale_fill_manual(values=c('white', 'seagreen'), guide='none')+
  theme_classic()+
  coord_cartesian(xlim=c(.5,2.5), expand = F, ylim = c(0,35))
ggsave('figures/degs.png', width = 2.5, height = 3.5, dpi = 300)

#volcano
keyvals <- ifelse(shr.res$padj >= 0.05, "gray",
                  ifelse(shr.res$log2FoldChange <= -1, cores[1],
                         ifelse(shr.res$log2FoldChange >= 1, cores[2], "gray")))
keyvals[is.na(keyvals)] <- "gray"
names(keyvals)[keyvals == cores[2]] <- "Up-regulated"
names(keyvals)[keyvals == "gray"] <- "Not significant"
names(keyvals)[keyvals == cores[1]] <- "Down-regulated"

EnhancedVolcano(shr.res,
                x= 'log2FoldChange',
                y='pvalue',
                max.overlaps = Inf,
                lab = shr.res$Gene, 
                labSize = 3, pointSize = 4,
                pCutoffCol = 'padj',
                colCustom = keyvals,
                subtitle = NULL, caption = NULL, title=NULL,
                FCcutoff = 1, pCutoff = 0.05)+
  theme_classic()+
  theme(legend.position = 'top',
        legend.title = element_blank())+
  coord_cartesian(xlim=c(-4,10), ylim=c(0,12.5))
ggsave('figures/volcano.png', width = 6, height = 4)
