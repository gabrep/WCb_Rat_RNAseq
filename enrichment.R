## Functions used to plot dendrogram and to run enrichment analysis can be pulled from https://github.com/gabrep/funcoes_R

#Enrichment analysis
source('../../../Doutorado/Bioinformática/funcoes_R/run_GSEA.R')
source('../../../Doutorado/Bioinformática/funcoes_R/run_ORA.R')
library(org.Rn.eg.db)
ora.shr <- run_ORA(shr.res, lfc_threshold = 1, 
                   OrgDb = org.Rn.eg.db, organismKEGG = 'rno', organismWP = 'Rattus norvegicus')

ora.shr$go_up@result %>% View
ora.shr$go_down@result %>% View


go.up <- ora.shr$go_up@result %>% mutate(reg='up')
go.down <- ora.shr$go_down@result %>% mutate(reg='down')
#order GO results by fold enrichment score
go.up <- go.up[order(go.up$FoldEnrichment, decreasing = T),]
go.down <- go.down[order(go.down$FoldEnrichment, decreasing = T),]

kegg.up <- ora.shr$kegg_up@result
wp.up <- ora.shr$wp_up@result

## GO Plots ##
rbind(go.up %>% slice_max(n = 20, order_by = FoldEnrichment),
      go.down %>% slice_max(n = 20, order_by = FoldEnrichment)) %>%
  filter(ONTOLOGY =='BP') %>% 
  ggplot(aes(reg, reorder(Description, -FoldEnrichment), fill=FoldEnrichment))+
  geom_tile(color='white')+
  geom_text(aes(label=round(FoldEnrichment,2)), color='white')

go.up %>% filter(ONTOLOGY == 'BP') %>% slice_max(n = 10, order_by = FoldEnrichment) %>% 
  ggplot(aes(log10(FoldEnrichment),reorder(Description, FoldEnrichment), fill=FoldEnrichment))+
  geom_col()+
  geom_text(aes(label=paste('p-adj =',round(p.adjust, 3))), color='white', hjust=1, nudge_x = -.1, size=3)+
  labs(y=NULL, title='GO:BP')+
  scale_fill_gradient(low = '#104020', high = 'springgreen3')+
  theme_classic()
ggsave('figures/GOup_bp.png', width = 8, height = 4)

go.up %>% filter(ONTOLOGY == 'MF') %>% slice_max(n = 10, order_by = FoldEnrichment) %>% 
  ggplot(aes(log10(FoldEnrichment),reorder(Description, FoldEnrichment), fill=FoldEnrichment))+
  geom_col()+
  geom_text(aes(label=paste('p-adj =',round(p.adjust, 3))), color='white', hjust=1, nudge_x = -.1, size=3)+
  labs(y=NULL, title='GO:MF')+
  scale_fill_gradient(low = '#104020', high = 'springgreen3')+
  theme_classic()
ggsave('figures/GOup_mf.png', width = 6, height = 4)

go.up %>% filter(ONTOLOGY == 'CC') %>% slice_max(n = 10, order_by = FoldEnrichment) %>% 
  ggplot(aes(log10(FoldEnrichment),reorder(Description, FoldEnrichment), fill=FoldEnrichment))+
  geom_col()+
  geom_text(aes(label=paste('p-adj =',round(p.adjust, 3))), color='white', hjust=1, nudge_x = -.1, size=3)+
  labs(y=NULL, title='GO:CC')+
  scale_fill_gradient(low = '#104020', high = 'springgreen3')+
  theme_classic()
ggsave('figures/GOup_cc.png', width = 6, height = 4)


enrichplot::pairwise_termsim(ora.shr$go_up) %>% enrichplot::treeplot()
ggsave('figures/treeplot go up.png', width = 18, height = 8)


#GSEA
gsea.shr <- run_GSEA(shr.res, org = 'Rattus norvegicus')
