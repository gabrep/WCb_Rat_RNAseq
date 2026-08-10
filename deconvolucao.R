#prep para deconvolution com LM22
library(ggpubr)

tpm_orth <- tpm %>% rownames_to_column('ENSEMBL') %>% 
  left_join(dplyr::select(annot, ENSEMBL, SYMBOL)) %>% 
  left_join(conv, by=c('SYMBOL'='rat_gene'))
conv[4624,]  

tpm_orth <- tpm_orth %>% mutate(GENE = ifelse(is.na(human_gene), toupper(SYMBOL), human_gene))

tpm_orth <- tpm_orth %>% dplyr::select(GENE, W1, W2, WCb1, WCb2, WCb3)
write_tsv(tpm_orth, file = 'Deconvolution/tpm_gene_orth_deconv.txt')

#========================================#
#DECONVOLUÇÃO COM LM22#

deconv <- read.csv('Deconvolution/CIBERSORTx_Job3_Results.csv')
deconv[1:23] %>% writexl::write_xlsx('Deconvolution/CYBERSORTx-RESULTADO.xlsx')

deconv$group <- c('W', 'W', 'WCb75', 'WCb75', 'WCb75')
deconv$group <- factor(deconv$group, levels=c('W', 'WCb75'))

pdf('Deconvolution/NK-activated.pdf', width = 4, height = 4)
deconv %>% 
  ggplot(aes(group, NK.cells.activated, fill=group))+
  stat_summary(fun.data = 'mean_sd', geom='errorbar', width=.2)+
  stat_summary(fun='mean', geom='bar', width=.75)+
  scale_fill_manual(values=cores, guide='none')+
  theme_classic()+
  labs(x='', y='NK Cells Activated (proportion)')+
  coord_cartesian(xlim=c(0.5,2.5), expand = F, ylim=c(0,0.05))+
  geom_bracket(xmin=1, xmax = 2, label='T-test\n p-value = 0.044\n*', inherit.aes = F, y.position = 0.04)+
  theme(axis.text = element_text(color='black'))
dev.off()
t.test(NK.cells.activated ~ group, deconv)

pdf('Deconvolution/NK-RESTING.pdf', width = 4, height = 4)
deconv %>% 
  ggplot(aes(group, NK.cells.resting, fill=group))+
  stat_summary(fun.data = 'mean_sd', geom='errorbar', width=.2)+
  stat_summary(fun='mean', geom='bar', width=.75)+
  scale_fill_manual(values=cores, guide='none')+
  theme_classic()+
  labs(x='', y='NK Cells Resting (proportion)')+
  coord_cartesian(xlim=c(0.5,2.5), expand = F, ylim=c(0,0.5))+
  theme(axis.text = element_text(color='black'))
dev.off()

t.test(NK.cells.resting ~ group, deconv)

deconv %>% 
  ggplot(aes(group, (NK.cells.activated+NK.cells.resting), fill=group))+
  stat_summary(fun.data = 'mean_sd', geom='errorbar', width=.2)+
  stat_summary(fun='mean', geom='bar', width=.75)+
  scale_fill_manual(values=cores)+
  theme_classic()+
  coord_cartesian(xlim=c(0.5,2.5), expand = F, ylim=c(0,0.7))

t.test(NK.cells.activated+NK.cells.resting ~ group, deconv)

pdf('Deconvolution/total_pop.pdf', width = 6, height = 8)
deconv[,c(1:23,27)] %>% 
  reshape2::melt(id.vars=c('Mixture', 'group')) %>% 
  
  ggplot(aes(value, variable, fill=group))+
  stat_summary(fun.data = 'mean_sd', geom='errorbar', width=.9, position=position_dodge2(reverse=T))+
  stat_summary(fun='mean', geom='bar', position=position_dodge2(reverse=T))+
  labs(x='Cell proportion', y='')+
  scale_fill_manual(values=cores, name='Group')+
  theme_classic()+
  coord_cartesian(expand = F, xlim=c(0,.5))+
  theme(axis.text = element_text(color='black'))
dev.off()

calcular_ttest_celulas <- function(df, coluna_grupo = "group", colunas_ignorar = c("Mixture", "P.value", "Correlation", "RMSE")) {
  
  # 1. Isolar apenas as colunas que correspondem às células
  colunas_celulas <- setdiff(names(df), c(coluna_grupo, colunas_ignorar))
  
  # 2. Verificar se há exatamente 2 grupos na coluna de comparação
  grupos <- na.omit(unique(df[[coluna_grupo]]))
  if(length(grupos) != 2) {
    stop(paste("O teste T exige exatamente 2 grupos. Grupos encontrados:", paste(grupos, collapse = ", ")))
  }
  
  # 3. Executar o Teste T em loop para cada tipo celular
  lista_resultados <- lapply(colunas_celulas, function(coluna) {
    
    # tryCatch protege o código de erros (ex: variância zero na coluna)
    tryCatch({
      # Executa o teste
      teste <- t.test(df[[coluna]] ~ df[[coluna_grupo]], data = df)
      
      # Retorna a linha de resultado
      data.frame(
        Tipo_Celular = coluna,
        Estatistica_T = round(unname(teste$statistic), 4),
        P_valor = teste$p.value
      )
    }, error = function(e) {
      # O que fazer se der erro (ex: B.cells.memory todo zerado)
      data.frame(
        Tipo_Celular = coluna,
        Estatistica_T = NA,
        P_valor = NA,
      )
    })
  })
  
  # 4. Juntar a lista em um único data.frame
  df_final <- do.call(rbind, lista_resultados)
  
  return(df_final)
}

# ---------------------------------------------------------
resultados_teste <- calcular_ttest_celulas(deconv)
resultados_teste
writexl::write_xlsx(resultados_teste,'Deconvolution/T_test_Deconv_LM22.xlsx')


#================================================#
#RESULTADOS NSCLS

deconv_nscls <- read.csv('Deconvolution/CIBERSORTx_Job4_Results_nscls.csv')
deconv_nscls$group <- c('W', 'W', 'WCb75', 'WCb75', 'WCb75')
result_nscls <- calcular_ttest_celulas(deconv_nscls)

