## Use tximport to import transcript abundances estimated by Salmon through
## nextflow nf-core/rnaseq pipeline ##

## Run DESeq2

## Rattus norvegicus genome, transcriptome and gene annotation ENSEMBL 
## version mRatBN7.2 ##

library(tidyverse)
library(tximport)
library(DESeq2)

#read tx2gene
#tx2gene.tsv created from Rattus_norvegicus.mRatBN7.2.113.gtf.gz
tx2gene <- read_tsv('../tx2gene.tsv', col_names = F)

#retrieve paths to quants.sf
dir <- '../Results/salmon/'
samples <- list.dirs(dir, recursive = F, full.names = T)

files <- file.path(samples, 'quant.sf')
names(files) <- basename(samples)
files

#Change files order to match the samples group
group.order <- data.frame(samples = c('rna_1', 'rna_11', 'rna_3', 'rna_8', 'rna_9', 'rna_10'),
                          group = c(rep('W', 3), rep('WCb75',3)))
files <- files[match(group.order$samples, names(files))]

txi <- tximport(files,
                type = 'salmon',
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

head(txi$counts)
head(txi$length)
txi$abundance
write.csv(txi$abundance, 'tpm.txt')
write.csv(group.order, 'groups.txt')



dds <- DESeqDataSetFromTximport(txi,
                                colData = group.order,
                                design = ~ group)

## sample filtering
## dendrogram clusterization reveals that sample rna_1 shows unusual pattern compared to its groups
dds <- dds[,-1]
counts <- counts(dds)
counts <- counts %>% as.data.frame %>% rownames_to_column('tx')
keep <- rowSums(counts(dds) >= 10) >= 2
sum(keep)
dds <- dds[keep, ]

dds <- DESeq(dds)
res <- results(dds)
res %>% as.data.frame %>%  View

resultsNames(dds)
shr.res <- lfcShrink(dds = dds, coef = "group_WCb75_vs_W")
shr.res <- as.data.frame(shr.res)
hist(res$pvalue)

table(shr.res$padj < 0.05)

##Gene annotation
library(AnnotationDbi)
library(org.Rn.eg.db)
keytypes(org.Rn.eg.db)
org.Rn.eg.db

annot <- AnnotationDbi::select(org.Rn.eg.db, 
                               keys = rownames(res),
                               keytype = 'ENSEMBL',
                               columns = c('GENENAME','SYMBOL','ENTREZID'))

res <- as.data.frame(res) %>% rownames_to_column('ENSEMBL') %>% 
  left_join(., annot)

res %>% filter(abs(log2FoldChange) > 1, padj <0.05) %>% View
