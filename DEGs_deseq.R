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

dds <- DESeqDataSetFromTximport(txi,
                                colData = group.order,
                                design = ~ group)
head(dds@assays@data$counts)
dds@assays@data$avgTxLength
dds
  