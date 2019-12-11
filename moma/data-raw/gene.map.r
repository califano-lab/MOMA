library(dplyr)
library(readr)

raw <- read_tsv("data-raw/hgnc.gene.information.txt", col_types = "cccc")
colnames(raw) <- c("Gene.Symbol", "Cytoband", "Entrez.IDs","Ensembl")

gene.map <- raw %>% select(Gene.Symbol, Entrez.IDs, Cytoband, Ensembl)
gene.map <- as.data.frame(gene.map)

usethis::use_data(gene.map)
