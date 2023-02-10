## ---- include = FALSE----------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup, include = FALSE----------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(dplyr)


## ---- eval=FALSE---------------------------------------------------------------------------------------------
## install.packages("vaRHC")


## ---- eval=FALSE---------------------------------------------------------------------------------------------
## if(!require("remotes", quietly = TRUE)) install.packages('remotes') ## Only the first time
## library(remotes)
## devtools::install_github("emunte/vaRHC")


## ------------------------------------------------------------------------------------------------------------
library("vaRHC")


## ---- eval=FALSE, cache=TRUE---------------------------------------------------------------------------------
## eg.gene <- "CDH1"
## eg.variant <- "c.1137+1G>A"
## 
## 
## var.information <- vaR(gene = eg.gene, variant = eg.variant)
## 


## ---- echo=FALSE---------------------------------------------------------------------------------------------
 gene.specific.df <- vaRHC:::connectionDB("SELECT * from  gene_specific;") %>% as.data.frame()
    knitr::kable(gene.specific.df %>% filter(genename %in% c("general", "ATM", "CDH1", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "PTEN", "TP53")) %>% dplyr::select(1:10))
     knitr::kable(gene.specific.df %>% filter(genename %in% c("general", "ATM", "CDH1", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "PTEN", "TP53")) %>% dplyr::select(11:20))
    knitr::kable(gene.specific.df %>% filter(genename %in% c("general", "ATM", "CDH1", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "PTEN", "TP53"))%>% dplyr::select(21:30))
       knitr::kable(gene.specific.df %>% filter(genename %in% c("general", "ATM", "CDH1", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "PTEN", "TP53"))%>% dplyr::select(31:40))
       knitr::kable(gene.specific.df %>% filter(genename %in% c("general", "ATM", "CDH1", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "PTEN", "TP53"))%>% dplyr::select(41:50))
       knitr::kable(gene.specific.df %>% filter(genename %in% c("general", "ATM", "CDH1", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "PTEN", "TP53"))%>% dplyr::select(51:58))
       knitr::kable(gene.specific.df %>% filter(genename %in% c("general", "ATM", "CDH1", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "PTEN", "TP53"))%>% dplyr::select(59:67))


## ---- eval=FALSE---------------------------------------------------------------------------------------------
## eg.gene <- "CDH1"
## eg.variant <- "c.1137+1G>A"
## data("gene_specific")
## eg.gene.specific <- gene_specific
## 
## var.information <- vaR(gene = eg.gene, variant = eg.variant, gene.specific.df = eg.gene.specific)


## ---- include=FALSE, cache=TRUE------------------------------------------------------------------------------
eg.gene <- "CDH1"
eg.variant <- "c.1137+1G>A"

#var.all <- vaR(assembly = eg.assembly, gene = eg.gene, variant = eg.variant, gene.specific.df = NULL)
var.information <- vaR(variant = eg.variant, gene= eg.gene, gene.specific.df = NULL)


## ------------------------------------------------------------------------------------------------------------
names(var.information$vaRinfo)



## ---- echo=F-------------------------------------------------------------------------------------------------
var.info <- var.information$vaRinfo$Variant.Info %>% as.data.frame()
str(var.info)


## ------------------------------------------------------------------------------------------------------------
#nomenclature
var.information$vaRinfo$gnomAD$nomenclature
#coverage
var.information$vaRinfo$gnomAD$coverage
#information from exomes non cancer separated by subpopulations
knitr::kable(var.information$vaRinfo$gnomAD$info$exomes$non.cancer$subpopulations)
#information from exomes + genomes non neuro overall frequency
knitr::kable(var.information$vaRinfo$gnomAD$info$exomes.genomes$non.neuro$overall)


## ---- echo=FALSE---------------------------------------------------------------------------------------------
knitr::kable(var.information$vaRinfo$predictors$predictor.table)


## ---- echo=FALSE---------------------------------------------------------------------------------------------
  knitr::kable(var.information$vaRinfo$codon.stop$variant.exon)


## ---- echo=FALSE---------------------------------------------------------------------------------------------
  knitr::kable(var.information$vaRinfo$codon.stop$canonical.skip.pred)


## ---- echo=FALSE---------------------------------------------------------------------------------------------
  knitr::kable(var.information$vaRinfo$codon.stop$exons)


## ------------------------------------------------------------------------------------------------------------
names(var.information$vaRclass)


## ----echo=FALSE----------------------------------------------------------------------------------------------
knitr::kable(var.information$vaRclass$final.criteria$criteria.res)
var.information$vaRclass$final.criteria[2:length(var.information$vaRclass$final.criteria)]


## ----eval=FALSE----------------------------------------------------------------------------------------------
## #Example 1
## data("ex_vaRbatch")
## ex_vaRbatch[] <- lapply(ex_vaRbatch, as.character) #convert to character
## all <- vaRbatch( all.variants = ex_vaRbatch, spliceai.program = TRUE, splieai.reference= "./hg19.fa", excel.results = TRUE)
## 
## #Example2
## eg.variants <- data.frame(gene=c("ATM", "MSH6", "BRCA1"), variants = c("c.8420A>T", "c.1559G>A", "c.211A>G"))
## batch.results <- vaRbatch( all.variants = eg.variants, spliceai.program =FALSE, print.data.frame = FALSE excel.results = TRUE)


## ----echo=FALSE----------------------------------------------------------------------------------------------
all.transcript <- vaRHC:::connectionDB("SELECT ensembltranscriptID, namegene, NM, NC, CCDS from transcript WHERE maintranscript = 'T';") %>% as.data.frame() %>% dplyr::arrange(namegene)
knitr::kable(all.transcript)


## ------------------------------------------------------------------------------------------------------------
sessionInfo()

