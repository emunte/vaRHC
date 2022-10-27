## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ---- eval=FALSE--------------------------------------------------------------
#  install.packages('remotes') ## solo primera vez
#  library(remotes)
#  devtools::install_github("emunte/vaRHC")

## -----------------------------------------------------------------------------
library("vaRHC")

## ---- eval=FALSE, cache=TRUE--------------------------------------------------
#  eg.gene <- "CDH1"
#  eg.variant <- "c.1137+1G>A"
#  eg.assembly <- "hg19"
#  
#  var.information <- vaR(assembly = eg.assembly, gene = eg.gene, variant = eg.variant)
#  

## ----echo=FALSE---------------------------------------------------------------
all.transcript <- vaRHC:::connectionDB("SELECT ensembltranscriptID, namegene, NM, NC, CCDS from transcript WHERE maintranscript = 'T';") %>% as.data.frame() %>% dplyr::arrange(namegene)
DT::datatable(all.transcript)

## ---- echo=FALSE--------------------------------------------------------------
 gene.specific.df <- vaRHC:::connectionDB("SELECT * from  gene_specific;") %>% as.data.frame()
    knitr::kable(gene.specific.df %>% filter(genename %in% c("general", "ATM", "CDH1", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "PTEN", "TP53")))

## ---- eval=FALSE--------------------------------------------------------------
#  eg.gene <- "CDH1"
#  eg.variant <- "c.1137+1G>A"
#  eg.assembly <- "hg19"
#  eg.gene.specific <- read.csv("./gene.specific.df")
#  
#  var.information <- vaR(assembly = eg.assembly, gene = eg.gene, variant = eg.variant, gene.specific.df = eg.gene.specific)

