################################################################################
## predictors
################################################################################

## General
predicInfo<- function(object, gene.specific, bbdd, gnomad, spliceai.program=FALSE, spliceai.reference=NULL, spliceai.annotation = system.file("docs", "gencode_spliceai_hg19.txt", package="vaRHC"), spliceai.distance=1000, spliceai.masked=1, provean.program=FALSE){
  httr::set_config(httr::config(ssl_verifypeer = 0L))
  ensembl.id <- ensemblTranscript(object$NM, object$gene) 
  #scores
  phylop <- ucscRPredictor(track = "phyloP100wayAll", object)
  phastcons <- ucscRPredictor(track = "phastCons100way",object)
  gerp <- ucscRPredictor(track = "allHg19RS_BW", object)
  dbnsfp <- bbdd$dbnsfp %>% 
                        dplyr::select(REVEL_score, VEST4_score, PROVEAN_score, BayesDel_noAF_score)
  if(nrow(dbnsfp)==0){
    dbnsfp <- dbnsfp %>% 
                     dplyr::add_row (REVEL_score=NA, VEST4_score=NA, PROVEAN_score=NA, BayesDel_noAF_score=NA)
  }
  
  if (object$most.severe.consequence %in% c("inframe_deletion", "inframe_insertion") && provean.program == TRUE){
      dbnsfp$provean <- proveanR(object)
  }
  align.gvgd <- alignGvgd(object, bbdd)
  spliceai.score <- spliceaiR(object = object, ext.spliceai = gnomad, bbdd = bbdd, spliceai.program = spliceai.program, reference.splice = spliceai.reference, annotation.splice = spliceai.annotation, distance = spliceai.distance, mask = spliceai.masked)
  # while(is.na(spliceai.score$Acceptor_Gain)[1]){
  #   spliceai.score <-  spliceaiR(object, ensembl.id)
  # }
  prior.utah <- bbdd$prior %>% 
                           dplyr::select(Polyphen, MAPP, prior)
  prior.utah.splice <- bbdd$prior %>% 
                                  dplyr::select(refsplice_prior, splice_severity, de_novo_prior, dn_severity, protein_prior, applicable_prior) %>%
                                  dplyr::mutate(splice_severity=ifelse(splice_severity=="", 
                                                                       "not applicable", 
                                                                       splice_severity), 
                                                dn_severity= ifelse(dn_severity=="", 
                                                                    "not applicable", 
                                                                    dn_severity))
  if(nrow(prior.utah)==0){
    prior.utah <-prior.utah %>% 
      dplyr::add_row (Polyphen=NA, MAPP=NA, prior=NA)
  }
  if(nrow(prior.utah.splice)==0){
    prior.utah.splice <- prior.utah.splice %>%
                                           dplyr::add_row(refsplice_prior = NA, splice_severity = NA, de_novo_prior = NA, dn_severity = NA, protein_prior = NA , applicable_prior = NA)
  }
  if(nrow(prior.utah.splice)>0){
    prior.utah.splice.denovo <- c(prior.utah.splice$dn_severity[1], prior.utah.splice$de_novo_prior[1])
    prior.utah.splice.reference <- c(prior.utah.splice$splice_severity[1], prior.utah.splice$refsplice_prior[1])
  } else{
    prior.utah.splice.denovo <- NA
    prior.utah.splice.reference <- NA
  }
  
  trap.score <- traPscore(gnomADnomen(object = object))
  trap.score <- ifelse(is.null(trap.score), NA, trap.score)
  
  #gene specific
  gene.sp <- gene.specific %>% row.names()
  
  #cut offs
  op.ben <- c(gene.specific %>% select(op_phylop_ben, op_phastcons_ben, op_gerp_ben, 
                                       op_revel_ben, op_VEST4_ben, op_provean_ben, op_bayesDel_noAF_ben,op_agvgd_ben, op_polyphen_ben, op_MAPP_ben, op_prior_utah_prot_ben) 
              %>% as.matrix() %>% c(), NA, NA, rep(gene.specific %>% select(op_spliceai_ben) %>% as.matrix() %>% c(), 4),  gene.specific %>% select(op_trap_ben) %>% as.matrix() %>% c())
  
  op.pat <- c(gene.specific %>% select(op_phylop_pat, op_phastcons_pat, op_gerp_pat, 
                                       op_revel_pat, op_VEST4_pat, op_provean_pat, op_bayesDel_noAF_pat,op_agvgd_pat, op_polyphen_pat, op_MAPP_pat, op_prior_utah_prot_pat_sup) 
              %>% as.matrix() %>%c(), NA, NA, rep(gene.specific %>% select(op_spliceai_pat) %>% as.matrix() %>% c(), 4),  gene.specific %>% select(op_trap_pat) %>% as.matrix() %>% c())
  
  cut.ben <- c(gene.specific %>% select(phylop_ben, phastcons_ben, gerp_ben, 
                                        revel_ben, VEST4_ben, provean_ben, bayesDel_noAF_ben,agvgd_ben, polyphen_ben, MAPP_ben, prior_utah_prot_ben) 
               %>% as.matrix() %>%c(), NA, NA, rep(gene.specific %>% select(spliceai_ben) %>% as.matrix() %>% c(), 4),  gene.specific %>% select(trap_ben) %>% as.matrix() %>% c())
  cut.pat <- c(gene.specific %>% select(phylop_pat, phastcons_pat, gerp_pat, 
                                        revel_pat, VEST4_pat, provean_pat, bayesDel_noAF_pat,agvgd_pat, polyphen_pat, MAPP_pat, prior_utah_prot_pat_sup) 
               %>% as.matrix() %>%c(), NA, NA, rep(gene.specific %>% select(spliceai_pat) %>% as.matrix() %>% c(), 4),   gene.specific %>% select(trap_pat) %>% as.matrix() %>% c())
  
  #table
  predictors.table <- data.frame(
    type = c(rep("Nucleotide conservation", 3), 
             rep("Protein effect", 8), 
             rep("Splicing Predictor",7)
    ),
    predictor = c("Phylop", "Phastcons", "Gerp", "Revel", 
                  "VEST4", "Provean", "BayesDel_noAF", "aGVGD_Zebrafish","PolyPhen", "MAPP", 
                  "Prior_utah(MAPP/PP2)", "Prior_utah_splicing_reference","Prior_utah_splicing_de_novo","SpliceAI-AcceptorGain", "SpliceAI-AcceptorLoss", 
                  "SpliceAI-DonorGain", "SpliceAI-DonorLoss", "TraP"
    ),
    classification= NA,
    use = NA,
    values = as.numeric(c(unlist(phylop), unlist(phastcons), 
                          unlist(gerp), dbnsfp, align.gvgd,
                          unlist(prior.utah), prior.utah.splice.reference[2], prior.utah.splice.denovo[2], unlist(purrr::map(spliceai.score,1)), trap.score
    )),
    position = c(rep(NA,13), unlist(purrr::map(spliceai.score,2)), NA
    ),
    operator.benign= op.ben,
    BE.cut.off = cut.ben,
    operator.pathogenic= op.pat,
    TD.cut.off = cut.pat, stringsAsFactors = FALSE
  )
  row.names(predictors.table) <- c("Phylop", "Phastcons", "Gerp", "Revel", 
                                   "VEST4", "Provean", "BayesDel_noAF","aGVGD_zebrafish", "PolyPhen", "MAPP", 
                                   "Prior_utah(MAPP/PP2)", "Prior_utah_splicing_reference","Prior_utah_splicing_de_novo", "SpliceAI-AcceptorGain", "SpliceAI-AcceptorLoss", 
                                   "SpliceAI-DonorGain", "SpliceAI-DonorLoss", "TraP")
  predictors.table2 <- predictorResult(table.predictor = predictors.table) %>% predictorsUse( object = object)
  predictors.table2["Prior_utah_splicing_de_novo", "classification"] <- prior.utah.splice.denovo[1]
  predictors.table2["Prior_utah_splicing_reference", "classification"] <- prior.utah.splice.reference[1]
  predictors.table2["Prior_utah_splicing_reference", "classification"] [predictors.table2["Prior_utah_splicing_reference", "use"]=="yes"& is.na(predictors.table2["Prior_utah_splicing_reference", "classification"])] <- "not applicable"
  predictors.table2["Prior_utah_splicing_de_novo", "classification"] [predictors.table2["Prior_utah_splicing_de_novo", "use"]=="yes"& is.na(predictors.table2["Prior_utah_splicing_de_novo", "classification"])] <- "not applicable"
  return(predictors.table2)
}

predictorResult <- function(table.predictor){
  for ( i in 1:nrow(table.predictor)){
    if (!anyNA(table.predictor[i,c("values", "operator.benign", "operator.pathogenic", "TD.cut.off")])){
      if((table.predictor$values[i] >= as.numeric(table.predictor$TD.cut.off[i]) & table.predictor$operator.pathogenic[i] == ">=")|
         (table.predictor$values[i] < as.numeric(table.predictor$TD.cut.off[i]) & table.predictor$operator.pathogenic[i] == "<")|
         (table.predictor$values[i] > as.numeric(table.predictor$TD.cut.off[i]) & table.predictor$operator.pathogenic[i] == ">")|
         (table.predictor$values[i] <= as.numeric(table.predictor$TD.cut.off[i]) & table.predictor$operator.pathogenic[i] == "<=")|
         (table.predictor$values[i] == as.numeric(table.predictor$TD.cut.off[i]) & table.predictor$operator.pathogenic[i] == "=")){
        table.predictor$classification[i] <- ifelse(table.predictor$type[i]=="Nucleotide conservation","Strongly conserved", "Pathogenic")
        
      }else if ((table.predictor$values[i] <= as.numeric(table.predictor$BE.cut.off[i])& table.predictor$operator.benign[i] =="<=")|
                (table.predictor$values[i] > as.numeric(table.predictor$BE.cut.off[i]) & table.predictor$operator.benign[i] == ">")|
                (table.predictor$values[i] < as.numeric(table.predictor$BE.cut.off[i]) & table.predictor$operator.benign[i] == "<")|
                (table.predictor$values[i] != as.numeric(table.predictor$BE.cut.off[i]) & table.predictor$operator.benign[i] == "!=")|
                (table.predictor$values[i] >= as.numeric(table.predictor$BE.cut.off[i]) & table.predictor$operator.benign[i] == ">=")){
        table.predictor$classification[i] <- ifelse(table.predictor$type[i]=="Nucleotide conservation","Not strongly conserved", "Benign")
      } else{
        table.predictor$classification[i] <-  "Grey"
      }
    }
  }
  return(table.predictor)
}

#' intronic
#' 

intronicCutoff <- function (object){
  intronic <- NA
  if(!(object$most.severe.consequence %in% c("5_prime_UTR_variant"))){
    pos.intronic.signe <-  stringr::str_extract(object$variant, "[-|+]")
    pos.intronic <- stringr::str_extract(object$variant, "[-|+]+[0-9]*") 
    pos.intronic <- as.numeric(stringr::str_extract(pos.intronic, "[0-9]+[0-9]*" ))
    exp.gene <- object$gene %in% c("ATM", "PALB2", "CHEK2", "TP53")
    intronic <- ifelse(exp.gene==FALSE & (pos.intronic.signe=="+"& pos.intronic >= 7 | pos.intronic.signe=="-"&pos.intronic >= 21) |
                         (object$gene=="ATM"& (pos.intronic.signe=="+"& pos.intronic>7 | pos.intronic.signe=="-"&pos.intronic >40))|
                         (object$gene=="PALB2"& (pos.intronic.signe=="+"& pos.intronic>20 | pos.intronic.signe=="-"&pos.intronic >40)),
                       "yes",
                       ifelse(pos.intronic.signe%in% c("+", "-"),
                              "no",
                              NA))
  }
  return (intronic)
} 

predictorsUse <- function (table.predictor, gene.specific, object){
  # Nucleotide conservation 
  #In general it is only taken into account in synonymous variants.
  no.conservation <- object$gene %in% c("ATM", "CDH1", "MLH1", "MSH2", "MSH6", "PM2")
  table.predictor[c("Phylop", "Phastcons"), "use"] <- ifelse(no.conservation, 
                                                             "no",
                                                             ifelse(object$most.severe.consequence=="synonymous_variant",
                                                                    "yes",
                                                                    "no"))
  #Exceptions: it is used as well for:
  #PTEN and MMR -> (A synonymous (silent) or intronic variant at or beyond +7/–21 (5′/3′ exonic) )  
  #PALB2 -> Can be used for deep intronic variants further than +20 (donor) and -40 (acceptor).
  exp.gene <- object$gene %in% c("PTEN")
  intronic <- intronicCutoff (object)
  table.predictor[c("Phylop", "Phastcons"), "use"][intronic=="yes" & (exp.gene|object$gene=="PALB2")] <-"yes"
  table.predictor[c("Gerp"), "use"]<- "no" 
  
  # REVEL is used for missense variants  in all genes except except MLH1, MSH2, MSH6 and PMS2 , PTEN, CDH1,PALB2 , TP53
  sel.gen <- object$gene %in% c("MLH1", "MSH2", "MSH6", "PMS2", "PTEN", "CDH1", "TP53")
  #sel.consequence <- object$most.severe.consequence %in% c( "synonymous_variant", "splice", "frameshift_variant")
  table.predictor["Revel", "use"] <- ifelse ( sel.gen|object$most.severe.consequence !="missense_variant",  "no", "yes")
  
  #Provean is used for inframe_deletions and inframe_insertions
  table.predictor["VEST4", "use"]<-"no" #we set to FALSE 
  table.predictor["Provean", "use"] <- ifelse(object$most.severe.consequence %in% c("inframe_deletion", "inframe_insertion"),
                                              "yes",
                                              "no")
  #polyphen and MAPP are set to no because we use the prior prob that mixes both and variant needs to be missense
  table.predictor[c("PolyPhen","MAPP"), "use"] <- "no"
  MMR <- object$gene %in% c("MLH1", "MSH2", "MSH6", "PMS2")
  table.predictor["Prior_utah(MAPP/PP2)", "use"] <- ifelse ( MMR&object$most.severe.consequence == "missense_variant", "yes", "no")
  
  #BayesDel and aGVGD_zebrafish are used in missense variants for TP53
  table.predictor[c("BayesDel_noAF", "aGVGD_zebrafish"), "use"] <- ifelse ( object$gene=="TP53" & object$most.severe.consequence == "missense_variant", c("yes", "yes"), c("no", "no"))
  
  #Prior_splicing is used for MMR
  table.predictor[c("Prior_utah_splicing_reference", "Prior_utah_splicing_de_novo"), "use"] <- ifelse( MMR, "yes", "no")
  
  #SpliceAI is used for all variants
  table.predictor[c("SpliceAI-AcceptorGain","SpliceAI-AcceptorLoss","SpliceAI-DonorGain","SpliceAI-DonorLoss"), "use"] <- "yes"
  
  #TraP is used for CDH1
  table.predictor["TraP", "use"] <- ifelse(object$gene=="CDH1"&object$most.severe.consequence=="missense_variant", "yes", "no")
  
  return(table.predictor)
}





#' UCSC predictors
#' @param track can be either phyloP100wayAll, phastCons100way or allHg19RS_BW
#' @param chr chromosome of interest
#' @param start start genomic position
#' @param end end genomic position
#' @return the value of the ucsc predictor 
#' @example
#' ucscRPredictor ("phyloP100wayAll", 17, 41258474, 41258474)
#' @author Elisabet Munté

ucscRPredictor <- function(track, object){
  ext.predictor <- ifelse(stringr::str_detect(object$variant, "dup")|stringr::str_detect(object$variant, "ins"), 
                          paste0("genome=hg19;track=", track, ";chrom=chr", object$chr, ";start=", object$end,";end=", object$start),
                          paste0("genome=hg19;track=", track, ";chrom=chr", object$chr, ";start=", object$start-1,";end=", object$end))
  server.ucsc <- "http://api.genome.ucsc.edu/getData/track?" #UCSC's REST API
  pred.score <- api2 (server.ucsc, ext.predictor)
  cromo <- paste0("chr", object$chr)
  if (track == "allHg19RS_BW"){
    pred.score2 <- mean(as.numeric(unlist(pred.score[[track]][[cromo]]["value"])))
  }else{
    pred.score2 <- mean(as.numeric(unlist(pred.score[[cromo]]["value"])))
  }
  pred.score2[is.null(pred.score2[[1]])]<- NA
  return(pred.score2)
}

#AlignGvgd
alignGvgd <- function(object, bbdd){
  if(object$most.severe.consequence=="missense_variant" & object$gene=="TP53"){
    prot <- toProtein(object$protein)
    alignGvgd <- bbdd$align.gvd %>% 
                                stringr::str_replace("Class C", "") %>% 
                                as.numeric()
  }else{
    alignGvgd <- NA
  }
  return(alignGvgd)
}

#' TrapScore
#' @references Gelfman, S., Q. Wang, K.M McSweeney, Z. Ren, F. La Carpia, M. Halvorsen, K. Schoch, F. Ratzon, E.L. Heinzen, M.J. Boland, S. Petrovski and D. B. Goldstein. 2017. Annotating Pathogenic Non-Coding Variants in Genic Regions. Accepted June 5th, 2017. 
traPscore <- function(gnom){
  url <- paste0("http://trap-score.org/Search?query=", gnom)
  traPscore <- readUrl(url) 
  if(!is.na(traPscore)){
    traPscore <- traPscore %>% 
                           rvest::html_table()  %>% 
                           purrr::map(8) %>% 
                           unlist()
  }
  return(traPscore)
}

dbnsfpDbQuery <- function(object, ensembl.id, bbdd){
  dbnsfp.variant <- bbdd$dbnsfp
  revel.score <- dbnsfp.variant[1,"REVEL_score"] %>%
                                                 as.numeric()#the value
  provean.score <- dbnsfp.variant[1,"PROVEAN_score"] %>%
                                                 as.numeric()#the value
  vest4.score <- dbnsfp.variant[1,"VEST4_score"] %>%
                                                 as.numeric()#the value
  bayesDel.score <-dbnsfp.variant[1,"BayesDel_noAF_score"] %>%
                                                  as_numeric()#the value
  
  dbnsfp.scores <- list(revel=revel.score, 
                        vest4=vest4.score, 
                        provean=provean.score, 
                        bayesDel=bayesDel.score)
  return(dbnsfp.scores)
}


#Prior Utah
priorUtahProb <- function(object, gene=NULL, variant =NULL){
  if(is.null(object)) {
    object <- data.frame(gene=gene, variant=variant)
  }
  query <- paste0("SELECT Polyphen, MAPP, prior  from prior_db WHERE gene= '", object$gene ,"'AND cdna='",object$variant,"';")
  prior <- connectionDB(query)[[1]]
  if (nrow(prior)==0){
    prior <- rep(NA,3)
  }
  return(list(prior[1], prior[2], prior[3]))
}

proveanR <- function(object){
  prot.cor <- object$protein
  score.provean <- NA
  if(object$most.severe.consequence %in% c("missense_variant", "synonymous_variant")){
    prot <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)"),2)  %>% 
                                                                   stringr::str_extract_all( "[A-z]+") %>% 
                                                                   unlist
    prot.num <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)"),2)  %>% 
                                                                       stringr::str_extract_all("[0-9]+") %>% 
                                                                       unlist
    protein <- paste0(aaShort(prot[1]), prot.num,  aaShort(prot[2]))
  }else if (object$most.severe.consequence %in% c("inframe_deletion")){
    prot1 <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|del"),2)  %>% 
                                                                          stringr::str_extract_all( "[A-z]+") %>% 
                                                                          unlist
    prot2 <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|del"),3)  %>% 
                                                                          stringr::str_extract_all( "[A-z]+") %>% 
                                                                          unlist
    
    prot.num1 <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|del"),2) %>% 
                                                                             stringr::str_extract_all("[0-9]+") %>% 
                                                                             unlist
    prot.num2 <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|del"),3) %>% 
                                                                             stringr::str_extract_all("[0-9]+") %>% 
                                                                             unlist
    protein <- ifelse(stringr::str_detect(prot.cor, "_"),
                      paste0(aaShort(prot1), prot.num1, "_", aaShort(prot2), prot.num2, "del"),
                      paste0(aaShort(prot1), prot.num1, "del"))
    
  }else if (object$most.severe.consequence %in% c("inframe_insertion")){
    if(stringr::str_detect(object$protein, "dup")){
      if(stringr::str_detect(object$protein, "_")){
        prot1 <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|dup"),2)  %>% 
                                                                              stringr::str_extract_all( "[A-z]+") %>% 
                                                                              unlist
        prot2 <-  purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|dup"),3)  %>% 
                                                                              stringr::str_extract_all( "[A-z]+") %>% 
                                                                              unlist
        num1 <-  purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|dup"),2)  %>% 
                                                                              stringr::str_extract_all( "[0-9]+") %>% 
                                                                              unlist
        num2 <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|dup"),3)  %>% 
                                                                             stringr::str_extract_all( "[0-9]+") %>% 
                                                                             unlist
        protein <- paste0(aaShort(prot1), num1, "_", aaShort(prot2), num2, "dup")
        
      }else{
        prot <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)"),2)  %>% 
                                                                       stringr::str_extract_all( "[A-z]+") %>% 
                                                                       unlist
        prot.num <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)"),2)  %>% 
                                                                           stringr::str_extract_all("[0-9]+") %>% 
                                                                           unlist
        protein <- paste0(aaShort(prot[1]), prot.num,  prot[2])
      }
      
    }else if(stringr::str_detect(object$protein, "ins")){
      prot1 <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|ins"),2)  %>% 
                                                                            stringr::str_extract_all( "[A-z]+") %>% 
                                                                            unlist
      prot2 <-  purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|ins"),3) %>% 
                                                                            stringr::str_extract_all( "[A-z]+") %>% 
                                                                            unlist
      prot3 <-  purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|ins"),4) %>% 
                                                                            unlist()
      prot.ins <- NULL
      for(i in 1:stringr::str_length(prot3)){
        proti <- aaShort(stringr::str_sub(prot3, i, i +2))
        prot.ins <- paste0(prot.ins, proti)
        i <- i+2
      }
      prot.num1 <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|ins"),2) %>% 
                                                                               stringr::str_extract_all("[0-9]+") %>% 
                                                                               unlist
      prot.num2 <- purrr::map(stringr::str_split(prot.cor, "\\(|\\)|_|ins"),3) %>% 
                                                                               stringr::str_extract_all("[0-9]+") %>% 
                                                                               unlist
      protein <- paste0(aaShort(prot1), prot.num1, "_", aaShort(prot2), prot.num2, "ins", prot.ins)
    }
  }else{
    protein <- NA
  }
  
  ##Web scrapping -> web-based version has been retired
  # if(!is.na(protein)){
  #   server.ensembl <- "http://grch37.rest.ensembl.org" 
  #   ext.ensembl <- paste0("/sequence/id/", ensemblTranscript(object$NM, object$gene)$id ,"?multiple_sequences=0;content-type=text/x-seqxml%2Bxml;type=protein")
  #   prot.seq <- api(server.ensembl, ext.ensembl)$seq
  #   rD <- RSelenium::rsDriver(browser=c("firefox"), port=4570L)
  #   driver <- rD[["client"]]
  #   url.provean <- "http://provean.jcvi.org/seq_submit.php"
  #   driver$navigate(url.provean)
  #   Sys.sleep(3)
  #   
  #   
  #   webElem <- driver$findElement(using = "css", "[name = 'query']")
  #   webElem$sendKeysToElement(list(prot.seq))
  #   webElem2 <- driver$findElement(using = "css", "[name = 'variant']")
  #   webElem2$sendKeysToElement(list(protein))
  #   webElem3 <- driver$findElement("xpath", '//input[@type="submit"]')
  #   webElem3$clickElement()
  #   
  #   windHand <- driver$getWindowHandles()
  #   driver$switchToWindow(windHand[[2]])
  #   information <- driver$getPageSource()[[1]] %>% read_html() %>% html_table()
  #   while(length(information)==0){
  #     Sys.sleep(5)
  #     information <- driver$getPageSource()[[1]] %>% read_html() %>% html_table()
  #   }
  #   score.provean <- information[[1]]$X2[4] %>% as.numeric()
  #   driver$close()
  #   rD[["server"]]$stop()
  # }
  # 
  # 
  
  
  return(score.provean)
}



spliceaiR <- function(object, ext.spliceai, genome = 37, distance = 1000, precomputed = 1, mask = 1, bbdd, spliceai.program = FALSE, .dat = "../docs", reference.splice = NULL, annotation.splice = system.file("docs", "gencode_spliceai_hg19.txt", package="vaRHC")) {
  assertthat::assert_that(is.logical(spliceai.program), msg="spliceai.program param must be TRUE or FALSE")
  cond <- precomputed==1 & genome==37 & distance==1000 & mask==1
  ensembl.id <- object$ensembl.id
  if(cond==T){
    spliceai.score <- bbdd$spliceai
    splice.return <- list(Acceptor_Gain = spliceai.score[1,c("AG", "AG_dis")] %>% as.matrix, 
                          Acceptor_Loss = spliceai.score[1, c("AL", "AL_dis")] %>% as.matrix,
                          Donor_Gain = spliceai.score[1,c("DG", "DG_dis")] %>% as.matrix,
                          Donor_Loss = spliceai.score[1,c("DL", "DL_dis")] %>% as.matrix)
    
  }
  
  if((cond==F || (cond==T && nrow(spliceai.score)==0)) && spliceai.program == TRUE){
    assertthat::assert_that(!is.null(reference.splice), msg="Reference file must be provided if spliceAI has to be computed")
    assertthat::assert_that(file.exists(reference.splice) & stringr::str_detect(reference.splice, ".fa"), msg="Please enter a valid reference file")
    assertthat::assert_that(!is.null(annotation.splice), msg="Reference file must be provided if spliceAI has to be computed")
    assertthat::assert_that(file.exists(annotation.splice), msg="Please enter a valid annotation file")
    
    #NM exeptions
    if(object$NM=="NM_007194.3"){#CHEK2
      ensembl.id <- api ("http://rest.ensembl.org", "/xrefs/symbol/homo_sapiens/NM_007194.4?content-type=application/json") %>% 
                                                                                                                            dplyr::filter(type=="transcript") %>%
                                                                                                                            dplyr::select(id) %>% 
                                                                                                                            unlist() %>%
                                                                                                                            as.character()
    } else if (object$gene=="MUTYH"){ #MUTYH
      ensembl.id <- api ("http://rest.ensembl.org", "/xrefs/symbol/homo_sapiens/NM_001128425.2?content-type=application/json") %>% 
                                                                                                                               dplyr::filter(type=="transcript") %>%
                                                                                                                               dplyr::select(id) %>%
                                                                                                                               unlist() %>%
                                                                                                                              as.character()
    }else if (object$gene=="ERCC5"){ #ERCC5 
      ensembl.id <- api ("http://rest.ensembl.org", "/xrefs/symbol/homo_sapiens/NM_000123.4?content-type=application/json") %>% 
                                                                                                                            dplyr::filter(type=="transcript") %>%
                                                                                                                            dplyr::select(id) %>%
                                                                                                                            as.character()
    }else if (object$NM=="NM_004448.2"){ #ERBB2
      ensembl.id <- api ("http://rest.ensembl.org", "/xrefs/symbol/homo_sapiens/NM_004448.4?content-type=application/json") %>% 
                                                                                                                            dplyr::filter(type=="transcript") %>%
                                                                                                                            dplyr::select(id) %>%
                                                                                                                            unlist() %>%
                                                                                                                            as.character()
    }
    
    #SPLICEAI PROGRAM
    #temporarly directory to store spliceAI output
    .tmp <- file.path(getwd(), ".tmp")
    dir.create(.tmp, showWarnings = FALSE)
    
    var.tros <- stringr::str_split(ext.spliceai, "-") %>% 
      unlist()
    var.tros2 <- c(var.tros[1:2], object$variant, var.tros[3:4], ".", ".", ".")  %>% 
      as.data.frame() %>% 
      t()
    time <-  Sys.time() %>% 
      stringr::str_replace_all("-|:| ", "_") 
                           
    write.table(var.tros2, 
                file  = file.path(.tmp, paste0("var_splice", time, ".txt")), 
                col.names=FALSE, 
                row.names=FALSE, 
                sep="\t", 
                quote=FALSE)
    vcf.header <- system.file("docs", "vcf_spliceAI.vcf", package="vaRHC")
    cmd <- paste0("cat ", vcf.header, " ", file.path(.tmp, paste0("var_splice", time, ".txt"))," > ", file.path(.tmp, paste0("spliceAI_R",time,".vcf")))
    print(cmd); try(system(cmd))
    inputVCF <- file.path(.tmp, paste0("spliceAI_R",time,".vcf"))
    outputFile <-file.path(.tmp, paste0("outputR",time,".vcf"))
    
    run.splice <- system.file("docs", "runSpliceAI.sh", package="vaRHC")
    cmd2 <- paste(run.splice, "spiceAI_env/bin/activate", inputVCF, outputFile, distance, mask, reference.splice, annotation.splice)
    print(cmd2); try(system(cmd2))
    vcf <- vcfR::read.vcfR(outputFile, "hg19")
    splice <- vcf@fix [,8] %>% 
      stringr::str_split (paste0(",[A-z]\\|", object$gene, "| [A-Z]\\|", object$gene, "|", object$gene)) %>% 
      unlist()
    #delete intermediate files and directory
    unlink(.tmp, recursive=TRUE)
    spliceai.score <- splice[stringr::str_detect(splice, ensembl.id)==T]
    spliceai.score <- stringr::str_split(spliceai.score, "[|]") %>% 
      unlist() 
    spliceai.score <- spliceai.score[2:9]
    spliceai.score <- matrix (spliceai.score, nrow=4)
   
    splice.return <- list(Acceptor_Gain = spliceai.score[1,], 
                          Acceptor_Loss = spliceai.score[2,],
                          Donor_Gain = spliceai.score[3,],
                          Donor_Loss = spliceai.score[4,])
  } else if (spliceai.program == FALSE && (cond==FALSE  || (cond==T &&nrow(spliceai.score)==0))){
    stop("This variant has not a precalculated score, spliceai program should be installed to calculate it")
  }
  return(splice.return)
}



insilico <- function (predictors.df, predictor.type, veredict ){
  in.silico <- predictors.df %>% 
    dplyr::filter(type %in% predictor.type, 
                  stringr::str_detect(classification, veredict), 
                  use %in% "yes")
  return (in.silico)
}
