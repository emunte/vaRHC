

################################################################################
## PVS1
################################################################################
PVS1 <- function (all.information, final.criteria){
  PVS1.message <- NULL
  spli.loss.pvs1.df <- all.information$predictors$predictor.table %>%
    dplyr::filter(grepl("Loss", .data$predictor),
                  .data$classification %in% c("Pathogenic" ,"Grey"))
  spli.gain.pvs1.df <- all.information$predictors$predictor.table %>%
    dplyr::filter(grepl("Gain", .data$predictor),
                  .data$classification == "Pathogenic")
  gene <- all.information$Variant.Info$gene
  #For frameshifts, nonsense and canonicals(except from CDH1 and ATM splice donor and acceptors)
  ###1. SpliceAi gains are pat -> supporting check manually
  sel.vars <- all.information$Variant.Info$most.severe.consequence %in% c("frameshift_variant", "stop_gained", "splice_donor_variant", "splice_acceptor_variant")
  final.criteria$criteria.res["PVS1", 1:4][nrow(spli.gain.pvs1.df)>0 & sel.vars ] <- c(0,0,0,1)
  PVS1.message[nrow(spli.gain.pvs1.df)>0 & sel.vars] <- paste("SpliceAI predicts a Gain, PVS1_supporting is assigned because it is a", all.information$Variant.Info$most.severe.consequence, "but please check PVS1 strength manually.")

  ##Last exon variants, only ATM, CDH1 and MMR
  if(gene %in% "CDH1" && !is.na(all.information$Variant.Info$most.severe.consequence.1) && all.information$Variant.Info$most.severe.consequence.1 %in% c("splice_region_variant")){
    pos <- stringr::str_extract(all.information$Variant.Info$variant, "[0-9]+")
    donor <- all.information$codon.stop$exons %>%
      dplyr::filter(.data$cStop==pos) %>%
      nrow()
    final.criteria$criteria.res["PVS1", 1:4][donor>0&all.information$Variant.Info$ref=="G"] <- c(0,0,1,0)
    PVS1.message[donor>0 & all.information$Variant.Info$ref=="G"] <- "PVS1_moderate is assigned because it is a G to non-G variant disrupting the last nucleotide of an exon."
  }else if (gene %in% c("ATM", "MLH1", "MSH2", "MSH6", "PMS2")&& !is.na(all.information$Variant.Info$most.severe.consequence.1) && all.information$Variant.Info$most.severe.consequence.1=="splice_region_variant" && !(all.information$Variant.Info$most.severe.consequence%in% c("intron_variant"))){
    last.nt <- all.information$class.info$last.nt
    final.criteria$criteria.res["PVS1", 1:4][last.nt$PVS1=="strong"] <- c(0,1,0,0)
    PVS1.message <- ifelse(nrow(last.nt)!=0&&last.nt$PVS1=="strong",
                           paste("PVS1_strong is assigned because the variant is located at the last nt of the exon and is a G>non-G and the first 6 bases of the intron are not GTRRGT (the bases are", last.nt$intron6, ")"),
                           ifelse(nrow(last.nt)!=0&&last.nt$PVS1=="no",
                                  paste(PVS1.message, "PVS1_strong is not assigned because although the variant is located at the last nt of the exon or the nucleotide is not a G (the nucleotide is", last.nt$nt, ") or because the first 6 bases of the intron are  GTRRGT (the bases are", last.nt$intron6, ")" ),
                                  NA))


  }


  #IF NO GAIN (benign or grey), FOR CANONICAL VARIANTS WE CHECK THAT THE SKIPPING IS AT THE FIRST OR LAST EXON (CODING OR NON-CODING)
  # A) First or last exon -> calculate manually
  # B) not located first or last : check if:
  # I. SpliceAI losses are pat or grey and is in-frame deletion
  #a) prot > 0.1 -> moderate
  #b) prot < 0.1 -> supporting+

  # II. Out of frame deletion : will enter in Tayoun WF for non-sense/frameshifts
  if (all.information$Variant.Info$most.severe.consequence %in% c("splice_donor_variant", "splice_acceptor_variant")){
    if(gene=="CDH1"){
      strength <- all.information$class.info$pvs1.cdh1 %>%
        as.character()
      PVS1.message <- paste0("PVS1_", strength, " is assigned according to site-specific recomendations in the splicing table (CDH1 Guidelines Version 3).")
      final.criteria$criteria.res["PVS1", 1:4] <- ifelse(rep(strength=="strong",4),
                                                         c(0,1,0,0),
                                                         c(0,0,1,0))
    }else if (gene=="ATM"){
      strength.query <- all.information$class.info$pvs1.atm
      PVS1.message <- paste0(strength.query$PVS1_strength, "is assigned because", strength.query$reasoning)
      final.criteria$criteria.res["PVS1", 1:4] <- ifelse(rep(strength.query$PVS1_strength=="PVS1",4),
                                                         c(1,0,0,0),
                                                         ifelse(rep(strength.query$PVS1_strength=="PVS1_strong",4),
                                                                c(0,1,0,0),
                                                                ifelse(rep(strength.query$PVS1_strength=="PVS1_supporting",4),
                                                                       c(0,0,0,1),
                                                                       c(0,0,0,0))))
    }else if (nrow(spli.gain.pvs1.df)==0){
      #A
      if ((all.information$codon.stop$variant.exon$exon ==1 |
           all.information$codon.stop$variant.exon$exon == nrow(all.information$codon.stop$exons)) |
          stringr::str_detect(all.information$codon.stop$exons$cStop[all.information$codon.stop$exons$exon==all.information$codon.stop$variant.exon$exon], "-") |
          stringr::str_detect(all.information$codon.stop$exons$cStart[all.information$codon.stop$exons$exon==all.information$codon.stop$variant.exon$exon], "-")){

        PVS1.message <-  "PVS1 must be calculated manually: Skipping of the first (coding) or last exon cannot be calculated automatically. "

      }else
      {  #B
        if (nrow(spli.loss.pvs1.df)>0 &
            all.information$codon.stop$canonical.skip.pred$most.severe.consequence == "inframe_deletion")
        {
          sel.prot <- all.information$codon.stop$canonical.skip.pred$porc.prot.splicing > 0.1
          final.criteria$criteria.res["PVS1", 1:4] <- ifelse (rep(sel.prot == TRUE,4) ,
                                                              c(0, 1, "NC", "NC"),
                                                              c(0, 0, 1, "NC"))
          PVS1.message <- ifelse (rep(sel.prot == TRUE),
                                  paste("PVS1_strong is met because the predicted exon skipping produced by the variant is inframe and the % of the protein lost is", round(as.numeric(all.information$codon.stop$canonical.skip.pred$porc.prot.splicing), 3), "which is > 0.1."),
                                  paste("PVS1_moderate is met because the predicted exon skipping produced by the variant is inframe and the % of the protein lost is",round( as.numeric(all.information$codon.stop$canonical.skip.pred$porc.prot.splicing),3), "which is <0.1. Check if it is a critical domain in order to achieve PVS1_strong."))
        }
      }
    }


  }

  if (is.null(PVS1.message)||is.na(PVS1.message)) PVS1.message <- "PVS1 does not apply for this variant."
  final.criteria[["PVS1.message"]] <- PVS1.message

  final.criteria <- NMD(all.information, final.criteria, spli.loss.pvs1.df, spli.gain.pvs1.df)
  final.criteria <- startVariants(all.information, final.criteria)

  return(final.criteria)
}
NMD <- function (all.information, final.criteria, spli.loss.pvs1.df,  spli.gain.pvs1.df){
  PVS1.message <- NA
  type <- c("frameshift_variant", "stop_gained")
  gene <- all.information$Variant.Info$gene
  #sel <- any(type %in% c(all.information$Variant.Info$most.severe.consequence, all.information$codon.stop$canonical.skip.pred))
  sel1 <- type %in% all.information$Variant.Info$most.severe.consequence
  sel2 <- type %in% all.information$codon.stop$canonical.skip.pred & nrow(spli.loss.pvs1.df) >0 & !(gene%in% c("CDH1"))  #for canonical spliceAI has to indicate a loss


  if ((any(sel1 == TRUE) | any(sel2 == TRUE)) & nrow(spli.gain.pvs1.df)==0 ) {
    codon <- ifelse(all.information$Variant.Info$most.severe.consequence %in% c("frameshift_variant", "stop_gained"),
                    ifelse(stringr::str_detect(all.information$Variant.Info$protein, "_"),
                           as.numeric(unlist(stringr::str_extract(all.information$Variant.Info$protein, "[0-9]+"))),
                           sum(as.numeric(unlist(stringr::str_extract_all(all.information$Variant.Info$protein, "[0-9]+"))))),
                    sum(as.numeric(unlist(stringr::str_extract_all(all.information$codon.stop$canonical.skip.pred$protein, "[0-9]+")))))

    if (!(gene %in% c("MLH1", "MSH2", "MSH6", "PMS2", "PTEN", "CDH1", "ATM"))){
      NMD <- ifelse(all.information$codon.stop$premature.ter.codon$exon<nrow(all.information$codon.stop$exons)-1|
                      (all.information$codon.stop$premature.ter.codon$exon==nrow(all.information$codon.stop$exons)-1&&all.information$codon.stop$premature.ter.codon$bp>50), "YES", "NO")
    }else {


      NMD <- ifelse((gene=="MLH1"&codon > 756)|
                      (gene=="MSH2"&codon > 933)|
                      (gene=="MSH6"&codon > 1360) |
                      (gene=="PMS2"&codon > 855) |
                      (gene=="PTEN"&codon > 373)|
                      (gene=="CDH1"&codon > 882) |
                      (gene=="ATM"&codon >3047),
                    "NO_effect",
                    ifelse (gene=="CDH1" &codon >796 & codon <= 836,
                            "PVS1_strong",
                            ifelse((gene=="MLH1" &codon>753 & codon <=756)|
                                     (gene=="MSH2" &codon>891 & codon <=933)|
                                     (gene=="MSH6" &codon>1341 & codon <=1360)|
                                     (gene=="PMS2" &codon>798 & codon <=855)|
                                     (gene=="CDH1" &codon>836 & codon <= 882),
                                   "PVS1_moderate",
                                   ifelse((gene=="ATM" & codon>2979 & codon <= 3047),
                                          "PVS1",
                                          "YES"))))
    }

    NMD[gene=="BRCA1"&(all.information$codon.stop$variant.exon$exon %in% c(8,9))] <- "NO_alternate_transcript"
    final.criteria$criteria.res["PVS1", 1:4] [NMD=="YES" | NMD=="PVS1"] <- c(1, 0,0,0)
    final.criteria$criteria.res["PVS1", 1:4] [NMD=="NO_effect"] <- c(0,0,0,0)
    final.criteria$criteria.res["PVS1", 1:4][NMD=="PVS1_moderate"] <- c(0, 0,1,0)
    final.criteria$criteria.res["PVS1", 1:4][NMD=="PVS1_strong"] <- c(0, 1,0,0)
    final.criteria$criteria.res["PVS1", 1:4][NMD=="No_alternate_transcript"] <- c(0, 0,0,0)

    if (NMD == "NO"){
      porc.prot <- ifelse (!is.na(all.information$codon.stop$porc.prot),
                           all.information$codon.stop$porc.prot,
                           all.information$codon.stop$canonical.skip.pred$porc.prot.splicing )

      sel.prot <- porc.prot > 0.1
      final.criteria$criteria.res["PVS1", 1:4] <- ifelse (rep(sel.prot == TRUE,4) ,
                                                          c(0, 1, "NC", "NC"),
                                                          c(0, 0, 1, "NC"))
      ##For ATM variants in the FATKIN domain (codons 1893-3056) its always PVS1_strong (when NMD==NO) no matter the % of protein
      final.criteria$criteria.res["PVS1", 1:4][gene=="ATM" & codon>=1893 & codon<=3056] <- c(0,1,0, "NC")
    }
    PVS1.message[NMD=="YES"] <-  paste("It exists NMD so PVS1 is met with a very strong strength. A stop codon is generated", all.information$codon.stop$premature.ter.codon$num.codon.stop, "codons since the variant, at the exon", all.information$codon.stop$premature.ter.codon$exon, "out of", nrow(all.information$codon.stop$exons),".")
    PVS1.message[NMD=="NO" && porc.prot > 0.1] <- paste("It does not exist NMD so PVS1 is not met with a very strong strength  because the stop codon is generated", all.information$codon.stop$premature.ter.codon$num.codon.stop, "codons since the variant at exon ", all.information$codon.stop$premature.ter.codon$exon, " out of", nrow(all.information$codon.stop$exons), ". PVS1_strong is met because variant removes >10% of protein.")
    PVS1.message[NMD=="NO" && porc.prot < 0.1] <- paste("It does not exist NMD so PVS1 is not met with a very strong strength  because the stop codon is generated", all.information$codon.stop$premature.ter.codon$num.codon.stop, "codons since the variant at exon ",all.information$codon.stop$premature.ter.codon$exon, " out of", nrow(all.information$codon.stop$exons), ". PVS1_moderate is met because variant removes <10% of protein.")
    PVS1.message[NMD== "No_alternate_transcript"] <- paste ("The stop codon is generated", all.information$codon.stop$premature.ter.codon$num.codon.stop, "codons since the variant at exon", all.information$codon.stop$premature.ter.codon$exon, "out of", nrow(all.information$codon.stop$exons), " There is an alternate transcript without exons 8 and 9 so PVS1 is not met with a very strong strength.")
    PVS1.message[NMD== "NO_effect"] <- paste("It does not exist NMD so PVS1 is not met because the codon stop is at position", codon, ".")
    PVS1.message[NMD== "PVS1_moderate"] <- paste("It  exist NMD but PVS1 ony mets a PVS1_moderate strength because the codon stop is at position", codon, ".")
    PVS1.message[NMD== "PVS1_strong"] <- paste("It  exist NMD but PVS1 ony mets a PVS1_strong strength because the codon stop is at position", codon, ".")
    final.criteria[["PVS1.message"]] <- PVS1.message
    final.criteria[["NMD"]] <- NMD

  }

  return(final.criteria)

}
startVariants <- function (all.information, final.criteria){
  message.second.met <- NULL
  pos.second.met <- all.information$second.met$start.lost.variants$pos.second.met
  if (!is.na(pos.second.met)){
    pathogenic.second.met <- all.information$second.met$start.lost.variants$clinvar.pvs1.info %>%
      dplyr::filter(.data$clinSign=="Pathogenic" & .data$reviewStatus=="reviewed by expert panel")
    #General rules
    gene <- all.information$Variant.Info$gene
    final.criteria$criteria.res["PVS1", 1:4] <- ifelse( rep(nrow(pathogenic.second.met)>0, 4),
                                                        c(0,0,1,"NC"), c(0,0, "NC", "NC"))
    message.second.met<- ifelse( nrow(pathogenic.second.met)>0,
                                 paste("Mutation at the start codon of gene", gene,". PVS1_moderate is met because there are", nrow(pathogenic.second.met),
                                       "variant(s) upstream of closest potential in-frame start codon classified in ClinVar by expert Panel ( next metionine is at codon",
                                       pos.second.met, ")"),
                                 paste("Mutation at the start codon of gene", gene,". PVS1_moderate is not met because there are", nrow(pathogenic.second.met),
                                       "variants usptream of closest potential in-frame start codon classified in ClinVar by expert Panel( next metionine is at codon",
                                       pos.second.met, "). The criteria is set to NC because you should check the variants not classified by expert panel and if different functional transcripts use alternative start codon."))

    # Specific rules
    very_strong <- gene %in% c("MLH1", "CDH1", "PTEN", "ATM")
    strong <- gene %in% c("PMS2", "MSH6", "BAP1")
    no.criteria <- gene %in% c("MSH2")
    final.criteria$criteria.res["PVS1", 1:4][very_strong] <- c(1,0, "NC", "NC")
    final.criteria$criteria.res["PVS1", 1:4][strong] <- c(0,1, "NC", "NC")
    final.criteria$criteria.res["PVS1", 1:4][no.criteria] <- rep(0,4)

    message.second.met[very_strong] <- paste("Mutation at the start codon of gene", gene, ", so PVS1 is met with a very strong strength.")
    message.second.met[strong] <- paste("Mutation at the start codon of gene", gene,". Further ATGs exist inframe in exon 1, so PVS1 is not met.")
    message.second.met[no.criteria] <- paste("Mutation at the start codon of gene", gene, ", so PVS1 is met with a strong strength.")
    final.criteria[["PVS1.message"]] <- message.second.met
  }


  return (final.criteria)
}

################################################################################
## PS1
################################################################################

PS1 <- function (all.information, final.criteria){
  PS1.message <- NULL
  gene <- all.information$Variant.Info$gene
  #server.mutalyzer <- "https://mutalyzer.nl/json/" #Mutalyzer's REST API
  final.criteria$criteria.res["PS1", c(1,3,4)]<- NA
  if (all.information$Variant.Info$most.severe.consequence=="missense_variant" & ! (gene %in% c("CDH1", "PALB2"))){
    #clinvar.evidence.calc <- calculate_clinvar_evidences(all.information=all.information)

    clinvar.evidence.calc <- clinVarEvidence(all.information)

    if(length(clinvar.evidence.calc)>0){
      rownames(clinvar.evidence.calc) <- c("ben_guideline", "pben_guideline", "unc_guideline",  "ppat_guideline", "pat_guideline",  "ben_expert",     "pben_expert",    "unc_expert",     "ppat_expert",    "pat_expert",     "ben_submitter",  "pben_submitter", "unc_submitter", "ppat_submitter", "pat_submitter")
      #server.mutalyzer <- "https://mutalyzer.nl/json/"
      ps1.df <- clinvar.evidence.calc %>%
        t() %>% as.data.frame %>%
        dplyr::mutate(variants.name=colnames(clinvar.evidence.calc)) %>%
        dplyr::filter (stringr::str_detect(.data$variants.name,all.information$Variant.Info$protein), !stringr::str_detect(.data$variants.name, all.information$Variant.Info$variant) & ((all.information$Variant.Info$gene!="TP53" & (.data$pat_guideline>0|.data$ppat_guideline>0|.data$pat_expert>0|.data$ppat_expert>0))| (all.information$Variant.Info$gene=="TP53" &(.data$pat_guideline>0|.data$pat_expert>0))))  %>%
        dplyr::mutate (ext_vep_ps1 = paste0("/vep/human/hgvs/", all.information$Variant.Info$NM,":", stringr::str_sub(purrr::map(stringr::str_split(.data$variants.name, "\\(|\\):"),3),1,-2)))


      if (nrow(ps1.df)>0){
        NC <- stringr::str_extract(all.information$Variant.Info$genomic, "NC_[0-9]+.[0-9]+")
        ps1.df <- ps1.df %>%
          dplyr::rowwise %>%
          dplyr::mutate(genomic=toGenomic(all.information$Variant.Info$NM, NC, stringr::str_sub(purrr::map(stringr::str_split(.data$variants.name, "\\(|\\):"),3),1,-2) ,all.information$Variant.Info$gene) ) %>%
          dplyr::mutate(ext.spliceai.ps1= paste0(all.information$Variant.Info$chr, "-", purrr::map(stringr::str_split(.data$genomic, "g\\.|[A-Z]"),4),"-", stringr::str_sub(.data$genomic, -3,-3), "-", stringr::str_sub(.data$genomic,-1,-1) )) %>%
          dplyr::mutate(scores.spliceai =  connectionDB( paste0("SELECT *  from spliceAI  WHERE var_chr= '", .data$ext.spliceai.ps1 ,"'AND max_dis= 1000 AND transcript='", all.information$Variant.Info$ensembl.id,"' AND masked='",TRUE,"';"))%>% purrr::map(.data, function(x) c(x[[8]], x[[10]],x[[12]],x[[14]]))) %>%
          dplyr::mutate(score.spliceai= max(unlist(.data$scores.spliceai), na.rm=FALSE))

        ps1.no.splicing <- ps1.df %>%
          dplyr::filter (as.numeric(.data$score.spliceai ) < 0.5)



        #gene_specificities
        #MMR
        if (gene %in% c("MLH1", "MSH2", "MSH6", "PMS2")){
          ps1.mmr.pat <- ps1.df %>%
            dplyr::filter (as.numeric(.data$score.spliceai)< 0.5, (.data$pat_guideline>0| .data$pat_expert>0))
          final.criteria$criteria.res["PS1", "strong"] <- ifelse(nrow(ps1.mmr.pat)>0,
                                                                 1,
                                                                 0)
          PS1.message[nrow(ps1.mmr.pat)>0] <-  paste("PS1 is assigned because it exists the variant", ps1.mmr.pat$variants.name, "in ClinVar classified as pathogenic by at least expert panel and it is not predicted to alter splicing.")
          if (nrow(ps1.mmr.pat)==0){
            ps1.mmr.ppat <- ps1.df %>%
              dplyr::filter (as.numeric(.data$score.spliceai)< 0.5, 
                             (.data$ppat_guideline>0|.data$ppat_expert>0))
            final.criteria$criteria.res["PS1", "moderate"] <- ifelse(nrow(ps1.mmr.ppat)>0,
                                                                     0,
                                                                     0)
            PS1.message <- ifelse(nrow(ps1.mmr.ppat)>0,
                                  paste("PS1_moderate should be considered because the variant ", ps1.mmr.ppat$variants.name, " is classified in ClinVar as ppat by expert panel. Your should check if the variant has normal RNA result and if is absent from appropiate population control references groups.  Moreover,  be careful because the variant cannot exceed Class 4 (if this criterion is necessary in order to achieve class 5)."),
                                  paste("PS1_moderate is not assigned because although it exists the variant", ps1.df$variants.name, "classified as PAT/pPAT, the spliceai value suggests a deffect on splicing value (", ps1.df$score.spliceai,">=0.5)."))
          }
        }else{
          if(nrow(ps1.no.splicing)>0){
            PS1.message <- ifelse(!(gene %in% c("TP53")),
                                  paste("PS1 is assigned because the variant", ps1.no.splicing$variants.name, "is classified in ClinVar as PAT or pPAT by at least expert panel and it is not predicted to alter splicing."),
                                  paste("PS1_moderate is assigned because it exists the variant", ps1.no.splicing$variants.name, "is classified in ClinVar as pathogenic  by at least expert panel and it is not predicted to alter splicing."))
            final.criteria$criteria.res["PS1", c("strong","moderate","supporting")] <- ifelse(rep(!(gene %in% c("TP53")),3),
                                                                                              c(1,0,0) ,
                                                                                              c(0,1,0))
            if(gene %in% c("ATM", "CHEK2")){
              our.score <- max(as.numeric(all.information$predictors$predictor.table[c("SpliceAI-AcceptorGain", "SpliceAI-AcceptorLoss", "SpliceAI-DonorGain", "SpliceAI-DonorLoss"), "values"]))
              sel <- our.score >= as.numeric(ps1.df$score.spliceai)
              final.criteria$criteria.res["PS1", c("strong", "moderate")] <- ifelse(rep(sel==TRUE,2),
                                                                                    c(0,1),
                                                                                    c(0,0))
              PS1.message <- ifelse(sel==TRUE,
                                    paste("PS1_moderate is assigned because although the variant reported in ClinVar is predicted to alter splicing, its max spliceai value is", as.numeric(ps1.df$score.spliceai), "which is less or equal than the max spliceai value of our variant of interest(", our.score, "). Please check if the variant we are comparing to has a confirmed splice defect, if not, it should not be applied."),
                                    paste ("PS1_moderate is not assigned because although there exist a variant in ClinVar, its max spliceai value is", as.numeric(ps1.df$score.spliceai), "which is greather than the score of our variant of interest (", our.score, ")."))
            }
          }else {
            final.criteria$criteria.res["PS1", c("strong","supporting")] <- c(0,0)
            PS1.message <- paste("PS1 is not assigned because although it exists the variant", ps1.df$variants.name, "classified as pat/ppat, the spliceai value suggests a deffect on splicing value (", ps1.df$score.spliceai,">=0.2) so it cannot be taken into account. ")
          }
        }
      } else{
        final.criteria$criteria.res["PS1", c("strong","supporting")] <- c(0,0)
        PS1.message <- "PS1 is not assigned because no pat/ppat variant with same amino acid change and at least reviewd by expert panel is found in ClinVar."
      }}else{
        PS1.message <- "PS1 is not assigned, it is not a missense variant or it is not used for this gene."
      }
  }else{
    final.criteria$criteria.res["PS1", 2]<- NA
    PS1.message <- "PS1 does not apply for this variant."
  }
  final.criteria[["PS1.message"]] <- PS1.message
  return (final.criteria)

}
clinVarEvidence <- function(all.information){
  score.total <- c()
  clinvar.evidence.calc <- sapply(all.information$clinVar$clinVar.info$same_codon,  function (x){
    score <- x[[8]]
    score.total <-rbind(score.total, as.matrix(score))
    return(score.total)
  })
  return(clinvar.evidence.calc)
}

################################################################################
## PM1
################################################################################
PM1 <- function (all.information, final.criteria){
  PM1.message <- NA
  final.criteria$criteria.res["PM1", 1:2]<- NA
  hot.spots <- hotSpList()

  if ( 1 %in% final.criteria$criteria.res["PVS1",]){
    final.criteria$criteria.res["PM1",c(3,4)] <- NA
    PM1.message <- "As PVS1 is assigned, PM1 is not calculated"
  }else if (all.information$Variant.Info$gene == "TP53" && all.information$Variant.Info$most.severe.consequence %in% c("missense_variant")&& nrow(all.information$cancer.hotspots$variant)>0 && as.numeric(all.information$cancer.hotspots$variant$count_alt) >= 10 ){
    if (1 %in% final.criteria$criteria.res["PS1",]){
      final.criteria$criteria.res["PM1",c(3,4)] <- c(0,1)
      PM1.message <- paste("PM1_supporting is assigned because there are ", all.information$cancer.hotspots$variant$count_alt, "somatic occurrences for the same amino acid reported in cancerhotspots.org, which is > or = than 10. It has been downgraded to supporting because PS1 is assigned.")
    }else{
      final.criteria$criteria.res["PM1",c(3,4)] <- c(1,0)
      PM1.message <- paste("PM1 is assigned because there are ", all.information$cancer.hotspots$variant$count_alt, "somatic occurrences for the same amino acid reported in cancerhotspots.org, which is > or = than 10.")
    }
  }else if(1 %in% final.criteria$criteria.res["PS1",]){
    final.criteria$criteria.res["PM1",c(3,4)] <- NA
    PM1.message <- "As PS1 is assigned by functional domain, PM1 is not calculated"
  }else{
    pos.aa <- stringr::str_extract(all.information$Variant.Info$prot, "[0-9]+") #we get the position of the aa

    if (!is.null(hot.spots[[all.information$Variant.Info$gene]])&
        (all.information$Variant.Info$most.severe.consequence %in% c("missense_variant"))){
      suma.mod <- hotSp("moderate", pos.aa, hot.spots, all.information) #we use the hot_sp function
      suma.sup <- hotSp("supporting", pos.aa, hot.spots, all.information) #we use the hot_sp function
      final.criteria$criteria.res["PM1",c(3,4)] <- c(suma.mod, suma.sup)
      PM1.message <- ifelse ((!is.na(suma.mod)&suma.mod ==1),
                             paste("PM1 is assigned because the variant is reported in a relevant functional domain (codon:", pos.aa,")."),
                             ifelse(!is.na(suma.sup) & suma.sup==1,
                                    paste("PM1_supporting is assigned because the variant is reported in a relevant functional domain (codon:", pos.aa,")."),
                                    "The variant is not located in a relevant functional domain."))
      PM1.message[all.information$Variant.Info$gene=="CHEK2"&((!is.na(suma.mod)&suma.mod ==1)|!is.na(suma.sup) & suma.sup==1)]<- paste(PM1.message , "Check if the aminoacid is highly conserved.")
    } else{ #meaning there is not
      final.criteria$criteria.res["PM1",c(3,4)] <- NA
      PM1.message <- "PM1 is not used for this variant or gene."
    }

  }
  final.criteria[["PM1.message"]] <- PM1.message
  return (final.criteria)
}
hotSp <- function(evidence, pos.aa, hot.spots, all.information){
  if (length(pos.aa)==1){
    hot.sp <- unlist(hot.spots[[all.information$Variant.Info$gene]][evidence]) == pos.aa
  }else {
    hot.sp <- unlist(hot.spots[[all.information$Variant.Info$gene]][evidence]) >= pos.aa[1] & unlist(hot.spots[[all.information$Variant.Info$gene]][evidence]) <= pos.aa[2]
  }
  sum <- sum(hot.sp==TRUE)
  return(sum)
}
hotSpList <- function(){
  hot.spots <- list(
    PTEN = list(moderate=c(90:94, 123:130, 166:168),
                supporting=c()),
    CHEK2 = list(moderate=c(116:118, 121,140,141,143:146, 148, 159, 162, 164, 166:171,173, 226,227,229,232,234,235,237,242,246,247,249,250,253, 256, 273, 275:278,282:284, 287,296, 298:302, 304:310, 319,321, 325:332,334,335,337:339,341:357,360,362,365:370,372,373,377,381,383:388,390,392:396,400,401,403,404,408,409,411:415,417,418,420,421,423:427,429,439,440,443,452,455,456,459,461:463,465:471,474,480,481,483:486),
                 supporting=c()), #revisar perque ha de ser highly conserverd
    TP53 = list(moderate=c(175,245, 248,249,273,282),
                supporting=c()), #tambe cal afegir variants with >OR=10 somatic observations in cancerspots.org (v)
    BAP1 = list(moderate=c(85,91,169,184),
                supporting=c()),
    BMPR1A = list (moderate=c(124,125,130),
                   supporting=c(59:123,126:132, 235:520)), #Revisar perque a més ha de ser considerat un highly conserved aa
    POLE = list(moderate=c(275,277,286,297,411,436,456,458,461),
                supporting=c( 276, 277, 278, 279, 280, 281, 284, 285, 286, 291, 292, 294, 363, 367, 368, 422, 423, 424, 425, 439, 440, 441, 442, 444, 445, 458, 462)),
    POLD1 = list(moderate=c(316,318),
                 supporting=c(317, 319, 320, 321, 325, 326, 327, 335, 396, 397, 400, 401, 402, 460, 471, 472, 473, 474, 475, 488, 489, 490, 491, 494, 515)))
  return(hot.spots)
}

################################################################################
## PM2
################################################################################
PM2 <- function (all.information, final.criteria){
  final.criteria$criteria.res["PM2", 1:2] <- NA
  PM2.message <- NULL
  gene <- all.information$Variant.Info$gene
  if (all.information$gnomAD$coverage$exomes >= 20){
    PM2.cutoff <- all.information$gene.specific.info$PM2
    freq.PM2 <- ifelse(all.information$gnomAD$coverage$exomes>=20 & all.information$gnomAD$coverage$genomes>=20,
                       all.information$gnomAD$info$exomes.genomes$non.cancer$overall$AF,
                       as.data.frame(all.information$gnomAD$info$exomes$non.cancer$overall)$AF)
    freq.PM2 <- ifelse(is.nan(freq.PM2),
                       0,
                       freq.PM2)
    gene <- all.information$Variant.Info$gene
    if( !is.na(PM2.cutoff)){
      sel.PM2 <- PM2.cutoff > freq.PM2
      subpopu.to.look <- c("non_cancer_nfe", "non_cancer_amr", "non_cancer_afr", "non_cancer_sas", "non_cancer_eas")

      if(all.information$gnomAD$coverage$exomes >= 20 & all.information$gnomAD$coverage$genomes >=20 ){
        selected <- all.information$gnomAD$info$exomes.genomes$non.cancer$subpopulations  %>%
          dplyr::filter(.data$rowname %in% subpopu.to.look,
                        .data$AC > 1,
                        .data$AF >= 0.00002)
      }else{
        selected <- all.information$gnomAD$info$exomes$non.cancer$subpopulations %>%
          tibble::rownames_to_column %>%
          dplyr::filter(.data$rowname %in% subpopu.to.look,
                        .data$AC > 1,
                        .data$AF >= 0.00002)
      }
      final.criteria$criteria.res["PM2",3] <- ifelse( sel.PM2==TRUE & nrow(selected) ==0,
                                                      1,
                                                      0)
      PM2.message[sel.PM2==TRUE & nrow(selected) ==0] <- paste ("PM2 is assigned because the variant is observed in ", round(freq.PM2, 6), " frequency which is < ", PM2.cutoff, " (PM2 cutoff) and no population with >=2 individuals has more than 1/50000 alleles.")
      PM2.message[sel.PM2==TRUE & nrow(selected) > 1] <- paste ("PM2 supporting is assigned because although", round(freq.PM2, 6), " frequency is < ", PM2.cutoff, "the subpopulation", paste(selected$rowname, collapse=","),"has/have", paste(selected$AC, collapse=","), "counts and an AF ", paste(selected$AF, collapse=","), "which is > than 0.00002.")


    }else{
      final.criteria$criteria.res["PM2",3] <- 0
    }


    if (final.criteria$criteria.res["PM2",3]==0 ){
      PM2.sup.cutoff <- all.information$gene.specific.info$PM2_sup
      sel.PM2.sup <- PM2.sup.cutoff >=  freq.PM2
      sel.PM2.sup[is.na(sel.PM2.sup)] <- FALSE
      if (gene %in% c("CDH1")){
        subpopu.to.look <- c("non_cancer_nfe", "non_cancer_amr", "non_cancer_afr", "non_cancer_sas", "non_cancer_eas")
        if (all.information$gnomAD$coverage$exomes >= 20 & all.information$gnomAD$coverage$genomes >= 20){
          selected.sup <- all.information$gnomAD$info$exomes.genomes$non.cancer$subpopulations  %>%
            dplyr::filter(.data$rowname %in% subpopu.to.look,
                          .data$AC > 1,
                          .data$AF >= 0.00002)

        }else{
          selected.sup <- all.information$gnomAD$info$exomes$non.cancer$subpopulations %>%
            tibble::rownames_to_column %>%
            dplyr::filter(.data$rowname %in% subpopu.to.look,
                          .data$AC > 1,
                          .data$AF >= 0.00002)

        }
        final.criteria$criteria.res["PM2",4] <- ifelse( !is.na(PM2.sup.cutoff)&&sel.PM2.sup==TRUE & nrow(selected.sup) ==0,
                                                        1,
                                                        0)
        PM2.message <- ifelse(sel.PM2.sup==TRUE& nrow(selected.sup) ==0,
                              paste("PM2_supporting is assigned because the variant is observed in", round(freq.PM2, 6), "frequency which is <", PM2.sup.cutoff, "(PM2_sup cutoff)."),
                              ifelse(sel.PM2.sup==TRUE& nrow(selected.sup) > 0,
                                     paste("PM2 supporting is denied because although", round(freq.PM2, 6), " frequency is < ", PM2.sup.cutoff, "the subpopulation", paste(selected.sup$rowname, collapse=","),"has/have", paste(selected.sup$AC, collapse=","), "counts and an AF ", paste(selected.sup$AF, collapse=","), "which is > than 0.00002."),
                                     paste( "PM2_supporting is denied because the variant is observed in", round(freq.PM2, 6), "frequency which is >", PM2.sup.cutoff, "(PM2 cutoff)." ))
        )
      }else{
        final.criteria$criteria.res["PM2",4] <- ifelse (!is.na(PM2.sup.cutoff) & sel.PM2.sup == TRUE,
                                                        1,
                                                        ifelse(!is.na(PM2.sup.cutoff),
                                                               0,
                                                               NA))
        if (is.null(PM2.message)||is.na(PM2.message)){
          PM2.message <- ifelse(sel.PM2.sup ==TRUE,
                                ifelse(gene=="TP53", "PM2_supporting is assigned because the variant is absent from gnomAD non-cancer",
                                       paste("PM2_supporting is assigned because the variant is observed in", round(freq.PM2, 6), "frequency which is <", PM2.sup.cutoff, "(PM2_sup cutoff).")),
                                paste( "PM2_supporting is denied because the variant is observed in", round(freq.PM2, 6), "frequency which is >", PM2.sup.cutoff, "(PM2 cutoff)." ))
          PM2.message[is.na(PM2.sup.cutoff)] <- paste("PM2 is denied because  either the variant is observed in a frequency higher than PM2 cut-off or some population with population with >=2 individuals has more than 1/50000 alleles (PM2 cut-off=  ", PM2.cutoff, "frequency: ", round(freq.PM2, 6),".")
        }}
    } else{
      final.criteria$criteria.res["PM2",4] <- NA
    }



  } else{
    final.criteria$criteria.res["PM2", c(3,4)] <- NA
    PM2.message <- ifelse(all.information$gnomAD$coverage$exomes < 20 & all.information$gnomAD$coverage$genomes < 20,
                          "PM2 cannot be calculated because there is not enough coverage in exomes neither genomes in gnomAD non_cancer v2.1 dataset.",
                          "PM2 cannot be calculated because there is not enough coverage in exomes and although genomes has >=20x there is only information about 15,708 whole genomes which is not representative enough.")
  }
  final.criteria[["PM2.message"]] <- PM2.message
  return (final.criteria)
}

################################################################################
## PM4
################################################################################
PM4 <- function (all.information, final.criteria){
  PM4.message <- NA
  final.criteria$criteria.res["PM4", 1:2] <- NA
  hot.spots <- hotSpList()
  gene <- all.information$Variant.Info$gene

  if ( 1 %in% final.criteria$criteria.res["PVS1",]){
    #| "error" %in% names(spliceai)){
    final.criteria$criteria.res["PM4",4] <- NA
    PM4.message <- "As PVS1 is assigned, PM4 is not calculated"
  }else if(gene %in% c("ATM", "CDH1")){
    final.criteria$criteria.res["PM4", 3] <- ifelse(all.information$Variant.Info$most.severe.consequence=="stop_lost",
                                                    1,
                                                    0)
    PM4.message <- ifelse(all.information$Variant.Info$most.severe.consequence=="stop_lost",
                          "PM4 is assigned because it is an stop-loss variant.",
                          "PM4 is not assigned because it is not an stop-loss variant.")
  }else if (gene %in% c("PTEN")){
    if(all.information$Variant.Info$most.severe.consequence=="stop_lost"){
      final.criteria$criteria.res["PM4", 3] <- 1
      PM4.message <- "PM4 is assigned because it is an stop-loss variant."
    }else if(all.information$Variant.Info$most.severe.consequence %in% c("frameshift_variant", "stop_gained") && final.criteria$NMD=="NO_effect"){
      final.criteria$criteria.res["PM4", 3] <- 1
      PM4.message <- "PM4 is assigned because the trucation is 3 prime to c.1121(NM_000314.6)."
    }}else if (!is.null(hot.spots[[all.information$Variant.Info$gene]])&
               (all.information$Variant.Info$most.severe.consequence %in% c("inframe_deletion"))){
      pos.aa.del <- stringr::str_extract_all(all.information$Variant.Info$prot, "[0-9]+") %>%
        unlist %>%
        as.numeric #we get the position of the aa
      suma.mod.del <- hotSp("moderate", pos.aa.del, hot.spots, all.information) #we use the hot_sp function
      suma.sup.del <- hotSp("supporting", pos.aa.del, hot.spots, all.information) #we use the hot_sp function
      final.criteria$criteria.res["PM4",3] <- ifelse(suma.mod.del>0, 1, 0)
      final.criteria$criteria.res["PM4",4] <- ifelse(suma.sup.del>0, 1, 0)

      PM4.message <- ifelse ((!is.na(suma.mod.del)&suma.mod.del > 0),
                             "PM4 is assigned because it is an inframe variant affecting a relevant functional domain.",
                             ifelse(!is.na(suma.sup.del) & suma.sup.del > 0,
                                    "PM4_supporting is assigned because it is an inframe variant affecting a relevant functional domain.",
                                    "The inframe variant is not located in a relevant functional domain."))

    }else{ #meaning there is not
      final.criteria$criteria.res["PM4",c(3,4)] <- NA
      PM4.message <- "PM4 is not used for this variant or gene."
    }

  final.criteria[["PM4.message"]] <- PM4.message
  return (final.criteria)
}

################################################################################
## PM5
################################################################################
PM5 <- function (all.information, final.criteria){
  PM5.message <- NULL
  server.mutalyzer <- "https://mutalyzer.nl/json/" #Mutalyzer's REST API
  gene <- all.information$Variant.Info$gene
  if(all.information$Variant.Info$most.severe.consequence=="missense_variant" & !(gene %in% c("PALB2", "CDH1", "ATM")) ){
    clinvar.evidence.calc <- clinVarEvidence(all.information)
    if(length(clinvar.evidence.calc)>0){
      rownames(clinvar.evidence.calc) <- c("ben_guideline", "pben_guideline", "unc_guideline",  "ppat_guideline", "pat_guideline",  "ben_expert",     "pben_expert",    "unc_expert",     "ppat_expert",    "pat_expert",     "ben_submitter",  "pben_submitter", "unc_submitter", "ppat_submitter", "pat_submitter")
    }
  }

  if( 1 %in% final.criteria$criteria.res["PS1",]){
    final.criteria$criteria.res["PM5",3] <- NA
    PM5.message <- "PM5 is not calculated because PS1 is assigned and co-existence is not allowed."
  }else if(1 %in% final.criteria$criteria.res["PM1",3] & !gene %in% c("TP53")){
    PM5.message <- "PM5 is not calculated because PM1 is assigned and co-existence is not allowed."
  }else if(all.information$Variant.Info$gene %in% c("PALB2") ){
    PM5.message <- "PM5 does not apply for this gene."
    final.criteria$criteria.res["PM5",3] <- NA
  }else if(all.information$Variant.Info$gene %in% c("CDH1")){
    if (all.information$Variant.Info$most.severe.consequence %in% c("frameshift_variant","stop_gained")){
      alt.splicing <- insilico (all.information$predictors$predictor.table, "Splicing Predictor", "Pathogenic")
      final.criteria$criteria.res["PM5",c(3,4)] <- ifelse(rep(nrow(alt.splicing)==0 && final.criteria$NMD=="YES", 2),
                                                          c(0,1),
                                                          c(0,0))
      PM5.message <- ifelse(nrow(alt.splicing)==0&&final.criteria$NMD=="YES",
                            paste("PM5_supporting is assigned because the variant is", all.information$Variant.Info$most.severe.consequence, "and is predicted to undergo NMD and SpliceAI suggests it has no impact."),
                            "PM5_supporting is not assigned because either there not exist NMD or SpliceAI suggests it has effect on splicing.")
    }else if (all.information$Variant.Info$most.severe.consequence %in% c("splice_donor_variant", "splice_acceptor_variant")){
      con <- DBI::dbConnect(RMySQL::MySQL(), user='userguest', dbname='class_variants', host='varhcdb001.cluster-ro-ca55bxrovxyt.eu-central-1.rds.amazonaws.com', password='jNU%cd%Xjw*tY*%')
      on.exit(DBI::dbDisconnect(con))
      pos <- stringr::str_extract(all.information$Variant.Info$variant, "[0-9]+")
      query <- paste0("SELECT PM5 FROM canonicals WHERE location='c.", pos ,"'")
      strength <- DBI::dbGetQuery(con, query) %>%
        as.character()

      final.criteria$criteria.res["PM5",c(3,4)] <- ifelse(rep(strength=="supporting",2),
                                                          c(0,1),
                                                          c(0,0))
      PM5.message <- ifelse(strength=="supporting",
                            "PM5_supporting is assigned according to site-specific recommendations for canonical splicing variants.",
                            "PM5_supporting is not assigned according to site-specific recommendations for canonical splicing variants.")

    }

  }else if (all.information$Variant.Info$gene %in% c("ATM")){
    if (all.information$Variant.Info$most.severe.consequence %in% c("frameshift_variant","stop_gained")){
      final.criteria$criteria.res["PM5",c(3,4)] <- ifelse(rep(final.criteria$NMD=="YES" | final.criteria$NMD=="PVS1",2),
                                                          c(0,1),
                                                          c(0,0))
      PM5.message <- ifelse(final.criteria$NMD=="YES"|final.criteria$NMD=="PVS1",
                            "PM5_supporting is assigned because the premature termination codon is upstream of p.Arg3047* which is the most C-terminal pathogenic variant, so it is expected to be more severe.",
                            "PM5_supporting is not assigned because the premature termination codon is downstream of p.Arg3047* which is the most C-terminal pathogenic variant, so it is expected to be less severe.")
    }
  }else if (all.information$Variant.Info$most.severe.consequence=="missense_variant"){
    if (length(clinvar.evidence.calc) > 0){
      data("BLOSUM62", package="Biostrings", envir = environment())
      prot.look <- stringr::str_sub(all.information$Variant.Info$protein, 4, -2)
      pm5.df <- clinvar.evidence.calc %>%
                                      t() %>%
                                      as.data.frame %>%
                                      dplyr::mutate(variants.name = colnames(clinvar.evidence.calc)) %>%
                                      dplyr::mutate (variant = unlist(purrr::map(stringr::str_split(.data$variants.name,":| \\("),2)),
                                                     protein = unlist(purrr::map(stringr::str_split(.data$variants.name,"\\(|\\)"),4))) %>%
                                      dplyr::filter ((.data$pat_guideline>0 | .data$ppat_guideline>0 | .data$pat_expert>0 | .data$ppat_expert>0)&!(stringr::str_detect(.data$protein, prot.look)))
      if(gene=="TP53"){
        pm5.df <- pm5.df %>%
                         dplyr::filter(.data$pat_guideline > 0 | .data$pat_expert > 0)
      }

      if (nrow(pm5.df) > 0){
        aa.our.variant <- aaShort(toProtein(all.information$Variant.Info$protein)$aa.alt)
        aa.ref <- aaShort(toProtein(all.information$Variant.Info$protein)$aa.ref)
        #server.spliceai <- "https://spliceailookup-api.broadinstitute.org/spliceai/?hg=37&distance=1000&precomputed=0&mask=1&variant="
        ensembl.id <- ensemblTranscript(all.information$Variant.Info$NM, gene)$id
        NC <- stringr::str_extract(all.information$Variant.Info$genomic, "NC_[0-9]+.[0-9]+")
        pm5.df <- pm5.df %>%
          dplyr::rowwise() %>%
          dplyr::mutate (a1 = toProtein(.data$protein)$aa.ref %>%
                           aaShort(), a2 = toProtein(.data$protein)$aa.alt %>%
                           aaShort()) %>%
          dplyr::mutate (blosum = BLOSUM62[.data$a1, .data$a2], prior = NA, grantham = calculateGrantham(.data$a1, .data$a2)) %>%
          dplyr::mutate(genomic = toGenomic(all.information$Variant.Info$NM, NC, stringr::str_sub(purrr::map(stringr::str_split(.data$variants.name, "\\(|\\):"),3),1,-2) ,all.information$Variant.Info$gene) ) %>%
          #dplyr::mutate(genomic=toGenomic(all.information$Variant.Info$NM, NC, stringr::str_sub(purrr::map(stringr::str_split(variants.name, "\\(|\\):"),3),1,-2) , all.information$Variant.Info$gene) ) %>%
          dplyr::mutate(ext.spliceai.pm5= paste0(all.information$Variant.Info$chr, "-", purrr::map(stringr::str_split(.data$genomic, "g\\.|[A-Z]"),4),"-", stringr::str_sub(.data$genomic, -3,-3), "-", stringr::str_sub(.data$genomic,-1,-1) )) %>%
          dplyr::mutate(scores.spliceai =  connectionDB( paste0("SELECT *  from spliceAI  WHERE var_chr= '", .data$ext.spliceai.pm5 ,"'AND max_dis= 1000 AND transcript='", all.information$Variant.Info$ensembl.id,"' AND masked='",TRUE,"';"))%>% purrr::map(.data, function(x) c(x[[8]], x[[10]],x[[12]],x[[14]]))) %>%
          dplyr::mutate(score.spliceai = max(unlist(.data$scores.spliceai), na.rm=FALSE))

        for (i in 1:nrow(pm5.df)){
          prior.pm5 <- priorUtahProb(object=NULL, gene=all.information$Variant.Info$gene, variant =pm5.df$variant[i])[[3]] %>%
            unlist
          pm5.df$prior[i] <- ifelse(length(prior.pm5)==0,
                                    NA,
                                    prior.pm5)
        }

        if (gene %in% c("MLH1", "MSH2", "MSH6", "PMS2")){
          if(all.information$predictors$predictor.table$values[11]>0.68){
            pm5.mmr.pat <- pm5.df %>%
              dplyr::filter(.data$pat_guideline>0 | .data$pat_expert>0)
            PM5.message[nrow(pm5.mmr.pat)>0] <-  paste("PM5 is assigned because because the prior value of the variant achieves PP3 and there are other variants at the same codon : ", pm5.mmr.pat$variants.name, "in ClinVar classified as pathogenic by at least expert panel (with a prior value", round(as.numeric(pm5.df$prior),6), " )", collapse=".")
            final.criteria$criteria.res["PM5","moderate"][nrow(pm5.mmr.pat)>0] <- 1

            if (nrow(pm5.mmr.pat)==0){
              pm5.mmr.ppat <- pm5.df %>%
                dplyr::filter(.data$ppat_guideline>0 | .data$ppat_expert>0)
              final.criteria$criteria.res["PM5", "supporting"] <- ifelse(nrow(pm5.mmr.ppat) > 0,
                                                                         1,
                                                                         0)
              PM5.message <- ifelse(nrow(pm5.mmr.ppat)>0,
                                    paste("PM5_supporting is assigned because the variant ", pm5.mmr.ppat$variants.name, " is classified in ClinVar as pPAT by expert panel. Be careful because the variant cannot exceed Class 4 (if this criterion is necessary in order to achieve class 5)."),
                                    paste("PM5_supporting is not assigned because there are not other variants in ClinVar classified as pPAT."))

            }
          }else{
            PM5.message<-  paste("PM5 is denied  because the prior value of the variant i below 0.68 (value=",all.information$predictors$predictor.table$values[11],")." )
            final.criteria$criteria.res["PM5","moderate"] <- 0
          }


        }else if(gene %in% c("TP53")){
          grantham.variant <- calculateGrantham(aa.ref, aa.our.variant)
          pm5.df.gratham <- pm5.df %>%
            dplyr::filter (.data$grantham <= grantham.variant && .data$score.spliceai<=0.2)
          final.criteria$criteria.res["PM5", c("moderate","supporting")] <- ifelse(rep(nrow(pm5.df.gratham)>1,2),
                                                                                   c(1,0),
                                                                                   ifelse(rep(nrow(pm5.df.gratham)==1,2),
                                                                                          c(0,1),
                                                                                          c(0,0)))
          PM5.message <- ifelse(nrow(pm5.df.gratham) > 1,
                                "PM5 is assigned because there are 2 or more PAT variants reported by the ClinGen TP53 VCE and the variant of interest has equal or higher Grantham score than the clinVar variants and no splicing alteration is observed.",
                                ifelse(nrow(pm5.df.gratham)==1,
                                       paste0("PM5 is assigned because there is 1 PAT variant (", pm5.df$variants.name, ") reported by the ClinGen TP53 VCE and the variant of interest has equal or higher Grantham score (", grantham.variant, ") than the ClinGen TP53 VCE variant (", pm5.df$grantham, ") and no splicing alteration is observed."),
                                       paste0("PM5 is not assigned because although the variant (", pm5.df$variants.name, ") is classified by ClinGen TP53 VCE, it does not meet one of the following items:the variant of interest has equal or higher Grantham score than the clinVar variants and no splicing alteration is observed ")))
          if(nrow(pm5.df.gratham)>=1 & 1 %in% final.criteria$criteria.res["PM1",3]){
            final.criteria$criteria.res["PM1",3] <- 0
            final.criteria$criteria.res["PM1",4] <- 1
            final.criteria$PM1.message <- paste(final.criteria$PM1.message, "However, it has been downgraded to supporting because PM5 is assigned.")
          }

        }else {

          blosum62.variant <- BLOSUM62[aa.ref, aa.our.variant]
          pm5.df.blosum <- pm5.df %>%
            dplyr::filter (.data$blosum >= blosum62.variant)
          final.criteria$criteria.res["PM5", "moderate"] <- ifelse(nrow(pm5.df.blosum)>0,
                                                                   1,
                                                                   0)
          PM5.message <- ifelse(nrow(pm5.df.blosum)>0,
                                paste("PM5 is assigned because the variant ", pm5.df.blosum$variants.name, "is classified as pPAT or PAT in ClinVar by Expert Panel and its BLOSUM is" , pm5.df.blosum$blosum, "and the BLOSUM of the queried variant is", blosum62.variant, "." ),
                                paste("PM5 is not assigned because although there are variants (", pm5.df$variants.name, ")  classified as pPAT or PAT in  ClinVar their BLOSUM is " , pm5.df$blosum, "and the BLOSUM of the queried variant is", blosum62.variant, "." ))
        }
      }else {
        PM5.message <- ifelse(gene=="TP53",
                              "PM5 is not assigned because there are not any pathogenic variants in ClinVar reviewed by expert panel.",
                              "PM5 is not assigned because there are not any pat/ppat variants in ClinVar reviewed by expert panel.")
        final.criteria$criteria.res["PM5",3] <- "NC"
      }

    }
  }else {
    PM5.message <- "PM5 does not apply for this type of variant."
    final.criteria$criteria.res["PM5",3] <- NA
  }

  final.criteria$criteria.res["PM5", 4][all.information$Variant.Info$gene=="BRCA1"&&
                                          (all.information$Variant.Info$most.severe.consequence %in% c("frameshift_variant","stop_gained"))&&
                                          !(all.information$codon.stop$variant.exon$exon %in% c(1,7,8,9))] <- 1
  if (all.information$Variant.Info$most.severe.consequence %in% c("frameshift_variant","stop_gained")){
    PM5.message[all.information$Variant.Info$gene=="BRCA1"&&
                  (all.information$Variant.Info$most.severe.consequence %in% c("frameshift_variant","stop_gained"))&&
                  !(all.information$codon.stop$variant.exon$exon %in% c(1,7,8,9))] <- paste("PM5 is assigned because the variant is a ", all.information$Variant.Info$most.severe.consequence, "variant located at exon", all.information$codon.stop$variant.exon$exon)
    final.criteria$criteria.res["PM5", 4][all.information$Variant.Info$gene=="BRCA2"&&
                                            (all.information$Variant.Info$most.severe.consequence %in% c("frameshift_variant", "stop_gained")) &&
                                            !(all.information$codon.stop$variant.exon$exon %in% c(1,4,6,12,21,24,26))] <- 1
    PM5.message[all.information$Variant.Info$gene=="BRCA2"&&
                  (all.information$Variant.Info$most.severe.consequence %in% c("frameshift_variant", "stop_gained")) &&
                  !(all.information$codon.stop$variant.exon$exon %in% c(1,4,6,12,21,24,26))] <- paste("PM5 is assigned because the variant is a ", all.information$Variant.Info$most.severe.consequence, "variant located at exon", all.information$codon.stop$variant.exon$exon)

  }
  if(is.null(PM5.message))PM5.message <- NA
  final.criteria[["PM5.message"]] <- PM5.message
  return (final.criteria)

}

################################################################################
## PP2
################################################################################
PP2 <- function(all.information, final.criteria){
  final.criteria$criteria.res["PP2", 1:3] <- NA
  if(all.information$Variant.Info$gene %in% c("PTEN")&&
     (all.information$Variant.Info$most.severe.consequence %in% c("missense_variant") |
      all.information$Variant.Info$most.severe.consequence.1 %in% c("missense_variant") |
      all.information$Variant.Info$most.severe.consequence.1 %in% c("missense_variant"))){
    final.criteria$criteria.res["PP2",4] <- 1
    PP2.message <- "PP2 is assigned because it is a missense variant in a gene that has a low rate of benign missense variation and where missense variants are a common mechanism of disease."
  }else{
    if(all.information$Variant.Info$most.severe.consequence !="missense_variant" | all.information$Variant.Info$gene %in% c("ATM", "CDH1", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "TP53" )){
      PP2.message <- "PP2 does not apply for this variant or gene."
      final.criteria$criteria.res["PP2", 4] <- NA
    }else{
      PP2.message <- "PP2 is not calculated. You should check if this gene has a low rate of benign missense variation and if missense variants are a common mechanism of disease to assign or deny PP2."
    }

  }
  final.criteria[["PP2.message"]] <- PP2.message

  return(final.criteria)
}

################################################################################
## PP3
################################################################################
PP3 <- function (all.information, final.criteria){
  PP3.message <- NULL
  final.criteria$criteria.res["PP3", 1:3] <- NA
  gene <- all.information$Variant.Info$gene
  variant.exon <- all.information$codon.stop$exons %>%
    dplyr::filter( .data$V1 <= all.information$Variant.Info$start,
                   .data$V2 >= all.information$Variant.Info$end)
  exceptions <- 1 %in% c(final.criteria$criteria.res["PVS1",],
                         final.criteria$criteria.res["PS1",],
                         final.criteria$criteria.res["PM4",])
  #Only allowed PM1 and PP3 when PM1 is for mutational hot-spots as is the case of TP53
  exceptions2 <- (1 %in% final.criteria$criteria.res["PM1",3]) & !(gene %in% c("TP53"))
  if (exceptions == TRUE | exceptions2 == TRUE){
    final.criteria$criteria.res["PP3",4] <- NA
    if (1 %in% final.criteria$criteria.res["PVS1",]){
      PP3.message <- "PP3 is not calculated because PVS1 is met and co-usage is not permitted"
    }else if (1 %in% final.criteria$criteria.res["PM4",]){
      PP3.message <- "PP3 is not calculated because PM4 is met and co-usage is not permitted"
    }else if (1 %in% final.criteria$criteria.res["PS1",]){
      PP3.message <- "PP3 is not calculated because PS1 is met and co-usage is not permitted"
    }else if (1 %in% final.criteria$criteria.res["PM1",3]){
      PP3.message <- "PP3 is not calculated because PM1_moderate is met and co-usage is not permitted (it can only be combined if PM1 is supporting)."
    }
  }else{
    if (all.information$Variant.Info$gene!="PTEN"|
        (all.information$Variant.Info$gene=="PTEN" && all.information$Variant.Info$most.severe.consequence %in% c("synonymous_variant", "intron_variant", "splice_donor_variant", "splice_acceptor_variant", "splice_donor_region_variant", "splice_acceptor_region_variant") & !stringr::str_detect(all.information$Variant.Info$variant, "c.79+1|c.79+2|c.80-2|c.80-1"))){
      # First SpliceAI is checked
      alt.splicing <- insilico(all.information$predictors$predictor.table, "Splicing Predictor", "Pathogenic")
      alt.splicing.prior <- insilico(all.information$predictors$predictor.table, "Splicing Predictor", "Increased")
      alt.splicing.prior2 <- insilico(all.information$predictors$predictor.table, "Splicing Predictor", "High")
      if(nrow(alt.splicing) >0 |nrow(alt.splicing.prior) >0 | nrow(alt.splicing.prior2) > 0){
        final.criteria$criteria.res["PP3", 4] <- 1
        predictor.altered <- ifelse(nrow(alt.splicing.prior)>0,
                                    paste0(rownames(alt.splicing.prior), ": ", alt.splicing.prior$values),
                                    ifelse(nrow(alt.splicing.prior2)>0,
                                           paste0(rownames(alt.splicing.prior2), ": ", alt.splicing.prior2$values),
                                           ifelse(nrow(alt.splicing)>0,
                                                  paste0(rownames(alt.splicing), ": ", alt.splicing$values))))
        PP3.message <- paste0(predictor.altered,", since splicing effect alteration is predicted PP3 is assigned.")

      }else if (all.information$Variant.Info$most.severe.consequence %in% c("inframe_deletion", "inframe_insertion")){

        alt.provean <- insilico(all.information$predictors$predictor.table, "Protein effect", "Pathogenic")
        final.criteria$criteria.res["PP3", 4] <- ifelse(nrow(alt.provean)>0,
                                                        1,
                                                        0)
        PP3.message <-  ifelse(is.na(all.information$predictors$predictor.table["Provean", "values"]),
                               "Please consider installing or checking Provean score manually to assign or deny PP3.",
                               ifelse(nrow(alt.provean)>0,
                               paste("PP3 is assigned because it is an", all.information$Variant.Info$most.severe.consequence, "and Provean score is", alt.provean$values, "that is =< than -2.5."),
                               paste("PP3 is denied because SpliceAi predicts no impact and it is an", all.information$Variant.Info$most.severe.consequence, "and Provean score is", all.information$predictors$predictor.table["Provean", "values"], "that is > than -2.5.")))

        }else if (nrow(alt.splicing)==0 & nrow(alt.splicing.prior)==0 & nrow(alt.splicing.prior2)==0 & all.information$Variant.Info$most.severe.consequence =="missense_variant" & !(all.information$Variant.Info$gene %in% c("CDH1", "BAP1", "PALB2"))){
        PP3.sentence <- NULL
        predictors.prot <- insilico(all.information$predictors$predictor.table, "Protein effect", "Pathogenic")
        final.criteria$criteria.res["PP3", 4][(nrow(predictors.prot)==1 & !(gene %in% c("TP53"))) |
                                                (nrow(predictors.prot)==2 &gene %in% c("TP53"))] <- 1
        PP3.sentence[(nrow(predictors.prot)==1 & !(gene %in% c("TP53"))) |
                       (nrow(predictors.prot)==2 &gene %in% c("TP53"))] <- paste0(rownames(predictors.prot), ": ", round(as.numeric(predictors.prot$values),4), " which is ", predictors.prot$operator.pathogenic, " than ", predictors.prot$TD.cut.off, collapse=" and ")


        #Exception : MMR genes -> PP3_Moderate:Missense with MAPP+PolyPhen-2 prior probability for pathogenicity >0.81
        final.criteria$criteria.res["PP3", c(3,4)][gene %in% c("MLH1", "MSH2", "MSH6", "PMS2") && !is.na(all.information$predictors$predictor.table$values[11]) & as.numeric(all.information$predictors$predictor.table$values[11])>0.81] <- c(1,0)
        final.criteria$criteria.res["PP3", 3][gene %in% c("MLH1", "MSH2", "MSH6", "PMS2") && !is.na(all.information$predictors$predictor.table$values[11]) & as.numeric(all.information$predictors$predictor.table$values[11])<=0.81] <- 0
        PP3.sentence[gene  %in% c("MLH1", "MSH2", "MSH6", "PMS2") && !is.na(all.information$predictors$predictor.table$values[11]) & as.numeric(all.information$predictors$predictor.table$values[11])>0.81] <- stringr::str_replace(PP3.sentence, "which is > than 0.68", "which is > than 0.81")

        #Exception: TP53 -> PP3_moderate: aGVGD Zebrafish Class C65 required and BayesDel score 0.16
        final.criteria$criteria.res["PP3", c(3,4)][gene  == "TP53"&&nrow(predictors.prot)==2 && predictors.prot["aGVGD_zebrafish", "values"] >= 65] <-  c(1,0)
        final.criteria$criteria.res["PP3", c(3)][gene  == "TP53"&&nrow(predictors.prot)==2 && predictors.prot["aGVGD_zebrafish", "values"] < 65] <-  0
        PP3.sentence[gene  == "TP53" && nrow(predictors.prot)==2 && predictors.prot["aGVGD_zebrafish", "values"]>=65] <- stringr::str_replace(PP3.sentence, "which is >= than 25", "which is >= than 65")

        #Exception: BRCA1 exons 8 and 9
        final.criteria$criteria.res["PP3",1:4][(variant.exon$exon %in% c(8, 9)) & gene=="BRCA1"] <- rep(NA,4)

        #final.criteria$criteria.res["PP3",4][is.na(final.criteria$criteria.res["PP3",4])] <- 0

        final.criteria$criteria.res["PP3",4][final.criteria$criteria.res["PP3",4] == "NC"] <- 0
        PP3.message[1 %in%  final.criteria$criteria.res["PP3",4]] <- paste(PP3.sentence, ".Protein effect is predicted to be altered so PP3 is met.")
        PP3.message[1 %in%  final.criteria$criteria.res["PP3",3]] <- paste(PP3.sentence, ".Protein effect is predicted to be altered so PP3_moderate is met.")
        PP3.message[variant.exon$exon %in% c(8, 9) &gene=="BRCA1"] <- "PP3 is not taken into account for variants located at exons 8 or 9 in BRCA1, alternate transcript exists."
      }else{
        PP3.message <- "PP3 is not assigned because splicing predictors are not altered."
      }
    }else{
      PP3.message <- ifelse(stringr::str_detect(all.information$Variant.Info$variant, "c.79+1|c.79+2|c.80-2|c.80-1")& gene %in% c("PTEN"),
                            "PP3 is denied because PTEN guidelines specify that it should not to be applied for variants which may impact the intron 1 splice donor or acceptor sites.",
                            "PP3 is not used with this type of variant.")
    }
  }
  PP3.message <- ifelse(gene=="PTEN" & stringr::str_detect(all.information$Variant.Info$variant, "c.635-2|c.635-1"),
                        paste(PP3.message, "Variant impacts the intron 6 splice acceptor, this criterion should be used cautiously according to PTEN guidelines."),
                        PP3.message)
  final.criteria[["PP3.message"]] <- PP3.message

  return (final.criteria)
}

################################################################################
## PS3/BS3
################################################################################
PS3_BS3 <- function (all.information, final.criteria){
  PS3 <- NULL
  BS3 <- NULL
  PS3.message <- NULL
  BS3.message <- NULL

  if(all.information$Variant.Info$gene %in% c("MLH1", "MSH2", "MSH6", "PMS2")){
    if(all.information$Variant.Info$most.severe.consequence %in% c("intron_variant", "synonymous_variant") | all.information$Variant.Info$most.severe.consequence.1 %in% c("intron_variant", "synonymous_variant") | all.information$Variant.Info$most.severe.consequence.2 %in% c("intron_variant", "synonymous_variant")){
      if(any(!is.na(all.information$insight.info$classification)) && length(all.information$insight.info$classification)>0){
        for (y in 1:length(all.information$insight.info$classification)){
          names <- names(all.information$insight.info$classification[[y]])
          true.names <- names %in% ("Splicing/transcript expression")
          if(any(true.names)){
            splice.insight <- all.information$insight.info$classification[[y]]$`Splicing/transcript expression`
            splice.insight[is.na(splice.insight)] <- ""
            final.criteria$criteria.res["BS3", 2] <- ifelse(stringr::str_detect(splice.insight, "No effect on splicing") && stringr::str_detect(splice.insight, "NMD inhibitor"),
                                                            1,
                                                            0)
            BS3.message <- ifelse(stringr::str_detect(splice.insight, "No effect on splicing") && stringr::str_detect(splice.insight, "NMD inhibitor"),
                                  paste("BS3 is assigned because the variant is inSIGHT with the explanation:", splice.insight),
                                  "BS3 cannot be assigned because we do not have information of splicing in inSIGHT. You should check more bibliography.")
          }
        }

      }}
    if (!is.na(all.information$functional.assays$cimra) && nrow(all.information$functional.assays$cimra)>0){
      cimra <- all.information$functional.assays$cimra
      final.criteria$criteria.res["PS3", 2] <- ifelse(cimra$Odds_pat > 18.7, 1, 0)
      final.criteria$criteria.res["PS3", 3] <- ifelse(cimra$Odds_pat > 4.3 & cimra$Odds_pat <= 18.7, 1, 0)
      final.criteria$criteria.res["PS3", 4] <- ifelse(cimra$Odds_pat > 2.08 & cimra$Odds_pat <= 4.3 , 1, 0)
      PS3.message <- ifelse(cimra$Odds_pat > 18.7,
                            paste("PS3 is assigned because the variant is found in CIMRA with an Odds_pat of", cimra$Odds_pat, "which is > 18.7 ."),
                            ifelse (cimra$Odds_pat > 4.3 & cimra$Odds_pat <= 18.7,
                                    paste ("PS3_mod is assigned because the variant is found in CIMRA with an Odds_pat of", cimra$Odds_pat, "which is > 4.3 and <=18.7 ."),
                                    ifelse ((cimra$Odds_pat > 2.08 & cimra$Odds_pat <= 4.3),
                                            paste("PS3_supporting is assigned because the variant is found in CIMRA with an Odds_pat of", cimra$Odds_pat, "which is > 2.08 and <= 4.3 ."),
                                            paste("PS3 is not met because although the variant is found in CIMRA, the odds_pat is", cimra$Odds_pat, "which is <= 2.08 ."))))

      final.criteria$criteria.res["BS3", 2] <- ifelse(cimra$Odds_pat <= 0.05 , 1, 0)
      final.criteria$criteria.res["BS3", 4] <- ifelse(cimra$Odds_pat <= 0.48 & cimra$Odds_pat > 0.05 , 1, 0)
      BS3.message <- ifelse ( cimra$Odds_pat <= 0.05,
                              paste("BS3 is assigned because the variant is found in CIMRA with an Odds_pat of", cimra$Odds_pat, "which is <= 0.05."),
                              ifelse(cimra<= 0.48 & cimra$Odds_pat > 0.05,
                                     paste ("BS3_supporting is assigned because the variant is found in CIMRA with an Odds_pat of", cimra$Odds_pat, "whivh is <=0.48 and > 0.05 ."),
                                     paste ("BS3 is not met because although the variant is found in CIMRA, the odds_pat is of", cimra$Odds_pat, "which is > 0.48 .")))
    }
  } else if(all.information$Variant.Info$gene %in% c("TP53")){
    functionals.tp53 <- all.information$functional.assays$tp53.functionals
    final.criteria$criteria.res["PS3", 2] <- ifelse(functionals.tp53$Kato == "non-functional" &
                                                      (functionals.tp53$Giacomelli=="DNE_LOF"|
                                                         as.numeric(functionals.tp53$Kotler) >= -1),
                                                    1,
                                                    0)


    final.criteria$criteria.res["PS3", 3] <- ifelse(final.criteria$criteria.res["PS3", 2]!=1 & ((functionals.tp53$Kato == "partially functional" &
                                                                                                   (functionals.tp53$Giacomelli=="DNE_LOF"| as.numeric(functionals.tp53$Kotler)>= -1))|
                                                                                                  (functionals.tp53$Kato=="NA" & functionals.tp53$Giacomelli=="DNE_LOF" & as.numeric(functionals.tp53$Kotler)>= -1)),
                                                    1,
                                                    0)
    PS3.message <- ifelse(final.criteria$criteria.res["PS3",2]==1,
                          paste("PS3 is assigned because the variant is classified in Kato et al.(2003) as non-functional", ifelse(functionals.tp53$Giacomelli=="DNE_LOF",
                                                                                                                                   "and by Giacomelli et al.(2018) as DNE+LOF.",
                                                                                                                                   "and by Kotler et al.,(2018) with a RFS >= -1.")),
                          ifelse(final.criteria$criteria.res["PS3",3]==1,
                                 paste("PS3_moderate is assigned because",
                                       ifelse(functionals.tp53$Kato == "partially functional",
                                              paste("the variant is classified in Kato et al.(2003) as partially functional and", ifelse(functionals.tp53$Giacomelli=="DNE_LOF",
                                                                                                                                         "by Giacomelli et al.(2018) as DNE+LOF",
                                                                                                                                         "by Kotler et al.(2018) with a RFS >= -1 ")),
                                              "the variant is not found in Kato et al.,(2003) but it is classified by Giacomelli et al.(2018) as DNE+LOF and by Kotler et al.(2018) with a RFS >= -1")),
                                 "PS3 is not assigned."))

    final.criteria$criteria.res["BS3", 2] <- ifelse((functionals.tp53$Kato == "supertrans"| functionals.tp53$Kato == "functional") &
                                                      (functionals.tp53$Giacomelli=="notDNE_notLOF"| as.numeric(functionals.tp53$Kotler) < -1),
                                                    1,
                                                    0)
    final.criteria$criteria.res["BS3", 4] <- ifelse(final.criteria$criteria.res["BS3", 2]!=1 & ((functionals.tp53$Kato == "partially functional" &
                                                                                                   (functionals.tp53$Giacomelli=="notDNE_notLOF"| as.numeric(functionals.tp53$Kotler)< -1))|
                                                                                                  (functionals.tp53$Kato=="NA" & functionals.tp53$Giacomelli=="notDNE_notLOF" & as.numeric(functionals.tp53$Kotler)< -1)),
                                                    1,
                                                    0)

    BS3.message <- ifelse(final.criteria$criteria.res["BS3",2]==1,
                          paste("BS3 is assigned because the variant is classified in Kato et al.(2003) as", functionals.tp53$kato, ifelse(functionals.tp53$Giacomelli == "notDNE_notLOF",
                                                                                                                                           "and by Giacomelli et al.(2018) as notDNE+notLOF.",
                                                                                                                                           "and by Kotler et al.,(2018) with a RFS < -1.")),
                          ifelse(final.criteria$criteria.res["BS3",4]==1,
                                 paste("BS3_supporting is assigned because",
                                       ifelse(functionals.tp53$Kato == "partially functional",
                                              paste("the variant is classified in Kato et al.(2003) as partially functional and", ifelse(functionals.tp53$Giacomelli=="notDNE_notLOF",
                                                                                                                                         "by Giacomelli et al.(2018) as DNE+LOF",
                                                                                                                                         "by Kotler et al.(2018) with a RFS < -1 ")),
                                              "the variant is not found in Kato et al.(2003) but it is classified by Giacomelli et al.(2018) as DNE+LOF and by Kotler et al.(2018) with a RFS < -1.")),
                                 "BS3 is not assigned."))

  } else if (all.information$Variant.Info$gene %in% c("PTEN")){
    final.criteria$criteria.res["PS3", 2] <- ifelse(!is.na(all.information$functional.assays$rna.pten) &&
                                                      nrow(all.information$functional.assays$rna.pten)>0 &&
                                                      all.information$functional.assays$rna.pten$criteria_given=="PS3",
                                                    1,
                                                    0)
    PS3.message <- ifelse(!is.na(all.information$functional.assays$rna.pten) &&
                            nrow(all.information$functional.assays$rna.pten)>0 &&
                            all.information$functional.assays$rna.pten$criteria_given=="PS3",
                          paste("PS3 is assigned because there is", all.information$functional.assays$rna.pten$RNA_change,  "(proven splicing impact functional) demonstrated by", all.information$functional.assays$rna.pten$functional_study_source, "."),
                          "We have not found information about functional studies of this variant in our database.")

  }else if (all.information$Variant.Info$gene %in% c("ATM")){
    atm.functionals <- all.information$functional.assays$atm.functionals
    if (!is.na(atm.functionals) && nrow(atm.functionals)>0){
      final.criteria$criteria.res["PS3", c(2:4)] <- ifelse(rep(atm.functionals$final_criteria == "PS3_mod",3),
                                                           c(0,1,0),
                                                           ifelse(rep(atm.functionals$final_criteria =="PS3_sup",3),
                                                                  c(0,0,1),
                                                                  c(0,0,0)))
      final.criteria$criteria.res["BS3", c(2:4)] <- ifelse(rep(atm.functionals$final_criteria == "BS3_mod",3),
                                                           c(0,1,0),
                                                           ifelse(rep(atm.functionals$final_criteria =="BS3_sup",3),
                                                                  c(0,0,1),
                                                                  c(0,0,0)))
      PS3.message <- atm.functionals$reasoning_PS3
      BS3.message <- atm.functionals$reasoning_BS3
    }}else if(all.information$Variant.Info$gene %in% c("CHEK2")){
      chek2.functionals <- all.information$functional.assays$chek2.functionals
      if (!is.na(chek2.functionals) && nrow(chek2.functionals)>0){
        final.criteria$criteria.res["PS3", 2] <- ifelse(chek2.functionals$PS3_criteria=="PS3", 1, 0)
        PS3.message <- ifelse(chek2.functionals$PS3_criteria=="PS3",
                              "PS3 is assigned because there are some functional assays that prove the pathogenic effect",
                              paste0(chek2.functionals$PS3_criteria, "PS3 is denied."))
        BS3.message <- "BS3 is denied, not enough functional studies in favour of benignity."
      }}


  PS3.message[is.null(PS3.message)] <- "PS3 is not assigned because variant is not found in the automated functional studies. Please check if there is more literature."
  BS3.message[is.null(BS3.message)] <- "BS3 is not assigned because variant is not found in the automated functional studies. Please check if there is more literature"
  final.criteria[["PS3.message"]] <- PS3.message
  final.criteria[["BS3.message"]] <- BS3.message
  return(final.criteria)
}

################################################################################
## BA1
################################################################################
ba1_bs1_df <- function (dataset, all.information, criterion){
  if (criterion %in%c("BA1", "BS1")){
    cutoff <- all.information$gene.specific.info[criterion]
  }else if(criterion %in% c("BS1_sup")& all.information$Variant.Info$gene %in% c("PTEN")){
    cutoff <- 0.000043
  }else{
    stop("No correct value entered.")
  }
  num.AN <- all.information$gene.specific.info$alleles
  ic <- all.information$gene.specific.info$IC
  subpopu.to.look <- c("non_cancer_nfe", "non_cancer_amr", "non_cancer_afr", "non_cancer_sas", "non_cancer_eas", ".", "1")
  if(!is.na(ic)){
    selected <- dataset %>%
      dplyr::filter(.data$rowname %in% subpopu.to.look,
                    CI >= as.numeric(cutoff))
  }else{
    selected <- dataset %>%
      dplyr::filter(.data$rowname %in% subpopu.to.look,
                    .data$AF >=as.numeric(cutoff),
                    .data$AN >= num.AN, 
                    .data$AC>= 5)
  }
  return(selected)
}
BA1 <- function (all.information, final.criteria){
  BA1.message <- NULL
  final.criteria$criteria.res["BA1", c(2:4)] <- NA
  ic <- all.information$gene.specific.info$IC
  ba1.cutoff <- all.information$gene.specific.info$BA1
  if (all.information$gnomAD$coverage$exomes >= 20) {
    if (all.information$gnomAD$coverage$exomes >= 20&all.information$gnomAD$coverage$genomes >= 20){
      ba1.all <- ba1_bs1_df(all.information$gnomAD$info$exomes.genomes$non.cancer$overall, all.information, "BA1")
      ba1.df <- ba1_bs1_df(all.information$gnomAD$info$exomes.genomes$non.cancer$subpopulations, all.information, "BA1")
    }else{
      ba1.all <- ba1_bs1_df(all.information$gnomAD$info$exomes$non.cancer$overall %>%
                                                                                  tibble::as_tibble() %>%
                                                                                  tibble::rownames_to_column(), all.information, "BA1")
      ba1.df <- ba1_bs1_df(all.information$gnomAD$info$exomes$non.cancer$subpopulations %>%
                                                                                        tibble::as_tibble() %>%
                                                                                        tibble::rownames_to_column(), all.information, "BA1")
    }

    final.criteria$criteria.res["BA1", 1] <- ifelse (nrow(ba1.df)==0,  0,  1)
    BA1.message <- ifelse (nrow(ba1.df)==0,
                           paste0("The location is well covered in gnomAD (at least in exomes) but the variant is found in AF <  ", ba1.cutoff, " in european(non-finish), lationo, African, South Asian and East Asian so BA1 is not assigned."),
                           paste0("The location is well covered in gnomAD (at least in exomes) and BA1 is assigned because", row.names(ba1.df), "subpopulation is ", ifelse(rep(!is.na(ic), nrow(ba1.df)), round(ba1.df$CI, 6) , round(ba1.df$AF, 6)), "that is > than", ba1.cutoff, " .", collapse=""))

  } else{
    final.criteria$criteria.res["BA1", 1] <- NA
    BA1.message <- "There is not enough coverage so criteria based on gnomAD cannot be calculated."
  }
  max.subpop <- all.information$gnomAD$info$exomes.genomes$non.cancer$subpopulations %>%
                                                                                     dplyr::filter(CI==max(all.information$gnomAD$info$exomes.genomes$non.cancer$subpopulations$CI))
  popu <- max.subpop %>%
                      dplyr::select("rowname") %>%
                      dplyr::mutate(rowname = stringr::str_sub(.data$rowname, -3,-1)) %>%
                      as.character()
  sentence <- paste0("Freq gnomAD all non-cancer v2.1.1 = ", round(as.numeric(all.information$gnomAD$info$exomes.genomes$non.cancer$overall$AF)*100,5), "% (MCAF ",  as.numeric(all.information$gene.specific.info$IC)*100,
                     "% = ", round(as.numeric(all.information$gnomAD$info$exomes.genomes$non.cancer$overall$AF)*100,5), "%).",
                     ifelse(nrow(max.subpop)==8, "",
                            paste0(".Max freq in ",  popu, " subpopulation ", round(as.numeric(max.subpop$AF)*100,5), "% (MCAF ", as.numeric(all.information$gene.specific.info$IC)*100, "% = ", round(as.numeric(max.subpop$CI)*100,5), "%).", collapse=" and ") ))
  BA1.message <- paste(sentence, BA1.message)

  final.criteria[["BA1.message"]] <- BA1.message

  return (final.criteria)
}

################################################################################
## BS1
################################################################################
BS1 <- function (all.information, final.criteria){
  final.criteria$criteria.res["BS1", c(1, 3,4)] <- NA
  ic <- all.information$gene.specific.info$IC
  bs1.cutoff <- all.information$gene.specific.info$BS1
  if (!(1 %in% final.criteria$criteria.res["BA1",])){
    if (all.information$gnomAD$coverage$exomes >= 20){
      if (all.information$gnomAD$coverage$exomes >= 20 & all.information$gnomAD$coverage$genomes >= 20){
        bs1.all <- ba1_bs1_df(all.information$gnomAD$info$exomes.genomes$non.cancer$overall, all.information, "BS1")
        bs1.df <- ba1_bs1_df(all.information$gnomAD$info$exomes.genomes$non.cancer$subpopulations, all.information, "BS1")
      }else{
        bs1.all <- ba1_bs1_df(all.information$gnomAD$info$exomes$non.cancer$overall %>% tibble::as_tibble() %>% tibble::rownames_to_column(), all.information , "BS1")
        bs1.df <- ba1_bs1_df(all.information$gnomAD$info$exomes$non.cancer$subpopulations %>% tibble::as_tibble() %>% tibble::rownames_to_column(), all.information , "BS1")
      }

      final.criteria$criteria.res["BS1", 2] <- ifelse (nrow(bs1.df)==0 & nrow(bs1.all)==0,  0,  1)
      BS1.message <- ifelse (nrow(bs1.all)>0,
                             paste0("The location is well covered in gnomADv2.1 (at least in exomes) and the overall AF in non_cancer v2.1.1 dataset is", bs1.all$AF, "which is >", bs1.cutoff),
                             ifelse(nrow(bs1.df)==0,
                                    paste0("The location is well covered in gnomAD (at least in exomes) but the variant is found in AF <  ", bs1.cutoff, " in european(non-finish), lationo, African, South Asian and East Asian so BS1 is not assigned.", collapse=""),
                                    paste0("The location is well covered in gnomAD (at least in exomes) and BS1 is assigned because", row.names(bs1.df), "subpopulation is ", ifelse(rep(!is.na(ic), nrow(bs1.df)), round(bs1.df$CI, 6) , round(bs1.df$AF, 6)), "that is > than", bs1.cutoff, " .", collapse="")))
      if (all.information$Variant.Info$gene %in% c("PTEN") & nrow(bs1.df)==0){
        bs1.sup.cutoff <- 0.000043
        if(all.information$gnomAD$coverage$exomes >= 20 & all.information$gnomAD$coverage$genomes >= 20){
          bs1.sup.df <- ba1_bs1_df(all.information$gnomAD$info$exomes.genomes$non.cancer$subpopulations, all.information, "BS1_sup")
          bs1.sup.all <- ba1_bs1_df(all.information$gnomAD$info$exomes.genomes$non.cancer$overall, all.information, "BS1_sup")
        }else{
          bs1.sup.df <- ba1_bs1_df(all.information$gnomAD$info$exomes$non.cancer$subpopulations %>% tibble::rownames_to_column() %>% tibble::as_tibble() , all.information, "BS1_sup")
          bs1.sup.all <- ba1_bs1_df(all.information$gnomAD$info$exomes$non.cancer$overall %>% tibble::as_tibble() %>% tibble::rownames_to_column() , all.information, "BS1_sup")
        }
        final.criteria$criteria.res["BS1", 4] <- ifelse (nrow(bs1.sup.df)==0 & nrow(bs1.sup.all)==0,  0,  1)
        BS1.message <- ifelse (nrow(bs1.sup.all)>0,
                               paste0("The location is well covered in gnomAD and the overall AF in non_cancer v2.1.1 dataset is", bs1.all$AF, "which is >", bs1.cutoff),
                               ifelse(nrow(bs1.sup.df)==0,
                                      paste0("The location is well covered in gnomAD but the variant is found in AF <  ", bs1.sup.cutoff, " in european(non-finish), lationo, African, South Asian and East Asian so BS1_supporting is not assigned.", collapse=""),
                                      paste0("The location is well covered and BS1_supporting is assigned because ",bs1.sup.df$rowname, " subpopulation is ", ifelse(rep(!is.na(ic), nrow(bs1.sup.df)), round(bs1.sup.df$CI, 6) , round(bs1.sup.df$AF, 6)), " that is > than 0.000043 but < 0.001.", collapse="")))

      }
    }else{
      BS1.message <- "There is not enough coverage so criteria based on gnomAD cannot be calculated."
    }

  }else{
    BS1.message <- "Coexistence with BA1 criterion is not possible."
  }
  final.criteria[["BS1.message"]] <- BS1.message
  return (final.criteria)
}

################################################################################
## BS2
################################################################################
BS2 <- function (all.information, final.criteria){
  BS2.message <- NULL
  #gene.specific.information
  BS2.strong <- all.information$gene.specific$BS2 %>%
    as.numeric()
  BS2.strong.db <- all.information$gene.specific$BS2_db
  BS2.sup <- all.information$gene.specific$BS2_sup %>% as.numeric()
  BS2.sup.db <- all.information$gene.specific$BS2_sup_db
  zigosity <- all.information$gene.specific$status

  if(!(BS2.strong.db %in% c("FLOSSIES", "GNOMAD_non_cancer", "GNOMAD_non_neuro", NA, "")))  stop("'BS2_db' cannot be accessed")
  if(!(BS2.sup.db %in% c("FLOSSIES", "GNOMAD_non_cancer", "GNOMAD_non_neuro", NA, "")))  stop("'BS2_sup_db' cannot be accessed")
  if (!(zigosity %in% c("homo_healthy", "hete_healthy", NA, "")))stop("'status' of gene.specific file is not correct please modify it to homo_healthy, hete_healthy or NA")

  if (all.information$Variant.Info$gene %in% c("MLH1", "MSH2", "MSH6", "PMS2")){
    BS2.message <- "BS2 cannot be calculated automatically, please check if there exists co-occurrence in trans with a known pathogenic sequence variant in the same gene in a patient with colorectal cancer after age 45 "
  }else{
    #information from db
    #healthy individuals from flossies
    homo.healthy <- all.information$flossies.db %>%
      dplyr::filter (.data$population == "all") %>%
      dplyr::select ("homo") %>%
      as.numeric()
    hete.healthy <- all.information$flossies.db %>%
      dplyr::filter (.data$population == "all") %>%
      dplyr::select ("hete", "homo") %>%
      as.numeric() %>%
      sum()
    #homo from gnomad
    homo.gnomad <- ifelse((is.na(BS2.sup.db)||BS2.sup.db=="FLOSSIES"),
                          NA,
                          ifelse(BS2.sup.db == "GNOMAD_non_cancer",
                                 all.information$gnomAD$info$exomes.genomes$non.cancer$overall %>%
                                   dplyr::select("nhomalt") %>%
                                   as.numeric(),
                                 all.information$gnomAD$info$exomes.genomes$non.neuro$overall %>%
                                   dplyr::select("nhomalt") %>%
                                   as.numeric()))

    BS2.strong.select <- (!is.na(BS2.strong.db)&& ((zigosity=="homo_healthy" & homo.healthy >= BS2.strong) |(zigosity=="hete_healthy" & hete.healthy >= BS2.strong)))
    final.criteria$criteria.res["BS2", "strong"] <- ifelse (BS2.strong.select==TRUE & all.information$Variant.Info$gene!="CDH1", 1, 0)
    if (BS2.strong.select == TRUE){
      BS2.message <- ifelse ( zigosity=="homo_healthy",
                              paste("BS2 is assigned because there are", homo.healthy, "homozygotes in flossies db which is > than", BS2.strong),
                              ifelse(all.information$Variant.Info$gene !="CDH1",
                                     paste("BS2 is assigned because there are", hete.healthy, "healthy heterozygotes in flossies db which is > than", BS2.strong),
                                     paste("BS2 should be considered because there are", hete.healthy, "healthy heterozygotes in flossies db which is > than", BS2.strong, "but we have not assigned because we  do not know if their family  suggests HDGC or no.")))

    }


    if (!(1 %in% final.criteria$criteria.res["BS2","strong"])){
      BS2.sup.select <- (!is.na(BS2.sup.db)&& ((zigosity=="homo_healthy" & homo.gnomad >= BS2.sup) |(zigosity=="hete_healthy" & hete.healthy >= BS2.sup & all.information$Variant.Info$gene !="CDH1")))
      final.criteria$criteria.res["BS2", "supporting"] <- ifelse (BS2.sup.select==TRUE, 1, 0)
      if (BS2.sup.select == TRUE){
        BS2.message <- ifelse ( zigosity=="hete_healthy"& !(all.information$Variant.Info$gene %in% c("CDH1")),
                                paste("BS2_supporting is assigned because there are", hete.healthy, "healthy heterozygotes in flossies db which is > than", BS2.sup),
                                ifelse(all.information$Variant.Info$gene %in% c("CDH1"),
                                       paste("BS2_supporting whould be considered because there are", hete.healthy, "healthy heterozygotes in flossies db which is > than", BS2.sup, "but we have not assigned because we  do not know if their family  suggests HDGC or no."),
                                       ifelse (all.information$Variant.Info$gene %in% c("MLH1", "MSH2", "MSH6", "PMS2"),
                                               paste("BS2_supporting is assigned because there are", homo.gnomad, "homozygotes in", BS2.sup.db, " db which is >= than", BS2.sup, " (BS2_sup cutoff). Please check if they are 45 years old or more"),
                                               paste("BS2_supporting is assigned because there are", homo.gnomad, "homozygotes in", BS2.sup.db, " db which is >= than", BS2.sup, "(BS2_sup cutoff)."))))
      }
    }
    #Exceptions
    ##criteria not used for that gene
    final.criteria$criteria.res["BS2", c("strong", "supporting")][is.na(BS2.strong)] <- c(NA, NA)
    BS2.message[is.na(BS2.strong)] <- "BS2 does not apply for this gene."
    ##if it meets BS1, BS2 can only be supporting
    BS2.message[final.criteria$criteria.res["BS1", "strong"]==1 & final.criteria$criteria.res["BS2", "strong"]==1] <- "BS2_supporting is assigned because as BS1 is assigned BS2 can only achieve a supporting strength."
    final.criteria$criteria.res["BS2", c("strong", "supporting")][!is.na(final.criteria$criteria.res["BS1", "strong"]) && final.criteria$criteria.res["BS1", "strong"]==1 &&!is.na(final.criteria$criteria.res["BS2", "strong"]) && final.criteria$criteria.res["BS2", "strong"]==1] <- c(0,1)
  }
  BS2.message[is.null(BS2.message) || is.na(BS2.message)] <- "BS2 is not assigned taking into account the db queried. Please check if there exist more."
  final.criteria[["BS2.message"]] <- BS2.message
  return (final.criteria)
}

################################################################################
## BP1
################################################################################
BP1 <- function (all.information, final.criteria){
  final.criteria$criteria.res["BP1", 1:3] <- NA
  if(all.information$Variant.Info$most.severe.consequence !="missense_variant" | all.information$Variant.Info$gene %in% c("ATM", "CDH1", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "TP53" )){
    BP1.message <- "BP1 does not apply for this variant or gene."
    final.criteria$criteria.res["BP1", 4] <- NA
  }else{
    BP1.message <- "BP1 is not calculated. You should check if the variant is a missense and if in this gene the most cause of disease are truncating variants."
  }
  final.criteria[["BP1.message"]] <- BP1.message

  return(final.criteria)
}

################################################################################
## BP2
################################################################################
BP2 <- function (all.information, final.criteria){
  BP2.message <- NULL
  gene <- all.information$Variant.Info$gene
  final.criteria$criteria.res["BP2", 1:3] <- NA
  if (gene == "CDH1"){
    homo.strong <- as.integer(all.information$flossies.db$homo[3])
    homoz <- as.integer(all.information$gnomAD$info$exomes.genomes$non.cancer$overall$nhomalt)
    final.criteria$criteria.res["BP2", 4] <- ifelse(homoz>0, 1, 0)
    BP2.message <- ifelse(homoz>0,
                          paste ("BP2 is assigned because the variant is observed", homoz, "times as homozygous state in gnomAD v2.1.1 non_cancer.",
                                 ifelse(homo.strong==0,
                                        "",
                                        paste("It has also been found in the homozygous state in", homo.strong, "healthy individuals over 70 (FLOSSIES db) but we do not know their fanily history, for these reason we do not give BP2_strong. Also it should be considered Variant observed in trans w/known pathogenic variant."))),
                          "BP2 is denied because the variant is not observed in homozygous state in gnomAD v2.1.1 non_cancer. Please check if the variant is observed in cis or trans with a pathogenic variant. ")

  }else{
    final.criteria$criteria.res["BP2", 4] <- NA
    BP2.message <- "BP2 has to be calculated manually."
  }
  final.criteria[["BP2.message"]] <- BP2.message
  return(final.criteria)
}

################################################################################
## BP4
################################################################################
BP4 <- function (all.information, final.criteria){
  BP4.message <- NULL
  gene <- all.information$Variant.Info$gene
  final.criteria$criteria.res["BP4", 1:3] <- NA
  if (1 %in% final.criteria$criteria.res["PVS1",]){
    final.criteria$criteria.res["BP4", 4] <- NA
    BP4.message <- "BP4 is not calculated because PVS1 is met and co-usage is not permitted."
  }else if (all.information$Variant.Info$most.severe.consequence=="missense_variant" & gene %in% c("PALB2", "CDH1", "PTEN")){
    final.criteria$criteria.res["BP4", 4] <- NA
    BP4.message <- paste("BP4 is not calculated because it is a missense variant in", gene,"gene.")
  }else if ((gene %in% "PTEN")& stringr::str_detect(all.information$Variant.Info$variant, "c.79+1|c.79+2|c.80-2|c.80-1")){
    final.criteria$criteria.res["BP4", 4] <- NA
    BP4.message <- "BP4 is denied because PTEN guidelines specify that it should not to be applied for variants which may impact the intron 1 splice donor or acceptor sites. "
  }else{
    ben.splicing <- insilico (all.information$predictors$predictor.table, "Splicing Predictor", "Benign")
    pat.splicing <- insilico (all.information$predictors$predictor.table, "Splicing Predictor", "Pathogenic")
    ben.splicing.prior <- insilico (all.information$predictors$predictor.table, "Splicing Predictor", "Low")
    ben.splicing.prior2 <- insilico (all.information$predictors$predictor.table, "Splicing Predictor", "Innocuous")
    ben.splicing.prior3 <- insilico (all.information$predictors$predictor.table, "Splicing Predictor", "Improved")
    ben.splicing.prior4 <- insilico (all.information$predictors$predictor.table, "Splicing Predictor", "Minimal")
    ben.splicing.notapplicable <- insilico (all.information$predictors$predictor.table, "Splicing Predictor", "not applicable")
    alt.splicing.prior <- insilico(all.information$predictors$predictor.table, "Splicing Predictor", "Increased")
    alt.splicing.prior2 <- insilico(all.information$predictors$predictor.table, "Splicing Predictor", "High")
    sum.prior.splice.pat <- sum(nrow(alt.splicing.prior), nrow(alt.splicing.prior2))
    MMR <- all.information$Variant.Info$gene %in% c("MLH1", "MSH2", "MSH6") &
      !(stringr::str_detect(all.information$Variant.Info$variant, "del")) &
      (!(all.information$Variant.Info$most.severe.consequence %in% c("intron_variant")) |
         (all.information$Variant.Info$most.severe.consequence %in% c("intron_variant") &&
            !(purrr::map(stringr::str_extract_all(all.information$Variant.Info$variant, "[0-9]+"),2) %>% unlist %>% as.numeric > 25)))
    sum.prior.splice <- sum(nrow(ben.splicing.prior), nrow(ben.splicing.prior2), nrow(ben.splicing.prior3), nrow(ben.splicing.prior4), nrow(ben.splicing.notapplicable))
    all.utah.splicing <- rbind(ben.splicing.prior, ben.splicing.prior2,ben.splicing.prior3, ben.splicing.prior4)
    final.criteria$criteria.res["BP4",4][nrow(ben.splicing) <4 & MMR==FALSE] <-0
    final.criteria$criteria.res["BP4",4][MMR &&
                                           ((nrow(pat.splicing)>1)| sum.prior.splice!=2)] <- 0
    pathogenic.all <- rbind(pat.splicing, alt.splicing.prior, alt.splicing.prior2)

    #intronic and silence: we apply BP4 if all splicing predictors are in the benign cut-off
    intronic <- (stringr::str_detect(all.information$Variant.Info$variant, "[+]|[-]")&
                   (!all.information$Variant.Info$most.severe.consequence %in% c("5_prime_UTR_variant", "3_prime_UTR_variant")))
    ben.splicing.pred <- (nrow(ben.splicing) ==4 &&MMR==FALSE)| (MMR && (nrow(pat.splicing)==0 & sum.prior.splice==2))
    final.criteria$criteria.res["BP4",4] <- ifelse (ben.splicing.pred & (intronic == TRUE | all.information$Variant.Info$most.severe.consequence == "synonymous_variant"), 1, 0)
    BP4.message <- ifelse(ben.splicing.pred & (intronic == TRUE | all.information$Variant.Info$most.severe.consequence == "synonymous_variant"),
                          ifelse(MMR==FALSE,
                                 paste("BP4 is assigned because splicing predictors do not predict any effect (max Splice-AI value=", max(ben.splicing$values),")."),
                                 paste("BP4 is assigned because splicing predictors do not predict any effect (max Splice-AI value=", max(ben.splicing$values)," and", paste0(rownames(all.utah.splicing), "=", all.utah.splicing$values,collapse=" and"),").")),
                          #paste ("BP4 is not met because " , ben.splicing$predictor, "value is", ben.splicing$values, "which is > than", ben.splicing$BE.cut.off, "."))
                          ifelse( MMR==FALSE,
                                  ifelse( nrow(ben.splicing)<4,
                                          paste0("BP4 is denied because at least one splicing predictor score is above the benign cut-off (max Splice-AI value =", max(pat.splicing$values), ")."),
                                          ""),
                                  paste0("BP4 is denied because at least one splicing predictor score is above the benign cut-off (", paste0(rownames(pathogenic.all), "=", pathogenic.all$values, "which is not", pathogenic.all$operator.benign, pathogenic.all$BE.cut.off, collapse = " and "))))

    if (ben.splicing.pred & all.information$Variant.Info$most.severe.consequence == "missense_variant"){
      ben.prot <- insilico (all.information$predictors$predictor.table, "Protein effect", "Benign")
      condition <- (nrow(ben.prot)==1 & !gene %in% c("TP53")) |(nrow(ben.prot)==2 &gene %in% c("TP53"))
      final.criteria$criteria.res["BP4", 4] <- ifelse (condition, 1, 0)
      BP4.message <- ifelse ( condition,
                              ifelse(MMR==FALSE,
                                     paste0("BP4 is assigned because SpliceAi does not predict any effect (max value SpliceAI ", max(ben.splicing$values), ") and protein predictors neither (", paste(rownames(ben.prot), "=", ben.prot$values, "which is", ben.prot$operator.benign, "than", ben.prot$BE.cut.off, collapse=" and "),")."),
                                     paste0("BP4 is assigned because SpliceAi does not predict any effect (max value SpliceAI ", max(ben.splicing$values), "and ", paste0(rownames(all.utah.splicing), " = ", all.utah.splicing$values,collapse=" and "),") and protein predictors neither (", paste(rownames(ben.prot), "=", ben.prot$values, "which is", ben.prot$operator.benign, "than", ben.prot$BE.cut.off, collapse=" and "),").")),
                              # paste("BP4 is not met because", ben.prot$predictor, " value is ", ben.prot$values , "which is not ", ben.prot$operator.benign, ben.prot$BE.cut.off, "."))
                              paste0("BP4 is denied because the protein predictor is in the grey or pathogenic zone (", paste(rownames(ben.prot), "=", ben.prot$values, "which is not ", ben.prot$operator.benign, "than", ben.prot$BE.cut.off, collapse=" and "), ")."))
    }else if( all.information$Variant.Info$most.severe.consequence %in% c("inframe_deletion", "inframe_insertion", "inframe_duplication") && nrow(ben.splicing) == 4){
      provean <- insilico(all.information$predictors$predictor.table, "Protein effect", "Benign")
      final.criteria$criteria.res["BP4",4] <- ifelse(is.na(all.information$predictors$predictor.table["Provean", "values"]) || nrow(provean)==0,
                                                     0,
                                                     1)
      BP4.message <- ifelse(is.na(all.information$predictors$predictor.table["Provean", "values"]),
                            "Please consider installing or checking Provean score manually to assign or deny PP3.",
                            ifelse(nrow(provean)==0,
                                   paste("BP4 is denied because although Splice AI does not predict splicing alteration  Provean predicts pathogenic protein alteration (score", all.information$predictors$predictor.table["Provean", "values"], "which is < than -2.5)."),
                                   paste("BP4 is assigned because Splice AI does not predict splicing alteration and Provean no pathogenic protein effect (score", all.information$predictors$predictor.table["Provean", "values"], "which is > than -2.5).")))

    }

      }
  BP4.message <- ifelse(gene=="PTEN" & stringr::str_detect(all.information$Variant.Info$variant, "c.635-2|c.635-1"),
                        paste(BP4.message, "Variant impacts the intron 6 splice acceptor, this criterion should be used cautiously according to PTEN guidelines."),
                        BP4.message)
  final.criteria[["BP4.message"]] <- BP4.message
  return ( final.criteria)
}

################################################################################
## BP7
################################################################################
BP7 <- function (all.information, final.criteria){
  gene.specific.df <- all.information$gene.specific.info
  BP7.message <- NULL
  gene <- all.information$Variant.Info$gene
  intronic <- intronicCutoff (all.information$Variant.Info)
  #MMR.CDH1 <- gene %in% c("MLH1", "MSH2", "MSH6", "PMS2", "CDH1", "ATM")
  final.criteria$criteria.res["BP7", 1:3] <- NA
  if(!is.na(intronic) && intronic=="no"){
    final.criteria$criteria.res["BP7", 4] <-0

    BP7.message <- ifelse(!gene %in% c( "ATM", "PALB2", "CHEK2", "TP53"),
                          "BP7 is denied bc.1399T>Cecause variants must be positioned at or beyond  +7/-21.",
                          ifelse(gene %in% c("ATM"),
                                 "BP7 is denied because variants must be positioned beyond (but not including)  +7/-40.",
                                 ifelse(gene %in% c("PALB2"),
                                        "BP7 is denied because variants must be positioned beyond (but not including)  +20/-40.",
                                        "BP7 is denied because it not applies for intronic variants."
                          )))

  }else if (all.information$gene.specific.info$BP7_splicing =="Independent" && (!is.na(intronic)&&intronic =="yes" ||
                                all.information$Variant.Info$most.severe.consequence == "synonymous_variant")){
    final.criteria$criteria.res["BP7", 4] <- 1
    BP7.message <- ifelse( !is.na(intronic) && intronic=="yes" && gene %in% c("ATM"),
                           "BP7 is assigned because it is an intronic variant  located beyond +7/-40 .",
                           ifelse(!is.na(intronic) && intronic=="yes" ,
                                  "BP7 is assigned because it is an intronic variant  located beyond +7/-21 .",
                                  "BP7 is assigned because it is a synonymous variant."))
  }else if ((!all.information$gene.specific.info$BP7_splicing =="Independent"  & !(all.information$Variant.Info$gene%in%c("TP53")) & all.information$Variant.Info$most.severe.consequence %in% c("synonymous_variant", "intron_variant","splice_region_variant", "splice_donor_region_variant", "splice_acceptor_region_variant")) |
             (all.information$Variant.Info$gene%in% c("TP53") & all.information$Variant.Info$most.severe.consequence %in% c("synonymous_variant"))){
    alt.splicing <- insilico (all.information$predictors$predictor.table, "Splicing Predictor", "Benign") %>%
                                                                                                          dplyr::filter(stringr::str_detect(type, "SpliceAI"))
    if(nrow(alt.splicing)==0){
      if (gene %in% c("CHEK2", "PTEN", "TP53", "PALB2")){
      predictors.conservation <- insilico (all.information$predictors$predictor.table, "Nucleotide conservation", "Not strongly conserved")
      condition <- (nrow(predictors.conservation)==1 & gene != "PTEN" )| (nrow(predictors.conservation)==2 & gene =="PTEN")
      final.criteria$criteria.res["BP7", 4] <- ifelse( condition, 1, 0)
      BP7.message <- ifelse(condition,
                            "BP7 is assigned because nucleotide predictors are not strongly conserved and splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site",
                            "BP7 is denied because nucleotide predictors are strongly conserved")
      }else{
        BP7.message <-  ifelse(!is.na(intronic) && intronic=="yes" ,
                               "BP7 is assigned because it is an intronic variant  located beyond +7/-21 and splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site.",
                               "BP7 is assigned because it is a synonymous variant and splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site.")
        final.criteria$criteria.res["BP7", 4] <- 1
        }
      }else{
      BP7.message <- "BP7 is denied because at least one splicing predictor is above the pathogenic theresehold."
      final.criteria$criteria.res["BP7", 4] <- 0
    }
  }else{
    BP7.message <- "BP7 does not apply in these kind of variants."
    final.criteria$criteria.res["BP7", 4] <- NA
  }

  final.criteria[["BP7.message"]] <- BP7.message
  return ( final.criteria)
}

################################################################################
## Other criteria (don't apply)
################################################################################
notAppCriteria <- function (all.information , final.criteria){
  final.criteria$criteria.res["PP5", c(1:4)] <- NA
  final.criteria$criteria.res["BP6", c(1:4)] <- NA
  PP5.message <- "According to ClinGen PP5 should not be applied since variants in reputable sources are not always linked to the evidence on which they were based (Biesecker et al., 2018) "
  BP6.message <- "According to ClinGen BP6 should not be applied since variants in reputable sources are not always linked to the evidence on which they were based (Biesecker et al., 2018) "
  final.criteria[["PP5.message"]] <- PP5.message
  final.criteria[["BP6.message"]] <- BP6.message

  not.auto <- paste0(c("PS2", "PS4", "PM3", "PM6", "PP1", "PP4", "BS4", "BP1", "BP3", "BP5"),".message")
  for( i in 1:length(not.auto)){
    final.criteria[[not.auto[i]]] <- "Not automated criteria"
  }
  return(final.criteria)
}

