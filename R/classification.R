#' vaRclass
#' @param all.information obtained with the getVaRinfo function
#' @return the classification of the variant  and the criteria given following ACMG updated rules and gene specific for BRCA1, BRCA2, PALB2, ATM, MMR
#' @author Elisabet Munté Roca
#' @examples
#' assembly <- "hg19"
#' gene <- "BRCA1"
#' variant <- "c.211A>G"
#' gene.specific.df <- read.csv ("./gen_specific.csv")
#' variant.info <- getVaRinfo(assembly, gene, variant, nm.info = NULL, gene.specific.df)
#' vaRclass (variant.info)
vaRclass <- function (all.information){
  criteria.res <- matrix ("NC",nrow=29, ncol=4)
  colnames(criteria.res) <- c("very_strong", "strong", "moderate", "supporting")
  criteria <- c("PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3",
                "PM4", "PM5", "PM6","PP1", "PP2", "PP3", "PP4", "PP5", "BA1",
                "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4",
                "BP5", "BP6", "BP7", "external")
  rownames(criteria.res) <- criteria
  final.criteria <- list (criteria.res = criteria.res)
  print("PVS1")
  final.criteria <- PVS1 (all.information = all.information, final.criteria = final.criteria)
  print("PS1")
  final.criteria <- PS1 (all.information = all.information, final.criteria =final.criteria)
  print("PS3_BS3")
  final.criteria <- PS3_BS3(all.information = all.information, final.criteria = final.criteria)
  print("PM1")
  final.criteria <- PM1 (all.information = all.information, final.criteria =final.criteria)
  print("PM2")
  final.criteria <- PM2 (all.information = all.information, final.criteria =final.criteria)
  print("PM4")
  final.criteria <- PM4 (all.information = all.information, final.criteria =final.criteria)
  print("PM5")
  final.criteria <- PM5 (all.information = all.information, final.criteria =final.criteria)
  print("PP2")
  final.criteria <- PP2 (all.information = all.information , final.criteria = final.criteria)
  print("PP3")
  final.criteria <- PP3 (all.information = all.information , final.criteria = final.criteria)
  print("BA1")
  final.criteria <- BA1 (all.information = all.information, final.criteria = final.criteria)
  print("BS1")

  final.criteria <- BS1 (all.information = all.information, final.criteria = final.criteria)
  print("BS2")
  final.criteria <- BS2 (all.information = all.information, final.criteria = final.criteria)
  print("BP1")
  final.criteria <- BP1 (all.information = all.information, final.criteria)
  print("BP2")
  final.criteria <- BP2 (all.information = all.information, final.criteria)
  print("BP4")
  final.criteria <- BP4 (all.information = all.information, final.criteria)
  print("BP7")
  final.criteria <- BP7 (all.information = all.information, final.criteria)
  print("other")
  final.criteria <- notAppCriteria(all.information = all.information, final.criteria = final.criteria)
  final.class <- finalClass(all.information = all.information, final.criteria=final.criteria)
  return (list(final.classification = final.class, final.criteria = final.criteria))
}

################################################################################
## Final Classification functions
################################################################################
finalClass <- function (all.information, final.criteria){
  discrep.reason <- NA
  results.pat <- final.criteria$criteria.res[1:16,]
  results.pat[results.pat=="NC"|is.na(results.pat)] <- as.numeric(0)
  sum.pat <- apply (results.pat, 2, function(x){
    as.numeric(x) %>%
                  sum()
  })

  total.pat <- sum(sum.pat)
  weigth.pat <- sum.pat * c(8,4,2,1)

  results.ben <- final.criteria$criteria.res[17:29,]
  results.ben[results.ben=="NC" | is.na(results.ben)] <- as.numeric(0)
  sum.ben <- apply (results.ben, 2, function(x){
    as.numeric(x) %>% sum()
  })
  total.ben <- sum(sum.ben[2:4])
  total.ben.ba1 <- sum.ben[1]
  weigth.ben <- sum.ben * c(-8,-4,-2,-1)

  total.pat.weigth <- sum(weigth.pat)
  total.ben.weigth <- sum(weigth.ben)
  #checking if at least are 2 or more criteria of the same nature
  if ( total.pat >1 | total.ben > 1 | total.ben.ba1 ==1){

    pm2.mod <-  final.criteria$criteria.res["PM2",3] == 1
    pm2.sup <- final.criteria$criteria.res["PM2",4] == 1

    discrep <- FALSE
    discrep.lb.lp <- FALSE
    #discrepancies
    if ((total.pat.weigth >1 & total.ben.weigth < (-1) & pm2.mod == FALSE)|
        (total.pat.weigth >2 & total.ben.weigth < (-1) & pm2.mod == TRUE)){
      discrep <- TRUE
      discrep.reason <- "There are discordant evidence elements and the classification cannot be achieved"
    }else if (total.pat.weigth==2 & total.ben.weigth < (-1) & pm2.mod == TRUE){
      discrep.lb.lp <- TRUE
      total.pat.weigth <- total.pat.weigth - 2
      discrep.reason <- "PM2_mod has not been taken into account because it is a discordant evidence. The variant can only achieve a LB classification."
    }else if (total.pat.weigth == 1 & total.ben.weigth < (-1) & pm2.sup == TRUE){
      discrep.lb.lp <- TRUE
      total.pat.weigth <- total.pat.weigth - 1
      discrep.reason <- "PM2_sup has not been taken into account in the sum because it is a discordant evidence. The variant can only achieve a LB classification."
    }else if ((total.pat.weigth == 1 & total.ben.weigth < (-1) & !(1 %in% final.criteria$criteria.res["PP3",]))|
              (total.pat.weigth > 1 & total.ben.weigth == (-1) & !(1 %in% final.criteria$criteria.res["BP4",]))){
      discrep.lb.lp <- TRUE
      discrep.reason <- "It exists a discordant evidence. The variant can only achieve a LB/LP classification."
    }else if ((total.pat.weigth == 1 & total.ben.weigth < (-1) & 1 %in% final.criteria$criteria.res["PP3",])|
              (total.pat.weigth > 1 & total.ben.weigth == (-1) & 1 %in% final.criteria$criteria.res["BP4",])){
      discrep.reason <- "It exists a discordant evidence. As it is an in silico evidence, it will be taken into account in the sum but the variant can achieve a pat or benign classification."
    }else{
      discrep.reason <- "There are no discordant evidences"
    }

    sum.total <- total.pat.weigth +  total.ben.weigth
    pat  <- discrep == FALSE & sum.total >= 10 & discrep.lb.lp == FALSE
    ppat <- discrep == FALSE & ((sum.total >= 6 & sum.total <= 9 ) | (sum.total>=10&discrep.lb.lp==TRUE))
    vus  <- discrep == FALSE & sum.total >= 0 & sum.total <= 5
    pben <- discrep == FALSE & ((sum.total >=- 5 & sum.total <= -1 ) | (sum.total <=-6 & discrep.lb.lp==TRUE))
    ben  <- discrep == FALSE & sum.total<=-6 & discrep.lb.lp==FALSE

    class.final <- ifelse (pat==TRUE,
                           "Pathogenic",
                           ifelse ( ppat == TRUE,
                                    "Likely Pathogenic",
                                    ifelse (vus ==TRUE,
                                            "VUS",
                                            ifelse (pben == TRUE,
                                                    "Likely Benign",
                                                    ifelse(ben == TRUE,
                                                           "Benign",
                                                           "VUS, discrepancies.")))))

  }else if ( total.pat==1 & total.ben ==1){
    sum.total <- total.pat.weigth +total.ben.weigth
    class.final <- "Cannot be calculated, one evidence supporting pathogenicity and one benignity."
    discrep.reason <- "One evidence supporting pathogenicity and one benignity."
  }else if (total.pat == 0 & total.ben == 0){
    sum.total <- 0
    class.final <-  "N/A: no evidences."
  } else{
    sum.total <- total.pat.weigth + total.ben.weigth
    class.final <- "N/A: single evidence element."
  }

  if (class.final == "N/A: single evidence element" & all.information$Variant.Info$gene=="CDH1" & final.criteria$criteria.res["BS1", 2]== 1){
    class.final <- "Likely Benign"
    discrep.reason <- "We allow a variant to reach a likely benign classification based on BS1 stand alone."
  }
  criteria.assigned <- criteriaAssigned (final.criteria)
  return (list (final.class = class.final, criteria.assigned = criteria.assigned, sum.criteria = sum.total, discrep.reason=discrep.reason))
}

criteriaAssigned <- function(final.criteria){
  vs <- which(final.criteria$criteria.res[,1]==1)
  s <- which(final.criteria$criteria.res[,2]==1)
  mod <- which(final.criteria$criteria.res[,3]==1)
  sup <- which(final.criteria$criteria.res[,4]==1)

  vs_names <- ifelse(rep(length(vs)>0, length(vs)),
                     paste0(names(vs), "_veryStrong"),
                     NA)
  s_names <- ifelse(rep(length(s)>0,length(s)),
                    paste0(names(s), "_strong"),
                    NA)
  mod_names <- ifelse(rep(length(mod)>0,length(mod)),
                      paste0(names(mod), "_moderate"),
                      NA)
  sup_names <- ifelse(rep(length(sup)>0, length(sup)),
                      paste0(names(sup), "_supporting"),
                      NA)
  criteria.given <- individualCriteria(c(vs_names, s_names, mod_names, sup_names))
  criteria.given <-criteria.given[!is.na(criteria.given)]

  return (criteria.given)
}

individualCriteria <- function(criteria.met.for){
  criteria.met.for <- stringr::str_replace_all ( criteria.met.for,
                                        c("PVS1_veryStrong"="PVS1",
                                          "PS1_strong"="PS1",
                                          "PS2_strong"="PS2",
                                          "PS3_strong"="PS3",
                                          "PS4_strong"="PS4",
                                          "PM1_moderate"="PM1",
                                          "PM2_moderate"="PM2",
                                          "PM3_moderate"="PM3",
                                          "PM4_moderate"="PM4",
                                          "PM5_moderate"="PM5",
                                          "PM6_moderate"="PM6",
                                          "PP1_supporting"="PP1",
                                          "PP2_supporting"="PP2",
                                          "PP3_supporting"="PP3",
                                          "PP4_supporting"="PP4",
                                          "PP5_supporting"="PP5",
                                          "BA1_veryStrong"="BA1",
                                          "BS1_strong"="BS1",
                                          "BS2_strong"="BS2",
                                          "BS3_strong"="BS3",
                                          "BS4_strong"="BS4",
                                          "BP1_supporting"="BP1",
                                          "BP2_supporting"="BP2",
                                          "BP3_supporting"="BP3",
                                          "BP4_supporting"="BP4",
                                          "BP5_supporting"="BP5",
                                          "BP6_supporting"="BP6",
                                          "BP7_supporting"="BP7"
                                        ))

  return (criteria.met.for)
}

sumCriteria <- function (vaRclass){
  results.pat <- vaRclass$final.criteria$criteria.res[1:16,]
  results.pat[results.pat=="NC" | is.na(results.pat)] <- as.numeric(0)
  results.ben <- vaRclass$final.criteria$criteria.res[17:29,]
  results.ben[results.ben=="NC" | is.na(results.ben)] <- 0


  sum.criteria <- list( verystrong.pat = sum(as.numeric(results.pat[,1])),
                        strong.pat = sum(as.numeric(results.pat[,2])),
                        moderate.pat = sum(as.numeric(results.pat[,3])),
                        supporting.pat = sum(as.numeric(results.pat[,4])),
                        verystrong.benign = sum(as.numeric(results.ben[,1])),
                        strong.benign <- sum(as.numeric(results.ben[,2])),
                        moderate.benign <- sum(as.numeric(results.ben[,3])),
                        support.benign <- sum(as.numeric(results.ben[,4])))

  return(sum.criteria)
}

completCriteria <- function(criteria.res){
  final.evidence <- c()
  crit <- c()
  for (i in 1:nrow(criteria.res)){
    crit <- criteriaMet(rownames(criteria.res)[i], criteria.res)
    final.evidence <- rbind(final.evidence, crit)
  }
  return(final.evidence)
}

criteriaMet <- function(criterion, criteria.res){
  h <- names(which(criteria.res[criterion,]==1))
  h0 <- names(which(criteria.res[criterion,]==0))
  hNA <- names(which(is.na(criteria.res[criterion,])))
  hNC <- names(which(is.na(criteria.res[criterion,])))
  final <- ifelse(rep(length(h)==1,2), c(1,h),
                  ifelse(rep(length(h0)>0, 2), c(0, NA),
                         ifelse(rep(length(hNA)>0,2), c("NA", NA),
                                c("NC", NA))))
  return(final)
}

