

vaRreport <-   function(vaRinfo, vaRclass, output.dir=NULL ){
  if(!requireNamespace("XLConnect", quietly=TRUE)){
    warning("Please install package 'XLConnect' when using 'excel.results = FALSE'")
  }else{
  path.original.file = system.file("extdata", "template.xlsx", package="vaRHC")
  #if(!(file.exists(path.original.file))) stop("The excel template path must be provided, you can download it from https://github.com/emunte/vaRHC/inst/extdata/template.xlsx")
  file.name <- paste0(vaRinfo$Variant.Info$gene,"_", stringr::str_replace(vaRinfo$Variant.Info$variant,">","-"),"_", Sys.Date(),".xlsx", "")  %>%
    stringr::str_replace_all("\\*", "asterisk")#name of the file, with the variant

  ## -------------------Check output directory and create excel subfolder
  output.dir.excels <- checkDir (output.dir, "excels")

  file.copy(from = path.original.file,
            to = file.path(output.dir.excels, file.name),
            recursive = FALSE,
            copy.mode = TRUE,
            copy.date = FALSE,
            overwrite = TRUE) #Plantilla cloned

  ###-----------------Load workbook; creating if not existing``
  wb <- XLConnect::loadWorkbook(file.path(output.dir.excels,file.name), create = FALSE)
  XLConnect::setStyleAction(wb,XLConnect::XLC$"STYLE_ACTION.NONE")

  ###-----------------Every sheet ----
  gene <- vaRinfo$Variant.Info$gene
  variant <- vaRinfo$Variant.Info$variant
  protein <- vaRinfo$Variant.Info$protein
  XLConnect::writeWorksheet(wb, gene, sheet=c(2:5,7:14),
                            startRow=1, startCol=2,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, variant, sheet=c(2:5,7:14),
                            startRow=2, startCol=2,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, protein, sheet=c(2:5,7:14),
                            startRow=3, startCol=2,
                            header=FALSE)

  XLConnect::writeWorksheet(wb, gene, sheet=6,
                            startRow=2, startCol=2,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, variant, sheet=6,
                            startRow=3, startCol=2,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, protein, sheet=6,
                            startRow=4, startCol=2,
                            header=FALSE)

  ###-----------------Classification Summary sheet ----
  XLConnect::writeWorksheet(wb,  vaRinfo$Variant.Info[,1:17], sheet="Classification Summary", #variant info
                            startRow=4, startCol=2,
                            header=FALSE)
  XLConnect::writeWorksheet(wb,  vaRinfo$Variant.Info[,18:23] %>% as.data.frame(), sheet="Classification Summary", #variant info
                            startRow=4, startCol=19,
                            header=FALSE)

  XLConnect::writeWorksheet(wb,  vaRinfo$Variant.Info$domain.info, sheet="Classification Summary", #variant info
                            startRow=4, startCol=25,
                            header=FALSE)

  if(length(vaRinfo$Variant.Info.other) > 0){
    for(i in 1:length(vaRinfo$Variant.Info.other)){
      j=1
      if(nrow(vaRinfo$Variant.Info.other[[i]])>0){
      XLConnect::writeWorksheet(wb,  vaRinfo$Variant.Info.other[[i]][,1:24], sheet="Classification Summary", #variant info
                                startRow=4+j, startCol=2,
                                header=FALSE)
      j=j+1
      # XLConnect::writeWorksheet(wb,  vaRinfo$Variant.Info.other[[i]][,18], sheet="Classification Summary", #variant info
      #                           startRow=4+i, startCol=20,
      #                           header=FALSE)
      }
    }
  }

  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$criteria.res, sheet="Classification Summary",
                            startRow=11, startCol=2,
                            header=FALSE)


  XLConnect::writeWorksheet(wb, as.matrix(sumCriteria(vaRclass), nrow=8), sheet="Classification Summary",
                            startRow=12, startCol=8,
                            header=FALSE)

  weigth <- c(8,4,2,1,-8,-4,-2,-1)
  sum.criteria2 <- unlist(sumCriteria(vaRclass))* weigth
  XLConnect::writeWorksheet(wb, as.matrix(sum.criteria2, nrow=8), sheet="Classification Summary",
                 startRow=12, startCol=9,
                 header=FALSE)
  XLConnect::writeWorksheet(wb, paste(vaRclass$final.classification$criteria.assigned, collapse=","), sheet="Classification Summary",
                            startRow=24, startCol=8,
                            header=FALSE)
  # XLConnect::writeWorksheet(wb, pasteCriteria(vaRclass$final.classification$criteria.assigned), sheet="Classification Summary",
  #                           startRow=24, startCol=8,
  #                           header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.classification$final.class, sheet="Classification Summary",
                            startRow=25, startCol=8,
                            header=FALSE)

  XLConnect::writeWorksheet(wb, as.matrix(vaRclass$final.classification$sum.criteria), sheet="Classification Summary",
                            startRow=17, startCol=10,
                            header=FALSE)

  XLConnect::writeWorksheet(wb, as.matrix(vaRinfo$variant.correction$warning), sheet="Classification Summary",
                            startRow=17, startCol=13,
                            header=FALSE)
  XLConnect:: writeWorksheet(wb,vaRclass$final.classification$discrep.reason , sheet="Classification Summary",
                             startRow=12, startCol=10,
                             header=FALSE)

  if(length(vaRinfo$clinVar$clinVar.info$variant) >0 &&length(vaRinfo$clinVar$clinVar.info$variant[[1]]) > 0){
    XLConnect:: writeWorksheet(wb,vaRinfo$clinVar$clinVar.info$variant[[1]]$Classification$status , sheet="Classification Summary",
                               startRow=9, startCol=13,
                               header=FALSE)
    XLConnect:: writeWorksheet(wb,vaRinfo$clinVar$clinVar.info$variant[[1]]$Classification$interpret , sheet="Classification Summary",
                               startRow=12, startCol=13,
                               header=FALSE)

  }

  functionals.messages <- ifelse(vaRclass$final.criteria$PS3.message==vaRclass$final.criteria$BS3.message,
                                 "BS3 and PS3 are not assigned because variant is not found in the automated literature. Please check for more. ",
                                 paste(vaRclass$final.criteria$PS3.message, vaRclass$final.criteria$BS3.message))

  all.messages <- paste0(vaRclass$final.criteria$PVS1.message, functionals.messages, vaRclass$final.criteria$PM1.message, vaRclass$final.criteria$PM2.message, vaRclass$final.criteria$PM4.message, vaRclass$final.criteria$PM5.message, vaRclass$final.criteria$PP3.message, vaRclass$final.criteria$BA1.message, vaRclass$final.criteria$BS1.message, vaRclass$final.criteria$BS2.message, vaRclass$final.criteria$BP2.message, vaRclass$final.criteria$BP4.message, vaRclass$final.criteria$BP7.message)
  XLConnect:: writeWorksheet(wb, all.messages , sheet="Classification Summary",
                             startRow=26, startCol=8,
                             header=FALSE)


  ###----------------Evidence sheet ----
  criteria <- completCriteria(vaRclass$final.criteria$criteria.res) %>%
                                                                    as.data.frame()
  XLConnect::writeWorksheet(wb, introMessage(vaRinfo), sheet="Evidence",
                            startRow=6, startCol=17,
                            header=FALSE)

  XLConnect::writeWorksheet(wb, criteria[1:16,], sheet="Evidence",
                            startRow=9, startCol=13,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, criteria[17:28,], sheet="Evidence",
                            startRow=28, startCol=13,
                            header=FALSE)

  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PVS1.message, sheet="Evidence",
                            startRow=9, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PS1.message, sheet="Evidence",
                            startRow=10, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PS2.message, sheet="Evidence",
                            startRow=11, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PS3.message, sheet="Evidence",
                            startRow=12, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PS4.message, sheet="Evidence",
                            startRow=13, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PM1.message, sheet="Evidence",
                            startRow=14, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PM2.message, sheet="Evidence",
                            startRow=15, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PM3.message, sheet="Evidence",
                            startRow=16, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PM4.message, sheet="Evidence",
                            startRow=17, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PM5.message, sheet="Evidence",
                            startRow=18, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PM6.message, sheet="Evidence",
                            startRow=19, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PP1.message, sheet="Evidence",
                            startRow=20, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PP2.message, sheet="Evidence",
                            startRow=21, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PP3.message, sheet="Evidence",
                            startRow=22, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PP4.message, sheet="Evidence",
                            startRow=23, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PP5.message, sheet="Evidence",
                            startRow=24, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BA1.message, sheet="Evidence",
                            startRow=28, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BS1.message, sheet="Evidence",
                            startRow=29, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BS2.message, sheet="Evidence",
                            startRow=30, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BS3.message, sheet="Evidence",
                            startRow=31, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BS4.message, sheet="Evidence",
                            startRow=32, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BP1.message, sheet="Evidence",
                            startRow=33, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BP2.message, sheet="Evidence",
                            startRow=34, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BP3.message, sheet="Evidence",
                            startRow=35, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BP4.message, sheet="Evidence",
                            startRow=36, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BP5.message, sheet="Evidence",
                            startRow=37, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BP6.message, sheet="Evidence",
                            startRow=38, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$BP7.message, sheet="Evidence",
                            startRow=39, startCol=17,
                            header=FALSE)




  ###----------------Control_freq_sheet ----

  ##Conclusion
  XLConnect::writeWorksheet(wb, paste0(vaRclass$final.criteria$BA1.message, vaRclass$final.criteria$BS1.message, vaRclass$final.criteria$BS2.message, vaRclass$final.criteria$PM2.message), sheet="Control_Freq", #exomes coverage
                            startRow=1, startCol=5,
                            header=FALSE)

  ## Variants name in gnomAD
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$nomenclature, sheet="Control_Freq", #exomes coverage
                            startRow=3, startCol=4,
                            header=FALSE)

  ##coverage

  message.cov <- ifelse ( vaRinfo$gnomAD$coverage$exomes <20 | vaRinfo$gnomAD$coverage$genomes <20,
                          "There is not enought coverage in gnomAD. Frequencies will not be taken into account.",
                          "The region is well covered.")
  # variant.found <- ifelse(sum(vaRinfo$gnomAD$info$exomes$non.cancer$overall$AC[1:2])==0 & sum(vaRinfo$gnomAD$info$geomes$non.cancer$overall[1:2])==0,
  #                         "The variant hasn't been found in gnomAD non-cancer",
  #                         "The variant has")


  XLConnect::writeWorksheet(wb, message.cov, sheet="Control_Freq", #exomes coverage
                            startRow=3, startCol=9,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$coverage$exomes, sheet="Control_Freq", #exomes coverage
                            startRow=3, startCol=13,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$coverage$genomes, sheet="Control_Freq", #exomes coverage
                            startRow=3, startCol=16,
                            header=FALSE)


  #non_cancer

  XLConnect::writeWorksheet(wb, paste("CI ", vaRinfo$gene.specific.info$IC) , sheet="Control_Freq", #non_cancer exomes
                            startRow=9, startCol=9,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, paste("CI", vaRinfo$gene.specific.info$IC) , sheet="Control_Freq", #non_cancer exomes
                            startRow=25, startCol=9,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, paste("CI", vaRinfo$gene.specific.info$IC) , sheet="Control_Freq", #non_cancer exomes
                            startRow=41, startCol=9,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, paste("CI",vaRinfo$gene.specific.info$IC) , sheet="Control_Freq", #non_cancer exomes
                            startRow=57, startCol=9,
                            header=FALSE)
  #·····Exomes
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$exomes$non.cancer$subpopulations, sheet="Control_Freq", #non_cancer exomes
                            startRow=10, startCol=5,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$exomes$non.cancer$gender, sheet="Control_Freq", #non_cancer exomes (male/female)
                            startRow=19, startCol=5,
                            header=FALSE)
  #·······Genomes
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$genomes$non.cancer$subpopulations, sheet="Control_Freq", #non_cancer exomes
                            startRow=26, startCol=5,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$genomes$non.cancer$gender, sheet="Control_Freq", #non_cancer exomes (male/female)
                            startRow=35, startCol=5,
                            header=FALSE)
  #·······exomes + genomes
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$exomes.genomes$non.cancer$subpopulations[,2:6], sheet="Control_Freq", #non_cancer exomes
                            startRow=42, startCol=5,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$exomes.genomes$non.cancer$gender[,2:6], sheet="Control_Freq", #non_cancer exomes (male/female)
                            startRow=51, startCol=5,
                            header=FALSE)

  #·······overall
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$exomes.genomes$non.cancer$overall[,2:6], sheet="Control_Freq", #non_cancer exomes
                            startRow=58, startCol=5,
                            header=FALSE)


  #non-neuro
  XLConnect::writeWorksheet(wb, paste("CI", vaRinfo$gene.specific.info$IC) , sheet="Control_Freq", #non_cancer exomes
                            startRow=9, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, paste("CI", vaRinfo$gene.specific.info$IC) , sheet="Control_Freq", #non_cancer exomes
                            startRow=25, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, paste("CI", vaRinfo$gene.specific.info$IC) , sheet="Control_Freq", #non_cancer exomes
                            startRow=41, startCol=17,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, paste("CI", vaRinfo$gene.specific.info$IC) , sheet="Control_Freq", #non_cancer exomes
                            startRow=57, startCol=17,
                            header=FALSE)
  #·······exomes
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$exomes$non.neuro$subpopulations, sheet="Control_Freq", #non_neuro exomes
                            startRow=10, startCol=13,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$exomes$non.neuro$gender, sheet="Control_Freq", #non_neuroexomes (male/female)
                            startRow=19, startCol=13,
                            header=FALSE)
  #·······genomes
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$genomes$non.neuro$subpopulations, sheet="Control_Freq", #non_cancer exomes
                            startRow=26, startCol=13,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$genomes$non.neuro$gender, sheet="Control_Freq", #non_cancer exomes (male/female)
                            startRow=35, startCol=13,
                            header=FALSE)
  #·······exomes+genomes
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$exomes.genomes$non.neuro$subpopulations[,2:6], sheet="Control_Freq", #non_cancer exomes
                            startRow=42, startCol=13,
                            header=FALSE)
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$exomes.genomes$non.neuro$gender[,2:6], sheet="Control_Freq", #non_cancer exomes (male/female)
                            startRow=51, startCol=13,
                            header=FALSE)
  #·······overall
  XLConnect::writeWorksheet(wb, vaRinfo$gnomAD$info$exomes.genomes$non.neuro$overall[,2:6], sheet="Control_Freq", #non_cancer exomes
                            startRow=58, startCol=13,
                            header=FALSE)


  #flossies db
  XLConnect::writeWorksheet(wb, vaRinfo$flossies.db, sheet="Control_Freq", #non_cancer exomes
                            startRow=66, startCol=4,
                            header=FALSE, rownames= TRUE)



  #ClinVar

  #··········VARIANT----

  if(nrow(vaRinfo$clinVar$clinVar.ids$table)>0){
    id.variant <- vaRinfo$clinVar$clinVar.ids$table %>%
                                                    dplyr::filter(variant==vaRinfo$Variant.Info$variant) %>%
                                                    dplyr::select("V1")
  }

  if(exists("id.variant")&&nrow(id.variant)>0 &&length(vaRinfo$clinVar$clinVar.info$variant)>0){
    XLConnect::writeWorksheet(wb, paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", id.variant), sheet="ClinVar",
                              startRow=13, startCol=2,
                              header=FALSE, rownames= FALSE)

    XLConnect::writeWorksheet(wb, purrr::map(vaRinfo$clinVar$clinVar.info$variant[[1]][8],2)%>%unlist, sheet="ClinVar", #non_cancer exomes
                              startRow=2, startCol=5,
                              header=FALSE, rownames= TRUE)
    XLConnect::writeWorksheet(wb, purrr::map(vaRinfo$clinVar$clinVar.info$variant[[1]][8],1)%>%unlist, sheet="ClinVar",
                              startRow=3, startCol=5,
                              header=FALSE, rownames= TRUE)
    XLConnect::writeWorksheet(wb, vaRinfo$clinVar$clinVar.info$variant[[1]]$Submitted.interpretations.and.evidence%>% as.data.frame(), sheet="ClinVar",
                              startRow=15, startCol=1,
                              header=FALSE, rownames= FALSE)
    XLConnect::writeWorksheet(wb, vaRinfo$clinVar$clinVar.info$variant[[1]]$Citations%>% as.data.frame(), sheet="ClinVar",
                              startRow=15, startCol=9,
                              header=FALSE, rownames= FALSE)

    XLConnect::writeWorksheet(wb, matrix(unlist(vaRinfo$clinVar$clinVar.info$variant[[1]][9]), nrow=5), sheet="ClinVar",
                              startRow=6, startCol=2,
                              header=FALSE, rownames= FALSE)
    #var.clin <- ifelse( exists("vaRinfo$clinVar$clinVar.ids$table$blosum"),
   var.clin <-    vaRinfo$clinVar$clinVar.ids$table %>%
     dplyr::mutate(url=paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", .data$V1))  %>%
     dplyr::filter(stringr::str_detect(.data$V2, vaRinfo$Variant.Info$variant)) %>%
     dplyr::select(-"V1"
                   #, -"variant", -"protein"
                   ) %>%
     dplyr::relocate(.data$url, .after= .data$V2)

    XLConnect::writeWorksheet(wb, var.clin , sheet="ClinVar Variants",
                              startRow=9, startCol=1,
                              header=FALSE, rownames= FALSE)
    XLConnect::writeWorksheet(wb, vaRinfo$clinVar$clinVar.info$variant[[1]]$Classification , sheet="ClinVar Variants",
                              startRow=9, startCol=8,
                              header=FALSE, rownames= FALSE)
  }else{
    XLConnect::writeWorksheet(wb, "The variant is not found in clinVar db.", sheet="ClinVar",
                              startRow=13, startCol=2,
                              header=FALSE, rownames= FALSE)
  }


  #·············OTHER VARIANTS----
  if (length(vaRinfo$clinVar$clinVar.info$same_codon)>0){
    for (j in 1:length(vaRinfo$clinVar$clinVar.info$same_codon)){
      variant.name <- purrr::map(stringr::str_split(names(vaRinfo$clinVar$clinVar.info$same_codon[j]),":"),2) %>%
                                                                                                              unlist
      var <- stringr::str_split(variant.name, " ") %>%
                                                   purrr::map(1) %>%
                                                   unlist
      sheet.name <- paste0("ClinVar ", variant.name)
      if(stringr::str_length(sheet.name) > 30){
        sheet.name <- stringr::str_sub(sheet.name,1,30)
      }
      XLConnect::cloneSheet(wb, 6, name=sheet.name)#we create sheets for each variant in the same codon found in clinvar
      XLConnect::writeWorksheet(wb, vaRinfo$clinVar$clinVar.ids$table %>%
                                                                      dplyr::filter(stringr::str_detect(.data$V2, vaRinfo$Variant.Info$variant)) %>%
                                                                      dplyr::select("blosum", "prior", "grantham") %>%
                                                                      t(),
                                sheet=sheet.name,
                                startRow=5, startCol=2,
                                header=FALSE, rownames= FALSE)
      XLConnect::setSheetPos(wb, sheet.name, 6+j)
      id.variant <- vaRinfo$clinVar$clinVar.ids$table %>%
                                                      dplyr::filter(.data$variant==var) %>%
                                                      dplyr::select("V1")
      num <- which(stringr::str_detect(as.character(vaRinfo$clinVar$clinVar.ids$table$V2), var))
      XLConnect::writeWorksheet(wb, vaRinfo$clinVar$clinVar.ids$table %>%
                                                                      dplyr::filter(stringr::str_detect(.data$V2, var)) %>%
                                                                      dplyr::mutate(url=paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", .data$V1))  %>%
                                                                      dplyr::select("V2", "url", "a1", "a2","blosum", "grantham", "prior") %>% as.data.frame(), sheet="ClinVar Variants",
                                startRow=9+j, startCol=1,
                                header=FALSE, rownames= FALSE)
      XLConnect::writeWorksheet(wb, vaRinfo$clinVar$clinVar.info$same_codon[[j]]$Classification , sheet="ClinVar Variants",
                                startRow=9+j, startCol=8,
                                header=FALSE, rownames= FALSE)
      XLConnect::writeWorksheet(wb, paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", id.variant), sheet=sheet.name,
                                startRow=15, startCol=2,
                                header=FALSE, rownames= FALSE)

      XLConnect::writeWorksheet(wb, variant.name, sheet=sheet.name,
                                startRow=3, startCol=5,
                                header=FALSE, rownames= FALSE)
      if(vaRinfo$Variant.Info$most.severe.consequence=="missense_variant"){
        XLConnect::writeWorksheet(wb, vaRinfo$clinVar$clinVar.ids$table %>%
                                                                        dplyr::filter(stringr::str_detect(.data$variant, vaRinfo$Variant.Info$variant)) %>%
                                                                        dplyr::select("blosum", "prior", "grantham") %>%
                                                                        t(),
                                  sheet=sheet.name,
                                  startRow=6, startCol=5,
                                  header=FALSE, rownames= FALSE)
        XLConnect::writeWorksheet(wb, t(vaRinfo$clinVar$clinVar.ids$table[num, c("blosum", "prior", "grantham")]), sheet=sheet.name,
                                  startRow=5, startCol=5,
                                  header=FALSE, rownames= FALSE)
      }
      XLConnect::writeWorksheet(wb, purrr::map(vaRinfo$clinVar$clinVar.info$same_codon[[j]][8],2) %>% unlist, sheet=sheet.name,
                                startRow=8, startCol=5,
                                header=FALSE, rownames= TRUE)
      XLConnect::writeWorksheet(wb, purrr::map(vaRinfo$clinVar$clinVar.info$same_codon[[j]][8],1) %>% unlist, sheet=sheet.name,
                                startRow=9, startCol=5,
                                header=FALSE, rownames= TRUE)
      XLConnect::writeWorksheet(wb, vaRinfo$clinVar$clinVar.info$same_codon[[j]][5], sheet=sheet.name,
                                startRow=20, startCol=1,
                                header=FALSE, rownames= FALSE)
      XLConnect::writeWorksheet(wb, vaRinfo$clinVar$clinVar.info$same_codon[[j]][7], sheet=sheet.name,
                                startRow=20, startCol=9,
                                header=FALSE, rownames= FALSE)
      XLConnect::writeWorksheet(wb, matrix(unlist(vaRinfo$clinVar$clinVar.info$same_codon[[j]][9]),nrow=5), sheet=sheet.name,
                                startRow=12, startCol=2,
                                header=FALSE, rownames= FALSE)
    }
  }
  XLConnect::removeSheet(wb, "ClinVar Extra")


  #············Predictors----
  XLConnect::writeWorksheet(wb, vaRinfo$predictors$predictor.table$predictors.table2, sheet="Predictors",
                            startRow=8, startCol=1,
                            header=FALSE, rownames= FALSE)

  XLConnect::writeWorksheet(wb, paste(vaRclass$final.criteria$PP3.message, vaRclass$final.criteria$BP4.message), sheet="Predictors",
                            startRow=2, startCol=5,
                            header=FALSE, rownames= FALSE)
  #···········NMD----
  if(!is.na(vaRinfo$codon.stop$variant.exon[1])){
    vari <- c(vaRinfo$codon.stop$variant.exon[,1], purrr::map(stringr::str_split(vaRinfo$Variant.Info$genomic, "g\\.|[A-z]"),5), purrr::map(stringr::str_split(vaRinfo$Variant.Info$variant, "c\\.|[A-z]"),2), vaRinfo$codon.stop$variant.exon[,2]) %>% unlist
    XLConnect::writeWorksheet(wb, t(vari), sheet="NMD",
                              startRow=7, startCol=1,
                              header=FALSE, rownames= FALSE)
    if(!is.na(vaRinfo$codon.stop$premature.ter.codon[1])&&nrow(vaRinfo$codon.stop$premature.ter.codon)!=0){
      XLConnect::writeWorksheet(wb, vaRinfo$codon.stop$premature.ter.codon[,c(1,3:4)], sheet="NMD",
                                startRow=7, startCol=5,
                                header=FALSE, rownames= FALSE)
    }
    XLConnect::writeWorksheet(wb,  nrow(vaRinfo$codon.stop$exons), sheet="NMD",
                              startRow=7, startCol=8,
                              header=FALSE, rownames= FALSE)

    XLConnect::writeWorksheet(wb,  vaRclass$final.criteria$NMD, sheet="NMD",
                              startRow=7, startCol=9,
                              header=FALSE, rownames= FALSE)


  }
  XLConnect::writeWorksheet(wb, vaRinfo$codon.stop$exons[,2:8], sheet="NMD",
                            startRow=21, startCol=8,
                            header=FALSE, rownames= FALSE)

  XLConnect::writeWorksheet(wb, vaRinfo$codon.stop$canonical.skip.pred %>% as.data.frame(), sheet="NMD",
                            startRow=7, startCol=10,
                            header=FALSE, rownames= FALSE)

  XLConnect::writeWorksheet(wb, vaRclass$final.criteria$PVS1, sheet="NMD",
                            startRow=1, startCol=6,
                            header=FALSE, rownames= FALSE)

if(nrow(vaRinfo$predictors$predictor.table$spliceai.10k)>0){
  XLConnect::writeWorksheet(wb, vaRinfo$predictors$predictor.table$spliceai.10k[1,c(1:5,10:19)], sheet="NMD",
                            startRow=12, startCol=1,
                            header=FALSE, rownames= FALSE)
  XLConnect::writeWorksheet(wb, vaRinfo$predictors$predictor.table$spliceai.10k[1,20:32], sheet="NMD",
                            startRow=14, startCol=1,
                            header=FALSE, rownames= FALSE)
  XLConnect::writeWorksheet(wb, vaRinfo$predictors$predictor.table$spliceai.10k[1,33:47], sheet="NMD",
                            startRow=16, startCol=1,
                            header=FALSE, rownames= FALSE)
}else{
  XLConnect::writeWorksheet(wb, "SpliceAI-10k not calculated.", sheet="NMD",
                            startRow=12, startCol=1,
                            header=FALSE, rownames= FALSE)
}


  #············Start Codon ----
  XLConnect::writeWorksheet(wb, vaRinfo$start.lost.variants$pos.second.met, sheet="Start Codon",
                            startRow=21, startCol=2,
                            header=FALSE, rownames= FALSE)

  XLConnect::writeWorksheet(wb, vaRinfo$start.lost.variants[[3]], sheet="Start Codon",
                            startRow=23, startCol=1,
                            header=FALSE, rownames= FALSE)

  #············Bibliography data vaRHC----

  #Parsons
  if(nrow(vaRinfo$functional.assays$parsons)==0){
    XLConnect::writeWorksheet(wb, "This variant is not listed in Parsons et al.",
                              sheet="Bibliography data vaRHC", startRow=8, startCol=1, header=FALSE)
    if (!(gene %in% c("BRCA1", "BRCA2"))) XLConnect::setRowHeight(wb, sheet="Bibliography data vaRHC", row=c(5:9), height = 0 )
  }else{
    XLConnect::writeWorksheet(wb,  vaRinfo$functional.assays$parsons[1,2:ncol(vaRinfo$functional.assays$parsons)],
                              sheet="Bibliography data vaRHC", startRow=8, startCol=1, header=FALSE)
  }

  #lyra
  if(is.na(vaRinfo$functional.assays$lyra[1])|| nrow(vaRinfo$functional.assays$lyra)==0){
    XLConnect::writeWorksheet(wb, "This variant is not listed in Lyra et al's Supplementary Table 9",
                              sheet="Bibliography data vaRHC", startRow=16, startCol=1, header=FALSE)
    if (!(gene %in% c("BRCA1"))) XLConnect::setRowHeight(wb, sheet="Bibliography data vaRHC", row=c(11:17), height = 0 )
  } else{
    XLConnect::writeWorksheet(wb, vaRinfo$functional.assays$lyra,
                              sheet="Bibliography data vaRHC", startRow=16, startCol=1, header=FALSE)
  }

  #Adamovich_Fayer

  if(is.na(vaRinfo$functional.assays$Adamovich_Fayer[1])|| nrow(vaRinfo$functional.assays$Adamovich_Fayer)==0){
    XLConnect::writeWorksheet(wb, "This variant is not listed in Fayer et al. (2021), nor in Adamovich et al. (2022)",
                              sheet="Bibliography data vaRHC", startRow=22, startCol=1, header=FALSE)
    if (!(gene %in% c("BRCA1"))) XLConnect::setRowHeight(wb, sheet="Bibliography data vaRHC", row=c(18:22), height = 0 )
  } else{
    XLConnect::writeWorksheet(wb, vaRinfo$functional.assays$Adamovich_Fayer[,2:12],
                              sheet="Bibliography data vaRHC", startRow=22, startCol=1, header=FALSE)
  }


  #jia
  if(is.na(vaRinfo$functional.assays$jia[1])||nrow(vaRinfo$functional.assays$jia)==0){
    XLConnect::writeWorksheet(wb,"This variant is not listed in Jia et al.",
                              sheet="Bibliography data vaRHC", startRow=27, startCol=1, header=FALSE)
    if (!(gene %in% c("MSH2"))) XLConnect::setRowHeight(wb, sheet="Bibliography data vaRHC", row=c(25:27), height = 0 )
  }else{
    XLConnect::writeWorksheet(wb, vaRinfo$functional.assays$jia[1, 2:ncol(vaRinfo$functional.assays$jia)],
                              sheet="Bibliography data vaRHC", startRow=27, startCol=1, header=FALSE)
  }

  #cimra
  if(is.na(vaRinfo$functional.assays$cimra[1])||nrow(vaRinfo$functional.assays$cimra)==0){
    XLConnect::writeWorksheet(wb, "This variant is not listed in Cimra et al.",
                              sheet="Bibliography data vaRHC", startRow=32, startCol=1, header=FALSE)
    if (!(gene %in% c("MLH1", "MSH2", "MSH6", "PMS2"))) XLConnect::setRowHeight(wb, sheet="Bibliography data vaRHC", row=c(30:32), height = 0 )

  }else{
    XLConnect::writeWorksheet(wb, vaRinfo$functional.assays$cimra[,2:19],
                              sheet="Bibliography data vaRHC", startRow=32, startCol=1, header=FALSE)
  }



  if((all(is.na(vaRinfo$insight.info$classification)) && all(is.na(vaRinfo$insight.info$multifactorial))) | (length(vaRinfo$insight.info$classification) == 0 & length(vaRinfo$insight.info$multifactorial) == 0 )){
    XLConnect::writeWorksheet(wb, "variant not found in InSight", sheet="Bibliography data vaRHC",
                              startRow=37, startCol=1,
                              header=TRUE)
    if (!(gene %in% c("APC", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "MUTYH", "CDH1", "GALNT12"))) XLConnect::setRowHeight(wb, sheet="Bibliography data vaRHC", row=c(35:43), height = 0 )
    #if (gene %in% c("APC", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "MUTYH", "CDH1", "GALNT12")) setRowHeight(wb, sheet="Bibliography data vaRHC", row=c(30:37), height = 0 )

  }else if (vaRinfo$insight.info$classification[1] == "Not working"){
    XLConnect::writeWorksheet(wb, "Please check manually, RSelenium is not working.", sheet="Bibliography data vaRHC",
                              startRow=37, startCol=1,
                              header=TRUE)
  }else if (length(vaRinfo$insight.info[[1]])>2){
    XLConnect::writeWorksheet(wb, vaRinfo$insight.info[[1]][1], sheet="Bibliography data vaRHC",
                              startRow=36, startCol=1,
                              header=TRUE)
    XLConnect::writeWorksheet(wb, vaRinfo$insight.info[[1]][2], sheet="Bibliography data vaRHC",
                              startRow=38, startCol=1,
                              header=TRUE)
    XLConnect::writeWorksheet(wb, vaRinfo$insight.info[[1]][3], sheet="Bibliography data vaRHC",
                              startRow=40, startCol=1,
                              header=TRUE)
    if (length(vaRinfo$insight.info[[2]])>0){
      XLConnect::writeWorksheet(wb, vaRinfo$insight.info[[2]], sheet="Bibliography data vaRHC",
                                startRow=42, startCol=1,
                                header=TRUE) }
  } else {
    XLConnect::writeWorksheet(wb, "variant not found in InSight", sheet="Bibliography data vaRHC",
                              startRow=37, startCol=1,
                              header=TRUE)
    #if (!(gene %in% c("APC", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "MUTYH", "CDH1", "GALNT12"))) setRowHeight(wb, sheet="Bibliography data vaRHC", row=c(27:37), height = 0 )
    #if (gene %in% c("APC", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "MUTYH", "CDH1", "GALNT12")) setRowHeight(wb, sheet="Bibliography data vaRHC", row=c(30:37), height = 0 )

  }

  #cancer_hotspots
  if (is.na(vaRinfo$cancer.hotspots$variant)[1] || nrow(vaRinfo$cancer.hotspots$variant)==0){
    XLConnect::writeWorksheet(wb,  "Not found in the db",
                              sheet="Bibliography data vaRHC", startRow=48, startCol=1, header=FALSE)
    XLConnect::setRowHeight(wb, sheet="Bibliography data vaRHC", row=c(46:50), height = 0 )

  }else{
    XLConnect::writeWorksheet(wb, unlist(vaRinfo$cancer.hotspots$variant),
                              sheet="Bibliography data vaRHC", startRow=48, startCol=1, header=FALSE)
    XLConnect::writeWorksheet(wb, vaRinfo$cancer.hotspots$cancer.type,
                              sheet="Bibliography data vaRHC", startRow=49, startCol=2, header=TRUE)
  }



  #tp53 functionals
  if(is.na(vaRinfo$functional.assays$tp53.functionals[1])){
    XLConnect::writeWorksheet(wb,  "Not found the variant in these functional studies",
                              sheet="Bibliography data vaRHC", startRow=56, startCol=1, header=FALSE)
    if (gene !="TP53") XLConnect::setRowHeight(wb, sheet="Bibliography data vaRHC", row=c(52:56), height = 0 )

  }else{
    XLConnect::writeWorksheet(wb, vaRinfo$functional.assays$tp53.functionals,
                              sheet="Bibliography data vaRHC", startRow=55, startCol=1, header=TRUE)
  }
  #chek2 functionals
  if (is.na(vaRinfo$functional.assays$chek2.functionals[1]) || nrow(vaRinfo$functional.assays$chek2.functionals)==0){
    XLConnect::writeWorksheet(wb, "Not found in the db",
                              sheet="Bibliography data vaRHC", startRow=62, startCol=1, header=FALSE)
  } else{
    XLConnect::writeWorksheet(wb, vaRinfo$functional.assays$chek2.functionals[,2:24],
                              sheet="Bibliography data vaRHC", startRow=62, startCol=1, header=FALSE)
  }


  #············Functional information----
  XLConnect::writeWorksheet(wb, paste(vaRclass$final.criteria$PS3.message,  vaRclass$final.criteria$BS3.message),
                            sheet="Functional information", startRow=4, startCol=4, header=TRUE)
  #CancerHotspots
  XLConnect::writeWorksheet(wb, vaRinfo$cancer.hotspots$variant, sheet="Bibliography data vaRHC", #non_cancer exomes
                            startRow=66, startCol=4,
                            header=FALSE, rownames= TRUE)
  XLConnect::writeWorksheet(wb, vaRinfo$cancer.hotspots$cancer.type, sheet="Bibliography data vaRHC", #non_cancer exomes
                            startRow=68, startCol=4,
                            header=TRUE, rownames= TRUE)

  #google references
  XLConnect::writeWorksheet(wb, vaRinfo$google.scholar.30.references$google.scholar.search,
                            sheet="Bibliography data vaRHC", startRow=71, startCol=3, header=FALSE)
  if(!is.na(vaRinfo$google.scholar.30.references$articles)[1]){
    if (stringr::str_detect(vaRinfo$google.scholar.30.references$articles[1],"Error 429:" )&& is.character(vaRinfo$google.scholar.30.references$articles)){
      XLConnect::writeWorksheet(wb, vaRinfo$google.scholar.30.references$articles,
                                sheet="Bibliography data vaRHC", startRow=73, startCol=1, header=FALSE)
    }else{
    XLConnect::writeWorksheet(wb, vaRinfo$google.scholar.30.references$articles[,1],
                              sheet="Bibliography data vaRHC", startRow=73, startCol=1, header=FALSE)
    XLConnect::writeWorksheet(wb, vaRinfo$google.scholar.30.references$articles[,2],
                              sheet="Bibliography data vaRHC", startRow=73, startCol=6, header=FALSE)
    XLConnect::writeWorksheet(wb, vaRinfo$google.scholar.30.references$articles[,3],
                              sheet="Bibliography data vaRHC", startRow=73, startCol=9, header=FALSE)
    XLConnect::writeWorksheet(wb, vaRinfo$google.scholar.30.references$articles[,4],
                              sheet="Bibliography data vaRHC", startRow=73, startCol=12, header=FALSE)

    }
  }

  if(gene=="CDH1"){
    XLConnect::setColumnWidth(wb,"Evidence",c(4,5:12),0)
  } else if (gene=="PTEN"){
    XLConnect::setColumnWidth(wb,"Evidence",c(3,4,6:12),0)
  } else if (gene=="ATM"){
    XLConnect::setColumnWidth(wb,"Evidence",c(3:5,7:12),0)
  } else if (gene=="CHEK2"){
    XLConnect::setColumnWidth(wb,"Evidence",c(3:6,8:12),0)
  } else if (gene=="TP53"){
    XLConnect::setColumnWidth(wb,"Evidence",c(3:7, 9:12),0)
  } else if (gene=="BRCA1"|gene=="BRCA2"){
    XLConnect::setColumnWidth(wb,"Evidence",c(3:8, 10:11),0)
  }else if (gene=="MLH1"|gene=="MSH2"|gene=="MSH6"|gene=="PMS2"){
    XLConnect::setColumnWidth(wb,"Evidence",c(3:9,11:12),0)
  }else if (gene=="PALB2"){
    XLConnect::setColumnWidth(wb,"Evidence",c(3:10,12),0)
  }else if (gene == "BRIP1"){
    XLConnect::setColumnWidth(wb,"Evidence",c(3:11),0)
  } else{
    XLConnect::setColumnWidth(wb,sheet="Evidence",column=c(4:12), width=0)
  }
  # Save workbook ----
  XLConnect::saveWorkbook(wb, )
}
}

################################################################################
## Introduction message
################################################################################
introMessage <- function(all.information){
  if(all.information$Variant.Info$most.severe.consequence == "synonymous_variant"){
    intro.message <- paste("The", all.information$Variant.Info$gene, all.information$Variant.Info$variant, "variant affects a not conserved nucleotide, resulting in no amino acid change.")
  }else if (all.information$Variant.Info$most.severe.consequence == "stop_gained"){ #nonesense
    intro.message <- paste0(all.information$Variant.Info$variant,
                            ", located in exon",
                            all.information$Variant.Info$exon_intron %>% as.character(),
                            "of the", all.information$Variant.Info$gene,
                            "gene, is a nonsense variant" )
  }else if (all.information$Variant.Info$most.severe.consequence == "frameshift_variant"){ #fer condicio que estigui a l'exó
    intro.message <- paste0(all.information$Variant.Info$variant,
                            ", located in exon ",
                            all.information$Variant.Info$exon_intron %>% as.character(),
                            " of the ",
                            all.information$Variant.Info$gene,
                            " gene, consistis in the deletion/insertion/duplication of",  nrExtract(all.information)$nr.nt ," nucleotide(s), causing a translational frameshift with a predicted alternate stop codon (",
                            all.information$Variant.Info$protein, ").")
  }else if (all.information$Variant.Info$most.severe.consequence %in% c("splice_donor_variant", "splice_acceptor_variant")){
    intro.message <- paste0(all.information$Variant.Info$variant,
                            ", located in a canonic splicing site of the ",
                            all.information$Variant.Info$gene,
                            " gene is predicted to alter splicing, probably causing the skiping of exon/s {skippedEXON(s)}. ")
  }else if(all.information$Variant.Info$most.severe.consequence == "missense_variant"){
    intro.message <- paste(all.information$Variant.Info$variant,
                           ", located in exon",
                           all.information$Variant.Info$exon_intron$exon,
                           "of the",
                           all.information$Variant.Info$gene,
                            "gene, is a missense variant predicted to result in a substitution",
                           all.information$Variant.Info$protein )
  }else{
    intro.message <- NA
  }
  return(intro.message)
}

showMessage <- function (all.information, crit.met, final.criteria){
  intro <- introMessage(all.information)
    max.subpop <- all.information$gnomAD$info$exomes.genomes$non.cancer$subpopulations %>%
    dplyr::filter(CI==max(all.information$gnomAD$info$exomes.genomes$non.cancer$subpopulations$CI))
  popu <- max.subpop %>%
    dplyr::select("rowname") %>%
    dplyr::mutate(rowname = stringr::str_sub(.data$rowname, -3,-1)) %>%
    as.character()
  freq.sum <- paste0("Freq gnomAD all non-cancer v2.1.1 = ", round(as.numeric(all.information$gnomAD$info$exomes.genomes$non.cancer$overall$AF)*100,5), "% (MCAF ",  as.numeric(all.information$gene.specific.info$IC)*100,
                     "% = ", round(as.numeric(all.information$gnomAD$info$exomes.genomes$non.cancer$overall$AF)*100,5), "%). ",
                     ifelse(nrow(max.subpop)==8, "",
                            paste0("Max freq in ",  popu, " subpopulation ", round(as.numeric(max.subpop$AF)*100,5), "% (MCAF ", as.numeric(all.information$gene.specific.info$IC)*100, "% = ", round(as.numeric(max.subpop$CI)*100,5), "%).", collapse=" and ") ))

  all.messages <- paste(intro,
                        ifelse(any(stringr::str_detect(crit.met, "PVS1")),
                               final.criteria$PVS1.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "PM1")),
                               final.criteria$PM1.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "PM4")),
                               final.criteria$PM1.message,
                               ""),
                        freq.sum,
                        ifelse(any(stringr::str_detect(crit.met, "BA1")),
                               final.criteria$BA1.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "BS1")),
                               final.criteria$BS1.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "PM2")),
                               final.criteria$PM2.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "BS2")),
                               final.criteria$BS2.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "PP3")),
                               final.criteria$PP3.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "BP4")),
                               final.criteria$BP4.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "BP7")),
                               final.criteria$BP7.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "BS3")),
                               final.criteria$BS4.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "PS3")),
                               final.criteria$BS4.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "PS1")),
                               final.criteria$PS1.message,
                               ""),
                        ifelse(any(stringr::str_detect(crit.met, "PM5")),
                               final.criteria$PM5.message,
                               "")

                        )
}
