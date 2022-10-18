#' @importFrom dplyr %>%
#' @import Biostrings
#' @import XLconnect
#' @import RMySQL
NULL

#' vaR()
#' @description  vaR function firstly collects different type of information related to the variant of interest (e.g. Population frequencies, in-silico predictions, some functional studies...). 
#' Secondly,  it mixes all this information and following modified ACMG rules it returns the final classification of the variant and it specifies which criteria are given. Optionally, all the information can be printed in an excel file.
#' @param assembly character hg19 or hg38
#' @param gene gene of interest
#' @param variant variant of interest in cdna
#' @param NM Accession number of the transcrit and mRNA from RefSeq. By default is NULL and vaRHC will consider the ones detailed in README file. Be careful if you use a different NM because the program has not been validated for it. If you provide a different NM,  NC and CCDS must also be provided.
#' @param NC Accession number of the chromosome RefSeq. It can be ommited if the variant is exonic. By default is NULL and vaRHC will consider the ones detailed in README file. Be careful if you use a different NM because the program has not been validated for it. If you provide a different NC, NM and CCDS must also be provided.
#' @param CCDS Consensus CDS id https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi. By default is NULL and vaRHC will consider the ones detailed in README file. Be careful if you use a different CCDS because the program has not been validated for it. If you provide a different CCDS, NM and NC must also be provided. Current version only works for hg19.
#' @param gene.specific.df By default is NULL, it uses the default parameters described in README. If you would like to change some defaults or include another gene, a template can be downloaded from Github: https://github.com/emunte/Class_variants/tree/main/documents/gen_especific.csv or in the package docs folder and some parameters can be modified taking into account your preferences
#' @browser Which browser to start Rselenium server. By default is "firefox" (the recommended). If you do not have firefox installed try either "chrome" or "phantomjs".
#' @spliceai.program Logical. By default is FALSE and it is assumed that SpliceAI program is not installed in your computer. If this parameter is FALSE, the program will only classify substitutions and simple deletion variants taking into account a spliceAI distance of 1000 and will show masked results. If you want to classify other variants please install SpliceAI (https://pypi.org/project/spliceai/) and set to TRUE the parameter.
#' @spliceai.reference Path to the Reference genome hg19 fasta file. Can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz . By default is NULL and it will only be taken into account if spliceai.program is set to TRUE.
#' @spliceai.annotation Path to gene annotation file. By default it uses the file stored in docs folder: "../docs/gencode_spliceai_hg19.txt"
#' @spliceai.distance  Integer. Maximum distance between the variant and gained/lost splice site (default: 1000)
#' @spliceai.masked Mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss (default: 1)
#' @provean.program Logical. By default is FALSE and it is assumed that provean program is not installed in your computer. 
#' @param excel.results Logical. By default is FALSE and no excel file would be produced. If TRUE and excel file will be saved
#' @param path.original.file If excel.results param is set to TRUE, path.original.file must contain the path to the excel template. By default is the template located in the package docs folder.
#' @param path.copy.file By default is NULL. If excel.results param is set to TRUE, path.copy.file must provide the path where the excel file must be saved. If not provided it will be saved in the working directory.
#' @author Elisabet Munté Roca
#' @examples
#' vaR (assembly= "hg19", gene= "BRCA1", variant= "c.211A>G",  excel.results=TRUE, path.copy.file="./excels/")
#' vaR.all <- vaR (assembly = "hg19", gene = "BRCA1", variant = "c.692C>T", spliceai.program = TRUE, spliceai.reference = "../docs/hg19.fa", excel.results = TRUE, path.copy.file="../excels" )
#' @references 
#' Richards, S., Aziz, N., Bale, S., Bick, D., Das, S., Gastier-Foster, J., Grody, W. W., Hegde, M., Lyon, E., Spector, E., Voelkerding, K., Rehm, H. L., & ACMG Laboratory Quality Assurance Committee (2015). Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genetics in medicine : official journal of the American College of Medical Genetics, 17(5), 405–424. https://doi.org/10.1038/gim.2015.30
#' Tavtigian, S. V., Harrison, S. M., Boucher, K. M., & Biesecker, L. G. (2020). Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. Human mutation, 41(10), 1734–1737. https://doi.org/10.1002/humu.24088
#' Garrett, A., Durkie, M., Callaway, A., Burghel, G. J., Robinson, R., Drummond, J., Torr, B., Cubuk, C., Berry, I. R., Wallace, A. J., Ellard, S., Eccles, D. M., Tischkowitz, M., Hanson, H., Turnbull, C., & CanVIG-UK (2021). Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. Journal of medical genetics, 58(5), 297–304. https://doi.org/10.1136/jmedgenet-2020-107248
#' ClinGen CDH1 Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 3: https://clinicalgenome.org/site/assets/files/7118/clingen_cdh1_acmg_specifications_v3.pdf
#' ClinGen PTEN Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 2: https://clinicalgenome.org/site/assets/files/4000/clingen_pten_acmg_specifications_v2.pdf
#' ClinGen InSiGHT Hereditary Colorectal Cancer/Polyposis Variant Curation Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 1 (draft): https://www.insight-group.org/content/uploads/2021/11/DRAFT_Nov_2021_TEMPLATE_SVI.ACMG_Specifications_InSiGHT_MMR_V1.pdf
#' Feliubadaló, L., Moles-Fernández, A., Santamariña-Pena, M., Sánchez, A. T., López-Novo, A., Porras, L. M., Blanco, A., Capellá, G., de la Hoya, M., Molina, I. J., Osorio, A., Pineda, M., Rueda, D., de la Cruz, X., Diez, O., Ruiz-Ponte, C., Gutiérrez-Enríquez, S., Vega, A., & Lázaro, C. (2021). A Collaborative Effort to Define Classification Criteria for ATM Variants in Hereditary Cancer Patients. Clinical chemistry, 67(3), 518–533. https://doi.org/10.1093/clinchem/hvaa250
vaR <- function(assembly, gene, variant, NM=NULL, NC = NULL, CCDS=NULL, gene.specific.df=NULL, browser="firefox", spliceai.program = FALSE, spliceai.reference = NULL, spliceai.annotation = system.file("docs", "gencode_spliceai_hg19.txt", package="vaRHC"), spliceai.distance = 1000, spliceai.masked = 1, provean.program = FALSE, excel.results = FALSE,  path.copy.file = NULL ){
  cat("looking for VariantInfo, please wait\n")
  info <- vaRinfo(assembly,  gene, variant, NM = NM, NC = NC, CCDS = CCDS, gene.specific.df = gene.specific.df, browser= browser, spliceai.program = spliceai.program, spliceai.reference = spliceai.reference, spliceai.annotation= spliceai.annotation, spliceai.distance= spliceai.distance, spliceai.masked = spliceai.masked, provean.program = provean.program)
  cat("calculating the final classification , please wait \n")
  class <- vaRclass(info)
  if (isTRUE(excel.results)){
    cat("printing information in an excel file,  please wait \n")
    vaRreport(info, class, path.copy.file)
  }
  return(list(vaRinfo = info, 
              vaRclass = class))
}



#' vaRbatch()
#' @description  vaRbatch function allows to perform vaR function in batch. It also returns a logfile.
#' @param assembly character hg19 or hg38
#' @param gene gene of interest
#' @param variant variant of interest in cdna
#' @param NM Accession number of the transcrit and mRNA from RefSeq. By default is NULL and vaRHC will consider the ones detailed in README file. Be careful if you use a different NM because the program has not been validated for it. If you provide a different NM,  NC and CCDS must also be provided.
#' @param NC Accession number of the chromosome RefSeq. It can be ommited if the variant is exonic. By default is NULL and vaRHC will consider the ones detailed in README file. Be careful if you use a different NM because the program has not been validated for it. If you provide a different NC, NM and CCDS must also be provided.
#' @param CCDS Consensus CD id https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi. By default is NULL and vaRHC will consider the ones detailed in README file. Be careful if you use a different CCDS because the program has not been validated for it. If you provide a different CCDS, NM and NC must also be provided. Current version only works for hg19.
#' @param gene.specific.df By default is NULL, it uses the default parameters described in README. If you would like to change some defaults or include another gene, a template can be downloaded from Github: https://github.com/emunte/Class_variants/tree/main/documents/gen_especific.csv or in the package docs folder and some parameters can be modified taking into account your preferences
#' @browser Which browser to start Rselenium server. By default is "firefox" (the recommended). If you do not have firefox installed try either "chrome" or "phantomjs".
#' @spliceai.program Logical. By default is FALSE and it is assumed that SpliceAI program is not installed in your computer. If this parameter is FALSE, the program will only classify substitutions and simple deletion variants taking into account a spliceAI distance of 1000 and will show masked results. If you want to classify other variants please install SpliceAI (https://pypi.org/project/spliceai/) and set to TRUE the parameter.
#' @spliceai.reference Path to the Reference genome hg19 fasta file. Can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz . By default is NULL and it will only be taken into account if spliceai.program is set to TRUE.
#' @spliceai.annotation Path to gene annotation file. By default it uses the file stored in docs folder: "../docs/gencode_spliceai_hg19.txt"
#' @spliceai.distance  Integer. Maximum distance between the variant and gained/lost splice site (default: 1000)
#' @spliceai.masked Mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss (default: 1)
#' @provean.program Logical. By default is FALSE and it is assumed that provean program is not installed in your computer. 
#' @param excel.results Logical. By default is FALSE and no excel file would be produced. If TRUE and excel file will be saved
#' @param path.copy.file By default is NULL. If excel.results param is set to TRUE, path.copy.file must provide the path where the excel file must be saved. If not provided it will be saved in the working directory.
#' @author Elisabet Munté Roca
#' @examples
# all <- vaRbatch(assembly = "hg19", all.variants, spliceai.program = TRUE, spliceai.reference= "/home/emunte/hg19.fa", excel.results = TRUE, path.copy.file="/media/emunte/Nuevo/Variants_diagnostic/NS49/excels")
# all <- vaRbatch(all.variants, spliceai.program = FALSE, excel.results = TRUE, path.original.file="../docs/template.xlsx", path.copy.file="/media/emunte/Nuevo/Variants_diagnostic/NS49/excels" )
#' @export
vaRbatch <- function (assembly = "hg19", all.variants, gene.specific.df=NULL, browser="firefox", spliceai.program = FALSE, spliceai.reference = NULL, spliceai.annotation = system.file("docs", "gencode.v38lift37.annotation.txt", package="vaRHC"), spliceai.distance = 1000, spliceai.masked = 1, provean.program = FALSE, print.data.frame = TRUE, excel.results = FALSE, path.copy.file = NULL){
  time <-  Sys.time() %>% 
    stringr::str_replace_all("-|:| ", "_")
  log.file <- file.path(getwd(), "log")
  dir.create(log.file, showWarnings = FALSE)
  file.create(paste0(log.file,"/", time, ".log"))
  all.variants.list <- list()
  df.variants2 <- data.frame(NM=NULL, gene=NULL, variant=NULL, prot=NULL, final_class=NULL, all_criteria=NULL, 
                             PVS1=NULL, PS1=NULL,  PS3=NULL, 
                             PM1=NULL, PM2=NULL, PM4=NULL, PM5=NULL,  PP3=NULL,  
                             BA1=NULL, BS1=NULL, BS2=NULL, BS3=NULL, BP4=NULL,BP7=NULL, 
                             clinvar_class=NULL, clinvar_review=NULL, 
                             spliceAI_AG_score=NULL, spliceAI_AG_dis=NULL,spliceAI_AL_score=NULL, spliceAI_AL_dis=NULL,spliceAI_DG_score=NULL, spliceAI_DG_dis=NULL, spliceAI_DL_score=NULL, spliceAI_DL_dis=NULL, 
                             revel=NULL )
  for (i in 1:nrow(all.variants)){
    cat(paste(i, "/", nrow(all.variants), "variants \n"))
    gene <- all.variants$gene[i]
    variant <- all.variants$variant[i]
    NM <- NULL
    NC <- NULL
    CCDS <- NULL
    tryCatch({
      time.start <- Sys.time()
      info.R <- vaR(assembly = assembly, 
                    gene = gene, 
                    variant = variant, 
                    NM = NM, 
                    NC = NC, 
                    CCDS = CCDS,  
                    gene.specific.df = gene.specific.df,  
                    browser = browser, 
                    spliceai.program = spliceai.program, 
                    spliceai.reference = spliceai.reference, 
                    spliceai.annotation = spliceai.annotation, 
                    spliceai.distance = spliceai.distance, 
                    spliceai.masked = spliceai.masked, 
                    provean.program = provean.program,  
                    excel.results = excel.results, 
                    path.copy.file= path.copy.file)
      crit <- paste0(info.R$vaRclass$final.classification$criteria.assigned  , collapse=";")
      if(print.data.frame == TRUE){
        df.variants2 <- rbind(df.variants2, 
                              c(info.R$vaRinfo$Variant.Info$NM,
                                info.R$vaRinfo$Variant.Info$gene,
                                info.R$vaRinfo$Variant.Info$variant, 
                                info.R$vaRinfo$Variant.Info$protein,
                                info.R$vaRclass$final.classification$final.class,
                                crit, 
                                info.R$vaRclass$final.criteria$PVS1.message,
                                info.R$vaRclass$final.criteria$PS1.message, 
                                info.R$vaRclass$final.criteria$PS3.message, 
                                info.R$vaRclass$final.criteria$PM1.message,
                                info.R$vaRclass$final.criteria$PM2.message,
                                info.R$vaRclass$final.criteria$PM4.message,
                                info.R$vaRclass$final.criteria$PM5.message, 
                                info.R$vaRclass$final.criteria$PP3.message,
                                info.R$vaRclass$final.criteria$BA1.message,
                                info.R$vaRclass$final.criteria$BS1.message,
                                info.R$vaRclass$final.criteria$BS2.message,
                                info.R$vaRclass$final.criteria$BS3.message,
                                info.R$vaRclass$final.criteria$BP4.message,
                                info.R$vaRclass$final.criteria$BP7.message,
                                ifelse(length(purrr::map(info.R$vaRinfo$clinVar$clinVar.info$variant,7) %>% purrr::map(1) %>% unlist() %>% as.character())==0,NA, purrr::map(info.R$vaRinfo$clinVar$clinVar.info$variant,7) %>% purrr::map(1) %>% unlist() %>% as.character()),
                                ifelse(length(purrr::map(info.R$vaRinfo$clinVar$clinVar.info$variant,7) %>% purrr::map(2) %>% unlist() %>% as.character())==0,NA, purrr::map(info.R$vaRinfo$clinVar$clinVar.info$variant,7) %>% purrr::map(2) %>% unlist() %>% as.character()),
                                info.R$vaRinfo$predictors$predictor.table[14,"values"],
                                info.R$vaRinfo$predictors$predictor.table[14,"position"],
                                info.R$vaRinfo$predictors$predictor.table[15,"values"],
                                info.R$vaRinfo$predictors$predictor.table[15,"position"],
                                info.R$vaRinfo$predictors$predictor.table[16,"values"],
                                info.R$vaRinfo$predictors$predictor.table[16,"position"],
                                info.R$vaRinfo$predictors$predictor.table[17,"values"],
                                info.R$vaRinfo$predictors$predictor.table[17,"position"],
                                info.R$vaRinfo$predictors$predictor.table[4,"values"]
                              ))
        names(df.variants2) <- c("NM", "gene", "variant", "prot", "final_class", "all_criteria", 
                                 "PVS1", "PS1", "PS3",  
                                 "PM1", "PM2",  "PM4", "PM5",  "PP3", 
                                 "BA1", "BS1", "BS2", "BS3", "BP4", "BP7", 
                                 "clinvar_class", "clinvar_review", 
                                 "spliceAI_AG_score", "spliceAI_AG_dis","spliceAI_AL_score", "spliceAI_AL_dis", "spliceAI_DG_score", "spliceAI_DG_dis",  "spliceAI_DL_score", "spliceAI_DL_dis", 
                                 "REVEL")
        output <- file.path(getwd(), "output")
        dir.create(output, showWarnings = FALSE)
        write.table(x = df.variants2, 
                    file = file.path(output, paste0(time, "variants.txt")), 
                    row.names = F)
      }
      
      time.end <- Sys.time()
      time.diff <- time.end - time.start
      write.table(x = paste( gene, variant, "executed in ", time.diff), 
                  file = file.path(log.file, paste0(time, ".log")), 
                  append=TRUE, 
                  col.names = FALSE, 
                  row.names = FALSE)
      name <-  paste0(gene,"_", variant)
      all.variants.list[[name]] <-  info.R
      
    }, error = function(e){
      write.table(x = paste(gene, variant, "ERROR:", conditionMessage(e)), 
                  file = file.path(log.file, paste0(time, ".log")), 
                  append=TRUE, 
                  col.names = FALSE, 
                  row.names = FALSE)
    })
  }
  cat (paste0("Find log file: ",  file.path(log.file, paste0(time, ".log"))))
  return(info.R)
} 
