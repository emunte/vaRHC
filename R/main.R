#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @import Biostrings
#' @import RMySQL
#' @import stringr
#' @import R.utils
#' @import testthat
#' @import XVector
#' @importFrom assertthat assert_that
#' @importFrom tibble as_tibble
#' @importFrom tibble rownames_to_column
#' @importFrom tibble tibble
#' @importFrom httr content
#' @importFrom httr set_config
#' @importFrom httr GET
#' @importFrom DBI dbConnect
#' @importFrom RCurl getURL
#' @importFrom jsonlite fromJSON
#' @importFrom jsonlite toJSON
#' @importFrom purrr map
#' @importFrom xml2 read_html
#' @importFrom XML readHTMLTable
#' @importFrom utils URLencode data write.table

NULL

# #' @import rtracklayer
# #' @importFrom GenomicRanges GRanges
# #' @importFrom IRanges IRanges
# #' @importFrom seqinr write.fasta

#' @name gencode_spliceaihg19_
#' @title gencode_spliceai_hg19
#' @docType data
#' @references adapted from SpliceAI-lookup \url{https://github.com/broadinstitute/SpliceAI-lookup/blob/master/annotations/gencode.v38lift37.annotation.txt.gz}
#' @keywords data
NULL

#' @name gene_specific
#' @title gene_specific
#' @docType data
#' @author Elisabet Munté Roca \email{emunte@@idibell.cat}
#' @references \url{https://github.com/emunte/vaRHC/blob/main/data/gene_specific.txt}
NULL

#' @name ex_vaRbatch
#' @title example_input_vaRbatch
#' @docType data
#' @author Elisabet Munté Roca \email{emunte@@idibell.cat}
#' @references \url{https://github.com/emunte/vaRHC/blob/main/data/ex_vaRbatch.txt}
#' @keywords data
NULL

#' @name BLOSUM62
#' @title BLOSUM62
#' @docType data
#' @author H. Pags and P. Aboyoun
#' @references K. Malde, The effect of sequence quality on sequence alignment, Bioinformatics, Feb 23, 2008. \url{https://academic.oup.com/bioinformatics/article/24/7/897/296249}
NULL


#' vaR()
#' @description  vaR function firstly collects different type of information related to the variant of interest (e.g. Population frequencies, in-silico predictions, some functional studies...).
#' Secondly,  it mixes all this information and following modified ACMG rules it returns the final classification of the variant and it specifies which criteria are given. Optionally, all the information can be printed in an excel file.
#' @param gene gene of interest
#' @param variant variant of interest in coding DNA nomenclature
#' @param NM Accession number of the transcrit and mRNA from RefSeq. By default is NULL and vaRHC will consider the ones detailed in README file. Be careful if you use a different NM because the program has not been validated for it.
#' @param CCDS Consensus CDS id https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi. By default is NULL and vaRHC will consider the ones detailed in README file. Be careful if you use a different CCDS because the program has not been validated for it.
#' @param gene.specific.df By default is NULL, it uses the default parameters described in README. If you would like to change some defaults or include another gene, a template can be downloaded from Github: https://github.com/emunte/Class_variants/tree/main/documents/gen_especific.csv or in the package extdata folder and some parameters can be modified taking into account your preferences
#' @param remote Logical. Connect remotely to RSelenium server? By default is TRUE and will start Rselenium server.If it is FALSE vaRHC will not connect to insight database.
#' @param browser Which browser to start Rselenium server. By default is "firefox" (the recommended). If you do not have firefox installed try either "chrome" or "phantomjs".
#' @param spliceai.program Logical. By default is FALSE and it is assumed that SpliceAI program is not installed in your computer. If this parameter is FALSE, the program will only classify substitutions and simple deletion variants taking into account a spliceAI distance of 1000 and will show masked results. If you want to classify other variants please install SpliceAI (https://pypi.org/project/spliceai/) and set to TRUE the parameter.
#' @param spliceai.genome Fasta file assembly provided. It can only be "hg19" or "hg38" assembly. By default it will be hg19.
#' @param spliceai.reference Path to the Reference genome hg19 hg38 fasta file. hg19 file can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz . By default is NULL and it will only be taken into account if spliceai.program is set to TRUE.
#' @param spliceai.annotation Path to gene annotation file. By default is null and it uses UCSC ncbiRefSeq table. This page is sometimes temporarily out of service. Then the file must be provided (see and example in package docs folder: "../docs/gencode.v38lift37.annotation.txt")
#' @param spliceai.distance  Integer. Maximum distance between the variant and gained/lost splice site (default: 1000)
#' @param spliceai.masked Mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss (default: 1)
#' @param provean.program Logical. By default is FALSE and it is assumed that provean program is not installed in your computer.
#' @param provean.sh Path to provean.sh. It will only be considered if provean.program is set to TRUE.
#' @param spliceai.10k Logical. By default is FALSE and SpliceAI-10k will not be computed. It can only be computed if provean.program is TRUE.
#' @param excel.results Logical. By default is FALSE and no excel file would be produced. If TRUE and excel file will be saved
#' @param output.dir  output.dir must provide the folder to store the results. By default is NULL.
#' @param google.search Logical. By default is FALSE and it will not look for variant results in google.
#' @param verbose logic. By default is FALSE. If TRUE, it includes status messages (if any).
#' @return A list of two items: vaRinfo and vaRclass.  The first one stores variant information collected from diverse databases and the second one it is related to the criteria assigned to the variant and its final classification.
#' @author Elisabet Munté Roca
#' @examples
#' \dontrun{
#' output.dir <- "/path/to/outputdir"
#' vaR (gene= "BRCA1", variant= "c.211A>G", output.dir = output.dir,  excel.results=FALSE)
#' }
#' @references
#' Richards, S., Aziz, N., Bale, S., Bick, D., Das, S., Gastier-Foster, J., Grody, W. W., Hegde, M., Lyon, E., Spector, E., Voelkerding, K., Rehm, H. L., & ACMG Laboratory Quality Assurance Committee (2015). Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genetics in medicine : official journal of the American College of Medical Genetics, 17(5), 405–424. https://doi.org/10.1038/gim.2015.30
#' Tavtigian, S. V., Harrison, S. M., Boucher, K. M., & Biesecker, L. G. (2020). Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. Human mutation, 41(10), 1734–1737. https://doi.org/10.1002/humu.24088
#' Garrett, A., Durkie, M., Callaway, A., Burghel, G. J., Robinson, R., Drummond, J., Torr, B., Cubuk, C., Berry, I. R., Wallace, A. J., Ellard, S., Eccles, D. M., Tischkowitz, M., Hanson, H., Turnbull, C., & CanVIG-UK (2021). Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. Journal of medical genetics, 58(5), 297–304. https://doi.org/10.1136/jmedgenet-2020-107248
#' ClinGen CDH1 Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 3: https://clinicalgenome.org/site/assets/files/7118/clingen_cdh1_acmg_specifications_v3.pdf
#' ClinGen PTEN Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 2: https://clinicalgenome.org/site/assets/files/4000/clingen_pten_acmg_specifications_v2.pdf
#' ClinGen InSiGHT Hereditary Colorectal Cancer/Polyposis Variant Curation Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 1 (draft): https://www.insight-group.org/content/uploads/2021/11/DRAFT_Nov_2021_TEMPLATE_SVI.ACMG_Specifications_InSiGHT_MMR_V1.pdf
#' Feliubadaló, L., Moles-Fernández, A., Santamariña-Pena, M., Sánchez, A. T., López-Novo, A., Porras, L. M., Blanco, A., Capellá, G., de la Hoya, M., Molina, I. J., Osorio, A., Pineda, M., Rueda, D., de la Cruz, X., Diez, O., Ruiz-Ponte, C., Gutiérrez-Enríquez, S., Vega, A., & Lázaro, C. (2021). A Collaborative Effort to Define Classification Criteria for ATM Variants in Hereditary Cancer Patients. Clinical chemistry, 67(3), 518–533. https://doi.org/10.1093/clinchem/hvaa250
#' @export

vaR <- function(gene, variant, NM=NULL, CCDS=NULL, gene.specific.df=NULL, remote=TRUE, browser="firefox", spliceai.program = FALSE, spliceai.genome="hg19", spliceai.reference = NULL, spliceai.annotation = NULL, spliceai.distance = 1000, spliceai.masked = 1, provean.program = FALSE, provean.sh = NULL, spliceai.10k= FALSE, excel.results = FALSE,  output.dir = NULL, google.search = FALSE, verbose = FALSE ){
  if(!is.null(output.dir))assertthat::assert_that(dir.exists(output.dir), msg = "Output directory does not exists, please enter a valid one.")
  cat("looking for VariantInfo, please wait\n")
  info <- vaRinfo(gene = gene,
                  variant = variant,
                  NM = NM,
                  CCDS = CCDS,
                  gene.specific.df = gene.specific.df,
                  output.dir = output.dir,
                  remote = remote,
                  browser = browser,
                  spliceai.program = spliceai.program,
                  spliceai.genome = spliceai.genome,
                  spliceai.reference = spliceai.reference,
                  spliceai.annotation = spliceai.annotation,
                  spliceai.distance = spliceai.distance,
                  spliceai.masked = spliceai.masked,
                  provean.program = provean.program,
                  provean.sh = provean.sh,
                  spliceai.10k = spliceai.10k,
                  google.search = google.search,
                  verbose = verbose)
  cat("calculating the final classification , please wait \n")
  class <- vaRclass(info)
  if (isTRUE(excel.results)){
    cat("printing information in an excel file,  please wait \n")
    vaRreport(info, class, output.dir)
  }
  return(list(vaRinfo = info,
              vaRclass = class))
}



#' vaRbatch()
#' @description  vaRbatch function allows to perform vaR function in batch. It also returns a logfile.
#' @param all.variants it requires a dataframe object or the path to a file vcf. For the dataframe object variants must be coding dna sequence and it may have two columns gene and variant (see vignette to know ho to prepare the dataframe)
#' @param assembly Assembly used for the vcf generation and the one that want to be used for SpliceAI calculation. It can only be "hg19" or "hg38" assembly. If is NULL it will be calculated in hg19.
#' @param annotation  character. Only needed if a vcf file is provided. It can only be LRG or MANE_select. Variants will be annotated considering LRG or MANE_select transcripts.
#' @param gene.specific.df By default is NULL, it uses the default parameters described in README. If you would like to change some defaults or include another gene, a template can be downloaded from Github: https://github.com/emunte/Class_variants/tree/main/documents/gen_especific.csv or in the package docs folder and some parameters can be modified taking into account your preferences
#' @param remote Logical. Connect remotely to RSelenium server? By default is FALSE and it will not start Rselenium server.If it is TRUE vaRHC will connect to insight database.
#' @param browser Which browser to start Rselenium server. By default is "firefox" (the recommended). If you do not have firefox installed try either "chrome" or "phantomjs".
#' @param spliceai.program Logical. By default is FALSE and it is assumed that SpliceAI program is not installed in your computer. If this parameter is FALSE, the program will only classify substitutions and simple deletion variants taking into account a spliceAI distance of 1000 and will show masked results. If you want to classify other variants please install SpliceAI (https://pypi.org/project/spliceai/) and set to TRUE the parameter.
#' @param spliceai.genome Fasta file assembly provided. It can only be "hg19" or "hg38" assembly. By default it will be hg19.
#' @param spliceai.reference Path to the Reference genome hg19 or hg38 fasta file. hg19 fasta file can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz . By default is NULL and it will only be taken into account if spliceai.program is set to TRUE.
#' @param spliceai.annotation Path to gene annotation file. By default is null and it uses UCSC ncbiRefSeq table. This page is sometimes temporarily out of service. Then the file must be provided (see and example in package docs folder: "../docs/gencode.v38lift37.annotation.txt")
#' @param spliceai.distance  Integer. Maximum distance between the variant and gained/lost splice site (default: 1000)
#' @param spliceai.masked Mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss (default: 1)
#' @param provean.program Logical. By default is FALSE and it is assumed that provean program is not installed in your computer.
#' @param provean.sh Path to provean.sh. It will only be considered if provean.program is set to TRUE.
#' @param spliceai.10k Logical. By default is FALSE and SpliceAI-10k will not be computed. It can only be computed if provean.program is TRUE
#' @param print.data.frame  Logical. By defaul is TRUE and the results will be stored in a txt file.
#' @param excel.results Logical. By default is FALSE and no excel file would be produced. If TRUE and excel file will be saved
#' @param google.search Logical. By default is FALSE and it will not look for variant results in google.
#' @param verbose logic. By default is FALSE. If TRUE, it includes status messages (if any).
#' @param output.dir output.dir must provide the folder to store the results. By default is NULL.
#' @return A list of n intems. Where n is the number of variants where vaR() was successfully.
#' @author Elisabet Munté Roca
#' @examples
#' \dontrun{
#' data("ex_vaRbatch")
#' output.dir <- "/path/to/output_dir"
#' ex_vaRbatch[] <- lapply(ex_vaRbatch, as.character) #convert to character
#' all <- vaRbatch( all.variants = ex_vaRbatch, spliceai.program = FALSE, output.dir = output.dir)
#' }
#' @export
vaRbatch <- function (all.variants, assembly = NULL, annotation = NULL, gene.specific.df=NULL, remote = FALSE, browser="firefox", spliceai.program = FALSE, spliceai.genome= "hg19", spliceai.reference = NULL, spliceai.annotation = NULL, spliceai.distance = 4999, spliceai.masked = 1, provean.program = FALSE, provean.sh = NULL, spliceai.10k= FALSE, print.data.frame = TRUE, excel.results = FALSE, output.dir = NULL, google.search = FALSE, verbose = FALSE){
  time <-  Sys.time() %>%
    stringr::str_replace_all("-|:| ", "_")
  log.folder<- checkDir (output.dir, "log")
  file.create(file.path(log.folder,paste0(time, ".log")))
  if(is.character(all.variants)){
  all.variants <- vcfToCodingDNA (file.vcf = all.variants , assembly = assembly , log.dir = output.dir )
  }
  all.variants.list <- list()
  df.variants2 <- data.frame(NM=NULL, gene=NULL, variant=NULL, prot=NULL, final_class=NULL, all_criteria=NULL,
                             PVS1=NULL, PS1=NULL,  PS3=NULL,
                             PM1=NULL, PM2=NULL, PM4=NULL, PM5=NULL,  PP3=NULL,
                             BA1=NULL, BS1=NULL, BS2=NULL, BS3=NULL, BP4=NULL,BP7=NULL,
                             clinvar_class=NULL, clinvar_review=NULL,
                             spliceAI_AG_score=NULL, spliceAI_AG_dis=NULL,spliceAI_AL_score=NULL, spliceAI_AL_dis=NULL,spliceAI_DG_score=NULL, spliceAI_DG_dis=NULL, spliceAI_DL_score=NULL, spliceAI_DL_dis=NULL,
                             revel=NULL )
  info.R <- list()

  for (i in 1:nrow(all.variants)){
    cat(paste(i, "/", nrow(all.variants), "variants \n"))
    gene <- all.variants$gene[i]
    variant <- all.variants$variant[i]
    if(!is.null(all.variants$NM[i])){
      NM <- all.variants$NM[i]
    }else{
      NM <- NULL
    }

    NC <- NULL
    CCDS <- NULL
    tryCatch({
      time.start <- Sys.time()
      info.R <- vaR(gene = gene,
                    variant = variant,
                    NM = NM,
                    CCDS = CCDS,
                    gene.specific.df = gene.specific.df,
                    remote = remote,
                    browser = browser,
                    spliceai.program = spliceai.program,
                    spliceai.genome = spliceai.genome,
                    spliceai.reference = spliceai.reference,
                    spliceai.annotation = spliceai.annotation,
                    spliceai.distance = spliceai.distance,
                    spliceai.masked = spliceai.masked,
                    provean.program = provean.program,
                    provean.sh = provean.sh,
                    spliceai.10k= spliceai.10k,
                    excel.results = excel.results,
                    google.search = google.search,
                    output.dir= output.dir,
                    verbose = verbose)
      crit <- paste0(info.R$vaRclass$final.classification$criteria.assigned  , collapse=";")
      if(print.data.frame == TRUE){
        var.pred <- ifelse(is.na(info.R$vaRinfo$codon.stop)[5],
               NA,
               info.R$vaRinfo$codon.stop$canonical.skip.pred$variant)

        prot.pred <- ifelse(is.na(info.R$vaRinfo$codon.stop)[5],
                            NA,
                            info.R$vaRinfo$codon.stop$canonical.skip.pred$protein)
        if (length(info.R$vaRinfo$predictors$predictor.table$spliceai.10k)==0){
          splice.10k <- rep(NA,47 ) %>% as.data.frame() %>% t()
        }else{
          splice.10k <- info.R$vaRinfo$predictors$predictor.table$spliceai.10k
        }


        df.var <- data.frame(info.R$vaRinfo$Variant.Info$NM,
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
                   info.R$vaRinfo$predictors$predictor.table$predictors.table2[14,"values"],
                   info.R$vaRinfo$predictors$predictor.table$predictors.table2[14,"position"],
                   info.R$vaRinfo$predictors$predictor.table$predictors.table2[15,"values"],
                   info.R$vaRinfo$predictors$predictor.table$predictors.table2[15,"position"],
                   info.R$vaRinfo$predictors$predictor.table$predictors.table2[16,"values"],
                   info.R$vaRinfo$predictors$predictor.table$predictors.table2[16,"position"],
                   info.R$vaRinfo$predictors$predictor.table$predictors.table2[17,"values"],
                   info.R$vaRinfo$predictors$predictor.table$predictors.table2[17,"position"],
                   info.R$vaRinfo$predictors$predictor.table$predictors.table2[4,"values"],
                   var.pred,
                   prot.pred,
                   splice.10k

        )
        names(df.var) <- c("NM", "gene", "variant", "prot", "final_class", "all_criteria",
                           "PVS1", "PS1", "PS3",
                           "PM1", "PM2",  "PM4", "PM5",  "PP3",
                           "BA1", "BS1", "BS2", "BS3", "BP4", "BP7",
                           "clinvar_class", "clinvar_review",
                           "spliceAI_AG_score", "spliceAI_AG_dis","spliceAI_AL_score", "spliceAI_AL_dis", "spliceAI_DG_score", "spliceAI_DG_dis",  "spliceAI_DL_score", "spliceAI_DL_dis",
                           "REVEL", "variant_pred", "prot.pred",
                           "CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",
                           "ALLELE",	"SYMBOL",	"DS_AG",	"DS_AL",	"DS_DG",
                           "DS_DL",	"DP_AG",	"DP_AL",	"DP_DG",	"DP_DL",
                           "Used_RefSeq_Transcript",	"strand",	"Cryptic_Acceptor_activation",
                           "Cryptic_Donor_activation",	"Any_splicing_aberration",	"bp_5prime",
                           "bp_3prime",	"Partial_intron_retention",	"Partial_exon_deletion",
                           "Partial_exon_start",	"Partial_exon_end",	"Partial_frameshift",
                           "Partial_intron_retention_aaseq",	"Partial_exon_deletion_aaseq",	"Gained_exon_size",
                           "Pseudoexon_activation",	"Pseudoexon_start",	"Pseudoexon_end",
                           "Pseudoexon_frameshift",	"Pseudoexon_intron",	"Pseudoexon_activation_aaseq",
                           "Exon_skipping",	"Lost_exons",	"Exon_skipping_frameshift",
                           "Exon_skipping_aaseq",	"Retained_intron_size",	"Intron_retention",	"Retained_intron",
                           "Intron_retention_frameshift",	"Intron_retention_aaseq"
)

        df.variants2 <- rbind(df.variants2,
                              df.var)

        dataframe.folder<- checkDir (output.dir, "dataframe")
        write.table(x = df.variants2,
                    file = file.path(dataframe.folder, paste0(time, "variants.txt")),
                    row.names = F, sep="\t")
      }

      time.end <- Sys.time()
      time.diff <- time.end - time.start
      write.table(x = paste( gene, variant, "executed in ", time.diff),
                  file = file.path(log.folder, paste0(time, ".log")),
                  append=TRUE,
                  col.names = FALSE,
                  row.names = FALSE)
      name <-  paste0(gene,"_", variant)
      all.variants.list[[name]] <-  info.R

    }, error = function(e){
      write.table(x = paste(gene, variant, "ERROR:", conditionMessage(e)),
                  file = file.path(log.folder, paste0(time, ".log")),
                  append=TRUE,
                  col.names = FALSE,
                  row.names = FALSE)
    })
  }
  message (paste0("Find log file: ",  file.path(log.folder, paste0(time, ".log"))))
  return(all.variants.list)
}



