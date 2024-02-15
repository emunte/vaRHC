#' @import assertthat
#' @import DBI
#' @import RMySQL


################################################################################
## Connection to IDIBELL DB
################################################################################
#' @noRd
connectionDB <- function(query){
  assertthat::assert_that(is.character(query)|is.list(query), msg="please enter a character vector or list")
  con <- DBI::dbConnect(RMySQL::MySQL(),
                        user='userguest',
                        dbname='class_variants',
                        port=3306,
                        host='idivarhcdb001-cluster.cluster-ro-c38m6cy2udz0.eu-central-1.rds.amazonaws.com',
                        password='jNU%cd%Xjw*tY*%')
  on.exit(DBI::dbDisconnect(con))
  results <- lapply(query, function(x){
    return(tryCatch( suppressWarnings(DBI::dbGetQuery(con, x)),
                     error=function(e) NULL))
  })
  return(results)
}


################################################################################
## APIS
################################################################################

api <- function(server, ext){
  r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
  result <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))
  return(result)
}

api2 <- function(server, ext){
  new_config <- httr::config(ssl_verifypeer = FALSE)
  httr::set_config(new_config, override = FALSE)
  result <- jsonlite::fromJSON(httr::content(httr::GET(paste0(server, ext)), "text", encoding="UTF-8"))
  #result <- jsonlite::fromJSON(paste0(server, ext))
  return(result)
}

# api2 <- function(server, ext){
#   result <- jsonlite::fromJSON(paste0(server, ext))
#   return(result)
# }



################################################################################
## Ensembl
################################################################################
ensemblTranscript <- function(NM, gene){
  #we need to know the ensembl ID
  ext.ensembl.id <- paste0("/xrefs/symbol/homo_sapiens/", NM, "?content-type=application/json") #ensembl id
  #CDH1 exception
  #if(NM == "NM_004360.5") ext.ensembl.id <- paste0("/xrefs/symbol/homo_sapiens/", "NM_004360.3", "?content-type=application/json")
  server.ensembl <- "http://grch37.rest.ensembl.org"
  ensembl.id <- api2(server.ensembl, ext.ensembl.id)
  if (length(ensembl.id)==0){
    NM <- stringr::str_extract(NM, "NM_[0-9]+")
    ext.ensembl.id <- paste0("/xrefs/symbol/homo_sapiens/", NM, "?content-type=application/json")
    ensembl.id <- api2(server.ensembl, ext.ensembl.id)
    if (length(ensembl.id)==0){
    ext.ensembl.id <- paste0("/xrefs/symbol/homo_sapiens/", gene, "?content-type=application/json")
    ensembl.id <- api2(server.ensembl, ext.ensembl.id)
    }
  }
  ensembl.id <- ensembl.id[stringr::str_detect(ensembl.id$id, "ENST")==T,]
  if(stringr::str_detect(NM, "NM_001127500")) ensembl.id <- data.frame(id="ENST00000318493", type = "transcript") #to solve ensembl error

  return(ensembl.id)
}


################################################################################
## Grantham distance
################################################################################

calculateGrantham <- function(a1, a2) {
  grantham <- readr::read_tsv("https://gist.githubusercontent.com/danielecook/501f03650bca6a3db31ff3af2d413d2a/raw/5583a134b36b60762be6cd54002a0f4044338cd0/grantham.tsv", show_col_types = FALSE) %>%
              tidyr:: gather("SECOND", "SCORE", -"FIRST") %>%
              dplyr::filter(.data$SCORE > 0)
  (grantham %>% dplyr::filter(.data$FIRST == a1|.data$FIRST==a2, .data$SECOND == a2 |.data$SECOND==a1))$SCORE
}


################################################################################
## Read web pages
################################################################################

## Without Java script
#' read url
#' @param description url to read
#' @param verbose logical. By default is FALSE``
readUrl <- function(description, verbose = FALSE){

  out <- tryCatch({
    if(verbose)message(paste0("Trying to enter the url: ", description))
    on.exit(close(getConnection(con)))
    con <- url(description, open = "rb")
    if(exists("con")){
      info <- rvest::read_html(x=con)
    }else{
      info <- NA
    }

    #on.exit(close(getConnection(url)))
    return(info)
  },
  error = function(cond){
    if(verbose)message(paste("URL does not seem to exist or it is not working: ", description))
    # Return value in case of error
    return(NA)
  },
  warning = function(cond){
    if(verbose)message(paste("URL caused a warning:", description))
    if(verbose)message("Here's the original warning message:")
    if(verbose)message(cond)
    return(NULL)
  }
)


  return (out)
}

## With Java script
#' readTableUrlJavascript
#' @param url url to query
#' @param browser Which browser to start Rselenium server. By default is "firefox" (the recommended). If you do not have firefox installed try either "chrome" or "phantomjs".
#' @param port port to use for the conneciton
#' @param verbose Logical If TRUE, include status messages (if any). By default is FALSE
#' @param chromever Firefox session
readTableUrlJavascript <- function (url, browser="firefox",  port=4568L, verbose = FALSE, chromever = "109.0.5414.25"){
assertthat::assert_that(browser %in% c("firefox", "chrome","phantomjs"), msg="Only supported for firefox, chrome or phantomjs browser")

if(!requireNamespace("RSelenium", quietly=TRUE)){
  warning("Please install package 'RSelenium' when using 'remote = TRUE'.")
  information <- "Not working"
}else{

if (browser=="firefox"){
  ex.cap <- list("moz:firefoxOptions" = list(args = list('--headless')))
}else if(browser == "chrome"){
  ex.cap <- list("chromeOptions" = list(args = list('--headless')))
}else{
  ex.cap <- NULL
}
  information <- "Not working"
  rD <- try(rD <- RSelenium::rsDriver(browser=browser, chromever=chromever, port=port, extraCapabilities = ex.cap, verbose = verbose ))

  if(class(rD)[1]=="try-error"){
    start.time <- Sys.time()
    end.time <- Sys.time()
    time.taken <- difftime(end.time,start.time, units = "secs")
    try(while(stringr::str_detect(rD, paste0("Selenium server signals port = [0-9]+")) && time.taken < 120){
      port <- port + 1 %>% as.integer()
      rD <- try(rD <- RSelenium::rsDriver(browser=browser, port=port, extraCapabilities = ex.cap  ))
    })
  }
  tryCatch(
    {
      driver <- rD[["client"]]
      driver$navigate(url)
      Sys.sleep(2)
      if(verbose)message("This is the 'try' part")
      try(driver$dismissAlert())
      information <- driver$getPageSource()[[1]] %>% rvest::read_html() %>% rvest::html_table()
      try(driver$close())
    },
    error=function(cond) {
      if(verbose)message("variant exists")
      if(verbose)message(cond)
      return(list())
    },
    warning=function(cond) {
      if(verbose)message(paste("URL caused a warning:", url))
      message(cond)
      return(list())
    },
    finally={
      if(verbose)message(paste("Processed URL:", url))

    }
  )
 try( rD[["server"]]$stop())
  }
  return(information)
}


################################################################################
## Check directory and create subfolder
################################################################################
checkDir <- function(output.dir, subfolder){
  assertthat::assert_that(!is.null(output.dir), msg = "Please set the output directory (output.dir)")
  if (!is.null(output.dir)) {
  assertthat::assert_that(dir.exists(output.dir), msg = "Output directory does not exists, please enter a valid one.")
  if (stringr::str_sub(output.dir, -1)=="/") output.dir <- stringr::str_sub(output.dir, 1, -2)
  }
  # output.dir.subfolder  <- ifelse(is.null(output.dir),
  #                            file.path(getwd(), subfolder),
  #                            file.path(output.dir, subfolder))

  output.dir.subfolder <- file.path(output.dir, subfolder)
  dir.create(output.dir.subfolder, showWarnings = FALSE)
  return(output.dir.subfolder)
}





#' vcfToCodingDNA()
#' @description This function allows to convert variant call format (vcf) file to a tibble with the variants annotated in coding DNA nomenclature considering MANE select transcript.
#' @param file.vcf A filename for a variant call format (vcf or compressed *.vcf.gz.format) file.
#' @param assembly character. It can only be "hg19" or "hg38" assembly.
#' @param annotation character.  It can only be LRG or MANE_select. Variants will be annotated considering LRG or MANE_select transcripts.
#' @param log.dir path where log folder will be created. By default log folder will be created in working diretory.
#' @param verbose logic. By default is FALSE. If TRUE, it includes status messages (if any).
#' @details  The function *vcfToCodingDNA* allows to annotate variants stored in a variant call format (vcf or compressed vcf) in coding DNA nomenclature.
#' To do so it uses Mutalyzer v3 API, thus connection to internet is needed.
#' This function is intended for small vcfs. Vcfs resulting from WES or WGS will take too long (although they work as well).
#' The output of this function can be used to run *vaRbatch* function.
#' @return vcfToHgvs returns a tibble  of four columns (gene, variant, NM, genomic). The number of rows will be equivalent to the number of variants annotated. A log file is generated and stored in log.dir path.
#' @author Elisabet Munté Roca
# AIXO NO FUNCIONA BE ARA MATEIX
# @references
# Morales, J., Pujar, S., Loveland, J.E. et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. Nature. 2022 Apr;604(7905):310-315. PubMed; PubMed Central; DOI: 10.1038/s41586-022-04558-8 . https://www.ncbi.nlm.nih.gov/refseq/MANE/
# Lefter M et al. (2021). Mutalyzer 2: Next Generation HGVS Nomenclature Checker. Bioinformatics, 2021 Sep 15; 37(28):2811-7 (version 3: https://mutalyzer.nl/ )
# Knaus, Brian J., and Niklaus J. Grünwald. "vcfr: a package to manipulate and visualize variant call format data in R." Molecular ecology resources 17.1 (2017): 44-53. (http://dx.doi.org/10.1101/041277)
#' @examples
#' \dontrun{
#' file <- "/path/to/vcf/file.vcf"
#' variants <- vcfToCodingDNA(file.vcf = file, assembly="hg19")
#' }
#' @export

vcfToCodingDNA <- function (file.vcf, assembly, annotation = "LRG", log.dir = NULL, verbose=FALSE){
  time <-  Sys.time() %>%
    stringr::str_replace_all("-|:| ", "_")
#Checks
  assembly <- tolower(assembly)
  assertthat::assert_that(assembly %in% c("hg19", "hg38"), msg = paste(assembly, "is not a valid value for assembly. It must be either hg19 or hg38."))
  assertthat::assert_that(annotation %in% c("LRG", "MANE_select"), msg = paste(annotation, "is not a valid value for annotation. It must be either LRG or MANE_select"))
  assertthat::assert_that(stringr::str_detect(file.vcf, ".vcf$") | stringr::str_detect(file.vcf, ".vcf.gz$"), msg = paste(file.vcf, "is not a vcf file."))
  cat("VCF file detected.")
  assertthat::assert_that(file.exists(file.vcf), msg = paste(file.vcf, "does not exist."))
  assertthat::assert_that(requireNamespace("vcfR", quietly=TRUE), msg = "Please install package 'vcfR' when using vcfToHgvs().")
  log.folder <- checkDir (log.dir, "log")

#Log file creation
  file.create(file.path(log.folder,paste0("vcfToHgvs",time, ".log")))
  write.table(x = "Not able to annotate variants:",
              file = file.path(log.folder,paste0("vcfToHgvs",time, ".log")),
              append=TRUE,
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)

#Process vcf file
  file.variants <- vcfR::read.vcfR(file.vcf, verbose = verbose) %>%
      vcfR::vcfR2tidy(single_frame = TRUE)
  all.variants <- file.variants$dat %>%
      dplyr::rowwise() %>%
      dplyr::mutate(genomic = ifelse(stringr::str_length(.data$REF)==1 && stringr::str_length(.data$ALT)==1,
                                 paste0("g.", .data$POS, .data$REF, ">", .data$ALT ),
                                 ifelse(stringr::str_length(.data$REF)>1 && stringr::str_length(.data$ALT)==1,
                                        paste0(as.numeric(.data$POS+1),"_",as.numeric(.data$POS+stringr::str_length(.data$REF)-1), "del", stringr::str_sub(.data$REF,2, -1)),
                                        ifelse(stringr::str_length(.data$REF)==1 && stringr::str_length(.data$ALT)>1,
                                               paste0(as.numeric(.data$POS), "_",as.numeric(.data$POS+1), "ins", stringr::str_sub(.data$ALT,2, -1)),
                                               ifelse(stringr::str_length(.data$ALT)>1 & stringr::str_length(.data$REF)>1,
                                                      paste0(as.numeric(.data$POS),"_", as.numeric(.data$POS+stringr::str_length(.data$REF)-1),"delins", .data$ALT ),
                                                      NA)))))

    query1 <- "SELECT * FROM NCs;"


    NCs <- connectionDB(query1)[[1]] %>%
      tibble::as_tibble()


    all.variants2 <- merge(all.variants, NCs, by.x="CHROM", by.y = "chr") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(genomic_final = ifelse(assembly=="hg19",
                                   paste0(.data$NC_hg19, ":", .data$genomic),
                                   paste0(.data$NC_hg38, ":", .data$genomic)))
    all.variants2$gene <- NA
    all.variants2$NM_ref <- NA
    all.variants2$cdna <- NA
    all.variants2$gene_2 <- NA
    all.variants2$NM_ref_2 <- NA
    all.variants2$cdna_2 <- NA
    all.variants2$warning <- NA

    if(annotation == "MANE_select"){
      query2 <- "SELECT * FROM MANE;"
      mane <- connectionDB(query2)[[1]] %>%
        tibble::as_tibble()
      anot <- mane %>%
        dplyr::select("Ensembl_Gene", "name", "NM", "chr", paste0("start_", assembly), paste0("end_", assembly))
      names(anot) <- c("id", "name", "NM", "chr", "start", "end")
    }else{
      query3 <- "SELECT * FROM LRG_genes;"
      lrg <- connectionDB(query3)[[1]] %>%
        tibble::as_tibble()
      anot  <- lrg %>%
        dplyr::select("LRG", "genename", "NM", "chr", paste0("start_", assembly), paste0("end_", assembly))
      names(anot) <- c("id", "name", "NM", "chr", "start", "end")
    }


    for (i in 1:nrow(all.variants2)){
      chr.var <- all.variants2$CHROM[i]
      pos.var <- all.variants$POS[i]
      info.var <- anot %>%
        dplyr::filter(.data$chr == chr.var,
                      .data$start <= pos.var ,
                      .data$end >= pos.var )  %>%
        dplyr::select(.data$name, .data$NM)

        tryCatch({
          all.variants2$gene[i] <- info.var$name[1]
          all.variants2$NM_ref[i] <- info.var$NM[1]
          server.mutalyzerv3 <- "https://mutalyzer.nl/api/"
          ext.mutalyzer.v3 <- paste0("normalize/", all.variants2$genomic_final[i])
          mutalyzerv3 <- NA
          mutalyzerv3 <- api2(server.mutalyzerv3, ext.mutalyzer.v3)$equivalent_descriptions$c %>%
            as.data.frame()
           mutalyzerv3.a<- mutalyzerv3 %>%
             dplyr::filter(stringr::str_detect(.data$V1, all.variants2$NM_ref[i] )) %>%
            dplyr::select("V1") %>%
            as.character()

           if(any(mutalyzerv3.a== "character(0)")){
             NM.search <-  stringr::str_extract(all.variants2$NM_ref[i], "NM_[0-9]++")
             mutalyzerv3.a<- mutalyzerv3 %>%
               dplyr::filter(stringr::str_detect(.data$V1, NM.search )) %>%
               dplyr::select("V1") %>%
               as.character()
             NM.sel <- stringr::str_split(mutalyzerv3.a, ":") %>%
               purrr::map(1) %>% as.character()
             all.variants2$warning[i] <- paste(all.variants2$NM_ref[i], "not found in mutalyzerv3", NM.sel, "has been used instead.")
           }


          all.variants2$cdna[i] <- stringr::str_split(mutalyzerv3.a, ":") %>%
            purrr::map(2) %>% as.character()
          if(nrow(info.var > 1)){
            all.variants2$gene_2[i] <- info.var$name[2]
            all.variants2$NM_ref_2[i] <- info.var$NM[2]
            server.mutalyzerv3 <- "https://mutalyzer.nl/api/"
            ext.mutalyzer.v3 <- paste0("normalize/", all.variants2$genomic_final[i])
            mutalyzerv3_2 <- NA
            mutalyzerv3_2 <- try(api2(server.mutalyzerv3, ext.mutalyzer.v3)$equivalent_descriptions$c %>%
                                   as.data.frame() %>%
                                   dplyr::filter(stringr::str_detect(.data$V1, all.variants2$NM_ref_2[i] )) %>%
                                   dplyr::select("V1") %>%
                                   as.character())
            if(!is.na(mutalyzerv3_2)){
              all.variants2$cdna_2[i] <- stringr::str_split(mutalyzerv3_2, ":") %>%
                purrr::map(2) %>% as.character()
            }}

        }, error = function(e){
          write.table(x = all.variants2$genomic_final[i],
                      file = file.path(log.folder,paste0("vcfToHgvs",time, ".log")),
                      append=TRUE,
                      col.names = FALSE,
                      row.names = FALSE)
        })
    }
  var1 <- all.variants2 %>%
    dplyr::select(.data$gene, .data$cdna, .data$NM_ref, .data$genomic_final, .data$warning) %>%
    dplyr::filter(!is.na(.data$cdna) && !is.null(.data$cdna))
  var2 <- all.variants2 %>%
    dplyr::select(.data$gene_2, .data$cdna_2, .data$NM_ref_2, .data$genomic_final, .data$warning) %>%
    dplyr::filter(!is.na(.data$gene_2)) %>%
    dplyr::filter(!is.null(.data$cdna_2))
  names(var1)<- c("gene", "variant", "NM", "genomic", "warning")
  names(var2)<- c("gene", "variant", "NM", "genomic", "warning")
  variants <- rbind(var1, var2)
  write.table(x = paste0("\n Summary: \n", nrow(variants), " variants were correctly annotated"),
              file = file.path(log.folder,paste0("vcfToHgvs",time, ".log")),
              append=TRUE,
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
  message(paste0("Find log file: ",  file.path(log.folder,paste0("vcfToHgvs",time, ".log"))))
  return(variants)
}


nrExtract <- function(all.information){
  nr <- ifelse(all.information$Variant.Info$most.severe.consequence %in% c("splice_donor_variant", "splice_acceptor_variant"),
               stringr::str_extract_all(all.information$codon.stop$canonical.skip.pred$variant, "[0-9]+"),
               ifelse(all.information$Variant.Info$most.severe.consequence %in% c("stop_gained","frameshift_variant" ),
                      stringr::str_extract_all(all.information$Variant.Info$variant, "[0-9]+"),
                      NA)) %>% unlist %>% as.numeric
nr.nt <- ifelse(length(nr)==2,
                nr[2]-nr[1],
                1)
nr.aa <- (nr.nt +1)/3

return(list(nr.nt = nr.nt, nr.aa = nr.aa))
}
