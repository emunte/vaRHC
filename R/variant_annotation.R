
#' Get transcript information
#'
#' @param gene Your gene of interest
#' If the gene doesn't exist in the table, an error message is obtained
#' @param NM Accession number of the transcrit and mRNA from RefSeq.
#' @param CCDS Consensus CDS id https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi.
#' @return The NM, NC and CCDS for that gene
#' @author Elisabet Munté Roca
NMparam <- function(gene , NM=NULL, CCDS=NULL){
  query <- paste0("SELECT * from  transcript WHERE namegene= '", gene ,"' AND maintranscript= 'T' ;")
  nm.nc <- connectionDB(query)[[1]] %>%
    tibble::as_tibble() %>%
    dplyr::select (NM, NC, CCDS)
  assertthat::assert_that(nrow(nm.nc)!=0, msg= "Gene not found in the database.")
  if (!is.null(NM)){
    assertthat::assert_that(is.character(NM) & stringr::str_detect(NM,"NM_[0-9]+\\.[0-9]"), msg="Invalid NM entered")
    NC <- nm.nc$NC
    CCDS <- NA
    #assertthat::assert_that(!is.null(CCDS), msg="You must provide the CCDS id.")
    #assertthat::assert_that(stringr::str_detect(CCDS,"CCDS[0-9]+"), msg="Invalid CCDS entered")
    nm.nc <- tibble::tibble(NM = NM,
                            NC = NC,
                            CCDS = CCDS)
  }

    #assertthat::assert_that(!is.null(NC), msg="You must provide a NC id.")


  return(nm.nc)
}


#' Correct variant hgvs nomenclature
#'
#' @param NM Accession number of the transcrit and mRNA from RefSeq
#' @param NC Accession number of the chromosome RefSeq.
#' @param gene Your gene of interest
#' @param variant Your cdna variant of interest
#' @param skip.pred By default is FALSE and will only be set to TRUE when calculin the predicted skipping effect
#' @return The correct hgvs nomenclature of your variant from Mutalyzer. If the variant is intronic,
#' Mutalyzer V3 is checked where as if the variant is not intronic Mutalyzer v2. A list is given, where the variant and the error message can be checked.
#' @author Elisabet Munté Roca
#' @references
#' "Lefter M et al. (2021). Mutalyzer 2: Next Generation HGVS Nomenclature Checker. Bioinformatics, btab051" (direct link)
#' @noRd

correctHgvsMutalyzerV1 <- function(NM, NC, gene, variant, skip.pred=FALSE){
  intronic <- stringr::str_detect(variant, "[0-9][+]|[0-9][-]")
  #utr <- stringr::str_detect(variant, "c.[+]|c.[-]")

  ###ext for rest apis
  #server.mutalyzer <- "https://v2.mutalyzer.nl/json/" #Mutalyzer's REST API
  server.mutalyzerv3 <- "https://mutalyzer.nl/api/"#Mutalyzer's V3 REST API

  ####we encode the URL
  variant.mutalyzer <- utils::URLencode(paste0(NC, "(",NM, "):",variant),reserved=TRUE)

  #Changing between equivalent versions, to avoid mutalyzer errors.

  message.mutalyzer <- NA

  ext.mutalyzer.v3 <- paste0("normalize/", variant.mutalyzer)
  mutalyzerv3 <- api2(server.mutalyzerv3, ext.mutalyzer.v3)
  #Checking possible mutalyzer errors
  assertthat::assert_that(!any(mutalyzerv3$custom$errors$code=="ERETR"), msg = paste(mutalyzerv3$custom$errors$details, ":try another NM Ref-Seq version"))
  assertthat::assert_that(!any(mutalyzerv3$custom$errors$code=="ESEQUENCEMISMATCH"), msg= mutalyzerv3$custom$errors$details)
  assertthat::assert_that(!any(mutalyzerv3$custom$errors$code=="EPARSE"), msg= paste(mutalyzerv3$custom$errors$details, ": Mutalyzer could not retrieve NM"))
  assertthat::assert_that(!any(mutalyzerv3$custom$errors$code=="EREF"), msg= mutalyzerv3$custom$errors$details)

  #checking problems with NM versions
  as <- "hg19"
  if( any(mutalyzerv3$custom$errors$code=="ENOSELECTORFOUND")){
   query<- paste0("SELECT NC_hg38 FROM NCs where NC_hg19 ='", NC,"';")
   NC.hg38 <- connectionDB(query)[[1]] %>%
     as.character
    variant.mutalyzer <- utils::URLencode(paste0(NC.hg38, "(",NM, "):",variant),reserved=TRUE)
    if(NM=="NM_000535.5")  variant.mutalyzer <- utils::URLencode(paste0(NC.hg38, "(","NM_000535.6", "):",variant),reserved=TRUE)
    if(NM=="NM_002528.5")  variant.mutalyzer <- utils::URLencode(paste0(NC.hg38, "(","NM_002528.6", "):",variant),reserved=TRUE)
    ext.mutalyzer.v3 <- paste0("normalize/", variant.mutalyzer)
    mutalyzerv3 <- api2(server.mutalyzerv3, ext.mutalyzer.v3)
    as <- "hg38"
  }

  #getting mutalyzer information
  cor.variant <- stringr::str_split(mutalyzerv3$normalized_description, ":")[[1]][2]
  html.prot <- mutalyzerv3$protein$description
  mutalyzer.prot.pred <- mutalyzerv3$protein$reference
  cor.prot <- stringr::str_split(html.prot,"\\:")
  if (length(cor.prot)==0){
    cor.prot <- "p.?"
  }else{
    cor.prot<-cor.prot[[1]][2]
    prot2<-cor.prot
  }

  mutalyzer.genomic <- mutalyzerv3$equivalent_descriptions$g
  if(as=="hg38"){
    #genomic nomenclature

    mutalyzer.genomic <- liftOverhg38_hg19(mutalyzer.genomic)

  }

  exons.mut <-   mutalyzerv3$selector_short$exon$c %>%
    tibble::as_tibble() %>% suppressWarnings()
  names(exons.mut) <- c("cStart", "cStop")



  message.mutalyzer2 <- ifelse(mutalyzerv3$corrected_description==mutalyzerv3$normalized_description,
                              "No errors found",
                              "The variant's nomenclature was not ok")

  message.mutalyzer.def <- data.frame(NM_nomenclature=message.mutalyzer, variant_nomenclature=message.mutalyzer2)



  #genomic nomenclature
  #mutalyzer.genomic <- mutalyzerv3$equivalent_descriptions$g
  mutalyzer.other.selected <- list()

  ###Other relevant transcripts
  if (gene=="TP53"){  #it has more important transcripts
    mutalyzer.other.selected <- mutalyzerv3$equivalent_descriptions$c[stringr::str_detect(mutalyzerv3$equivalent_descriptions$c, "NM_000546.5|NM_001126114.2|NM_001126113.2")==TRUE & stringr::str_detect(mutalyzerv3$equivalent_descriptions$c, NM)==FALSE]
  }else if( gene== "CDKN2A"){
    mutalyzer.other.selected <- mutalyzerv3$equivalent_descriptions$c[stringr::str_detect(mutalyzerv3$equivalent_descriptions$c, "NM_000077.4|NM_058195.3")==TRUE & stringr::str_detect(mutalyzerv3$equivalent_descriptions$c, NM)==FALSE]
  }


  if(skip.pred ==TRUE){
    variant.mutalyzer <- paste0(NM,":", variant)
    ext.mutalyzer.v3 <- paste0("normalize/", variant.mutalyzer)
    mutalyzerv3 <- api2(server.mutalyzerv3, ext.mutalyzer.v3)
    cor.variant <- stringr::str_split(mutalyzerv3$normalized_description, ":")[[1]][2]
    html.prot <- mutalyzerv3$protein$description
    cor.prot <- stringr::str_split(html.prot,"\\:")
    if (length(cor.prot)==0){
      cor.prot <- "p.?"
    }else{
      cor.prot<-cor.prot[[1]][2]
      prot2<-cor.prot
    }
  }


  ##final output
  correct.variant <- list(initial.variant = variant,
                          NM = NM,
                          gene = gene,
                          variant = cor.variant,
                          protein = cor.prot,
                          warning = message.mutalyzer.def,
                          genomic = as.character(mutalyzer.genomic),
                          exons = exons.mut,
                          protein_predicted = mutalyzer.prot.pred,
                          other.important.transcripts = mutalyzer.other.selected)
  return (correct.variant)
}


correctHgvsMutalyzer <- function(NM, NC, gene, variant, skip.pred=FALSE){
  intronic <- stringr::str_detect(variant, "[0-9][+]|[0-9][-]")
  #utr <- stringr::str_detect(variant, "c.[+]|c.[-]")

  ###ext for rest apis
  #server.mutalyzer <- "https://v2.mutalyzer.nl/json/" #Mutalyzer's REST API
  server.mutalyzerv3 <- "https://mutalyzer.nl/api/"#Mutalyzer's V3 REST API
  query<- paste0("SELECT NC_hg38 FROM NCs where NC_hg19 ='", NC,"';")
  NC.hg38 <- connectionDB(query)[[1]] %>%
    as.character
  variant.mutalyzer <- utils::URLencode(paste0(NC.hg38, "(",NM, "):",variant),reserved=TRUE)
  #if(NM=="NM_000535.5")  variant.mutalyzer <- URLencode(paste0(NC.hg38, "(","NM_000535.6", "):",variant),reserved=TRUE)
  #if(NM=="NM_002528.5")  variant.mutalyzer <- URLencode(paste0(NC.hg38, "(","NM_002528.6", "):",variant),reserved=TRUE)
  ext.mutalyzer.v3 <- paste0("normalize/", variant.mutalyzer)
  mutalyzerv3 <- api2(server.mutalyzerv3, ext.mutalyzer.v3)


  #Changing between equivalent versions, to avoid mutalyzer errors.
  message.mutalyzer <- NA
  i =0
  NM.old <- NM
 while (any(mutalyzerv3$custom$errors$code=="ENOSELECTORFOUND") && i <5){
   NM.split <- stringr::str_split(NM, "\\.") %>% unlist()
   NM <- paste0(NM.split[1], ".", NM.split[2] %>% as.numeric()+1)
   variant.mutalyzer <- utils::URLencode(paste0(NC.hg38, "(",NM, "):",variant),reserved=TRUE)
   ext.mutalyzer.v3 <- paste0("normalize/", variant.mutalyzer)
   mutalyzerv3 <- api2(server.mutalyzerv3, ext.mutalyzer.v3)
   message.mutalyzer <- paste("Not", NM.old, "selector found in reference, to continue", NM, "has been used instead")
    i = i +1
   }


  #Checking possible mutalyzer errors
  assertthat::assert_that(!any(mutalyzerv3$custom$errors$code=="ERETR"), msg = paste(mutalyzerv3$custom$errors$details, ":try another NM Ref-Seq version"))
  assertthat::assert_that(!any(mutalyzerv3$custom$errors$code=="ESEQUENCEMISMATCH"), msg= mutalyzerv3$custom$errors$details)
  assertthat::assert_that(!any(mutalyzerv3$custom$errors$code=="EPARSE"), msg= paste(mutalyzerv3$custom$errors$details, ": Mutalyzer could not retrieve NM"))
  assertthat::assert_that(!any(mutalyzerv3$custom$errors$code=="EREF"), msg= mutalyzerv3$custom$errors$details)
  assertthat::assert_that(!any(mutalyzerv3$custom$errors$code=="ENOSELECTORFOUND"), msg = paste(mutalyzerv3$custom$errors$details, ": Mutalyzer could not retrieve NM, try another version."))
  #checking problems with NM versions
  as <- "hg38"

  #getting mutalyzer information
  cor.variant <- stringr::str_split(mutalyzerv3$normalized_description, ":")[[1]][2]
  html.prot <- mutalyzerv3$protein$description
  mutalyzer.prot.pred <- mutalyzerv3$protein$reference
  cor.prot <- stringr::str_split(html.prot,"\\:")

  if (length(cor.prot)==0){
    cor.prot <- "p.?"
  }else{
    cor.prot<-cor.prot[[1]][2]
    prot2<-cor.prot
  }

  mutalyzer.genomic_hg38 <- mutalyzerv3$equivalent_descriptions$g
  # if(as=="hg38"){
  #   #genomic nomenclature
  #   mutalyzer.genomic_hg37 <- liftOverhg38_hg19(mutalyzer.genomic_hg38)
  #
  # }

  exons.mut <-   mutalyzerv3$selector_short$exon$c %>%
    tibble::as_tibble() %>% suppressWarnings()
  names(exons.mut) <- c("cStart", "cStop")

  exons.mut.genomic.hg38 <-   mutalyzerv3$selector_short$exon$g %>%
    tibble::as_tibble() %>% suppressWarnings()


  cds.pos.38 <- mutalyzerv3$selector_short$cds$g
  cds.c <- mutalyzerv3$selector_short$cds$c

  variant.mutalyzer.hg37 <- utils::URLencode(paste0(NC, "(",NM, "):",variant),reserved=TRUE)
  ext.mutalyzer.v3.hg37 <- paste0("normalize/", variant.mutalyzer.hg37)

    mutalyzerv3.hg37 <- api2(server.mutalyzerv3, ext.mutalyzer.v3.hg37 )
    if(!is.null(mutalyzerv3.hg37$selector_short)){
      exons.mut.genomic.hg37 <-   mutalyzerv3.hg37$selector_short$exon$g %>%
        tibble::as_tibble() %>%
        suppressWarnings()
      mutalyzer.genomic_hg37 <- mutalyzerv3.hg37$equivalent_descriptions$g
      cds.pos.37 <- mutalyzerv3.hg37$selector_short$cds$g
    }else{
      mutalyzer.genomic_hg37 <- liftOverhg38_hg19(mutalyzer.genomic_hg38)
      chr <- stringr::str_extract(mutalyzer.genomic_hg38, "[0-9]++") %>% as.integer()
      v1 <-liftOverPos(exons.mut.genomic.hg38$V1, chr)
      v2 <-liftOverPos(exons.mut.genomic.hg38$V2, chr)

      exons.mut.genomic.hg37 <- tibble::tibble(V1=v1, V2 = v2)

      v1.cds <- liftOverPos(cds.pos.38[1], chr)
      v2.cds <- liftOverPos(cds.pos.38[2], chr)
      cds.pos.37 <- cbind(v1.cds, v2.cds)
    }

    cds <- cbind(cds.pos.37, cds.pos.38, cds.c) %>% data.frame()
    names(cds) <- c("v1.cds.37", "v2.cd.37", "v1.cds.38", "v2.cds.38", "v1.cds.c", "v2.cds.c")


  message.mutalyzer2 <- ifelse(mutalyzerv3$corrected_description==mutalyzerv3$normalized_description,
                               "No errors found",
                               "The variant's nomenclature was not ok")

  message.mutalyzer.def <- data.frame(NM_nomenclature=message.mutalyzer, variant_nomenclature=message.mutalyzer2)



  #genomic nomenclature
  #mutalyzer.genomic <- mutalyzerv3$equivalent_descriptions$g
  mutalyzer.other.selected <- list()

  ###Other relevant transcripts
  if (gene=="TP53"){  #it has more important transcripts
    mutalyzer.other.selected <- mutalyzerv3$equivalent_descriptions$c[stringr::str_detect(mutalyzerv3$equivalent_descriptions$c, "NM_000546.5|NM_001126114.2|NM_001126113.2")==TRUE & stringr::str_detect(mutalyzerv3$equivalent_descriptions$c, NM)==FALSE]
  }else if( gene== "CDKN2A"){
    mutalyzer.other.selected <- mutalyzerv3$equivalent_descriptions$c[stringr::str_detect(mutalyzerv3$equivalent_descriptions$c, "NM_000077.4|NM_058195.3")==TRUE & stringr::str_detect(mutalyzerv3$equivalent_descriptions$c, NM)==FALSE]
  }


  if(skip.pred ==TRUE){
    variant.mutalyzer <- paste0(NM,":", variant)
    ext.mutalyzer.v3 <- paste0("normalize/", variant.mutalyzer)
    mutalyzerv3 <- api2(server.mutalyzerv3, ext.mutalyzer.v3)
    cor.variant <- stringr::str_split(mutalyzerv3$normalized_description, ":")[[1]][2]
    html.prot <- mutalyzerv3$protein$description
    cor.prot <- stringr::str_split(html.prot,"\\:")
    if (length(cor.prot)==0){
      cor.prot <- "p.?"
    }else{
      cor.prot<-cor.prot[[1]][2]
      prot2<-cor.prot
    }
  }

if(stringr::str_detect(cor.variant, "\\[[0-9+]\\]")){

server <- "https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg19/"
ext <- paste0(NM,"%3A", variant, "/", NM, "?content-type=application%2Fjson")
variant.validator <- api2(server, ext)
var.validator <- variant.validator[[1]]$hgvs_refseqgene_variant
cor.variant <- variant.validator[[1]]$hgvs_transcript_variant %>% stringr::str_split(":") %>% purrr::map(2) %>% unlist()
mutalyzer.genomic_hg37 <- variant.validator[[1]]$primary_assembly_loci$grch37$hgvs_genomic_description
mutalyzer.genomic_hg38 <- variant.validator[[1]]$primary_assembly_loci$grch38$hgvs_genomic_description
}

  ##final output
  correct.variant <- list(initial.variant = variant,
                          NM = NM,
                          gene = gene,
                          variant = cor.variant,
                          protein = cor.prot,
                          warning = message.mutalyzer.def,
                          genomic_hg37 = as.character(mutalyzer.genomic_hg37),
                          genomic_hg38 = as.character(mutalyzer.genomic_hg38),
                          exons = exons.mut,
                          exons_g_37 = exons.mut.genomic.hg37,
                          exons_g_38 = exons.mut.genomic.hg38,
                          cds = cds,
                          protein_predicted = mutalyzer.prot.pred,
                          other.important.transcripts = mutalyzer.other.selected)
  return (correct.variant)
}

#' Variant info from Ensembl VEP
#'
#' @param NM Accession number of the transcrit and mRNA from RefSeq
#' @param NC Accession number of the chromosome RefSeq. It can be ommited if the variant is exonic.
#' @param CCDS The Consensus CDS code
#' @param gene The gene where the variant is located
#' @param variant The variant in coding nomencature
#' @param variant.mutalyzer By default is null. Obtained with the correctHgvsMutalyzer function contains all information of the variant from mutalyzer.
#' @param skip.pred By default is false. If you are predicting the consequence of a skipping then it should be set to TRUE
#' @return The chromosome, the initial genomic coordinate, the final genomic
#' coordinate, the reference allele, the alternative allele, the most important
#' variant consequence and if the variant is located in a protein domain
#' @author Elisabet Munté Roca
#' @references McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. Genome Biology Jun 6;17(1):122. (2016) doi:10.1186/s13059-016-0974-4
#' @noRd

varDetails <- function (NM, NC=NULL, CCDS, gene, variant, variant.mutalyzer=NULL, skip.pred=FALSE){
  if(is.null(variant.mutalyzer)){
    variant.mutalyzer <- correctHgvsMutalyzer(NM,NC, gene, variant)
  }
  if(is.null(variant)) variant<- variant.mutalyzer$initial.variant
  if(is.null(gene)) gene <- variant.mutalyzer$gene

  #get ensembl id
  ensembl.id <- ensemblTranscript(NM, gene)$id

  #we construct vep extension
  server.ensembl <- "http://grch37.rest.ensembl.org" #Ensembl's REST API
  ext.vep <-  paste0("/vep/human/hgvs/",ensembl.id,":",variant.mutalyzer$variant)
  # ext.vep <-  paste0("/vep/human/hgvs/",NM,":",variant.mutalyzer$variant)
  # if(NM=="NM_000314.6")ext.vep <- paste0("/vep/human/hgvs/NM_000314.8:",variant.mutalyzer$variant)
  # if(NM=="NM_005228.4")ext.vep <- paste0("/vep/human/hgvs/NM_005228.5:",variant.mutalyzer$variant) #100% of identity
  # if(NM=="NM_030930.3")ext.vep <- paste0("/vep/human/hgvs/NM_030930.4:",variant.mutalyzer$variant) #100% of identity
  coordinates <- api2(server.ensembl, ext.vep)
  assertthat::assert_that(length(coordinates$error)==0, msg=coordinates$error)
  chr<-unlist(coordinates$seq_region_name)
  strand <- unlist(coordinates$strand)

  ##to get the correct start and end
  ext.vep2 <- paste0("/vep/human/hgvs/",chr,":",unlist(stringr::str_split(variant.mutalyzer$genomic_hg37, ":"))[2])

  #we get the information
  coordinates2 <- api2(server.ensembl, ext.vep2)
  assertthat::assert_that(length(coordinates2$error)==0, msg=coordinates2$error)
  end <-coordinates2$end[[1]]#coordenada final
  start <-coordinates2$start[[1]]#coordenada inicio
  allele.string <- coordinates2$allele_string[[1]] %>%
    stringr::str_split("/") %>% unlist
  ref <- allele.string[1]#reference
  alt <- allele.string[2]#alternative
  most.severe.consequence <- tibble::tibble(coordinates$transcript_consequences[[1]])
  if(ncol(most.severe.consequence)==0) most.severe.consequence <- tibble::tibble(coordinates2$transcript_consequences[[1]])
  most.severe.consequence <- most.severe.consequence %>%
    dplyr::filter (.data$transcript_id==as.character(ensembl.id)) %>%
    dplyr::select("consequence_terms") %>%
    unlist
  if(length(most.severe.consequence)>1){
    if(isFALSE(skip.pred)){
      num <- which(!stringr::str_detect("splice_region_variant", most.severe.consequence) & !stringr::str_detect("splice_polypyrimidine_tract_variant", most.severe.consequence))
      cons1 <- most.severe.consequence[num][1]
      num2 <- which(stringr::str_detect("splice_region_variant", most.severe.consequence))
      cons2 <- most.severe.consequence[num2]
      num3 <-  which(stringr::str_detect("splice_polypyrimidine_tract_variant", most.severe.consequence))
      cons3 <- most.severe.consequence[num3]
      most.severe.consequence <- c(cons1, cons2, cons3) %>% as.character()
    }else {
      num <- which(stringr::str_detect("inframe_deletion", most.severe.consequence) | stringr::str_detect("frameshift_variant", most.severe.consequence))
      if(length(num)!=0){
        most.severe.consequence <-most.severe.consequence[num]
      }else{
        conse <- stringr::str_split(variant.mutalyzer$genomic_hg37, ":g.") %>%
          purrr::map(2) %>%
          unlist() %>%
          stringr::str_split(pattern = "del") %>%
          purrr::map(1) %>% unlist()
          pos.g <- stringr::str_extract_all(conse, "[0-9]+") %>%
            unlist() %>%
            as.numeric()
          division <- (pos.g[2]-pos.g[1]+1)/3
          division.exact <- (pos.g[2]-pos.g[1]+1) %/% 3
          most.severe.consequence <- ifelse(division == division.exact,
                                            "inframe_deletion",
                                            "frameshift_variant")
      }
    }
  }
  object <- data.frame(NM=variant.mutalyzer$NM, variant=variant, protein=variant.mutalyzer$protein, strand=strand, gene=gene, most.severe.consequence=most.severe.consequence[1])
  exons <- coordNonCoding(variant.mutalyzer, object)
  #%>%
   # dplyr::mutate(cStop)

  if (!(most.severe.consequence[1] %in% c("3_prime_UTR_variant", "5_prime_UTR_variant"))){
    if(stringr::str_detect(variant, "-")){
      pos.var <- as.numeric(stringr::str_extract(variant, "[0-9]+"))-1
    }else{
      pos.var <- stringr::str_extract(variant, "[0-9]+")
    }
    variant.exon <- exons %>%
      dplyr::mutate (cStop = stringr::str_replace(.data$cStop, "\\*[0-9]*", "1000000000"),
                     cStart = stringr::str_replace(.data$cStart, "\\*[0-9]*", "NA")) %>%
      dplyr::filter( as.numeric(.data$cStart) <= as.numeric(pos.var),
                     as.numeric(.data$cStop) >= as.numeric(pos.var)) %>%
      dplyr::select("exon")
  }else{
    variant.exon <- NA
  }
  #adjust synonymous nomenclature
  protein <- protsyn2(object, variant.mutalyzer)
  if(skip.pred == FALSE && any( stringr::str_detect(most.severe.consequence, "frameshift_variant"))){
    most.severe.consequence <- ifelse(rep(most.severe.consequence=="frameshift_variant" && length(stringr::str_extract_all(variant.mutalyzer$protein, "[0-9]+"))==1,3),
                                      "stop_gained",
                                      most.severe.consequence) #sometimes frameshift variants are nonsense
    most.severe.consequence <- ifelse(rep(stringr::str_detect(variant.mutalyzer$protein, "fs"),3),
                                      "frameshift_variant",
                                      most.severe.consequence)
  }
  domain.info <- domainProt(chr, start, end) %>%
                 unique
  uvs.df <- tibble::tibble (gene = variant.mutalyzer$gene,
                            NM= variant.mutalyzer$NM,
                            initial.var = variant,
                            variant = variant.mutalyzer$variant,
                            protein = protein,
                            genomic = variant.mutalyzer$genomic_hg37,
                            genomic_hg38 = variant.mutalyzer$genomic_hg38,
                            chr = chr,
                            start = start,
                            end = end,
                            ref = ref,
                            alt = alt ,
                            strand = strand,
                            exon_intron = variant.exon,
                            ensembl.id = ensembl.id,
                            most.severe.consequence = most.severe.consequence[1],
                            most.severe.consequence.1 = most.severe.consequence[2],
                            most.severe.consequence.2= most.severe.consequence[3],
                            CCDS = CCDS,
                            domain.info = t(domain.info)) %>%
    tibble::remove_rownames()

  return(uvs.df)
}



geneLrgCoord <- function (object){
  query<- paste0("SELECT  l.transcript, l.namegene,l.coordinates, l.transcript2, l.cds_start, l.cds_end, l.strand FROM LRG l LEFT JOIN  transcript t ON t.ensembltranscriptID=l.transcript_id WHERE l.namegene= '",object$gene ,"' AND t.NM='", object$NM, "'; ")
  lrg.gene <- connectionDB(query)[[1]]
  # if(nrow(lrg.gene)==0){
  #   query <- paste0("SELECT l.transcript_id, l.coordinates, l.strand FROM transcript t INNER JOIN noLRG l ON t.ensembltranscriptID=l.transcript_id WHERE l.namegene= '",object$gene ,"' AND t.NM='", object$NM, "'; ")
  #   lrg.gene <- connectionDB(query)[[1]]
  #   names(lrg.gene)[1]<-"transcript"
  # }

  if(nrow(lrg.gene)>0){
  query.hg38 <- paste0("SELECT transcript, coordinates_hg38 FROM LRG_hg38 WHERE transcript ='", lrg.gene$transcript, "';")
  lrg.gene.hg38 <- connectionDB(query.hg38)[[1]]
  exon.nocodi <- unlist(stringr::str_split(lrg.gene["coordinates"], ",")) %>%
                 stringr::str_split(pattern="-")
  exon.nocodi.hg38 <- unlist(stringr::str_split(lrg.gene.hg38["coordinates_hg38"], ",")) %>%
    stringr::str_split(pattern="-")
  num.exons <- length(exon.nocodi)
  exons.ord <- ifelse(rep(object$strand ==1, num.exons), 1:num.exons, num.exons:1)

  coordinates.exon.nocodi <- data.frame(exon = exons.ord,
                                        V1 = unlist(purrr::map(exon.nocodi,1)),
                                        V2=unlist(purrr::map(exon.nocodi, 2)),
                                        V1_hg38 = unlist(purrr::map(exon.nocodi.hg38,1)),
                                        V2_hg38 = unlist(purrr::map(exon.nocodi.hg38,1)))  %>%
                                        dplyr::mutate(transcript = lrg.gene$transcript) %>%
                                        dplyr::arrange(.data$exon) %>%
                                        dplyr::relocate ("transcript")

  return (coordinates.exon.nocodi)
  } else{ return (NA)}

}

coordNonCoding <- function (variant.mutalyzer, object){
  coordinates.exon.nocodi <- try(geneLrgCoord(object))
  if(any(is.na(coordinates.exon.nocodi))){
    if(object$strand==1){
    coordinates.exon.nocodi <- data.frame(exon = 1:nrow(variant.mutalyzer$exons_g_37),
                                          V1 = variant.mutalyzer$exons_g_37$V1,
                                          V2=variant.mutalyzer$exons_g_37$V2,
                                          V1_hg38 = variant.mutalyzer$exons_g_38$V1,
                                          V2_hg38 = variant.mutalyzer$exons_g_38$V2)  %>%
      dplyr::mutate(transcript = object$NM) %>%
      dplyr::arrange(.data$exon) %>%
      dplyr::relocate ("transcript")
    }else{
      coordinates.exon.nocodi <- data.frame(exon = 1:nrow(variant.mutalyzer$exons_g_37),
                                            V1 = variant.mutalyzer$exons_g_37$V2,
                                            V2=variant.mutalyzer$exons_g_37$V1,
                                            V1_hg38 = variant.mutalyzer$exons_g_38$V2,
                                            V2_hg38 = variant.mutalyzer$exons_g_38$V1)  %>%
        dplyr::mutate(transcript = object$NM) %>%
        dplyr::arrange(.data$exon) %>%
        dplyr::relocate ("transcript")
    }
  }
  query <- paste0("SELECT exon, cStart, cStop FROM LRG_cds WHERE LRG_id= '", coordinates.exon.nocodi$transcript[1] ,"'; ")
  lrg.cds <- connectionDB(query)[[1]]

  if(nrow(lrg.cds)>0){
    all.exons <- merge(coordinates.exon.nocodi, lrg.cds, by="exon")
  }else{
    all.exons <- cbind( coordinates.exon.nocodi, variant.mutalyzer$exons[,c("cStart", "cStop")])
  }

  return(all.exons)
}


#' Variant's domain from UCSC
#' @param chr Chromosome where the variant is located
#' @param start genomic start of the variant
#' @param end genomic end of the variant
#' @return The protein domain where the variant is located (Uniprot)
#' @author Elisabet Munté Roca
#' @noRd
domainProt <- function(chr, start, end){
  server.ucsc <- "http://api.genome.ucsc.edu/getData/track?" #UCSC's REST API
  strandi <- ifelse(end >= start,
                    paste0(";start=",start-1, ";end=", end),
                    paste0(";start=",end, ";end=", start))
  ext.domain.prot <- paste0("genome=hg19;track=unipDomain;chrom=chr", chr, strandi)
  domain.ucsc <- api2(server.ucsc, ext.domain.prot)
  domain.info <- domain.ucsc$unipDomain$name
  domain.info[is.null(domain.info)]<-NA
  return(domain.info)
}

toProtein <- function(prot){
  aa.pos <- stringr::str_extract(prot, "[0-9]+")
  aa.ref <- purrr::map(stringr::str_extract_all(prot, "[A-z]+"),2) %>%
            unlist
  aa.alt <- purrr::map(stringr::str_extract_all(prot, "[A-z]+"),3) %>%
            unlist
  aa <- list(aa.pos=aa.pos, aa.ref=aa.ref, aa.alt = aa.alt)
  return(aa)
}

aaShort <- function (aa){
  aa.total <- names(Biostrings::AMINO_ACID_CODE)[Biostrings::AMINO_ACID_CODE==aa]
  return(aa.total)
}


gnomADnomen <- function(object){
  del.ins.sel <- stringr::str_detect(object$variant,"dup")|stringr::str_detect(object$variant,"del")|stringr::str_detect(object$variant,"ins")
  if (del.ins.sel== TRUE){
    web.mutalyzer <- xml2::read_html(utils::URLencode(paste0("https://v2.mutalyzer.nl/name-checker?description=",object$genomic)))
    sequence <- web.mutalyzer %>% rvest::html_nodes("pre")
    sequence1 <- sequence[1]
    seq1 <- stringr::str_split (as.character(sequence1), "\\>|\\<")[[1]][3]
    seq2 <- stringr::str_split (as.character(sequence1), "\\>|\\<")[[1]][5]

    #we split the seq1 for deletions and seq2 for dups and ins
    seq.split <- (ifelse(stringr::str_detect(object$variant, "del"), stringr::str_split(seq1, " "), stringr::str_split(seq2, " ")) %>%
                    unlist())
    n <- stringr::str_length(seq.split[1]) #the number of characters given before the pattern (the mutation)
    seq.pattern <- seq.split[2]
    length.pattern <- stringr::str_length(seq.pattern)

    #we get the wt sequence and the mutate sequence
    sequence.wt <- ifelse(stringr::str_detect(object$variant, "del"), gsub(" ", "", seq1), stringr::str_replace_all(seq1,"-| ","")) %>%
                   unlist()
    sequence.mutate <- ifelse(stringr::str_detect(object$variant, "del"), stringr::str_replace_all(seq2,"-| ",""), gsub(" ", "", seq2) ) %>%
                    unlist()

    #matrix for each position of this string, we will save if the mutation is in that position if we obtain the same result (true), or not (FALSE)
    y=1:n
    str.look<- ifelse(stringr::str_detect(object$variant, "del"), sequence.wt, sequence.mutate)
    str.compare <- ifelse(stringr::str_detect(object$variant, "del"),  sequence.mutate, sequence.wt)
    part1 <- stringr::str_sub(str.look, 1, y)
    part2 <- stringr::str_sub(str.look, y+length.pattern+1, stringr::str_length(str.look))
    string <- paste0(part1, part2)
    position<- string==str.compare #matrix for each position of this string, we will save if the mutation is in that position if we obtain the same result (true), or not (FALSE)
    values <- grep(TRUE, position) %>% as.numeric()
    num.rest<-n-values[1] #we obtain the value where is the mutation located to the most leff position, and we do n-that number
    start2<-object$start-num.rest-1 #-1 as gnomAD gives the position of one base before the change
    ref <-ifelse(stringr::str_detect(object$variant, "del"),
                 stringr::str_sub(str.look, values[1], values[1]+length.pattern),
                 stringr::str_sub(str.look, values[1],values[1]))
    alt<-ifelse(stringr::str_detect(object$variant, "del"),
                stringr::str_sub(str.compare, values[1],values[1]),
                stringr::str_sub(sequence.mutate, values[1], values[1]+length.pattern)) #we obtain the alt bp
    #if (object$strand==1){
      if(stringr::str_detect(object$variant, "del")){
        gnomAD<-c(object$chr,start2,ref, alt)#we obtain the gnomAD
      }else{
        gnomAD<-c(object$chr,start2,ref,alt)#we obtain the gnomAD
      }
    # }else if( object$strand == -1){
    #   seq1 <- DescTools::StrRev(seq1) #we need the reverse sequence
    #   if(stringr::str_detect(object$variant, "del")){
    #     gnomAD<- c(object$chr,object$start-1, ref, alt)
    #   }else{
    #     gnomAD<- c(object$chr,object$end-stringr::str_length(object$alt),ref, alt)
    #   }
    # }
  }else{
    gnomAD<-c(object$chr, object$start, object$ref, object$alt)
  }

  gnomAD.variant<-paste0(object$chr, "-", gnomAD[2],"-",gnomAD[3], "-", gnomAD[4])#variants name in gnomAD

  return(gnomAD.variant)
}

protsyn <- function(object, var.mut){
  if (object$most.severe.consequence %in% c("synonymous_variant")){
    pos.dna<-stringr::str_extract(object$variant, "[0-9]+")
    pos.aa<-as.integer(as.numeric(pos.dna)/3)
    if (as.numeric(pos.dna) %% 3 >0){ #if the division is not exact, we need to add 1 position
      pos.aa <- pos.aa + 1
    }
    seq.prot <- stringr::str_sub(var.mut$protein_predicted, pos.aa, pos.aa)
    aa <- Biostrings::AMINO_ACID_CODE[seq.prot]
    protein <- paste0("p.(", aa, pos.aa, aa,")")
  }
  prot <- ifelse(object$most.severe.consequence %in% c("synonymous_variant"), protein, object$protein)
  return(prot)
}

protsyn2 <- function(object, var.mut){
  if (object$most.severe.consequence %in% c("synonymous_variant")){
    pos.dna <- stringr::str_extract(object$variant, "[0-9]+")
    pos.aa <- as.integer(as.numeric(pos.dna)/3)
    if (as.numeric(pos.dna) %% 3 >0){ #if the division is not exact, we need to add 1 position
      pos.aa <- pos.aa + 1
    }
    seq.prot <- stringr::str_sub(var.mut$protein_predicted, pos.aa, pos.aa)
    aa <- Biostrings::AMINO_ACID_CODE[seq.prot]
    protein <- paste0("p.(", aa, pos.aa, "=)")
  }
  prot <- ifelse(object$most.severe.consequence %in% c("synonymous_variant"), protein, object$protein)
  return(prot)
}

toGenomic <- function(NM, NC, cor.variant, gene){
  mutalyzer <- correctHgvsMutalyzer(NM = NM, NC = NC, gene = gene, variant = cor.variant)
  return(mutalyzer$genomic_hg37)
}

CodingTranscriptCds <- function(object){
  server.ucsc <- "http://api.genome.ucsc.edu/getData/track?" #UCSC's REST API
  ext.ucsc <- paste0("genome=hg19;track=ccdsGene;chrom=chr", unique(object$chr))
  ucsc <- api2(server.ucsc, ext.ucsc)
  exo <- ucsc$ccdsGene[ucsc$ccdsGene$name==object$CCDS,]
  if(nrow(exo)!=0){
  coordinates.exon <- data.frame(start =unlist(stringr::str_split(exo$exonStarts, ",")), end=unlist(stringr::str_split(exo$exonEnds, ","))) %>%
    dplyr::filter (start!="")}else{
      coordinates.exon <- exo
    }
  return(coordinates.exon)
}

liftOverhg38_hg19  <- function (genomic){
  assertthat::assert_that(requireNamespace("rtracklayer", quietly=TRUE), msg = "Please install  'rtracklayer' package to be able to compute this variant'")
  assertthat::assert_that(requireNamespace("GenomicRanges", quietly=TRUE), msg = "Please install  'GenomicRanges' package to be able to compute this variant'")
  assertthat::assert_that(requireNamespace("IRanges", quietly=TRUE), msg = "Please install  'IRanges' package to be able to compute this variant'")

  chr <- stringr::str_split(genomic, "\\.") %>%
    purrr::map(1) %>%
    stringr::str_extract("[0-9]+") %>%
    as.integer()

  pos1 <- stringr::str_split(genomic, "g.") %>%
    purrr::map(2) %>%
    stringr::str_split("_") %>%
    stringr::str_extract("[0-9+]+") %>%
    as.integer()

  coord <- stringr::str_split(genomic, "g.") %>%
    purrr::map(2) %>%
    stringr::str_extract_all("[0-9+]+") %>%
    unlist() %>%
    as.integer()
  large <- ifelse(length(coord) == 1,
         1,
         coord[2] - coord[1] +1)

  # assertthat::assert_that(file.exists(system.file(package="vaRHC", "extdata", "hg38ToHg19.over.chain")),
  # msg = "Please download and unzip 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' in vaRHC/extdata folder to be able to compute this variant")
  #   path = system.file(package="vaRHC", "extdata", "hg38ToHg19.over.chain")

  if (!file.exists(system.file(package="vaRHC", "extdata", "hg38ToHg19.over.chain"))){
    #message("Trying to download ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz file")
    #try(utils::download.file(url = "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
    #                     destfile = file.path(.tmp, "hg38ToHg19.over.chain.gz")))
    assertthat::assert_that(file.exists(system.file(package="vaRHC", "extdata", "hg38ToHg19.over.chain.gz")))
    R.utils::gunzip(system.file(package="vaRHC", "extdata", "hg38ToHg19.over.chain.gz")) #unzip the file
  }

  path = system.file(package="vaRHC", "extdata", "hg38ToHg19.over.chain")
  ch = rtracklayer::import.chain(path)
  variantGR <- GenomicRanges::GRanges(seqnames = paste0("chr", chr),
                       IRanges::IRanges(pos1, width = large)
                       )
  lift <- rtracklayer::liftOver(variantGR, ch)
  start <- lift[[1]]@ranges@start

  query <- paste0("SELECT * FROM NCs WHERE chr='", chr, "'")
  NC.table <- connectionDB(query)[[1]]


 genomichg19 <- stringr::str_replace(genomic, pos1 %>% as.character(), start %>% as.character())
 genomichg19 <- stringr::str_replace(genomichg19, NC.table$NC_hg38 %>% as.character(), NC.table$NC_hg19 %>% as.character())
 if(large>=2) genomichg19 <-  stringr::str_replace(genomichg19, (pos1 +large-1 ) %>% as.character(), (start + large -1) %>% as.character())
  return (genomichg19)
}

liftOverPos <- function(vector.pos, chr){
  lapply(vector.pos, function(x){
    x <- as.numeric(x)
    path = system.file(package="vaRHC", "extdata", "hg38ToHg19.over.chain")
    ch = rtracklayer::import.chain(path)
    variantGR1 <- GenomicRanges::GRanges(seqnames = paste0("chr", chr),
                                         IRanges::IRanges(x, width = 1))
    lift1 <- rtracklayer::liftOver(variantGR1, ch)
    lift1[[1]]@ranges@start %>% as.character()
  }) %>% unlist()
}
#
# anotate_vcf <- function(vcf){
#   library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#   vcf <- VariantAnnotation::readVcf(file.vcf, genome = "hg19")
#   seqlevels(vcf) <- paste0("chr", seqlevels(vcf))
#   loc <- locateVariants(vcf, txdb, VariantAnnotation::CodingVariants())
#
# }





