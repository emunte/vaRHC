
#' Get transcript information
#'
#' @param gene Your gene of interest
#' If the gene doesn't exist in the table, an error message is obtained
#' @param nm.info By default NULL. If not a valid NM (Ref-seq Accession number) must be provided. Be careful because only default NM have been tested.
#' @return The NM, NC and CCDS for that gene
#' @author Elisabet Munté Roca
#' @examples
#' 
#' NMparam(gene="BRCA1")
#' NMparam(gene="TP53")
#' NMparam(gene="TP53", NM = "NM_001126113.2", NC = "NC_000017.10", CCDS= "CCDS45606")

NMparam <- function(gene , NM=NULL, NC= NULL, CCDS=NULL){
  query <- paste0("SELECT * from  transcript WHERE namegene= '", gene ,"' AND maintranscript= 'T' ;")
  if (is.null(NM)){
    nm.nc <- connectionDB(query)[[1]] %>% 
      tibble::as_tibble() %>% 
      dplyr::select (NM, NC, CCDS)
  }else{
    assertthat::assert_that(is.character(nm.info) & stringr::str_detect(nm.info,"NM_[0-9]+\\.[0-9]"), msg="Invalid NM entered")
    assertthat::assert_that(!is.null(NC), msg="You must provide a NC id.")
    assertthat::assert_that(stringr::str_detect(NC,"NC_[0-9]+\\.[0-9]"), msg="Invalid NC entered")
    assertthat::assert_that(!is.null(CCDS), msg="You must provide the CCDS id.")
    assertthat::assert_that(stringr::str_detect(CCDS,"CCDS[0-9]+"), msg="Invalid CCDS entered")
    nm.nc <- tibble::tibble(NM = NM, 
                            NC = NC, 
                            CCDS = CCDS) 
  }
  return(nm.nc)
}


#' Correct variant hgvs nomenclature
#'
#' @param NM Accession number of the transcrit and mRNA from RefSeq
#' @param NC Accession number of the chromosome RefSeq. It can be ommited if the variant is exonic.
#' @param gene Your gene of interest
#' @param variant Your cdna variant of interest
#' @return The correct hgvs nomenclature of your variant from Mutalyzer. If the variant is intronic,
#' Mutalyzer V3 is checked where as if the variant is not intronic Mutalyzer v2. A list is given, where the variant and the error message can be checked.
#' @author Elisabet Munté Roca
#' @references
#' "Lefter M et al. (2021). Mutalyzer 2: Next Generation HGVS Nomenclature Checker. Bioinformatics, btab051" (direct link) 
#' @examples
#' correctHgvsMutalyzer(NM="NM_000051.3", gene="ATM", variant="c.3436G>A")
#' correctHgvsMutalyzer(NM="NM_007294.3", NC="NC_000017.10", gene="BRCA1", variant="c.-19-85C>G")


correctHgvsMutalyzer <- function(NM, NC=NULL, gene, variant){
  intronic <- stringr::str_detect(variant, "[0-9][+]|[0-9][-]")
  utr <- stringr::str_detect(variant, "c.[+]|c.[-]")
  
  ###ext for rest apis
  server.mutalyzer <- "https://mutalyzer.nl/json/" #Mutalyzer's REST API
  server.mutalyzerv3 <- "https://v3.mutalyzer.nl/api/"#Mutalyzer's V3 REST API
  
  ####we encode the URL
  variant.mutalyzer <- ifelse(intronic, 
                              URLencode(paste0(NC, "(",NM, "):",variant),reserved=TRUE),
                              paste0(NM, ":",variant))
  
  
  # getting information about the variant from mutalyzer
  #Non intronic variants
  if (intronic == FALSE){
    ext.mutalyzer.run <-paste0("runMutalyzer?variant=", variant.mutalyzer)
    mutalyzer <- api2(server.mutalyzer,ext.mutalyzer.run)
    
    
    assertthat::assert_that(!any(mutalyzer$messages$errorcode=="EPARSE"), msg="Mutalyzer could not retrieve NM")
    assertthat::assert_that(!any(mutalyzer$messages$errorcode=="ERETR"), msg="try another NM Ref-Seq version")
    assertthat::assert_that(!any(mutalyzer$messages$errorcode=="EREF"), msg = mutalyzer$messages$message)    
    
    warnings <- mutalyzer$warnings
    errors <- mutalyzer$errors
    cor.variant <- stringr::str_split(mutalyzer$transcriptDescriptions,":")[[1]][2]
    html.prot <-  mutalyzer$proteinDescriptions[1]
    message.mutalyzer <- ifelse(warnings+errors>0,
                                paste(unlist(mutalyzer$messages[1]),unlist(mutalyzer$messages[2]), sep=":" ),
                                "No warnings found")
    exons.mut <- mutalyzer$exons
    mutalyzer.prot.pred <- mutalyzer$origProtein
    
  }else{ #Intronic variants are searched in mutalyzer v3 (alpha) because the other doesn't work properly for these variants
    assertthat::assert_that(!is.null(NC), msg="'NC' argument must be given")
    
    #Changing between equivalent versions, to avoid mutalyzer errors.
    #MSH2
    if(NM=="NM_000251.2")variant.mutalyzer <- ifelse(intronic, 
                                                     URLencode(paste0(NC, "(NM_000251.3):",variant),reserved=TRUE),
                                                     paste0(NM, ":",variant))
    #NTHL1
    if(NM=="NM_002528.5")variant.mutalyzer <- ifelse(intronic,
                                                     URLencode(paste0(NC, "(NM_002528.3):",variant),reserved=TRUE),
                                                     paste0(NM, ":",variant))
    #BRCA1
    if(NM=="NM_007294.3")variant.mutalyzer <- ifelse(intronic, 
                                                     URLencode(paste0(NC, "(NM_007294.4):",variant),reserved=TRUE),
                                                     paste0(NM, ":",variant))
    #PMS2
    if(NM=="NM_000535.5")variant.mutalyzer <- ifelse(intronic, 
                                                     URLencode(paste0(NC, "(NM_000535.7):",variant),reserved=TRUE),
                                                     paste0(NM, ":",variant))
    #CDH1
    if(NM=="NM_004360.3")variant.mutalyzer <- ifelse(intronic, 
                                                     URLencode(paste0(NC, "(NM_004360.5):",variant),reserved=TRUE),
                                                     paste0(NM, ":",variant))
   #BRIP1
     if(NM=="NM_032043.2")variant.mutalyzer <- ifelse(intronic, 
                                                     URLencode(paste0(NC, "(NM_032043.3):",variant),reserved=TRUE),
                                                     paste0(NM, ":",variant))
     #UNC93B1
     if(NM=="NM_030930.2")variant.mutalyzer <- ifelse(intronic, 
                                                    URLencode(paste0("NG_007581.1", "(", NM,"):",variant),reserved=TRUE),
                                                    paste0(NM, ":",variant))
     #TRAF3
     if(NM=="NM_145725.2")variant.mutalyzer <- ifelse(intronic, 
                                                      URLencode(paste0(NC, "(NM_145725.3):",variant),reserved=TRUE),
                                                      paste0(NM, ":",variant))
     #TYK2
     if(NM=="NM_003331.4")variant.mutalyzer <- ifelse(intronic, 
                                                      URLencode(paste0(NC, "(NM_003331.5):",variant),reserved=TRUE),
                                                      paste0(NM, ":",variant))
     
    #ext.mutalyzer.v3 <- paste0("name_check/", variant.mutalyzer)
    ext.mutalyzer.v3 <- paste0("normalize/", variant.mutalyzer)
    mutalyzerv3 <- api2(server.mutalyzerv3, ext.mutalyzer.v3)
    assertthat::assert_that(!any(mutalyzerv3$errors$code=="ERETR"), msg = mutalyzerv3$errors$details)
    assertthat::assert_that(!any(mutalyzerv3$errors$code=="ESEQUENCEMISMATCH"), msg= mutalyzerv3$errors$details)
    
    cor.variant <- stringr::str_split(mutalyzerv3$normalized_description, ":")[[1]][2]
    html.prot <- mutalyzerv3$protein$description
    mutalyzer.prot.pred <- mutalyzerv3$protein$predicted
    exons.mut <- cbind(tibble::as_tibble(mutalyzerv3$selector_short$exon$g),tibble::as_tibble(mutalyzerv3$selector_short$exon$c))
    names(exons.mut) <- c("gStart", "gStop", "cStart", "cStop")
    message.mutalyzer <- ifelse(mutalyzerv3$corrected_description==mutalyzerv3$normalized_description,
                                "No errors found",
                                "The variant's nomenclature was not ok")
  }
  cor.prot <- stringr::str_split(html.prot,"\\:")
  if (length(cor.prot)==0){
    cor.prot <- "p.?"
  }else{
    cor.prot<-cor.prot[[1]][2]
    prot2<-cor.prot
  }
  
  
  ###Genomic coordinates
  ext.mut.convert <- paste0("numberConversion?build=hg19&variant=", NM,":", cor.variant, "&gene=", gene)
  #TLR7 gene equivalent versions
  if(NM=="NM_016562.4")ext.mut.convert <- paste0("numberConversion?build=hg19&variant=NM_016562.3:", cor.variant, "&gene=", gene)
  if(NM=="NM_001572.5")ext.mut.convert <- paste0("numberConversion?build=hg19&variant=NM_001572.3:", cor.variant, "&gene=", gene)
  # STAT2
  if(NM=="NM_005419.4")ext.mut.convert <- paste0("numberConversion?build=hg19&variant=NM_005419.3:", cor.variant, "&gene=", gene)
  #TBK1
  if(NM=="NM_013254.4")ext.mut.convert <- paste0("numberConversion?build=hg19&variant=NM_013254.3:", cor.variant, "&gene=", gene)
  
 
  ext.mut.convert <- gsub("\\+", "%2B", ext.mut.convert)
  mutalyzer.genomic <- api2(server.mutalyzer, ext.mut.convert )
  mutalyzer.other.selected <- list()
  
  
  ###Other relevant transcripts
  if (gene=="TP53"){  #it has more important transcripts
    ext.mut.convert2 <- paste0("numberConversion?build=hg19&variant=", mutalyzer.genomic ,"&gene=", gene)
    mutalyzer.other.trans <- api2(server.mutalyzer, ext.mut.convert2 ) 
    mutalyzer.other.selected <- mutalyzer.other.trans[stringr::str_detect(mutalyzer.other.trans, "NM_000546.5|NM_001126114.2|NM_001126113.2")==TRUE & stringr::str_detect(mutalyzer.other.trans, NM)==FALSE]
  }else if( gene== "CDKN2A"){
    ext.mut.convert2 <- paste0("numberConversion?build=hg19&variant=", mutalyzer.genomic ,"&gene=", gene)
    mutalyzer.other.trans <- api2(server.mutalyzer, ext.mut.convert2 ) 
    mutalyzer.other.selected <- mutalyzer.other.trans[stringr::str_detect(mutalyzer.other.trans, "NM_000077.4|NM_058195.3")==TRUE & stringr::str_detect(mutalyzer.other.trans, NM)==FALSE]
  }
  
  ##final output
  correct.variant <- list(initial.variant = variant,
                          NM = NM, 
                          gene = gene, 
                          variant = cor.variant, 
                          protein = cor.prot, 
                          warning = message.mutalyzer, 
                          genomic = as.character(mutalyzer.genomic), 
                          exons = exons.mut, 
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
#' @param 
#' @param variant.mutalyzer By default is null. Obtained with the correctHgvsMutalyzer function contains all information of the variant from mutalyzer.
#' @param skip.pred By default is false. If you are predicting the consequence of a skipping then it should be set to TRUE
#' @return The chromosome, the initial genomic coordinate, the final genomic
#' coordinate, the reference allele, the alternative allele, the most important
#' variant consequence and if the variant is located in a protein domain
#' @author Elisabet Munté Roca
#' @references McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. Genome Biology Jun 6;17(1):122. (2016) doi:10.1186/s13059-016-0974-4
#' @examples
#' NM <- "NM_007294.3"
#' NC <- "NC_000017.10"
#' gene <- "BRCA1"
#' CCDS <- "CCDS11453.1"
#' variant <- "c.211A>G"
#' mutalyzer.info <- correctHgvsMutalyzer(NM, NC,  gene, variant)
#' varDetails(NM,NC, CCDS, gene, mutalyzer.info, skip.pred=FALSE)


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
  ext.vep <-  paste0("/vep/human/hgvs/",NM,":",variant.mutalyzer$variant)
  coordinates <- api2(server.ensembl, ext.vep)
  assertthat::assert_that(length(coordinates$error)==0, msg=coordinates$error)
  chr<-unlist(coordinates$seq_region_name)
  strand <- unlist(coordinates$strand)
  
  ##to get the correct start and end
  ext.vep2 <- paste0("/vep/human/hgvs/",chr,":",unlist(stringr::str_split(variant.mutalyzer$genomic, ":"))[2])
  
  #we get the information
  coordinates2 <- api2(server.ensembl, ext.vep2)
  assertthat::assert_that(length(coordinates2$error)==0, msg=coordinates2$error)
  end <-coordinates2$end[[1]]#coordenada final
  start <-coordinates2$start[[1]]#coordenada inicio
  allele.string <- coordinates2$allele_string[[1]] %>% 
    stringr::str_split("/") %>% unlist
  ref <- allele.string[1]#reference
  alt <- allele.string[2]#alternative
  most.severe.consequence <- tibble::tibble(coordinates2$transcript_consequences[[1]]) %>%
    dplyr::filter (transcript_id==as.character(ensembl.id)) %>% 
    dplyr::select(consequence_terms) %>% 
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
      most.severe.consequence <- most.severe.consequence[num]
    }
  }
  object <- data.frame(NM=NM, variant=variant, protein=variant.mutalyzer$protein, strand=strand, gene=gene, most.severe.consequence=most.severe.consequence[1])
  exons <- coordNonCoding(variant.mutalyzer, object) %>% 
    dplyr::mutate(cStop)
  
  if (!(most.severe.consequence[1] %in% c("3_prime_UTR_variant", "5_prime_UTR_variant"))){
    if(stringr::str_detect(variant, "-")){ 
      pos.var <- as.numeric(stringr::str_extract(variant, "[0-9]+"))-1
    }else{
      pos.var <- stringr::str_extract(variant, "[0-9]+")
    }
    variant.exon <- exons %>% 
      dplyr::mutate (cStop = stringr::str_replace(cStop, "\\*[0-9]*", "1000000000"), cStart = stringr::str_replace(cStart, "\\*[0-9]*", "NA")) %>% 
      dplyr::filter( as.numeric(cStart) <= as.numeric(pos.var), as.numeric(cStop) >= as.numeric(pos.var)) %>% 
      dplyr::select(exon) 
  }else{
    variant.exon <- NA
  }
  #adjust synonymous nomenclature
  protein <- protsyn2(object, variant.mutalyzer)
  if(any( stringr::str_detect(most.severe.consequence, "frameshift_variant"))){
    most.severe.consequence <- ifelse(rep(most.severe.consequence=="frameshift_variant" && length(stringr::str_extract_all(variant.mutalyzer$protein, "[0-9]+"))==1,3), "stop_gained", most.severe.consequence) #sometimes frameshift variants are nonsense
    most.severe.consequence <- ifelse(rep(stringr::str_detect(variant.mutalyzer$protein, "fs"),3), "frameshift_variant", most.severe.consequence)
  }
  domain.info <- domainProt(chr, start, end) %>% 
                 unique
  uvs.df <- tibble::tibble (gene = variant.mutalyzer$gene,
                            NM= NM,  
                            initial.var = variant,
                            variant = variant.mutalyzer$variant, 
                            protein = protein,
                            genomic = variant.mutalyzer$genomic, 
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
  if(nrow(lrg.gene)==0){
    query <- paste0("SELECT l.transcript_id, l.coordinates, l.strand  
                  FROM transcript t INNER JOIN noLRG l ON t.ensembltranscriptID=l.transcript_id 
                  WHERE l.namegene= '",object$gene ,"' AND t.NM='", object$NM, "'; ")  
    lrg.gene <- connectionDB(query)[[1]] 
    names(lrg.gene)[1]<-"transcript"
  }
  exon.nocodi <- unlist(stringr::str_split(lrg.gene["coordinates"], ",")) %>% 
                 stringr::str_split(pattern="-") 
  num.exons <- length(exon.nocodi)
  exons.ord <- ifelse(rep(object$strand ==1, num.exons), 1:num.exons, num.exons:1)
  coordinates.exon.nocodi <- data.frame(exon = exons.ord, 
                                        V1 = unlist(purrr::map(exon.nocodi,1)), 
                                        V2=unlist(purrr::map(exon.nocodi, 2)))  %>%
                                        dplyr::mutate(transcript = lrg.gene$transcript) %>% 
                                        dplyr::arrange(exon) %>% 
                                        dplyr::relocate (transcript)
  return (coordinates.exon.nocodi)
}

coordNonCoding <- function (variant.mutalyzer, object){
  coordinates.exon.nocodi<- geneLrgCoord(object)
  query <- paste0("SELECT exon, cStart, cStop FROM LRG_cds WHERE LRG_id= '", coordinates.exon.nocodi$transcript[1] ,"'; ")
  lrg.cds <- connectionDB(query)[[1]]
  
  if(nrow(lrg.cds)>0){
    all.exons <- merge(coordinates.exon.nocodi, lrg.cds, by="exon")
  }else if (object$gene=="APC"){
    all.exons <- cbind(coordinates.exon.nocodi, rbind(c("",""), variant.mutalyzer$exons[,c("cStart", "cStop")]))
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
#' @example
#' varDetails(17, 41258473,41258474)
#' @author Elisabet Munté Roca
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
    web.mutalyzer <- xml2::read_html(utils::URLencode(paste0("https://mutalyzer.nl/name-checker?description=",object$genomic)))
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
    if (object$strand==1){ 
      if(stringr::str_detect(object$variant, "del")){
        gnomAD<-c(object$chr,start2,ref, alt)#we obtain the gnomAD
      }else{
        gnomAD<-c(object$chr,start2,ref,alt)#we obtain the gnomAD
      }
    }else if( object$strand == -1){
      seq1 <- DescTools::StrRev(seq1) #we need the reverse sequence
      if(stringr::str_detect(object$variant, "del")){
        gnomAD<- c(object$chr,object$start-1, ref, alt)
      }else{
        gnomAD<- c(object$chr,object$end-stringr::str_length(object$alt),ref, alt)
      }
    }
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
    pos.dna<-stringr::str_extract(object$variant, "[0-9]+")
    pos.aa<-as.integer(as.numeric(pos.dna)/3)
    if (as.numeric(pos.dna) %% 3 >0){ #if the division is not exact, we need to add 1 position
      pos.aa <- pos.aa + 1  
    }
    seq.prot<-stringr::str_sub(var.mut$protein_predicted, pos.aa, pos.aa)
    aa <- Biostrings::AMINO_ACID_CODE[seq.prot]
    protein <- paste0("p.(", aa, pos.aa, "=)")
  }
  prot <- ifelse(object$most.severe.consequence %in% c("synonymous_variant"), protein, object$protein)
  return(prot)
}

toGenomic <- function(NM, cor.variant, gene){
  server.mutalyzer <- "https://mutalyzer.nl/json/"
  ext.mut.convert <- paste0("numberConversion?build=hg19&variant=", NM, ":", cor.variant, "&gene=", gene)
  ext.mut.convert <- gsub("\\+", "%2B", ext.mut.convert)
  mutalyzer.genomic <- api2(server.mutalyzer, ext.mut.convert )
  return(mutalyzer.genomic)
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