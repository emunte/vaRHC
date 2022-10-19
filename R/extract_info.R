
#' vaRinfo ()
#' @description  Collect different types of criteria used by ACMG to classify variants.
#' @param assembly character hg19 or hg38
#' @param gene gene of interest
#' @param variant variant of interest in cdna
#' @param NM refSeq nomenclature. By default is NULL and vaRHC will consider the ones detailed in README file. Be careful if you use a different NM because the program has not been validated for it. If you provide a different NM,  NC and CCDS must also be provided.
#' @param NC refSeq nomenclature. By default is NULL and vaRHC will consider the ones detailed in README file. Be careful if you use a different NM because the program has not been validated for it. If you provide a different NC, NM and CCDS must also be provided.
#' @param CCDS Consensus CD id https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi. By default is NULL and vaRHC will consider the ones detailed in README file. Be careful if you use a different CCDS because the program has not been validated for it. If you provide a different CCDS, NM and NC must also be provided. Current version only works for hg19.
#' @param gene.specific.df By default is NULL, it uses the default parameters described in README. If you would like to change some defaults or include another gene, a template can be downloaded from Github: https://github.com/emunte/Class_variants/tree/main/documents/gen_especific.csv and some parameters can be modified taking into account your preferences
#' @param browser Which browser to start Rselenium server. By default is "firefox" (the recommended). If you do not have firefox installed try either "chrome" or "phantomjs".
#' @param spliceai.program Logical. By default is FALSE and it is assumed that SpliceAI program is not installed in your computer. If this parameter is FALSE, the program will only classify substitutions and simple deletion variants taking into account a spliceAI distance of 1000 and will show masked results. If you want to classify other variants please install SpliceAI (https://pypi.org/project/spliceai/) and set to TRUE the parameter.
#' @param spliceai.reference Path to the Reference genome hg19 fasta file. Can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz . By default is NULL and it will only be taken into account if spliceai.program is set to TRUE.
#' @param spliceai.annotation Path to gene annotation file. By default it uses the file stored in extdata folder: "../extdata/gencode_spliceai_hg19.txt"
#' @param spliceai.distance  Integer. Maximum distance between the variant and gained/lost splice site (default: 1000)
#' @param spliceai.masked Mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss (default: 1)
#' @param provean.program Logical. By default is FALSE and it is assumed that provean program is not installed in your computer.
#' @return information about the variant.
#' @author Elisabet Munté Roca
#' @examples
#' vaRinfo("hg19", "MLH1", "c.1AG>")
#' vaRinfo("hg19", "MSH6", "c.211A>G", spliceai.program = TRUE, spliceai.reference = "./hg19.fa")
#' @references
#' Richards, S., Aziz, N., Bale, S., Bick, D., Das, S., Gastier-Foster, J., Grody, W. W., Hegde, M., Lyon, E., Spector, E., Voelkerding, K., Rehm, H. L., & ACMG Laboratory Quality Assurance Committee (2015). Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genetics in medicine : official journal of the American College of Medical Genetics, 17, 405–424. https://doi.org/10.1038/gim.2015.30
vaRinfo <- function(assembly,  gene, variant, NM=NULL, NC = NULL, CCDS=NULL, gene.specific.df=NULL, browser="firefox",  spliceai.program=FALSE, spliceai.reference=NULL, spliceai.annotation =  system.file("extdata", "gencode_spliceai_hg19.txt", package="vaRHC"), spliceai.distance=1000, spliceai.masked=1, provean.program=FALSE){
  nm.nc <- NMparam(gene, NM = NM, NC = NC, CCDS = CCDS)
  cat("0% completed... correcting variant nomenclature \n")
  variant.mutalyzer <- correctHgvsMutalyzer (NM = nm.nc$NM, NC = nm.nc$NC, gene = gene, variant = variant)
  cat("10% completed... getting variant coordinates\n")
  variant.info <- varDetails(NM=nm.nc$NM, NC=nm.nc$NC, CCDS=nm.nc$CCDS, gene = gene, variant = variant, variant.mutalyzer, skip.pred=FALSE)
  variant.info.other <- list()
  if (!is.na(variant.mutalyzer$other.important.transcripts[1]) && length(variant.mutalyzer$other.important.transcripts)>0){
    NM.other <- stringr::str_split(variant.mutalyzer$other.important.transcripts, ":") %>% purrr::map(1) %>% unlist()
    variant.other <- stringr::str_split(variant.mutalyzer$other.important.transcripts, ":") %>% purrr::map(2) %>% unlist()
    for (i in 1:length(NM.other)){
      cat(paste("getting information from other transcripts:", NM.other[i], "\n"))
      query <- paste0("SELECT * from  transcript WHERE namegene= '", gene ,"' AND NM = '", NM.other[i], "' ;")
      CCDS.other <- connectionDB(query)[[1]] %>%
                    tibble::as_tibble() %>%
                    dplyr::select (CCDS)
      variant.info.other[[NM.other[i]]]<- varDetails(NM.other[i], NC=nm.nc$NC, CCDS=CCDS.other, gene, variant.other[i])
      }
  }

  cat("15% completed ... extracting bbdd information\n")
  gnom <- gnomADnomen(object = variant.info)
  bbdd.info <- extractBBDD(mutalyzer = variant.mutalyzer, object = variant.info, gnom = gnom) %>% connectionDB()
  if(is.null(gene.specific.df)) gene.specific.df <- bbdd.info$gene.specific
  gene.specific <- geneSpecific(gene = gene, gene.specific.df= gene.specific.df)

  cat("20% completed ... getting information from gnomAD\n")

  exomes <- gnomADcov(assembly=assembly, track="Exomes", gnomAD.nomenclature=gnom)
  genomes <- gnomADcov(assembly=assembly, track="Genomes", gnomAD.nomenclature=gnom)
  gnomad.info <- gnomADmerge(gnom, gene.specific, bbdd= bbdd.info)

  cat ("30% completed ... getting information from flossiesdb\n")
  flossies <- flossiesInfo(nomen.flossies= NULL, nomen.gnomAD=gnom, gene=variant.info$gene)
  cat ("40% completed ... getting clinVar information\n")
  clinvar.ids <- clinVarIds(object = variant.info)
  clinvar.info <- clinVarInfo(clinvar.ids, object= variant.info)
  cat("50% completed ...  getting predictors information\n")
  predictor.table <- predicInfo(object = variant.info,
                                gene.specific =  gene.specific,
                                bbdd = bbdd.info,
                                gnomad = gnom,
                                spliceai.program=spliceai.program,
                                spliceai.reference = spliceai.reference,
                                spliceai.annotation =  spliceai.annotation,
                                spliceai.distance = spliceai.distance,
                                spliceai.masked = spliceai.masked,
                                provean.program = provean.program
                                )
  cat("60% completed ... getting loss of function information\n")
  codon.stop <- stopCodon(object = variant.info, variant.mutalyzer = variant.mutalyzer)

  cat("70% completed ... getting start codon information\n")
  second.met <- secondMet(variant.mutalyzer = variant.mutalyzer, object = variant.info, assembly=assembly)
  cat("80% completed ... getting functional studies information\n")
  insight.info <- insightInfo(variant.info, browser)
  articles <- articlesInfo(object = variant.info, variant.mutalyzer, bbdd = bbdd.info)
  references.google.30 = biblioScholar(variant.info, c(0,10,20))
  cat("90% completed ... getting cancerhotspots information\n")
  cancer.hotspots <- cancerHotspots(object = variant.info, variant.cor = variant.mutalyzer, bbdd = bbdd.info)
  cat("100% completed!!\n")
  out<- list(Variant.Info = variant.info, Variant.Info.other= variant.info.other,
             variant.correction = variant.mutalyzer,
             gene.specific.info = gene.specific,
             gnomAD = list(
               nomenclature=gnom,
               coverage=list(exomes=exomes, genomes=genomes),
               info= gnomad.info),
             flossies.db=flossies,
             clinVar = list(clinVar.ids=clinvar.ids, clinVar.info=clinvar.info),
             predictors= list(predictor.table=predictor.table) ,
             codon.stop=codon.stop,
             second.met = second.met,
             insight.info = insight.info,
             functional.assays = articles,
             google.scholar.30.references = references.google.30,
             cancer.hotspots= cancer.hotspots,
             class.info = list(last.nt = bbdd.info$last.nt,
                               pvs1.cdh1 = bbdd.info$pvs1.cdh1,
                               pvs1.atm = bbdd.info$pvs1.atm
             ))
  return(out)
}



#' @noRd
extractBBDD <- function(mutalyzer, object, gnom){
  #query IDIBELL database
  #prepare input
  prot <- toProtein(object$protein)
  pos.nt <- stringr::str_extract(object$variant, "[0-9]+")

  #gene specific
  gene.specific = paste0("SELECT * from  gene_specific;")

  #predictors

  ###prior
  prior.db <- paste0("SELECT p.refsplice_prior, p.splice_severity, p.de_novo_prior, p.dn_severity, p.protein_prior, p.applicable_prior, b.Polyphen, b.MAPP, b.prior from prior_db_all p LEFT JOIN prior_db b  ON (p.Gene=b.gene AND p.cdna=b.cdna) WHERE p.Gene= '", object$gene ,"' AND p.cdna='", object$variant, "';")
  ###alignGvgd (only for TP53)
  align.gvgd <- paste0("SELECT Prediction from align_gvgd WHERE Substitution= '", paste0(aaShort(unlist(prot$aa.ref)), unlist(prot$aa.pos),aaShort(unlist(prot$aa.alt)) ),"' AND namegene='", object$gene,"';")
  ###dbnsfp
  dbnsfp <- paste0("SELECT chr, start, ref , alt, PROVEAN_score, VEST4_score, REVEL_score, BayesDel_noAF_score from dbnsfp_all WHERE start= '", object$start ,"'AND chr='", object$chr,"'AND ref='",object$ref,"'AND alt='", object$alt,"';")
  ####spliceAI
  spliceai <- paste0("SELECT *  from spliceAI  WHERE var_chr= '", gnom ,"'AND max_dis= 1000 AND transcript='", object$ensembl.id,"' AND masked='",TRUE,"';")

  #geneLrgCoord
  #gene.LRG <- paste0("SELECT  l.transcript, l.namegene,l.coordinates, l.transcript2, l.cds_start, l.cds_end, l.strand,  c.exon, c.cStart, c.cStop FROM LRG l LEFT JOIN  transcript t ON t.ensembltranscriptID=l.transcript_id LEFT JOIN LRG_cds c ON l.transcript = c.LRG_id WHERE l.namegene= '",object$gene ,"' AND t.NM='", object$NM, "'; ")

  #no lrg
  # genenoLRG<- paste0("SELECT l.transcript_id, l.coordinates, l.strand
  #                 FROM transcript t INNER JOIN noLRG l ON t.ensembltranscriptID=l.transcript_id
  #                 WHERE l.namegene= '",object$gene ,"' AND t.NM='", object$NM, "'; ")

  #cancer hotspots
  cancer.hotspots <- paste0("SELECT * from cancer_hotspotsv2 WHERE genename= '", object$gene, "';")

  #lyra
  lyra <- "SELECT * from lyra_db;"

  #Adamovich
  short.prot <- paste0(aaShort(prot$aa.ref),prot$aa.pos,aaShort(prot$aa.alt))
  adamovich <- paste0("SELECT * from Adamovich_Fayer_db WHERE variantID='", short.prot , "' OR transcript_variant='", object$variant, "';")

  #Cimra

  cimra1 <- paste0("SELECT * from cimra_db WHERE gene= '", object$gene,  "' AND DNA='", object$variant, "';")
  cimra2 <- paste0("SELECT * from cimra_db WHERE gene= '", object$gene,  "' AND Protein='", paste0("p.", short.prot), "';")
  #jia
  jia <- "SELECT * from jia_db;"
  #parsons
  parsons <- paste0( "SELECT * from parsons_db WHERE Gene= '", object$gene, "' AND HGVS_Nucleotide ='", object$variant, "';")

  #funcionals chek2
  chek2 <- paste0("SELECT * from functionals_CHEK2 WHERE cDNA_variant='", object$variant, "';")

  #TP53 funcionals
  if (object$gene =="TP53"){
  tp53.1 <- paste0("SELECT * from functionals_TP53 WHERE c_description='", object$variant, "';")
  pos.dna<-stringr::str_extract(object$variant, "[0-9]+") %>% as.numeric/3 %>% as.integer()
  if (!(object$most.severe.consequence %in% c("frameshift_variant", "inframe_deletion", ""))){
    tp53.2 <- paste0("SELECT * from functionals_TP53_kotler WHERE  AA_WT='", aaShort(prot$aa.ref), "' AND AA_mutant='", aaShort(prot$aa.alt), "' AND Codon_num='", prot$aa.pos, "';")
    }else if(object$most.severe.consequence %in% c("frameshift_variant")){
    tp53.2 <- paste0("SELECT * from functionals_TP53_kotler WHERE AA_mutant= 'fs' AND Codon_num='", pos.dna, "';")
    }
  }else{
    tp53.1 <- NULL
    tp53.2 <- NULL
  }


  #pten
  pten <- paste0("SELECT * from rna_PTEN WHERE variant='", object$variant, "';")

  #ATM
  atm <- paste0("SELECT * from functionals_ATM WHERE variant='", object$variant, "';")

  #last nt
  last.nt <- paste0("SELECT * from last_nt_exons WHERE gene='", object$gene,"' AND pos=", as.integer(pos.nt))

  #PM5
  pos <- stringr::str_extract(object$variant, "[0-9]+")
  #PVS1 cdh1
  pvs1.cdh1 <- paste0("SELECT PVS1 FROM canonicals WHERE location='c.", pos.nt ,"'")
  #PVS1 ATM
  pvs1.atm <- paste0("SELECT PVS1_strength, reasoning FROM canonicals_ATM WHERE variant='", object$variant,"'")

  #gnomad
  exomes.gnomad <- queriesGnomad("exomes", gnom)
  genomes.gnomad <- queriesGnomad("genomes", gnom)

  queries <- list(gene.specific = gene.specific,
                  prior = prior.db,
                  align.gvd = align.gvgd,
                  dbnsfp = dbnsfp,
                  spliceai = spliceai,
                  cancer.hotspots = cancer.hotspots,
                  lyra = lyra,
                  adamovich = adamovich,
                  cimra1 = cimra1,
                  cimra2 = cimra2,
                  jia = jia,
                  parsons = parsons,
                  chek2 = chek2,
                  tp53.1 = tp53.1,
                  tp53.2 = tp53.2,
                  pten = pten,
                  atm = atm,
                  last.nt = last.nt,
                  pvs1.cdh1 = pvs1.cdh1,
                  pvs1.atm = pvs1.atm,
                  gnomad1.exomes = exomes.gnomad[1],
                  gnomad2.exomes = exomes.gnomad[2],
                  gnomad.up.exomes = exomes.gnomad[3],
                  gnomad.down.exomes = exomes.gnomad[4],
                  gnomad1.genomes = genomes.gnomad[1],
                  gnomad2.genomes = genomes.gnomad[2],
                  gnomad.up.genomes = genomes.gnomad[3],
                  gnomad.down.genomes = genomes.gnomad[4]
  )
  return(queries)
}

#' @noRd
queriesGnomad <- function (track, gnom){
  gnomad <- stringr::str_split(gnom, "-")[[1]]
  gnomad1 <- paste0("SELECT ", track, "_id from ", track, "_gnomad WHERE POS= '", gnomad[2] ,"' AND CHROM='", gnomad[1],"' AND ref='",gnomad[3],"' AND alt='", gnomad[4],"' ;")
  gnomad2 <- paste0("SELECT ", track, "_id from ", track, "_gnomad WHERE POS= '", gnomad[2] ,"' AND CHROM='", gnomad[1],"';")
  gnomad.up <- paste0("SELECT ", track, "_id from ", track, "_gnomad WHERE POS> '", gnomad[2] ,"'AND POS<'", as.numeric(gnomad[2])+50, "'AND CHROM='", gnomad[1],"';")
  gnomad.down <- paste0("SELECT ", track, "_id from ", track, "_gnomad WHERE POS> '", as.numeric(gnomad[2])-50 ,"'AND POS<'", gnomad[2], "'AND CHROM='", gnomad[1],"';")
  gnomad1.name = paste0("gnomad1.", track)
  gnomad2.name = paste0("gnomad2.", track)
  gnomad.up.name = paste0("gnomad.up.", track)
  gnomad.down.name = paste0("gnomap.down.", track)
  gnomad.all <- c(gnomad1,gnomad2,gnomad.up,gnomad.down)
  return(gnomad.all)
}

#' Know if the gene has special ACMG rules
#'
#' @param gene the gene of interest
#' @param gene.specific.df By default is NULL, it uses the default parameters described in README. If you would like to change some defaults or include another gene, a template can be downloaded from Github: https://github.com/emunte/Class_variants/tree/main/documents/gen_especific.csv and some parameters can be modified taking into account your preferences
#' @return The name of the gene if there are special ACMG rules for it and if not it returns general ACMG rules.
#' if it follows general ACMG rules (only for hereditary cancer related genes)
#' @author Elisabet Munté Roca
#' @examples
#' #using the default
#' geneSpecific("BRCA1")
#' geneSpecific("RAD51C")
#' geneSpecific("ERCC5")
#'
#' #changing the csv file
#' file1 <- file.path("./gene_specific.csv")
#' gene.specific.df <- read.csv(file1, sep=",", header=TRUE,row.names=1, stringsAsFactors = FALSE)
#' geneSpecific("BRCA1", gene.specific.df)
#' geneSpecific("RAD51C", gene.specific.df)
#' geneSpecific("ERCC5", gene.specific.df)
#' @references
#' Richards, S., Aziz, N., Bale, S., Bick, D., Das, S., Gastier-Foster, J., Grody, W. W., Hegde, M., Lyon, E., Spector, E., Voelkerding, K., Rehm, H. L., & ACMG Laboratory Quality Assurance Committee (2015). Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genetics in medicine : official journal of the American College of Medical Genetics, 17(5), 405–424. https://doi.org/10.1038/gim.2015.30
#' Tavtigian, S. V., Harrison, S. M., Boucher, K. M., & Biesecker, L. G. (2020). Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. Human mutation, 41(10), 1734–1737. https://doi.org/10.1002/humu.24088
#' Garrett, A., Durkie, M., Callaway, A., Burghel, G. J., Robinson, R., Drummond, J., Torr, B., Cubuk, C., Berry, I. R., Wallace, A. J., Ellard, S., Eccles, D. M., Tischkowitz, M., Hanson, H., Turnbull, C., & CanVIG-UK (2021). Combining evidence for and against pathogenicity for variants in cancer susceptibility genes: CanVIG-UK consensus recommendations. Journal of medical genetics, 58(5), 297–304. https://doi.org/10.1136/jmedgenet-2020-107248
#' ClinGen CDH1 Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 3: https://clinicalgenome.org/site/assets/files/7118/clingen_cdh1_acmg_specifications_v3.pdf
#' ClinGen PTEN Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 2: https://clinicalgenome.org/site/assets/files/4000/clingen_pten_acmg_specifications_v2.pdf
#' ClinGen InSiGHT Hereditary Colorectal Cancer/Polyposis Variant Curation Expert Panel Specifications to the ACMG/AMP Variant Interpretation Guidelines Version 1 (draft): https://www.insight-group.org/content/uploads/2021/11/DRAFT_Nov_2021_TEMPLATE_SVI.ACMG_Specifications_InSiGHT_MMR_V1.pdf

geneSpecific <- function (gene, gene.specific.df) {
  assertthat::assert_that(!is.null(gene.specific.df), msg="You must provide either the information collected with extractBBDD()$gene.specific or your own specific file")
  rownames(gene.specific.df) <- gene.specific.df$genename
  gene.specific.df <- gene.specific.df[2:ncol(gene.specific.df)]
  gn <- gene %in% row.names(gene.specific.df)
  gene.sp <- ifelse(gn,
                    gene,
                    "general")
  gene.sp.info <- gene.specific.df[gene.sp, ]
  return(gene.sp.info)
}

################################################################################
## GnomAD
################################################################################

#' gnomADcov
#' @param assembly character hg19 or hg38
#' @param track character.It can be either Exomes or Genomes. GnomAD
#' @param gnomAD.nomenclature character. GnomAD variant nomenclature. If ID is not known it can be set to NA. Chr, start and end must be given then
#' @param chr numeric. Chromosome of interest. By default is NULL
#' @param start numeric.Start of the region of interest. By default is NULL
#' @param end numeric. End of the region of interest. By default is NULL
#' @param  mean.value logical. By default is TRUE. If the region is more than 1 bp long it returns the mean coverage or if FALSE it returns the coverage of each position
#' @return Numeric value of the coverage region
#' @author Elisabet Munté Roca
#' @references
#' Karczewski, K.J., Francioli, L.C., Tiao, G. et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 581, 434–443 (2020). https://doi.org/10.1038/s41586-020-2308-7
#' @examples
#' gene <- "BRCA1"
#' variant <- "c.211A>G"
#'nm.nc <- NMparam(gene, nm.info)
#'variant.correction <- correctHgvsMutalyzer (nm.nc$NM, nm.nc$NC, gene, variant)
#'variant.info<- varDetails(NM=nm.nc$NM, NC=nm.nc$NC, CCDS=nm.nc$CCDS, gene, variant.correction, skip.pred=FALSE)
#'gnom<- gnomADnomen(object = variant.info)
#'exomes <- gnomADcov(assembly=assembly, track="Exomes", gnomAD.nomenclature=gnom)
gnomADcov <- function(assembly, track, gnomAD.nomenclature, chr=NULL, start=NULL, end=NULL, mean.value=TRUE){
  if (!(assembly%in% c("hg19","hg38"))) stop ("this assembly is not supported, please enter hg19 or hg38")
  if(assembly== "hg38") stop("We are sorry but right now only hg19 works, we will incorporate hg38 soon")
  if (!is.null(gnomAD.nomenclature)&!is.na(gnomAD.nomenclature)){
    gnomad <- stringr::str_split(gnomAD.nomenclature, "-")[[1]]
    if (length(gnomad)!=4) stop("GnomAD ID nomenclature is not correct, please check")
    chr <- gnomad[1]
    start <- gnomad[2] %>% as.numeric()
    diference <- abs(stringr::str_length(gnomad[3])-abs(stringr::str_length(gnomad[4])))
    if (diference==0){
      end <- start
      start <- start-1
    }else{
      end <- start+diference
    }
  } else {
    if (is.null(chr)|is.null(start)|is.null(end)) stop("No GnomADID and no chr, start and end has been provided, please check")
    if(start==end){
      start <- start-1
    }
  }
  chr.gnom<- paste0("chr",chr)
  ext.coverage<-paste0("genome=", assembly, ";track=gnomad",track,  "MeanCoverage;chrom=", chr.gnom, ";start=",start,";end=",end)
  server.ucsc <- "http://api.genome.ucsc.edu/getData/track?" #UCSC's REST API
  coverage <- api2(server.ucsc, ext.coverage)
  ge <- paste0("gnomad", track, "MeanCoverage")
  cover <- coverage[[ge]][chr.gnom] #we check the coverage for the specific genome and chromosome
  value.cov <- cover[[chr.gnom]]["value"]
  if (mean.value==TRUE &!is.null(value.cov)){
    mean.cov <- mean(unlist(value.cov$value))
    return(mean.cov)
  }else{
    value.cov[is.null(value.cov)]<- 0
    return(value.cov)
  }
}

#' @noRd
gnomAD <- function(track, gnomAD.ID, bbdd){
  id.same.position <- NA
  if (track!="exomes" & track!="genomes") stop ("'track' argument invalid. It must be exomes or genomes ")
  gnomad <- stringr::str_split(gnomAD.ID, "-")[[1]]
  id.gnomad <- bbdd[paste0("gnomad1.", track)] %>% unlist() %>% as.numeric()
  id.gnomad[length(id.gnomad)==0] <- NA
  if (is.na(id.gnomad)){
    id.same.position <- bbdd[paste0("gnomad2.", track)] %>% unlist() %>% as.numeric()
    if(length(id.same.position)==0){
      id.same.position <- NA
    }
  }
  id.up.down <- upDownPos (gnomAD.ID, track, id.gnomad, id.same.position)
  return(list=c(id.variant =id.gnomad, id.variant.same.pos=id.same.position, id.down=id.up.down[1], id.up=id.up.down[2]))
}

#' @noRd
upDownPos<-function(gnomAD.nomen, track, id.gnomad, id.same.position){ #x means location of the variant and table, exomes or genomes
  gnomad <- stringr::str_split(gnomAD.nomen, "-")[[1]]
  up.down <-c(NA,NA)
  if (is.na(id.gnomad) && is.na(id.same.position)){
    query.down <- paste0("SELECT ", track, "_id from ", track, "_gnomad WHERE POS> '", as.numeric(gnomad[2])-50 ,"'AND POS<'", gnomad[2], "'AND CHROM='", gnomad[1],"';")
    chrom1 <- connectionDB(query.down) %>%
              unlist() %>%
              as.data.frame()
    down<- chrom1[nrow(chrom1),]
    query.up <- paste0("SELECT ", track, "_id from ", track, "_gnomad WHERE POS> '", gnomad[2] ,"'AND POS<'", as.numeric(gnomad[2])+50, "'AND CHROM='", gnomad[1],"';")
    chrom2 <- connectionDB(query.up) %>%
              unlist() %>%
              as.data.frame()
    up <- chrom2[1,]
    up.down <- c(down, up)
  }
  return(up.down)
}

#'gnomADmerge
#' @param gnom character. Variant in gnomAD nomenclature. Use gnomADnomen() function if you don't have gnomAD nomenclature but yes cdna annotation.
#' @param specific.info data obtained with the geneSpecific () function
#' @return variant's information located in gnomAD in the gnomAD v2.1,1 (non-cancer) and gnomAD v2.1.1 (non-neuro). It uses connection to mySQL server.
#' @author Elisabet Munté Roca
#' @references
#' Karczewski, K.J., Francioli, L.C., Tiao, G. et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 581, 434–443 (2020). https://doi.org/10.1038/s41586-020-2308-7
#' @noRd
#' @examples
#' #####Knowing gnomAD nomenclature
#' nomenclature <- "17-41258474-T-C"
#' gene <- "BRCA1"
#' specific <- geneSpecific(nomenclature, geneSpecific(gene, gene.specific.df))
#' gnomADmerge(nomenclature, specific)
#' #####Not knowing gnomAD nomenclature
#' #' NM <- "NM_007294.3"
#' NC <- "NC_000017.10"
#' gene <- "BRCA1"
#' CCDS <- "CCDS11453.1"
#' variant <- "c.211A>G"
#' mutalyzer.info <- correctHgvsMutalyzer(NM, NC,  gene, variant)
#' variant.info <- varDetails(NM,NC, CCDS, gene, mutalyzer.info, skip.pred=FALSE)
#' nomenclature <- gnomADnomen(variant.info)
#' specific <- geneSpecific(nomenclature, geneSpecific(gene, gene.specific.df))
#' gnomADmerge(nomenclature, specific)
gnomADmerge <- function (gnom, specific.info, bbdd ){
  porc <- as.numeric(specific.info$IC)
  exomes.info.gnom <- gnomADinfo(track="exomes", gnomAD.ID = gnom, porc = porc, bbdd=bbdd)
  genomes.info.gnom <- gnomADinfo(track="genomes", gnomAD.ID = gnom, porc = porc, bbdd)

  total.non.cancer.subpopulations <- gnomADtotal ( exomes.info.gnom$non.cancer$subpopulations,genomes.info.gnom$non.cancer$subpopulations, porc)
  total.non.cancer.gender <-  gnomADtotal (exomes.info.gnom$non.cancer$gender, genomes.info.gnom$non.cancer$gender, porc)
  total.non.cancer.overall <- gnomADtotal (exomes.info.gnom$non.cancer$overall, genomes.info.gnom$non.cancer$overall, porc)

  total.non.neuro.subpopulations <- gnomADtotal(exomes.info.gnom$non.neuro$subpopulations,genomes.info.gnom$non.neuro$subpopulations, porc)
  total.non.neuro.gender <- gnomADtotal(exomes.info.gnom$non.neuro$gender, genomes.info.gnom$non.neuro$gender, porc)
  total.non.neuro.overall <- gnomADtotal(exomes.info.gnom$non.neuro$overall, genomes.info.gnom$non.neuro$overall, porc)

  return(list(exomes= exomes.info.gnom,
              genomes = genomes.info.gnom,
              exomes.genomes = list(non.cancer=list(subpopulations=total.non.cancer.subpopulations, gender= total.non.cancer.gender, overall= total.non.cancer.overall),
                                    non.neuro=list(subpopulations=total.non.neuro.subpopulations, gender= total.non.neuro.gender, overall= total.non.neuro.overall))))
}

#' @noRd
gnomADinfo <- function(track, gnomAD.ID, porc, bbdd){
  positions <- positionToLook(track, gnomAD.ID, bbdd)
  AC_columns <- c("non_cancer_AC_nfe", "non_cancer_AC_fin", "non_cancer_AC_amr", "non_cancer_AC_afr","non_cancer_AC_sas", "non_cancer_AC_eas",
                  "non_cancer_AC_asj",  "non_cancer_AC_oth", "non_cancer_AC_male", "non_cancer_AC_female","non_neuro_AC_nfe", "non_neuro_AC_fin", "non_neuro_AC_amr",
                  "non_neuro_AC_afr", "non_neuro_AC_sas", "non_neuro_AC_eas", "non_neuro_AC_asj", "non_neuro_AC_oth", "non_neuro_AC_male","non_neuro_AC_female")
  AN_columns <- c("non_cancer_AN_nfe", "non_cancer_AN_fin", "non_cancer_AN_amr", "non_cancer_AN_afr", "non_cancer_AN_sas", "non_cancer_AN_eas",
                  "non_cancer_AN_asj", "non_cancer_AN_oth", "non_cancer_AN_male", "non_cancer_AN_female", "non_neuro_AN_nfe",  "non_neuro_AN_fin", "non_neuro_AN_amr",
                  "non_neuro_AN_afr", "non_neuro_AN_sas", "non_neuro_AN_eas", "non_neuro_AN_asj", "non_neuro_AN_oth", "non_neuro_AN_male", "non_neuro_AN_female")
  nhomalt_columns<- c("non_cancer_nhomalt_nfe", "non_cancer_nhomalt_fin", "non_cancer_nhomalt_amr", "non_cancer_nhomalt_afr", "non_cancer_nhomalt_sas", "non_cancer_nhomalt_eas",
                      "non_cancer_nhomalt_asj",  "non_cancer_nhomalt_oth",  "non_cancer_nhomalt_male", "non_cancer_nhomalt_female",  "non_neuro_nhomalt_nfe", "non_neuro_nhomalt_fin",
                      "non_neuro_nhomalt_amr",  "non_neuro_nhomalt_afr", "non_neuro_nhomalt_sas",  "non_neuro_nhomalt_eas", "non_neuro_nhomalt_asj","non_neuro_nhomalt_oth",
                      "non_neuro_nhomalt_male","non_neuro_nhomalt_female" )
  #AC and nhomalt
  AC <- rep(0,20)
  nhomalt <- rep(0,20)
  if (!is.na(positions[["positions.gnomad"]]["id.variant"])){
    AC <- positions[["info.gnomad"]][,AC_columns] %>%
          as.numeric()
    nhomalt <- positions[["info.gnomad"]][,nhomalt_columns] %>%
          as.numeric()
  }
  #AN
  AN <- rep(0,20)
  if (all(!is.na(positions[["position.to.check"]]))){
    AN <- positions[["info.gnomad"]] %>%
          dplyr::select (all_of(AN_columns))
    AN<-apply(AN, 2, mean)}
  AF <- AC/AN
  #IC interval
  gnom.merge <- data.frame(AC=matrix(unlist(AC), nrow=20), AN=matrix(unlist(AN), nrow=20), nhomalt=matrix(unlist(nhomalt), nrow=20), AF=matrix(unlist(AF), nrow=20))
  gnom.merge <- gnom.merge %>%
                dplyr::rowwise() %>%
                dplyr::mutate (CI= ifelse(is.na(porc),NA,  CI(AC, AN, porc))) %>%
                as.data.frame()
  row.names(gnom.merge) <- gsub("_AC_", "_", AC_columns)

  non.cancer.subpopulations <- gnom.merge[1:8,]
  non.cancer.gender <- gnom.merge[9:10,]
  non.neuro.subpopulations <- gnom.merge[11:18,]
  non.neuro.gender <- gnom.merge[19:20,]

  non.cancer.overall <- apply(non.cancer.gender, 2, sum) %>%
                        as.data.frame %>%
                        t()
  non.neuro.overall <- apply(non.neuro.gender,2, sum) %>%
                        as.data.frame %>%
                        t()
  return (list(non.cancer=list(subpopulations=non.cancer.subpopulations,
                               gender= non.cancer.gender,
                               overall= non.cancer.overall),
               non.neuro=list(subpopulations=non.neuro.subpopulations,
                              gender= non.neuro.gender,
                              overall= non.neuro.overall)))
}

#' @noRd
positionToLook <- function (track, gnomAD.ID, bbdd){
  positions <- unlist(gnomAD(track, gnomAD.ID, bbdd))
  if (all(is.na(positions)==TRUE)) {
    position <- NA} else{
      pos.na <- !is.na(positions)
      i <- 1
      while (pos.na[i]==FALSE){
        i <- i+1
        if (i == 5){  #if i becomes 5 means that no position has been found
          break }
      }
      if(i==2 & stringr::str_detect(names(positions)[3], "id.variant.same.pos")){
        i <- c(i,3)
      }else if (i == 3 & !is.na(pos.na[4])){
        i <- c(i,4)
      } #if the value is 3, we need to know if value 4 is not NA in order to do the median
      i[i==5] <- NA
      position <- ifelse(!is.na(i), unlist(positions[i]), NA)
    }
  if (length(position)>1){
    position.query <- stringr::str_c(position[1], position[2], sep=paste0("' OR ", track,"_id ='"))}else{
      position.query <- position
    }
  con <- DBI::dbConnect(RMySQL::MySQL(), user='userguest', dbname='class_variants', host='varhcdb001.cluster-ro-ca55bxrovxyt.eu-central-1.rds.amazonaws.com', password='jNU%cd%Xjw*tY*%')
  query.info.gnomad <- paste0("SELECT * from ", track, "_gnomad WHERE ", track,"_id = '", position.query ,"';")
  on.exit(DBI::dbDisconnect(con))
  info <- DBI::dbGetQuery(con, query.info.gnomad)
  if (track=="genomes" & nrow(info > 0)) info[c("non_cancer_AC_sas", "non_neuro_AC_sas", "non_cancer_AN_sas","non_neuro_AN_sas", "non_cancer_nhomalt_sas", "non_neuro_nhomalt_sas")]<-0 #these two variables don't exist, we create them but give 0
  if (track=="genomes" & nrow(info == 0)) info[c("non_cancer_AC_sas", "non_neuro_AC_sas", "non_cancer_AN_sas","non_neuro_AN_sas", "non_cancer_nhomalt_sas", "non_neuro_nhomalt_sas")] #these two variables don't exist, we create them but give 0
  return(list(positions.gnomad=positions, position.to.check=position, info.gnomad=info))
}

#' @noRd
CI <- function (AC, AN, porc){
  CI.total <- c()
  for (i in 1:length(AC)){
    CI.values <- DescTools::PoissonCI(as.numeric(AC[i]), n = 1, conf.level = porc, sides = "left",
                           method = "exact")[2]/AN[i]#non-cancer
    CI.total <- rbind(CI.total, CI.values, deparse.level=1)
  }
  return(CI.total)
}

#' @noRd
gnomADtotal <- function (datasetA, datasetB, porc){
  total <- datasetA + datasetB %>% as.data.frame()
  total <- tibble::rownames_to_column(total) %>%
    dplyr::rowwise() %>%
    dplyr::mutate (AF= AC/AN) %>%
    dplyr::mutate (CI=ifelse(is.na(porc), NA, CI(AC,AN,porc)))
  return(total)
}


#' Information from FLOSSIES db
#' @param nomen.flossies  character. FlOSSIES Id of the variant of interest. If you don't know it it can be set to NULL
#' @param nomen.gnomAD character. By default is null. If  nomen.flossies is not known you can use gnomAD nomenclature. gnomAD nomenclature can be infered with gnomADnomen() function.
#' @param gene character. The gene of interest
#' @return A dataframe with the information of the variant listen in FLOSSIES db (Fabulous Ladies Over Seventy).
#' @author Elisabet Munté
#' @references
#' https://whi.color.com/
#' @examples
#' ####Knowing FLOSSIES nomenclature
#' nomenclature.flossies <- "17-41245700-AGA-"
#' flossiesInfo (nomen.flossies=nomenclature.flossies)
#' @noRd
flossiesInfo <- function(nomen.flossies, nomen.gnomAD=NULL, gene){
  if(is.null(nomen.flossies)){
    gnomad <- stringr::str_split(nomen.gnomAD, "-") %>%
              unlist()
    dup <- stringr::str_length(gnomad[4]) > 1
    del <- stringr::str_length(gnomad[3]) > 1
    nomen.flossies <- nomen.gnomAD
    nomen.flossies <- ifelse(del,
                             paste0(gnomad[1], "-", as.numeric(gnomad[2])+1, "-", stringr::str_sub(gnomad[3],2, stringr::str_length(gnomad[3])), "-"),
                             ifelse(dup,
                                    paste0(gnomad[1], "-", as.numeric(gnomad[2])+1, "--", stringr::str_sub(gnomad[4],2, stringr::str_length(gnomad[4]))),
                                    nomen.gnomAD))
  }
  url.flossies <- RCurl::getURL(paste0("https://whi.color.com/variant/", nomen.flossies))
  if (gene %in% c( "BRCA1", "BRCA2", "ATM", "ATR", "BAP1", "BARD1", "BRIP1", "CDH1", "CHEK1", "CHEK2", "CTNNA1", "FAM175A", "FANCM", "GEN1", "MRE11A", "NBN", "PALB2", "PTEN", "RAD51B", "RAD51C", "RAD51D", "RECQL", "RINT1", "SLX4", "STK11", "TP53", "XRCC2")
      && stringr::str_detect(url.flossies, "was not found in the FLOSSIES data set") == FALSE & stringr::str_detect(url.flossies, "have not sequenced and analyzed") == FALSE ){
    flossies <-XML::readHTMLTable(url.flossies, as.data.frame=T, which = 1, stringsAsFactors=FALSE )#we want the frequency_table of flossies
    flossies2 <- data.frame(population=c("European American", "African American"),
                            carrier=as.numeric(flossies$`Carrier Count`),
                            homo=as.numeric(flossies$Homozygotes),
                            hete=as.numeric(flossies$Heterozygotes),
                            total=as.numeric(flossies$`Total Subjects`),
                            freq=as.numeric(flossies$`Carrier Frequency`))
    all_flossies <- cbind(population="all", dplyr::summarise(flossies2, carrier=sum(carrier), homo=sum(homo), hete=sum(hete), total=sum(total), freq=mean(freq)))
  } else {
    flossies2 <- data.frame(population = c("European American","African American"),carrier=c(0,0), homo=c(0,0), hete=c(0,0), total=c(0,0), freq=c(0,0))
    all_flossies <- data.frame(population="all", carrier=0, homo=0, hete=0, total=0, freq=0)
  }
  return(rbind(flossies2, all_flossies))
}


################################################################################
## ClinVar
################################################################################
#' clinVarIds variant
#' @param object obtained with the function varDetails()
#' @return A list containing ClinVar id of the variant of interest and for missense variants it also returns the variants located at the same codon.
#' @author Elisabet Munté Roca
#' @references
#' https://clinicaltables.nlm.nih.gov/apidoc/variants/v4/doc.html
#' @noRd
#' @examples
#' NM <- "NM_007294.3"
#' NC <- "NC_000017.10"
#' gene <- "BRCA1"
#' CCDS <- "CCDS11453.1"
#' variant <- "c.211A>G"
#' mutalyzer.info <- correctHgvsMutalyzer(NM, NC,  gene, variant)
#' variant.info <- varDetails(NM,NC, CCDS, gene, mutalyzer.info, skip.pred=FALSE)
#' clinVarIds (variant.info)
clinVarIds <- function(object){
  clinvar.info <- data.frame()
  server.clinvar <- "https://clinicaltables.nlm.nih.gov/api/variants/v3/search?"
  if (object$most.severe.consequence == "missense_variant"){
    data(BLOSUM62, package="Biostrings")
    aa.search <- paste0("p.", toProtein(object$protein)$aa.ref, toProtein(object$protein)$aa.pos)
    variant.clinvar <- stringr::str_replace(object$variant, "\\+", "\\\\+" )
    ext.clinvar.all <- paste0("terms=(", object$gene,")", aa.search, "{1}")
    clinvar.info <- api2(server.clinvar, ext.clinvar.all)[[4]] %>% tibble::as_tibble()
    if(nrow(clinvar.info)>0){
      clinvar.info <- clinvar.info %>%
                      dplyr::filter((stringr::str_detect(V2, paste0(aa.search, "[A-Z]+")) | stringr::str_detect(V2, paste0(aa.search, "[a-z]+"))) &
                                      stringr::str_detect(V2,object$gene) & !stringr::str_detect(V2, "del")& !stringr::str_detect(V2, "Ter") &!stringr::str_detect(V2, "fs"))
    }
    if(nrow(clinvar.info)>0){
      clinvar.info <- clinvar.info %>%
                      dplyr::rowwise() %>%
                      dplyr::mutate(variant =as.character(purrr::map(stringr::str_split(V2, " \\(|\\):"),2)), protein=as.character(purrr::map(stringr::str_split(V2, "\\(|\\)"),4))) %>%
                      dplyr::mutate(a1=aaShort(toProtein(protein)$aa.ref), a2=aaShort(toProtein(protein)$aa.alt )) %>%
                      dplyr::mutate (blosum = BLOSUM62[a1, a2], grantham=calculateGrantham(a1,a2))
      for (i in 1:nrow(clinvar.info)){
        prior.pm5 <- priorUtahProb(object=NULL, gene=object$gene, variant =clinvar.info$variant[i])[[3]] %>% unlist
        clinvar.info$prior[i] <- ifelse(length(prior.pm5)==0, NA, prior.pm5)
      }
    }
  }else{
    ext.clinvar.all <- paste0("terms=(", object$gene,")",object$variant,"&ef=HGVS_exprs")
    clinvar.info <- api2(server.clinvar, ext.clinvar.all)
    if(length(clinvar.info[[3]]$HGVS_exprs)!=0){
      number <- stringr::str_detect(clinvar.info[[3]]$HGVS_exprs, object$variant)
      number.variant <- which(number)
      if(length(number.variant)>0){
        clinvar.info <- clinvar.info[[4]][number.variant] %>%
                        tibble::as_tibble() %>%
                        dplyr::mutate(V2=paste0(object$NM, "(", object$gene, ")", object$variant," ",object$protein))
      }else{
        clinvar.info <- tibble::tibble()
      }
    }else{
      clinvar.info <- clinvar.info[[4]] %>% as.data.frame()
    }
    if (nrow(clinvar.info) ==0){ #sometimes the nomenclature used in clinvar defers from the nomenclature of our variant. We will search for protein then
      ext.clinvar.all <- paste0("terms=(", object$gene,")",object$protein, "&ef=HGVS_exprs")
      clinvar.info <- api2(server.clinvar, ext.clinvar.all)
      if(length(clinvar.info[[3]]$HGVS_exprs)!=0){
        number <- stringr::str_detect(clinvar.info[[3]]$HGVS_exprs, object$variant)
        number.variant <- which(number)
        if(length(number.variant)>0){
          clinvar.info <- clinvar.info[[4]][number.variant] %>%
                          tibble::as_tibble() %>%
                          dplyr::mutate(V2=paste0(object$NM, "(", object$gene, ")", object$variant," ",object$protein))
        }else{
          clinvar.info <- clinvar.info[[4]] %>% as_tibble()
        }
      }else{
        clinvar.info <- clinvar.info[[4]] %>% as_tibble()
      }
    }
    if(nrow(clinvar.info)>0){
      clinvar.info <- clinvar.info%>%
                      dplyr::filter(stringr::str_detect(V2,object$variant)&stringr::str_detect(V2,object$gene))
      names(clinvar.info) <- c("V1", "V2")
    }
  }
  clinvar.info.message <- "No warnings"
  if (nrow(clinvar.info)==0){
    clinvar.info.message <- "Warning: There are not variants reported in clinvar at this codon"
  }
  return(list(message=clinvar.info.message,
              table=clinvar.info))
}

#' ClinVar info
#' @param clinvar.id a list containing clinvarIDs, obtained with function clinVarIds()
#' @return All the information contained in Clinvar for this variant
#' @author Elisabet Munté Roca
#' @references
#' Landrum MJ, Lee JM, Benson M, Brown GR, Chao C, Chitipiralla S, Gu B, Hart J, Hoffman D, Jang W, Karapetyan K, Katz K, Liu C, Maddipatla Z, Malheiro A, McDaniel K, Ovetsky M, Riley G, Zhou G, Holmes JB, Kattman BL, Maglott DR. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Res . 2018 Jan 4. PubMed PMID: 29165669 .
#' @noRd
#' @examples
#' NM <- "NM_007294.3"
#' NC <- "NC_000017.10"
#' gene <- "BRCA1"
#' CCDS <- "CCDS11453.1"
#' variant <- "c.211A>G"
#' mutalyzer.info <- correctHgvsMutalyzer(NM, NC,  gene, variant)
#' variant.info <- varDetails(NM,NC, CCDS, gene, mutalyzer.info, skip.pred=FALSE)
#' clinvar.ids <- clinVarIds (variant.info)
#' clinvarDetails (clinvar.id= clinvar.ids)

clinVarInfo <- function (clinvar.id, object){
  if (clinvar.id$message!="Warning: It is not a missense variant."& clinvar.id$message!="Warning: There are not variants reported in clinvar at this codon"){
    clinvar.ids <- clinvar.id$table %>%
      tibble::as_tibble() %>%
      dplyr::mutate(url = paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", V1))
    b <- clinVarTable(clinvar.ids)
    names(b) <- clinvar.ids$V2
    for (i in 1:nrow(clinvar.ids)){
      if(!is.null(b[[i]])){
        names(b[[i]]) <- c("Nomenclature", "Condition", "Gene(s)", "Submitted.interpretations.and.evidence", "Functional", "Citations", "Classification")
      }
    }
    b <- Filter(Negate(is.null), b)
    return(list(variant = b[stringr::str_detect(names(b), object$variant)], same_codon = b[!stringr::str_detect(names(b), object$variant)]))
  }


}

#' @noRd
reviewStatus <- function(table.clinvar){
  if(exists("table.clinvar")){
    i<-1
    ben.guideline<-0
    ben.expert<-0
    ben.submitter<-0
    pben.guideline<-0
    pben.expert<-0
    pben.submitter<-0
    unc.guideline<-0
    unc.expert<-0
    unc.submitter<-0
    ppat.guideline<-0
    ppat.expert<-0
    ppat.submitter<-0
    pat.guideline<-0
    pat.expert<-0
    pat.submitter<-0
    while (i<=nrow(table.clinvar)) {
      # likely benign
      if(stringr::str_detect(table.clinvar[i,1],"Likely benign")){
        if(stringr::str_detect(table.clinvar[i,2],"practice guideline")){
          pben.guideline<-pben.guideline+1
        } else if(stringr::str_detect(table.clinvar[i,2],"reviewed by expert panel")){
          pben.expert<-pben.expert+1
        } else if(stringr::str_detect(table.clinvar[i,2],"criteria provided, single submitter|no assertion criteria provided")){
          pben.submitter<-pben.submitter+1
        }
        # benign
      } else if(stringr::str_detect(table.clinvar[i,1],"Benign")){
        if(stringr::str_detect(table.clinvar[i,2],"practice guideline")){
          ben.guideline<-ben.guideline+1
        } else if(stringr::str_detect(table.clinvar[i,2],"reviewed by expert panel")){
          ben.expert<-ben.expert+1
        } else if(stringr::str_detect(table.clinvar[i,2],"criteria provided, single submitter|no assertion criteria provided")){
          ben.submitter<-ben.submitter+1
        }
        # uncertain significance
      } else if(stringr::str_detect(table.clinvar[i,1],"Uncertain significance")){
        if(stringr::str_detect(table.clinvar[i,2],"practice guideline")){
          unc.guideline<-unc.guideline+1
        } else if(stringr::str_detect(table.clinvar[i,2],"reviewed by expert panel")){
          unc.expert<-unc.expert+1
        } else if(stringr::str_detect(table.clinvar[i,2],"criteria provided, single submitter|no assertion criteria provided")){
          unc.submitter<-unc.submitter+1
        }
        # likely pathogenic
      } else if(stringr::str_detect(table.clinvar[i,1],"Likely pathogenic")){
        if(stringr::str_detect(table.clinvar[i,2],"practice guideline")){
          ppat.guideline<-ppat.guideline+as.numeric(table.clinvar[i,3])
        } else if(stringr::str_detect(table.clinvar[i,2],"reviewed by expert panel")){
          ppat.expert<-ppat.expert+1
        } else if(stringr::str_detect(table.clinvar[i,2],"criteria provided, single submitter|no assertion criteria provided")){
          ppat.submitter<-ppat.submitter+1
        }
        # pathogenic
      } else if(stringr::str_detect(table.clinvar[i,1],"Pathogenic")){
        if(stringr::str_detect(table.clinvar[i,2],"practice guideline")){
          pat.guideline<-pat.guideline+1
        } else if(stringr::str_detect(table.clinvar[i,2],"reviewed by expert panel")){
          pat.expert<-pat.expert+1
        } else if(stringr::str_detect(table.clinvar[i,2],"criteria provided, single submitter|no assertion criteria provided")){
          pat.submitter<-pat.submitter+1
        }
      }
      i<-i+1
    }
  }
  total <- data.frame(ben_guideline=ben.guideline, pben_guideline=pben.guideline,unc_guideline=unc.guideline,ppat_guideline=ppat.guideline,pat_guideline=pat.guideline,
                    ben_expert=ben.expert, pben_expert=pben.expert, unc_expert=unc.expert, ppat_expert=ppat.expert, pat_expert=pat.expert,
                    ben_submitter=ben.submitter,   pben_submitter=pben.submitter,   unc_submitter=unc.submitter,  ppat_submitter=ppat.submitter,   pat_submitter=pat.submitter)
  return(total)}

#' @noRd
clinVarTable <- function(table){
  a <- lapply(table$url, function(g){
    page <- readUrl(g)
    if (!is.na(page)){
      table.ex <- rvest::html_table(page, fill=TRUE)
      colnames(table.ex[[4]]) <- c("Interpretation", "Review_status",  "Condition", "Submitter", "More_information", "more")
      interpret <- table.ex[[4]][1] %>%
                   as.data.frame() %>%
                   dplyr::mutate (Interpretation= purrr::map(stringr::str_split(Interpretation, "\\n"),1))
      date.clin <- table.ex[[4]][1] %>%
                   as.data.frame() %>%
                   dplyr::mutate (date= purrr::map(stringr::str_split(Interpretation, "\\(|\\)"),2))
      review <- table.ex[[4]][2] %>%
                as.data.frame() %>%
                dplyr::mutate (review= purrr::map(stringr::str_split(Review_status, "\\n"),1)) %>%
                dplyr::select(review)
      review.mode <- table.ex[[4]][2] %>%
                     as.data.frame() %>%
                     dplyr::mutate (review2= purrr::map(stringr::str_split(Review_status, "\\n"),3)) %>%
                     dplyr::select(review2)
      table.ex[[4]][1] <- unlist(interpret$Interpretation)
      table.ex[[4]][2] <- unlist(review)
      table.ex[[4]][6] <- unlist(date.clin$date)
      table.ex[[4]]$review_mode <- unlist(review.mode)
      table.ex[[4]] <- table.ex[[4]] %>%
                       dplyr::rowwise() %>%
                       dplyr::mutate(More_information=stringr::str_replace_all(More_information, "\n", "") %>%
                       stringr::str_replace_all("  ", ""), Condition=stringr::str_replace_all(Condition, "\n", ""))
      summary <- reviewStatus(as.data.frame(table.ex[[4]]))
      text.ex <- rvest::html_nodes(page, "dd")
      interpret <- stringr::str_replace_all(rvest::html_text(text.ex)[1], "  |\\n","")
      status <- stringr::str_replace_all(rvest::html_text(text.ex)[2], "  |\\n","")
      submissions <- stringr::str_replace_all(rvest::html_text(text.ex)[3], "  |\\n","")
      mix <- data.frame(interpret=interpret, status=status, submissions=submissions)
      table.ex$mix <- mix
      table.ex[["summary"]] <- summary
      return(table.ex)
    }
  })
  return (a)
}


#' @noRd
secondMet <- function(variant.mutalyzer, object, assembly){
  if(!(assembly %in% c("hg19", "hg38"))) stop("Assembly not supported, please enter hg19 or hg38")
  if(assembly== "hg38") stop("We are sorry but right now only hg19 works, we will incorporate hg38 soon")
  clinvar.pvs1.info <- NA
  pvs1.start.significance <- NA
  pos.second.met <- NA
  if(object$most.severe.consequence == "start_lost"){
    pos.second.met <- unlist(stringr::str_locate_all(variant.mutalyzer$protein_predicted, "M"))[2]
    coordinates.exon <- CodingTranscriptCds(object)

    start.sec <- ifelse (object$strand == 1,
                         coordinates.exon[1,1],
                         as.numeric(coordinates.exon[nrow(coordinates.exon),2])-(pos.second.met-1)*3)
    end.sec   <- ifelse (object$strand == 1,
                         as.numeric(coordinates.exon[1,1])+(pos.second.met-1)*3,
                         as.numeric(coordinates.exon[nrow(coordinates.exon),2]))

    ext.clinvar.pvs1 <- paste0("track=clinvarMain;chrom=chr", object$chr,";start=",start.sec, ";end=", end.sec)
    server.ucsc <- "http://api.genome.ucsc.edu/getData/track?" #UCSC's REST API
    clinvar.pvs1.info <- api2(paste0(server.ucsc, "genome=", assembly,";"), ext.clinvar.pvs1)[["clinvarMain"]] %>%
      as_tibble() %>%
                  dplyr::mutate (reviewStatus = purrr::map(stringr::str_split(reviewStatus, ":"),2),
                                 url = paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", purrr::map(stringr::str_split(origName, "\\|"),1)),
                                 origNames = purrr::map(stringr::str_split(origName, "\\|"),2),
                                 clinSign = unlist(clinSign)) %>%
                  dplyr::select(origName,clinSign, reviewStatus, type, molConseq, '_hgvsProt', chromStart, chromEnd,phenotype, url ) %>%
                  dplyr::mutate(reviewStatus = stringr::str_sub(reviewStatus, 2, -9)) %>%
                  dplyr::group_by(clinSign) %>%
                  dplyr::arrange(dplyr::desc(clinSign))

    pvs1.start.significance <- clinvar.pvs1.info %>%
                                                 dplyr::count(clinSign, reviewStatus) %>%
                                                 dplyr::arrange(clinSign, reviewStatus)
    table(unlist(clinvar.pvs1.info$clinSign), clinvar.pvs1.info$reviewStatus)
  }
  return (list(start.lost.variants= list(pos.second.met = pos.second.met,
                                         clinvar.pvs1.info=clinvar.pvs1.info,
                                         pvs1.start.significance)))
}




################################################################################
## Stop Codon
################################################################################
#' @noRd
stopCodon <- function(object, variant.mutalyzer){
  fs.window <- fsWindow (object)
  variant.exon <- NA
  porc.prot <- NA
  porc.prot.splicing <- NA
  stop.pos <- NA
  exons <- coordNonCoding(variant.mutalyzer , object)
  canonical.skipping <- NA
  stop.pos.all <- NA
  length.transcript <- tryCatch({
    CodingTranscriptCds(object) %>% lenCodingTranscriptCds
  },
  error=function(cond){
    return(NA)
  })

  if(object$most.severe.consequence =="frameshift_variant"|
     object$most.severe.consequence =="stop_gained"){
    variant.exon <- exons %>%
                    dplyr::filter( V1 <= object$start, V2 >= object$end) %>%
                    dplyr::mutate( cdna.var.pos = fs.window$dntp.pos)
  }else if (object$most.severe.consequence =="splice_donor_variant"|
            object$most.severe.consequence =="splice_acceptor_variant"){
    pos.sum <- ifelse(stringr::str_detect(object$variant,"\\+1")| stringr::str_detect(object$variant,"\\-1"),1,2) #if it is +1 o -1 => 1 and if it is +2 o -2 => 2
    if(object$most.severe.consequence =="splice_donor_variant"){
      posi <- ifelse(object$start==object$end,
                     object$start,
                     object$start)
      variant.exon <- exons %>%
                            dplyr::filter((V1 == (object$start+pos.sum)&object$strand==-1) | (V2 == (object$start-pos.sum)&object$strand==1)) %>%
                            dplyr::mutate( cdna.var.pos=fs.window$dntp.pos)
    }else if (object$most.severe.consequence == "splice_acceptor_variant"){
      posi <- ifelse(object$start==object$end,
                     object$start,
                     object$end)
      variant.exon <- exons %>%
                            dplyr::filter((V2 == (posi-pos.sum) & object$strand == -1) | (V1 == (posi+pos.sum)&object$strand==1)) %>%
                            dplyr::mutate( cdna.var.pos=fs.window$dntp.pos)
    }
    if (nrow(variant.exon)==0){
      value.close <- DescTools::Closest(c(as.numeric(exons$V1),as.numeric(exons$V2)), as.numeric(posi))
      variant.exon <- exons %>%
                            dplyr::filter (V1==value.close|V2==value.close)
    }
    if (variant.exon$exon !=1 &
        variant.exon$exon != nrow(exons) &
        !stringr::str_detect(exons$cStop[exons$exon==variant.exon$exon], "-") &
        !stringr::str_detect(exons$cStart[exons$exon==variant.exon$exon], "-")){
      pvs1.mutalyzer.cdna <- paste0("c.", variant.exon$cStart, "_", variant.exon$cStop, "del")
      correct.skipping <- correctHgvsMutalyzer(NM=object$NM, gene=object$gene, variant=pvs1.mutalyzer.cdna )
      skipping.info<- varDetails(NM=object$NM, NC=object$NC, CCDS=object$CCDS, gene=NULL, variant=NULL, variant.mutalyzer=correct.skipping, skip.pred=TRUE)
      fs.window <- fsWindow (skipping.info)
      canonical.skipping <- skipping.info %>%
                                          dplyr::select (variant, protein, most.severe.consequence) %>%
                                          tibble::remove_rownames()
      if (skipping.info$most.severe.consequence=="inframe_deletion"){
        porc.prot.splicing <- abs(as.numeric(variant.exon$V1)-as.numeric(variant.exon$V2))/length.transcript
      }
      canonical.skipping <- skipping.info %>%
                                          dplyr::select (variant, protein, most.severe.consequence) %>%
                                          dplyr::mutate (porc.prot.splicing=porc.prot.splicing) %>%
                                          tibble::remove_rownames()
    }
  }
  if (!is.na(fs.window$dntp.stop)){
    stop.pos.all <- exons %>%
                          dplyr::mutate (cStop = stringr::str_replace(cStop, "\\*[0-9]*", "1000000000"), cStart = stringr::str_replace(cStart, "\\*[0-9]*", "NA")) %>%
                          dplyr::mutate (cStart =as.numeric(cStart), cStop=as.numeric(cStop)) %>%
                          dplyr::filter(cStop >= fs.window$dntp.stop, cStart <= fs.window$dntp.stop) %>%
                          dplyr::mutate(bp=ifelse(cStop==1000000000, 0, as.numeric(cStop)-as.numeric(fs.window$dntp.stop))) %>%
                          dplyr::select (exon, bp) %>%
                          dplyr::mutate (num.codon.stop = fs.window$salt.fs, cDNA.stop.pos = fs.window$dntp.stop)
    porc.prot<-(length.transcript-as.numeric(fs.window$dntp.pos))/length.transcript
  }
  return(list(variant.exon = variant.exon ,
              premature.ter.codon = stop.pos.all,
              length.transcript = length.transcript,
              porc.prot = porc.prot,
              canonical.skip.pred= canonical.skipping,
              exons = exons))
}

#' @noRd
fsWindow <- function(object){
  dntp.var<-stringr::str_extract_all(object$variant,"[0-9]+") %>%
            unlist %>%
            as.numeric()
  if (length(dntp.var)==2){
    dntp.var <- dntp.var[2]
  }else if (length(dntp.var)>2){
    dntp.var <- NA
  }
  salt.fs <- NA
  salt.fs[object$most.severe.consequence=="frameshift_variant"] <- stringr::str_extract_all(object$protein, "[0-9]+")  #obtain number of codon until stop
  salt.fs[object$most.severe.consequence=="frameshift_variant"]<-salt.fs[[1]][length(salt.fs[[1]])] %>%
                                                                 as.numeric()
  salt.fs[object$most.severe.consequence=="stop_gained"] <- 0
  dntp.stop <- dntp.var + as.numeric(salt.fs)*3
  return(list(salt.fs = unlist(salt.fs),
              dntp.pos =  dntp.var,
              dntp.stop = dntp.stop))
}

lenCodingTranscriptCds <- function(coordinates.exon){
  coordinates.cdna <- as.numeric(coordinates.exon[1:nrow(coordinates.exon)-1,] %>%
                                   dplyr::mutate(exon=as.numeric(end)-as.numeric(start)) %>%
                                   dplyr::summarise(sum(exon)))
  return(coordinates.cdna)

}

#' @noRd
CodingTranscriptCds <- function(object){
  server.ucsc <- "http://api.genome.ucsc.edu/getData/track?" #UCSC's REST API
  ext.ucsc <- paste0("genome=hg19;track=ccdsGene;chrom=chr", unique(object$chr))
  ucsc <- api2(server.ucsc, ext.ucsc)
  exo <- ucsc$ccdsGene[ucsc$ccdsGene$name==object$CCDS,]
  coordinates.exon <- data.frame(start =unlist(stringr::str_split(exo$exonStarts, ",")), end=unlist(stringr::str_split(exo$exonEnds, ","))) %>%
                      dplyr::filter (start!="")
  return(coordinates.exon)
}


################################################################################
## insight
################################################################################
#' @noRd
insightUrl <- function (object, list.genes, database, browser){
  insight <- NA
  if(object$gene %in% list.genes)insight <- readTableUrlJavascript(paste0("http://www.insight-database.org/classifications/",database, ".html?gene=", object$gene, "&variant=", URLencode(object$variant), "&protein="), browser = browser)
  return (insight)
}

#' @noRd
insightInfo <- function (object, browser){
  index <- NA
  functional <- NA
  list.genes1 <- c("APC", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "MUTYH", "CDH1", "GALNT12")
  list.genes2 <- c("MLH1", "MSH2", "MSH6", "PMS2")
  index <- insightUrl (object, list.genes1, "index", browser)
  multifactorial <- insightUrl (object, list.genes2, "mmr_integrative_eval", browser)

  return (list (classification = index,
                multifactorial = multifactorial))
}


################################################################################
## Articles info
################################################################################

#' Functional studies
#' @description articlesInfo() retrieves information from different multifactorial studies or functional assays
#' @param object  obtained with varDetails() function
#' @param variant.cor obtained with correctHgvsMutalyzer() function
#' @param bbdd obtained with extractBBDD and connectionDB functions
#' @returns   If the variants has undergone some of the following  multifactorial studies or functional assays, it returns the information: Parsons et al.(2019); Lyra et al. (2021); Jia et al.(2021), Drost et al. (2018).
#' @author Elisabet Munté Roca
#' @examples
#' #' #' NM <- "NM_007294.3"
#' NC <- "NC_000017.10"
#' gene <- "BRCA1"
#' CCDS <- "CCDS11453.1"
#' variant <- "c.211A>G"
#'
#' #Correct variant nomenclature
#' mutalyzer.info <- correctHgvsMutalyzer(NM, NC,  gene, variant)
#'
#'#variant information
#' variant.info <- varDetails(NM,NC, CCDS, gene, mutalyzer.info, skip.pred=FALSE)
#'
#' gnom <- gnomADnomen(object = variant.info)
#' bbdd.info <- extractBBDD(mutalyzer = mutalyzer.info, object = variant.info, gnom = gnom) %>% connectionDB()
#' articlesInfo(object = variant.info, variant.cor = mutalyzer.info, ddbb = bbdd.info)
#' @references
#' Parsons, M. T., Tudini, E., Li, H., Hahnen, E., Wappenschmidt, B., Feliubadaló, L., ... & Pohl‐Rescigno, E. (2019). Large scale multifactorial likelihood quantitative analysis of BRCA1 and BRCA2 variants: An ENIGMA resource to support clinical variant classification. Human mutation, 40(9), 1557-1578.
#' Lyra, P. C., Nepomuceno, T. C., de Souza, M. L., Machado, G. F., Veloso, M. F., Henriques, T. B., ... & Monteiro, A. N. (2021). Integration of functional assay data results provides strong evidence for classification of hundreds of BRCA1 variants of uncertain significance. Genetics in Medicine, 23(2), 306-315.
#' Jia, X., Burugula, B. B., Chen, V., Lemons, R. M., Jayakody, S., Maksutova, M., & Kitzman, J. O. (2021). Massively parallel functional testing of MSH2 missense variants conferring Lynch syndrome risk. The American Journal of Human Genetics, 108(1), 163-175.
#' Drost, M., Tiersma, Y., Thompson, B. A., Frederiksen, J. H., Keijzers, G., Glubb, D., ... & Tavtigian, S. V. (2019). A functional assay–based procedure to classify mismatch repair gene variants in Lynch syndrome. Genetics in Medicine, 21(7), 1486-1496.
#' Kato, S., Han, S. Y., Liu, W., Otsuka, K., Shibata, H., Kanamaru, R., & Ishioka, C. (2003). Understanding the function-structure and function-mutation relationships of p53 tumor suppressor protein by high-resolution missense mutation analysis. Proceedings of the National Academy of Sciences of the United States of America, 100(14), 8424–8429. https://doi.org/10.1073/pnas.1431692100
#' Giacomelli, A. O., Yang, X., Lintner, R. E., McFarland, J. M., Duby, M., Kim, J., Howard, T. P., Takeda, D. Y., Ly, S. H., Kim, E., Gannon, H. S., Hurhula, B., Sharpe, T., Goodale, A., Fritchman, B., Steelman, S., Vazquez, F., Tsherniak, A., Aguirre, A. J., Doench, J. G., … Hahn, W. C. (2018). Mutational processes shape the landscape of TP53 mutations in human cancer. Nature genetics, 50(10), 1381–1387. https://doi.org/10.1038/s41588-018-0204-y
#' Kotler, E., Segal, E., & Oren, M. (2018). Functional characterization of the p53 "mutome". Molecular & cellular oncology, 5(6), e1511207. https://doi.org/10.1080/23723556.2018.1511207
#' @noRd
articlesInfo <- function (object, variant.cor, bbdd){
  lyra <- NA
  cimra <- NA
  adam.fayer <- NA
  jia <- NA
  parsons <- NA
  chek2.functionals<- NA
  tp53.functionals <- NA
  rna.pten <- NA
  atm.functionals <- NA
  prot <- toProtein(object$protein)
  prot.cor <- protsyn(object, variant.cor)
  if (object$gene == "BRCA1"){
    query <- "SELECT * from lyra_db;"
    lyra <- bbdd$lyra %>%
                      dplyr::mutate (Wild_type = Biostrings::AMINO_ACID_CODE[Wild_type],
                            Variant = Biostrings::AMINO_ACID_CODE[Variant],
                            prote = paste0("p.", Wild_type, Codon,Variant)) %>%
                      dplyr::filter(prote== prot.cor)
    adam.fayer <- bbdd$adamovich
  }

  if (object$gene %in% c("MLH1", "MSH2", "MSH6", "PMS2")){
    cimra <-bbdd$cimra1
    if(is.null(cimra) || nrow(cimra)==0){
      cimra <- bbdd$cimra2
    }
  }

  if (object$gene== "MSH2" & object$most.severe.consequence=="missense_variant"){
    query3 <- "SELECT * from jia_db;"
    jia <- bbdd$jia %>%
                    dplyr::filter(Variant == paste0(aaShort(prot$aa.ref),prot$aa.pos,aaShort(prot$aa.alt)))
  }


  parsons <- bbdd$parsons


  if (object$gene == "CHEK2"){
    query5 <- paste0("SELECT * from functionals_CHEK2 WHERE cDNA_variant='", object$variant, "';")
    chek2.functionals <- bbdd$chek2
  }

  if (object$gene == "TP53"){
    tp53.functionals.Kato <- bbdd$tp53.1 %>%
                                         dplyr::select(c_description, Kato_TransactivationClass, StructureFunctionClass, DNEclass, Giacomelli_DNE_LOFclass)
    tp53.functionals.kotler <- bbdd$tp53.2
    tp53.functionals <- data.frame(variant=object$variant,
                                   Kato=ifelse(nrow(tp53.functionals.Kato)==0,
                                               NA,
                                               as.character(tp53.functionals.Kato["Kato_TransactivationClass"])),
                                   Giacomelli= ifelse(nrow(tp53.functionals.Kato)==0,
                                                      NA,
                                                      as.character(tp53.functionals.Kato["Giacomelli_DNE_LOFclass"])),
                                   Kotler= ifelse(nrow(tp53.functionals.kotler)==0,
                                                  NA,
                                                  as.character(tp53.functionals.kotler["RFS_H1299"])),
                                   IARC_occurrences= ifelse(nrow(tp53.functionals.kotler)==0,
                                                            NA,
                                                            as.character(tp53.functionals.kotler["IARC_occurrences"])) )

  }

  if (object$gene == "PTEN"){
    rna.pten <- bbdd$pten
  }

  if (object$gene == "ATM"){
    atm.functionals<- bbdd$atm
  }

  return (list (lyra = lyra,
                Adamovich_Fayer = adam.fayer,
                cimra = cimra,
                jia = jia,
                parsons = parsons,
                chek2.functionals = chek2.functionals,
                tp53.functionals = tp53.functionals,
                rna.pten = rna.pten,
                atm.functionals = atm.functionals))
}



################################################################################
## google scholar
################################################################################
#' @noRd
biblioScholar <- function(object, first.article){
  articles.summary <- tibble()
  nomenclature.google <-  nomenclatureScholaR(object)
  articles.summary <- "Error 429: too many requests, please do the search manually."

for (i in first.article){
  url <- paste0("https://scholar.google.com/scholar?hl=en&start=", i, "&as_vis=1&as_sdt=0%252C5&q=", URLencode(nomenclature.google))
  information <- NA
  try(information <- readUrl(url))
  if (!is.na(information)){
    articles.page <- information %>%
                                rvest::html_nodes("h3") %>%
                                rvest::html_nodes("a") %>%
                                rvest::html_text() %>%
                                tibble::as_tibble(.name_repair=~c("Article"))
    articles.link <- information %>%
                                rvest::html_nodes("h3") %>%
                                rvest::html_nodes("a") %>%
                                rvest::html_attrs() %>%
                                purrr::map(2) %>%
                                unlist() %>%
                                tibble::as_tibble(.name_repair=~c("Link"))
    if(ncol(articles.link)==0){articles.link <- tibble::tibble("Link"=character(0), .name_repair=~c("Link"))}
    articles.author <- information %>%
                                   rvest::html_elements(".gs_a") %>%
                                   rvest::html_text() %>%
                                   tibble::as_tibble(.name_repair=~c("Authors_date"))
    #articles.cited <- information %>% html_elements(".gs_fl") %>% html_nodes("a") %>% html_text() %>% as_tibble() %>% filter(str_detect(value,"Cited by"))
    articles.cited <- information %>%
                                  rvest::html_elements (".gs_ri") %>%
                                  rvest::html_text() %>%
                                  unlist() %>%
                                  stringr::str_extract ("Cited by [0-9]*")   %>%
                                  tibble::as_tibble(.name_repair=~c("Cited_by"))
    articles.all <- cbind(articles.page, articles.author, articles.link, articles.cited) %>%
                                                                                         tibble::as_tibble()
    articles.summary <- rbind(articles.summary, articles.all)
  }else{
    #articles.summary= NA
    articles.summary= information
  }
}

  return(list(google.scholar.search=nomenclature.google, articles=articles.summary))
}

#' Google scholar search pattern
#' @description This function creates a pattern to search for the variant in google scholar
#' @param object  obtained with varDetails() function
#' @noRd
nomenclatureScholaR <- function(object){

  #####bibliographic search

  no.c <- stringr::str_sub(object$variant, 3)
  no.c.1 <- stringr::str_replace(no.c, ">", "->")
  no.c.2 <- stringr::str_replace(no.c, ">", "-->")
  no.c.3 <- stringr::str_replace(no.c, ">", "/")

  resuweb <- readUrl(URLencode(paste0("https://www.ncbi.nlm.nih.gov/snp/?term=", object$genomic))) %>%
                                                                                                   rvest::html_nodes("span") %>%
                                                                                                   rvest::html_nodes("a") %>%
                                                                                                   rvest::html_text()
  rs <- ifelse(stringr::str_detect(resuweb[1], "rs[0-9]+"),resuweb[1], "")
  if (object$most.severe.consequence %in% c("synonymous_variant", "missense_variant", "stop_retained_variant")){
    protein <- ifelse(object$most.severe.consequence=="missense_variant",
                      stringr::str_extract(object$protein, "[A-z]+[0-9]+[A-z]+" ),
                      paste0(stringr::str_extract(object$protein, "[A-z]+[0-9]+"), stringr::str_sub(stringr::str_extract(object$protein, "[A-z]+[0-9]"),1,3)))

    no.p.1 <- paste0(aaShort(stringr::str_sub(protein,1,3)),
                     stringr::str_extract(protein,"[0-9]+"),
                     aaShort(stringr::str_sub(protein,-3,-1)))

    biblio.variant <- paste0("'", object$gene, "'('", no.c, "' | '", no.c.1,"' | '", no.c.2,"' | '",
                             no.c.3, "' | '",   protein, "' | '",  no.p.1,"' | '",rs, "')")

  }else if (object$most.severe.consequence%in% c("intron_variant", "splice_donor_variant", "splice_acceptor_variant", "splice_donor_region_variant", "splice_acceptor_region_variant")){
    no.c.1.b <- paste0("IVS", object$exon, stringr::str_extract(object$variant, "[-|+][0-9]+"))
    no.c.1.c <- paste0("IVS ", object$exon, stringr::str_extract(object$variant, "[-|+][0-9]+"))

    if(stringr::str_detect(object$variant, "del|dup")){
      # web_mutalyzer <- xml2::read_html(URLencode(paste0("https://mutalyzer.nl/name-checker?description=",object$genomic)))
      # bases <- web_mutalyzer %>%
      #                        rvest::html_nodes("pre") %>%
      #                        stringr::str_extract(" [A-Z]+ ") %>%
      #                        stringr::str_extract("[A-Z]+")
      bases <- ifelse(stringr::str_detect(object$variant, "del"),
                      object$ref,
                      object$alt)
      if (object$strand == "-1")bases <- Biostrings::reverseComplement(Biostrings::DNAString(bases)) %>% as.character()
      no.c.1 <- paste0(no.c.1, bases)
      no.c.2 <- stringr::str_replace(no.c, "_", "-")
      no.c.3 <- paste0(no.c.2, bases)

    }

    biblio.variant <- paste0("'", object$gene, "'('", no.c, "' | '", no.c.1.b,"' | '", no.c.1.c,"' | '",no.c.1,"' | '", no.c.2,"' | '",
                             no.c.3, "' | '", rs,"')")

  }else if(object$most.severe.consequence=="stop_gained"){
    protein <- stringr::str_extract(object$protein, "[A-z]+[0-9]+[*]" )
    no.p.1 <- stringr::str_replace_all(protein, "\\*", "X")
    no.p.2 <- stringr::str_replace_all(protein, "\\*", "Ter")
    no.p.3 <- paste0(aaShort(stringr::str_sub(protein,1,3)),
                     stringr::str_extract(protein,"[0-9]+"), "*")
    no.p.4 <- stringr::str_replace_all(no.p.3, "\\*", "X")

    biblio.variant <- paste0("'", object$gene, "'('", no.c, "' | '",no.c.1,"' | '", no.c.2,"' | '",
                             no.c.3, "' | '", protein, "' | '",no.p.1, "' | '", no.p.2, "' | '",
                             no.p.3, "' | '", no.p.4, "' | '", rs, "'")

  }else if(object$most.severe.consequence=="frameshift_variant"){

    web.mutalyzer <- xml2::read_html(URLencode(paste0("https://mutalyzer.nl/name-checker?description=",object$genomic)))
    bases <- web.mutalyzer %>%
                           rvest::html_nodes("pre") %>%
                           stringr::str_extract(" [A-Z]+ ") %>%
                           stringr::str_extract("[A-Z]+")
    no.c.4 <- stringr::str_replace(no.c, "_", "-")
    no.c.5 <- paste0(no.c, bases)
    no.c.6 <- paste0(no.c.4, bases)
    protein <- stringr::str_extract(object$protein, "[A-z]+[0-9]+[A-z]+[*][0-9]+" )
    no.p.1 <- stringr::str_replace_all(protein, "\\*", "X")
    no.p.2 <- stringr::str_replace_all(protein, "\\*", "Ter")
    no.p.3  <- stringr::str_replace_all(protein, "fs\\*[0-9]+", "fs")
    split.prot <- stringr::str_split(protein, "[0-9]+") %>%
                                                        unlist
    no.p.6 <- paste0(aaShort(split.prot[1]),
                     stringr::str_extract(protein, "[0-9]+"),
                     aaShort(stringr::str_sub(split.prot[2],1,3)), "fs")

    no.p.4 <- paste0(no.p.6, stringr::str_extract(protein, "\\*[0-9]+"))
    no.p.5 <- stringr::str_replace_all(no.p.4, "\\*", "X")
    biblio.variant <- paste0("'", object$gene, "'('", no.c, "' | '", no.c.5, ifelse(stringr::str_detect(no.c, "_"), paste0("' | '", no.c.4, "' | '", no.c.6),""), "' | '",
                             protein, "' | '",no.p.1, "' | '", no.p.2, "' | '",
                             no.p.3, "' | '", no.p.4,"' | '", no.p.5,"' | '",
                             no.p.6, "' | '", rs, "')")

  }else if(object$most.severe.consequence=="stop_lost"){
    protein <- stringr::str_extract(object$protein, "[*][0-9]+[A-z]+[*][0-9]+" )
    no.p.1 <- stringr::str_replace_all(protein, "\\*", "X")
    no.p.2 <- stringr::str_replace_all(protein, "\\*", "Ter")
    no.p.3 <- paste0(stringr::str_extract(protein, "[*][0-9]+"),
                     aaShort(stringr::str_sub(stringr::str_extract(protein, "[A-z]+"),1,3)),
                     str_sub(stringr::str_extract(protein, "[A-z]+[*]+[0-9]+"),4,-1))
    no.p.4 <- stringr::str_replace_all(no.p.3, "\\*", "X")
    biblio.variant <- paste0("'", object$gene, "'('", no.c, "' | '",no.c.1,"' | '", no.c.2,"' | '",
                             no.c.3, "' | '", protein, "' | '",no.p.1, "' | '", no.p.2, "' | '",
                             no.p.3, "' | '", no.p.4, "' | '", rs, "')")
  }else if (object$most.severe.consequence %in% c("5_prime_UTR_variant", "3_prime_UTR_variant", "start_lost")){
    biblio.variant <- paste("'",object$gene, "'('", no.c, "' | '", rs, "')")
  }else if(object$most.severe.consequence %in% c("inframe_deletion", "inframe_insertion")){
    if(stringr::str_detect(object$variant, "del|dup")){
      # web_mutalyzer <- xml2::read_html(URLencode(paste0("https://mutalyzer.nl/name-checker?description=",object$genomic)))
      # bases <- web_mutalyzer %>%
      #                        rvest::html_nodes("pre") %>%
      #                        stringr::str_extract(" [A-Z]+ ") %>%
      #                        stringr::str_extract("[A-Z]+")
      bases <- ifelse(stringr::str_detect(object$variant, "del"),
                      object$ref,
                      object$alt)
      if (object$strand == "-1")bases <- Biostrings::reverseComplement(Biostrings::DNAString(bases)) %>% as.character()
      no.c.1 <- paste0(no.c.1, bases)
      no.c.2 <- stringr::str_replace(no.c, "_", "-")
      no.c.3 <- paste0(no.c.2, bases)
      biblio.variant <- paste0("'", object$gene, "'('", no.c, "' | '",no.c.1,"' | '", no.c.2,"' | '",
                               no.c.3, "' | '", rs,"')")
    }
  }

  return(biblio.variant)

}

################################################################################
## Cancer hotspots
################################################################################
#' Cancer hotspots databse
#' @description This function extracts information from https://www.cancerhotspots.org/ stored in IDIBELL database
#' @param object  obtained with varDetails() function
#' @param variant.cor obtained with correctHgvsMutalyzer() function
#' @param bbdd obtained with extractBBDD and connectionDB functions
#' @references Chang et al., Accelerating discovery of functional mutant alleles in cancer. Cancer Discovery, 10.1158/2159-8290.CD-17-0321 (2017).
#' Chang et al., Identifying recurrent mutations in cancer reveals widespread lineage diversity and mutational specificity. Nature Biotechnology 34, 155–163 (2016)
#' @noRd
cancerHotspots <- function (object, variant.cor, bbdd){
  prot.cor <- protsyn(object, variant.cor)
  if(!(prot.cor %in% c("p.(=)","p.?")) & !(stringr::str_detect(prot.cor, "ext")) & !(object$most.severe.consequence %in% c("frameshift_variant"))){
    aa.pos <- stringr::str_extract(prot.cor, "[0-9]+")
    aa.ref <- purrr::map(stringr::str_extract_all(prot.cor, "[A-z]+"),2) %>%
                                                                         unlist
    aa.alt <- ifelse(object$most.severe.consequence=="stop_gained", "*", purrr::map(stringr::str_extract_all(prot.cor, "[A-z]+"),3)) %>%
                                                                                                                                      unlist
   cancer.hotspots <- bbdd$cancer.hotspots %>%
                                            dplyr::mutate (Wild_type = Biostrings::AMINO_ACID_CODE[REF],
              Variant =ifelse(ALT=="*", ALT, Biostrings::AMINO_ACID_CODE[ALT])) %>%
      dplyr::filter(Wild_type== aa.ref, Variant==aa.alt, Amino_Acid_Position==as.numeric(aa.pos) )
    cancer.hotspots.final <- list(variant = cancer.hotspots[,2:11], cancer.type=NA)
    if (nrow(cancer.hotspots)>0){
      cancer.hotspots.type <- cancer.hotspots$Samples %>%
                                                        stringr::str_split("\\|") %>%
                                                        unlist() %>%
                                                        stringr::str_split(":")
      cancer.type <- data.frame(cancer_type=as.numeric(unlist(purrr::map(cancer.hotspots.type,2))), row.names = unlist(purrr::map(cancer.hotspots.type,1))) %>% t()
      cancer.hotspots.final[["cancer.type"]] <- cancer.type
    }
  }else{
    cancer.hotspots.final <- list(variant = NA, cancer.type=NA)
  }

  return(cancer.hotspots.final)
}
