#' @import assertthat
#' @import DBI
#' @import RMySQL


################################################################################
## Connection to IDIBELL DB
################################################################################
connectionDB <- function(query){
  assertthat::assert_that(is.character(query)|is.list(query), msg="please enter a character vector or list")
  con <- DBI::dbConnect(RMySQL::MySQL(),
                        user='userguest',
                        dbname='class_variants',
                        host='varhcdb001.cluster-ro-ca55bxrovxyt.eu-central-1.rds.amazonaws.com',
                        password='jNU%cd%Xjw*tY*%')
  on.exit(DBI::dbDisconnect(con))
  results <- lapply(query, function(x){
    return(tryCatch( DBI::dbGetQuery(con, x), error=function(e) NULL))
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
  if(NM == "NM_004360.5") ext.ensembl.id <- paste0("/xrefs/symbol/homo_sapiens/", "NM_004360.3", "?content-type=application/json")
  server.ensembl <- "http://grch37.rest.ensembl.org"
  ensembl.id <- api2(server.ensembl, ext.ensembl.id)
  if (length(ensembl.id)==0){
    ext.ensembl.id <- paste0("/xrefs/symbol/homo_sapiens/", gene, "?content-type=application/json")
    ensembl.id <- api2(server.ensembl, ext.ensembl.id)
  }
  ensembl.id <- ensembl.id[stringr::str_detect(ensembl.id$id, "ENST")==T,]
  return(ensembl.id)
}


################################################################################
## Grantham distance
################################################################################

calculateGrantham <- function(a1, a2) {
  grantham <- readr::read_tsv("https://gist.githubusercontent.com/danielecook/501f03650bca6a3db31ff3af2d413d2a/raw/5583a134b36b60762be6cd54002a0f4044338cd0/grantham.tsv", show_col_types = FALSE) %>%
              tidyr:: gather(SECOND,SCORE, -FIRST) %>%
              dplyr::filter(SCORE > 0)
  (grantham %>% dplyr::filter(FIRST == a1|FIRST==a2, SECOND == a2 |SECOND==a1))$SCORE
}


################################################################################
## Read web pages
################################################################################

## Without Java script
#' read url
#' @param url url to read
readUrl <- function(url){
  out <- tryCatch({
    message(paste0("Trying to enter the url: ", url))
    rvest::read_html(url)
  },
  error = function(cond){
    message(paste("URL does not seem to exist or it is not working: ", url))
    print(cond)
    # Return value in case of error
    return(NA)
  },
  warning = function(cond){
    message(paste("URL caused a warning:", url))
    message("Here's the original warning message:")
    print(cond)
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

readTableUrlJavascript <- function (url, browser="firefox",  port=4568L){
assertthat::assert_that(browser %in% c("firefox", "chrome","phantomjs"), msg="Only supported for firefix, chrome or phantomjs browser")

if (browser=="firefox"){
  ex.cap <- list("moz:firefoxOptions" = list(args = list('--headless')))
}else if(browser == "chrome"){
  ex.cap <- list("chromeOptions" = list(args = list('--headless')))
}else{
  ex.cap <- NULL
}
  information <- "Not working"
  rD <- try(rD <- RSelenium::rsDriver(browser=browser, port=port, extraCapabilities = ex.cap ))

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
      message("This is the 'try' part")
      try(driver$dismissAlert())
      information <- driver$getPageSource()[[1]] %>% rvest::read_html() %>% rvest::html_table()
      driver$close()
    },
    error=function(cond) {
      message("variant exists")
      message(cond)
      return(list())
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", url))
      message(cond)
      return(list())
    },
    finally={
      message(paste("Processed URL:", url))

    }
  )
  rD[["server"]]$stop()
  return(information)
}


