# loading
library(tidyverse)
library(dummies)
library(readxl)


## clean Y colnames
clean_name <- function(s) {
  if (str_detect(s, "[a-z]\\d$")) {
    num <- str_extract(s, "[a-z]") %>% utf8ToInt() - 96
    paste0(str_remove(s, "[a-z]\\d$"), ".", num, ".", str_extract(s, "\\d$"))
  } else if (str_detect(s, "[a-z]$")) {
    num <- str_extract(s, "[a-z]$") %>% utf8ToInt() - 96
    paste0(str_remove(s, "[a-z]$"), ".", num)
  } else {
    paste0(s, ".1")
  }
}


## processing function
process_Y_M <- function(round, country, country_abbr) {
  
  # files paths
  M_file <- paste0("Data/M_matrix/", country, "/", country, "_", round, "_M.csv" ) 
  Y_file <- paste0("Data/Geo_AB/", "afb_full_", tolower(round), ".csv" ) 
  C_file <- paste0("Data/Conversion/", country, ".csv")
  
  
  # read files
  Y          <- read.csv(Y_file, stringsAsFactors = FALSE)
  conversion <- read.csv(C_file, stringsAsFactors = FALSE, na.strings = '')
  M          <- read.csv(M_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  
  
  # select M variables
  colnames(M)
  M <- M %>% dplyr::select("Survey.Question", "Theta.1", "Theta.2", "Theta.3")
  
  
  # country labels
  lbl <- Y %>% 
    dplyr::select(country, respno) %>% 
    mutate(label = substr(respno, 1, 3)) %>% 
    dplyr::select(-respno) %>% 
    distinct()
  country_label <- lbl$country[lbl$label == country_abbr]
  cat("country label :\n", country_label)
  
  
  # select country
  Y %>% group_by(country) %>% summarise( n = n() )
  Y <- Y %>% filter( country == country_label )
  N <- nrow(Y)
  cat("\n\nN :\n", N)
  
  
  # handle with colnames formatting
  name0 <- paste0(country, "_", round, "_conv")
  name1 <- paste0(country, "_", round %>% tolower(), "_conv")
  name2 <- paste0(country, "_", round, "_Conv")
  name  <-c(name0, name1, name2)[c(name0, name1, name2) %in% colnames(conversion)]
  # name <- "SaoTome_R6_conv"
  
  # conversion
  conv <- conversion[[ name ]] %>% na.omit() %>% unique()
  conv <- sapply(conv, function(s) {ifelse(str_detect(s, "^(q|Q)"), s, paste0("q", s))})
  conv <- conv %>% tolower() %>% str_remove_all("_.*") %>% unique()
  cat("\n\nconversion names : \n", conv)
  
  
  # Y match with conversion
  Y <- Y %>% dplyr::select( one_of(conv) )
  colnames(Y) <- colnames(Y) %>% str_remove_all("^q") %>% sapply(clean_name)
  cat("\n\nY colnames : \n", colnames(Y))
  
  
  # dummy Y
  Y <- Y %>% dummy.data.frame(dummy.classes = class(Y[,1]), sep = ".")
  Y[1:10, 1:10]
  
  
  # combine Y and M
  combine <- M %>% 
    inner_join(
      Y %>% 
        t() %>% 
        as.data.frame(stringsAsFactors = FALSE) %>% 
        type_convert() %>%
        rownames_to_column(var = "question"),
      by = c( "Survey.Question" = "question"  )
    ) 
  
  # rebuild M
  M <- array(NA, c(3, 3, nrow(combine)))
  for (i in 1:nrow(combine)) {
    M[,,i] <- diag(combine[i, 2:4])
  }
  M[1:3, 1:3, 1:5]
  cat("\n\nK :\n", dim(M)[3])
  
  
  # rebuild Y
  Y <- combine[, 5:ncol(combine)] %>% t() %>% as.data.frame()
  Y <- as.data.frame(sapply(Y, as.numeric))
  question <- combine[,1] %>% as.data.frame()
  colnames(Y) <- question[,1]
  Y[1:10, 1:10]
  class(Y[2,2])
  cat("\n\nDimension of Y :\n", dim(Y))
  cat("\n\nConsistent with original data?\n", nrow(Y) == N)
  
  
  # save data
  M_rds_file <- paste0("Data/Clean_Data/", country, "/", round, "_", "M_", country, ".rds")
  Y_rds_file <- paste0("Data/Clean_Data/", country, "/", round, "_", "Y_", country, ".rds")
  
  saveRDS(M, M_rds_file)
  saveRDS(Y, Y_rds_file)
  
  cat("\n\n File Saved!\n\n")
}



## ************ SET VARIABLES HERE ************** #

country_abbr  <-  "STP"
country       <-  "Sao_Tome_and_Principe"
round         <-  "R6"

## ************ SET VARIABLES HERE ************** #


process_Y_M(round, country, country_abbr)




