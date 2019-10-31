
# loading
library(tidyverse)
library(dummies)


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


missing_question <- function(country, country_abbr, round) {
  
  # file paths
  M_file <- paste0("Data/M_matrix/", country, "/", country, "_", round, "_M.csv" ) 
  Y_file <- paste0("Data/Geo_AB/", "afb_full_", tolower(round), ".csv" ) 
  C_file <- paste0("Data/Conversion/", country, ".csv")

  # read files
  Y <- read.csv(Y_file, stringsAsFactors = FALSE, na.strings = c("NA", ""))
  M <- read.csv(M_file, stringsAsFactors = FALSE, na.strings = c("NA", ""), fileEncoding = "UTF-8-BOM")
  
  # preprocess M
  M <- M[rowSums(is.na(M)) != ncol(M),] 
  M <- M %>% dplyr::select("Survey.Question", "Theta.1", "Theta.2", "Theta.3")
  
  # questions in M table
  M_question <- str_extract(M[,1], "[0-9]*\\.[0-9]*") %>% unique()
  
  # country labels
  lbl <- Y %>% 
    dplyr::select(country, respno) %>% 
    mutate(label = substr(respno, 1, 3)) %>% 
    dplyr::select(-respno) %>% 
    distinct()
  country_label <- lbl$country[lbl$label == country_abbr]
  
  # select country
  Y <- Y %>% filter( country == country_label )
  
  # question in Y
  Y <- Y %>% select( matches("q[0-9]") )
  # colnames(Y)
  Y_question <- colnames(Y) %>% 
    str_remove_all("^q") %>% 
    sapply(clean_name) %>%
    unname()
  
  # keep only questions appered in M
  Y <- Y[ , Y_question %in% M_question]
  colnames(Y) <- colnames(Y) %>% 
    str_remove_all("^q") %>% 
    sapply(clean_name) %>%
    unname()
  
  # dummy Y
  Y <- Y %>% dummy.data.frame(dummy.classes = class(Y[,1]), sep = ".")
  
  # check Y colnames with M first column
  missing <- colnames(Y)[!(colnames(Y) %in% M[,1])]
  return(missing)
}


country_abbr  <-  "GHA"
country       <-  "Ghana"
round         <-  "R5"

missing_question(country, country_abbr, round)




