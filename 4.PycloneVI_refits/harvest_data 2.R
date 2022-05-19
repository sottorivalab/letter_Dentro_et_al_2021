library(tidyverse)

# Here we just load the output of the pyclone refitting script

read_pyclone_BB_res <- function(x) 
{
  readr::read_tsv(file.path(x, "pyclone_output.tsv"))
}

#loading "names"
load("samples.rda")

all_pyclone_BB <-  easypar::run(lapply(names, list),FUN =  read_pyclone_BB_res, export = NULL, filter_errors = F )


saveRDS(all_pyclone_BB, file = "../letter_pyclone_BB_all.rds")
