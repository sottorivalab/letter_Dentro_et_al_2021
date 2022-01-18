# letter_Dentro_et_al_2021

This repository contains the code used to produce all the analyses and figures in the letter. The code has been organized into articles that can be accessed at the following [link](https://github.com/sottorivalab/letter_Dentro_et_al_2021). 

Each article provides the explanation for reproducing a particular point in the analysis. They cover the following points:

* Data generation and PCAWG CCF and clustering plots
* Germline overdispersion analysis
* WeMe consensus for beta-binomial vs binomial models
* WeMe consensus for a subset of deconvolution methods used in the paper
* Pyclone reanalysis
* Figure generation

We also release some processed dataset obtained from the analysis of overdispesion in germline SNPs (everything can be reproduced by following), in particular:

* `rho_selected_cases.tsv` a text file with the rho values for 41 selected sample
* `samples_letter.rda` which contains the tumour aliquot ids of the 41 sample
* `germline_analysis.rds` an rds file with `BMix` fits and common information criteria scores for each of the 41 sample
* `letter_SNPs.rds` an RDS file with the AD field for germline SNPs piled-up from tumour BAMs (download link is in the vignette).



We provide a docker image with everything installed to follow the vignettes in the website. It provides a convenient RStudio server interface to run the analysis (built upon `rocker/tidyverse`):

```{bash}

# user is by default rstudio
# you can then connect to localhost:8787 on your machine 
# and have an Rstudio server with all the packages ready  
docker run -p 8787:8787 -v$(pwd):/home/rstudio/workspace:rw -e PASSWORD=yourpasswordhere  smilite/letter_cell

```
