# Overdispersion in WGS influences subclonal reconstruction.

This repository contains the code used to produce all the analyses and figures in the letter. The code has been organized into articles that can be accessed at the following [link](https://sottorivalab.github.io/letter_Dentro_et_al_2021/). 

Each article provides the explanation for reproducing a particular point in the analysis. They cover the following points:

* Data generation and PCAWG CCF and clustering plots (1.data_fits_pcawg/1.1.data.Rmd, 1.data_fits_pcawg/1.2.fits.Rmd)
* Germline overdispersion analysis (2.germline_overdispersion/2.germline_overdispersion_analysis.Rmd)
* WeMe consensus for beta-binomial vs binomial models (3.Weme/3.1Binomila_BetaBinomial/3.1.CC_overdispersion.Rmd)
* WeMe consensus with the deconvolution methods used in the paper (3.Weme/3.2.PCAWG_pipeline/3.2.Example_simulated.Rmd, we also added a install_info.Rmd with some info about installing the required software)
* Pyclone re-analysis (4.PycloneVI_refits/4.pyclone_vi_fits.Rmd) 
* Figure generation (5.Figures/figure_*.Rmd)

We also release some processed dataset obtained from the analysis of overdispesion in germline SNPs in particular:

* `rho_selected_cases.tsv` a text file with the rho values for 41 selected sample
* `samples_letter.rda` which contains the tumour aliquot ids of the 41 sample
* `germline_analysis.rds` an rds file with `BMix` fits and common information criteria scores for each of the 41 sample
* `letter_SNPs.rds` an RDS file with the AD field for germline SNPs piled-up from tumour BAMs (download link is in the vignette).


