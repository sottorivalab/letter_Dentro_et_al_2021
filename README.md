# letter_Dentro_et_al_2021

### General content

Code used to generate figures and results for the response letter to Dentro et al. 2021

The repo contains the following files:
* `rho_selected_cases.csv`: rho values estimated for the 41 manually picked samples in the letter
* `germline_analysis.R`: this file contains all the auxiliary functions to make the supllementary figure S3, i.e. the overdispersion reports. It generates the six panels and returns the best binomial and beta binomial fits with the corresponding parameters.
* `germline_supp_figures.R`: here there is the code to generate supplementary figures 1 and 2
* `figure3.R`: Code to generate figure 3
* `example_simulated`: a directory with all the code to run the simulated weme example using some tools from the original PCWAG consensus pipeline, all the information to run the example are in the `Example_simulated.Rmd` file
* `weme_test_BB_B`: in this directory you can find the code to reproduce the example and the simulation done in figure 2 of the letter, all the information and the code to run a minimal example are contained in the `CC_overdispersion.Rmd` Rmarkdown. The file `runner.R` actually runs the whole simulation and plots the output (Figure 2), the only input needed id the histogram of rho values estimated from the data that is automatically loaded from the parent directory.
* `pyclone_refits`: a direectory with the scripts to rerun the pyclone inference as done in the letter, more information below.

### Germline Analysis

For the germline analysis the function `generate_report_overdispersion` provides a complete report for overdispersion, with the 6 plots included in the Supplementary Data and binomial and beta-binomial MLE paramters. The input requires 5 objects that can be obtained from the PCWAG data release:

* *snps*: germinal SNPs, they can be obtained from the normal and tumour BAMS (restricted access), the input table should have the 5 following field "chr", "ref", "alt", "AD", "DP". Where "ref" is the reference allele, “alt“ the alternative and "DP" and "AD" can be obtained from the VCF.
*  *cna*: PCWAG consensus CN files, can be downloaded from https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_cnv/consensus.20170119.somatic.cna.annotated.tar.gz
*  *meta*: a table with at least a coloumn named `purity` with the actual purity values for the sample. Can be downloaded from https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_cnv/consensus.20170119.somatic.cna.annotated.tar.gz
*  *snvs*: cna be obtained by inner joining all the SNVs files released on the `subclonal_reconstruction` folder of ICGC data portal
*  *name*: the sample name (just used to generate plot labels)

Once you have all the data you can simply do:

```{R}

require(tidyverse)
require(ggrepel)

source(germline_analysis.R)

res <- generate_report_overdispersion(snps, snvs, cna, name, meta)

# overdispersion report
res$plts

# information criteria and lk
res$scores

# Fitted objects more information on the fields here (https://github.com/caravagnalab/BMix)
res$BBMIX
res$BMIX

```
### Pyclone fit

Pyclone fits have the same inputs as the germline overdispersion analysis, in particular for each sample a *snvs* and *meta* file is needed as an input.
First we need to generate a valid input file for pyclone

```{R}
require(dplyr)

source("pyclone_refits/prepare_pyclone.R")

inp <- format_to_pyclone(snvs, meta)

inp %>% head()

```

Then we have to create a suitable directory structure, we assume that the varibale `name` stores the id of the sample

```{R}

dir.create(name, showWarnings = F)
write.table(x, file = paste0("./", name, "/pyclone_input.tsv"), sep = "\t", row.names = F, quote = F)

```

The script to run the pyclone analysis are then stored in the `run_pyclone.R` file

```{R}

source("pyclone_refits/run_pyclone.R")

run_pyclone(name) # generate pyclone_outputplot_pyclone.tsv file 
plot_pyclone(name) # generate pyclone_plot.pdf file 


```
