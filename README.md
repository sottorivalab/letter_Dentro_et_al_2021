# letter_Dentro_et_al_2021
Code used to generate figures and results for the response letter to Dentro et al. 2021

The repo contains the following files:
* `germline_analysis.R`: this file contains all the auxiliary functions to make the supllementary figure S3, i.e. the overdispersion reports. It generates the six panels and returns the best binomial and beta binomial fits with the corresponding parameters.
* `germline_supp_figures.R`: here there is the code to generate supplementary figures 1 and 2
* `figure3.R`: 


For the germline analysis the function `generate_report_overdispersion` provides a complete report for overdispersion, with the 6 plots included in the Supplementary Data and binomial and beta-binomial MLE paramters. The input requires 5 objects that can be obtained from the PCWAG data release:

* *snps*: germinal SNPs, they can be obtained from the normal and tumour BAMS (restricted access), the input table should have the 5 following field "chr", "ref", "alt", "AD", "DP". Where "ref" is the reference allele, “alt“ the alternative and "DP" and "AD" can be obtained from the VCF.
*  *cna*: PCWAG consensus CN files, can be downloaded from https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_cnv/consensus.20170119.somatic.cna.annotated.tar.gz
*  *meta*: a table with at least a coloumn named `purity` with the actual purity values for the sample. Can be downloaded from consensus.20170217.purity.ploidy.txt.gz
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


