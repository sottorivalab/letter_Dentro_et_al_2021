
# Used packages
require(dplyr)
require(ggplot2)
require(cowplot)
require(VGAM)
require(ccube)
require(mclust)
require(fpc)
require(matrixStats)
require(readr)


source("../weme_test_BB_B/sampler.R")

seed = 3
COVERAGE = 45
PURITY = 0.65
REPS = 20

set.seed(seed)

generate_data <-  function(i){
  
  overdispersion = runif(1, 0.001, 0.02)
  N = runif(1, 500, 5000) %>% round()
  coverage = rpois(1,lambda = COVERAGE) %>% round()
  
  input = simulate_data(overdispersion = overdispersion,
                N = N,
                coverage = coverage, seed = seed, purity = PURITY)
  input$chr = 1
  input$position = sample.int(n = 1000000, replace = T, size = nrow(input))
  dir.create(paste0("sim_",i), showWarnings = F)
  saveRDS(input, file.path(paste0("sim_",i), "input_data_sim.rds"))
  

}


run_ccube <-  function(dir) {
  

  setwd(dir)
  library(ccube)
  
  numOfClusterPool = 1:6
  numOfRepeat = 5
  
  input <-  readRDS("input_data_sim.rds")
  
  ccube_inp <- data.frame(major_cn = 1, minor_cn = 1, total_cn = 2, 
                          purity = PURITY, normal_cn = 1,
                          mutation_id = paste(input$chr, input$position, sep = "_"),
                          var_counts = input$NV,
                          ref_counts = input$DP - input$NV,
                          total_counts = input$DP) %>% as_tibble()
  
  
  results_ccube <- RunCcubePipeline(ssm = ccube_inp, 
                                    numOfClusterPool = numOfClusterPool, 
                                    numOfRepeat = numOfRepeat,
                                    runAnalysis = T, 
                                    runQC = T)
  
  
  save(results_ccube, file = "results_ccube.rda")
  
  setwd("..")
}


data <-  lapply(1:REPS, function(x) generate_data(x))

ccube <-  lapply(paste0("sim_", 1:N), function(x) run_ccube(x))






