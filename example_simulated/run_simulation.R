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
REPS = 10

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

run_pyclone_vi <-  function(dir) {
  
  
  setwd(dir)

  input = readRDS("input_data_sim.rds")
  
  data.frame(mutation_id = paste0("S_",1:length(input$DP)), sample_id = "example",
             alt_counts	 = input$NV, ref_counts = input$DP - input$NV,
             normal_cn	 = 2, major_cn	 = 1,
             minor_cn	= 1, tumour_content = PURITY) %>% 
    readr::write_tsv(., "simulation_input.tsv")
  system("bash ../run_pyclone_sim.sh")
  
  setwd("..")
}


run_dpclust <-  function(dir) {
  library(dplyr)
  
  setwd(dir)
  
  input = readRDS("input_data_sim.rds")
  
  data.frame(chr = input$chr, end = input$position,
             Reference_Allele = "A", Tumor_Seq_Allele2 = "T",
             mut.count	 = input$NV, WT.count = input$DP - input$NV,
             subclonal.CN	 = 2, mutation.copy.number	 = 1,
             subclonal.fraction	 = 1, no.chrs.bearing.mut	=1) %>% 
    readr::write_tsv(., "simulation_input.txt")
  
  data.frame(sample = "simulation", subsample = "01",
             datafile = "./simulation_input.txt", cellularity = PURITY
  ) %>%
    readr::write_tsv(., "simulation.txt")
  
  system("bash ../run_dpclust.sh")
  

  setwd("..")
}



format_for_weme <- function(dir) {
  library(tidyverse)
  library(readr)
  setwd(dir)
  dir.create("example_simulated", showWarnings = F)
  dir.create("./example_simulated/method1", showWarnings = F)
  
  load("results_ccube.rda", verbose = T)
  
  results_ccube$ssm$cluster_id <- results_ccube$res$label
  
  ccube_sc = results_ccube$ssm %>% group_by(cluster_id) %>% 
    summarise(n_ssms = n(), proportion = ccube_ccf_mean %>%  unique()) %>%  
    rename(cluster = cluster_id) %>%  select(cluster, n_ssms, proportion) %>%  
    mutate(proportion = proportion * PURITY)
  
  readr::write_tsv("./example_simulated/method1/example_subclonal_structure.txt", 
                   x = ccube_sc %>%  filter(n_ssms > 50) %>%  
                     arrange(-proportion)%>% mutate(cluster = 1:length(cluster)))
  
  dir.create("./example_simulated/method2", showWarnings = F)
  
  cluster_composition_dpclust <- readr::read_tsv("simulation_DPoutput_2000iters_1000burnin_seed123/simulation_2000iters_1000burnin_bestClusterInfo.txt")
  dpclust_sc <-  cluster_composition_dpclust %>%  rename(cluster = cluster.no, proportion = location, n_ssms = no.of.mutations) %>%  
    mutate(proportion= proportion * PURITY) %>% select(cluster, n_ssms, proportion)
  
  readr::write_tsv("./example_simulated/method2/example_subclonal_structure.txt", 
                   x = dpclust_sc %>%  filter(n_ssms > 50) %>% 
                     arrange(-proportion)%>% mutate(cluster = 1:length(cluster)))
  
  
  dir.create("./example_simulated/method3", showWarnings = F)
  
  pyclone_out <- readr::read_tsv("simulation.tsv", col_types = cols()) %>% group_by(cluster_id, cellular_prevalence) %>%  summarize(n_ssms = n())
  
  pyclone_sc <-  pyclone_out %>%  rename(cluster = cluster_id, proportion = cellular_prevalence) %>%  select(cluster, n_ssms, proportion) %>% 
    mutate(proportion = proportion * PURITY)
  
  readr::write_tsv("./example_simulated/method3/example_subclonal_structure.txt", x = pyclone_sc %>%  filter(n_ssms > 50) %>%  arrange(-proportion) %>% mutate(cluster = 1:length(cluster)))
  
  setwd("..")
}


calculate_consensus <- function(dir) {
  source("weme/weme.R")
  
  setwd(dir)
  setwd("example_simulated/")
  sids = find_sids()
  genconsensus(sids,rounddown=FALSE)
  
  
  setwd("../..")
}

harvest_data <- function(dir) {
  library(tidyverse)
  setwd(dir)
  setwd("example_simulated/")
  consensus <- readr::read_tsv("example_subclonal_structure.txt") %>%  
    mutate(tool = "consensus")
  ccube <- readr::read_tsv("method1/example_subclonal_structure.txt")%>%  
    mutate(tool = "ccube")
  dpclust <- readr::read_tsv("method2/example_subclonal_structure.txt")%>%  
    mutate(tool = "dpclust")
  pyclone <- readr::read_tsv("method3/example_subclonal_structure.txt")%>%  
    mutate(tool = "pyclone")
  setwd("../..")
  
  return(list(consensus, ccube, dpclust, pyclone) %>%  do.call(rbind,.) )
}


data <-  lapply(1:REPS, function(x) generate_data(x))

ccube <-  lapply(paste0("sim_", 1:REPS), function(x) run_ccube(x))

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_pyclone_vi, filter_errors = FALSE, export = "PURITY")

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_dpclust, filter_errors = FALSE, export = "PURITY")

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list),  
             FUN = format_for_weme, filter_errors = FALSE, export = "PURITY")

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list),  
             FUN = calculate_consensus, filter_errors = FALSE)

res <- easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list),  
             FUN = harvest_data, filter_errors = FALSE)


res_bind <- do.call(rbind, lapply(1:REPS, FUN = function(i) res[[i]] %>%  mutate(rep = i)))


p_sim <-  ggplot(res_bind %>% group_by(tool, rep) %>% summarise(nclust = max(cluster))
                   , aes(nclust %>%  as.character(), fill = nclust %>%  as.character())) + geom_histogram(alpha = 0.85, stat = "count") + 
  facet_grid(.~tool) + scale_fill_brewer(palette = "Blues", direction = -1) + theme_bw() + xlab("#cluster") +
  ylab("") + ggtitle("#clusters on simulated datasets") + theme(legend.position = "None")

saveRDS(p_sim, "p_sim.rds")

