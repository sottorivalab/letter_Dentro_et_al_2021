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


source("../BB_vs_B_weme/sampler.R")

seed = 3
COVERAGE = 45
PURITY = 0.65
REPS = 50

set.seed(seed)

generate_data <-  function(i){
  
  rhos <- readr::read_tsv("../rho_selected_cases.tsv") %>% filter(rho > 1e-12, rho < 0.025)
  
  X <- rhos$rho
  
  LL <- function(shape, rate) {
     R = dgamma(X, shape, rate)
    -sum(log(R))
     }
  
  mle_fit <- mle(LL, start = list(shape = 3.5, rate=194))
  
  overdispersion = rgamma(1, shape = mle_fit@coef[1], mle_fit@coef[2])
  N = runif(1, 1000, 5000) %>% round()
  coverage = rpois(N,lambda = COVERAGE) %>% round()
  
  input = simulate_data(overdispersion = overdispersion,
                        N = N,
                        coverage = coverage, seed = seed, purity = PURITY)
  input$chr = 1
  input$position = sample.int(n = 1000000, replace = T, size = nrow(input))
  input$overdispersion = overdispersion
  dir.create(paste0("sim_",i), showWarnings = F)
  saveRDS(input, file.path(paste0("sim_",i), "input_data_sim.rds"))
  
  
}


run_ccube <-  function(dir) {
  
  
  setwd(dir)
  library(ccube)
  library(tidyverse)
  
  numOfClusterPool = 1:6
  numOfRepeat = 5
  
  input <-  readRDS("input_data_sim.rds")
  
  ccube_inp <- data.frame(major_cn = 1, minor_cn = 1, total_cn = 2, 
                          purity = PURITY, normal_cn = 2,
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




run_mclust <-  function(dir) {
  
  library(dplyr)
  
  setwd(dir)
  
  inp = readRDS("input_data_sim.rds")
  
 
  ccf <- inp$VAF * 2
  
  library(mclust)
  
  res_mclust <- Mclust(ccf, G = 1:7)
  

  library(fpc)
  
  dir.create("mclust_out")
  
  if(res_mclust$G > 1){
    res_mclust_merged <- mergenormals(ccf,res_mclust, method = "ridge.uni")
    
    sm <-  summary(res_mclust_merged)
    
    
    data.frame(id = 1:length(sm$clusternumbers), n_ssm = sm$clustering %>%  table() %>%  as.numeric(), loc = sm$muarray) %>% readr::write_tsv("mclust_out/mclust_subclonal_structure.tsv")
    
    data.frame(clusters = sm$clustering) %>% readr::write_tsv("mclust_out/mclust_clusters.tsv")
    
  } else {
    sm <- summary(res_mclust)
    
    
    data.frame(id = 1, n_ssm = length(ccf), loc = sm$mean) %>% readr::write_tsv("mclust_out/mclust_subclonal_structure.tsv")
    
    data.frame(clusters = 1) %>% readr::write_tsv("mclust_out/mclust_clusters.tsv")
    
  }
  

  
  
  setwd("..")
  
  return(NULL)
  
} 


calculate_BIC_clip <- function(df, inp){
  library(matrixStats)
  
  pi <- df$num_SNV / sum(df$num_SNV)
  locs <- pmin(df$cellular_prevalence, 1)
  
  npar <- length(locs) + length(pi) - 1
  
  lk <- sapply(locs, function(loc) dbinom(inp$NV,prob = loc/2, size = inp$DP, log = T)) + log(pi)
  
  lk <-  apply(lk, 1, logSumExp) %>%  sum()
  
  NLL2 <- -lk * 2
  PEN <- log(nrow(df)) * npar
  
  return(NLL2 + PEN)
  
}


collapse_clusters_clip <- function(x, e = 0.01) {
  
  copy <- x
  pi <- x$num_SNV / sum(x$num_SNV)
  idx_to_collapse <- which(pi < e)
  idx_to_retain <- which(pi >= e)
  for (idx in idx_to_collapse){
    nearest_cluster <- which.min(abs(x$cellular_prevalence[idx] - x$cellular_prevalence[idx_to_retain]))
    copy$num_SNV[idx_to_retain[nearest_cluster]] <- copy$num_SNV[idx_to_retain[nearest_cluster]] + x$num_SNV[idx]
  }
  
  return(copy[idx_to_retain,])
  
}

run_clip <-  function(dir) {
  library(dplyr)
  
  setwd(dir)
  
  inp = readRDS("input_data_sim.rds")
  
  dir.create("clip_input", showWarnings = FALSE)
  
  data.frame(
    chromosome_index = inp$chr,
    position = inp$position,
    alt_count	 = inp$NV,
    ref_count = inp$DP - inp$NV
  ) %>%
    readr::write_tsv(., "clip_input/sample.snv.txt")
  
  data.frame(
    chromosome_index = 1,
    start_position = 1,
    end_position = max(inp$position) + 10000,
    major_cn = 1,
    minor_cn = 1,
    total_cn = 2
  ) %>%
    readr::write_tsv(., "clip_input/sample.cna.txt")
  
  as.data.frame(PURITY) %>%
    readr::write_tsv(., "clip_input/sample.purity.txt", col_names = FALSE)
  
  system(paste0("bash ../run_clip.sh ", dir))
  
  
  setwd("..")
}


best_cluster_clip <- function(dir) {
  
  library(dplyr)
  
  setwd(dir)
  
  inp = readRDS("input_data_sim.rds")
  
  
  res_dir <- paste0("../CliP/", dir, "/final_result/")
  
  dir.create("clip_output")
  
  files <- list.files(res_dir, full.names = T) %>% grep(pattern = "subclonal_structure", value = T)
  
  dfs <- lapply(files, readr::read_tsv)
  
  dfs <- lapply(dfs, function(x) collapse_clusters_clip(x, e = 0.01) )
  
  BICs <- lapply(dfs, function(df) calculate_BIC_clip(df, inp))
  
  best_model_idx <- unlist(BICs) %>%  which.min()
  
  (best_model_clip <- dfs[[best_model_idx]])
  
  best_model_clip %>% readr::write_tsv("clip_output/clip_subclonal_structure.tsv")
  
  best_lamda <- files[best_model_idx] %>%  stringi::stri_match(., regex = "lam.*") 
  
  best_ass <- readr::read_tsv(paste0(res_dir, "mutation_assignments_", best_lamda))
  
  best_ass %>% readr::write_tsv("clip_output/clip_clusters.tsv")
  
  setwd("..")
  
  
}


run_ctpsingle <- function(dir) {
  
  library(dplyr)
  
  setwd(dir)
  
  inp = readRDS("input_data_sim.rds")
  
  dir.create("ctp_input", showWarnings = FALSE)
  
  data.frame(
    Chromosome = inp$chr,
    Position = inp$position,
    Mutant = "A",
    Reference = "T",
    Mcount	 = inp$NV,
    Rcount = inp$DP - inp$NV,
    Multiplier = 2,
    Gender = "Female"
  ) %>%
    readr::write_tsv(., "ctp_input/sample.snv.txt")
  
  system("bash ../run_ctpsingle.sh ")
  
  
  setwd("..")
  
}


run_phylogic <- function(dir) {
  
  library(dplyr)
  
  setwd(dir)
  
  inp = readRDS("input_data_sim.rds")
  
  
  dir.create("phylogic_input", showWarnings = FALSE)
  
  data.frame(
    Hugo_Symbol = "",
    Chromosome = inp$chr,
    Start_position = inp$position,
    Reference_Allele = "A", # Made up, do not affect results
    Tumor_Seq_Allele2 = "T",
    t_alt_count	 = inp$NV,
    t_ref_count = inp$DP - inp$NV,
    local_cn_a1	 = 1,
    local_cn_a2	 = 1
  ) %>%
    readr::write_tsv(., "phylogic_input/sample0_input.txt")
  
  data.frame(
    sample = "sample0",
    maf_fn = "./phylogic_input/sample0_input.txt",
    seg_fn = "",
    cellularity = PURITY,
    timepoint = 0
  ) %>%
    readr::write_tsv(., "phylogic_input/sample0.sif")
  
  system("bash ../run_phylogic.sh ")
  setwd("..")
  
}

run_phylowgs <- function(dir) {
  
  library(dplyr)
  
  setwd(dir)
  
  inp = readRDS("input_data_sim.rds")
  
  
  dir.create("wgs_input", showWarnings = FALSE)
  
  data.frame(
    id = paste0("s", 0:(nrow(inp) - 1)),
    gene = paste0("gene",  0:(nrow(inp) - 1)),
    a = inp$DP - inp$NV,
    d	 = inp$DP,
    mu_r = 0.999,
    mu_v = 0.499
  ) %>% filter(row_number() <= 500) %>% 
    readr::write_tsv(., "wgs_input/ssm_data.txt")
  
  
  system("bash ../run_phylowgs.sh ")
  setwd("..")
  
}


filter_phylowgs <- function(dir) {
  
  library(dplyr)
  
  setwd(dir)
  
  inp = readRDS("input_data_sim.rds")
  
  system("bash ../decompress_output_phylowgs.sh ")
  
  library(R.utils)
  tryCatch(
  gunzip("wgs_output/example.summ.json.gz"), error = function(e) return(setwd("..")))
  
  library(jsonlite)
  
  wgs_inp <- readr::read_tsv("wgs_input/ssm_data.txt")
  
  lks <- jsonlite::read_json("wgs_output/example.summ.json")
  
  summary_trees <- lapply(lks$trees, function(s) {
    res <- s$populations %>% as.data.frame() %>% matrix(ncol = 3, byrow = T) %>%  as.data.frame()
    res <- res[-1,-2]
    names(res) <- c("n_ssm", "loc")
    res$id <- paste0(1:nrow(res))
    return(res %>% mutate(n_ssm = as.integer(n_ssm), loc = as.numeric(loc)))
  })
  
  has_supercluster <-  sapply(summary_trees, function(sm) {
    if(nrow(sm) > 1) {
    sm$n_ssm[2] > sm$n_ssm[1] * 3
      } else {FALSE}
    }
    )
  
  summary_trees <- summary_trees[!has_supercluster]
  
  nK <-  table(sapply(summary_trees, function(sm) nrow(sm))) 
  
  MAP_K <- names(nK)[which.max(nK)]
  
  which_K <- sapply(summary_trees, function(sm) nrow(sm) == MAP_K)
  
  summary_trees_K <- summary_trees[which_K]
  
  final_struct_wgs <- summary_trees_K %>%  do.call(rbind,.) %>%  group_by(id) %>%  summarize(n_ssm = mean(n_ssm) %>%  round(), loc = mean(loc))
  
  final_struct_wgs %>%  readr::write_tsv(., file = "wgs_output/final_struct_wgs.tsv")
  
  setwd("..")
  
}

run_sclust <- function(dir) {
  
  library(dplyr)
  
  setwd(dir)
  
  inp = readRDS("input_data_sim.rds")
  
  NUC <- c("A", "C", "G", "T")
  
  comp <- function(n) {
    if(n == "A") return("T")
    if(n == "T") return("A")
    if(n == "C") return("G")
    if(n == "G") return("C")
  }
  
  comp <- Vectorize(comp)
  
  dir.create("sclust_input", showWarnings = FALSE)
  
  data.frame(Mut_ID = paste0("example_",paste0("chr",inp$chr), "_SNM:", inp$position)) %>% mutate(
    
    Chr = paste0("chr",inp$chr),
    Position = inp$position,
    Wt = sample(NUC, replace = T, size = nrow(inp)), # Made up, do not affect results
    Mut = comp(Wt),
    Af_obs = inp$VAF, 
    Coverage	 = inp$DP,
    AF_exp = PURITY * 0.5,
    Mut_Copies = 1,
    Mut_Copies_Raw = rnorm(nrow(inp), mean = 1, sd = 0.05),
    Is_Subclonal_CN = 0,
    iCN = 2,
    # good approximation to the real one
    P_Is_Clonal = apply(inp, 1,  function(i) binom.test(i["NV"], n =  i["DP"], p = PURITY * 0.5)$p.value)
    
    
  ) %>%
    readr::write_tsv(., "sclust_input/example_muts_expAF.txt")
  
  system("bash ../run_sclust.sh ")

  
  setwd("..")
  
}

run_svclone <- function(dir) {
  
  library(dplyr)
  
  setwd(dir)
  
  inp = readRDS("input_data_sim.rds")
  
 
  dir.create("svclone_input", showWarnings = FALSE)
  
  data.frame(
    chrom = inp$chr,
    pos = inp$position,
    gtype = "1,1,1.0",
    var = inp$NV,
    ref = inp$DP - inp$NV
  ) %>% 
    readr::write_tsv(., "svclone_input/example_filtered_snvs.tsv")
  
  
  
  data.frame(
    sample = "example",
    purity = PURITY,
    ploidy = 2
  ) %>% 
    readr::write_tsv(., "svclone_input/pp_file.txt")

  
  system("bash ../run_svclone.sh ")
  
  
  setwd("..")
  
}


run_clonehd <- function(dir) {
  
  
  library(dplyr)
  
  setwd(dir)
  
  inp = readRDS("input_data_sim.rds")
  
  
  dir.create("clonehd_input", showWarnings = FALSE)
  
  data.frame(
    chrom = inp$chr,
    pos = inp$position,
    var = inp$NV,
    ref = inp$DP 
  ) %>% 
    readr::write_tsv(., "clonehd_input/clonehd_snvs.txt", col_names = FALSE)
  
  data.frame(
    chrom = inp$chr,
    pos = inp$position,
    ploidy = 2 + rnorm(nrow(inp), sd = 0.1)
  ) %>% 
    readr::write_tsv(., "clonehd_input/mean.tcn.txt", col_names = FALSE)
  
  data.frame(
    purity = PURITY
  ) %>% 
    readr::write_tsv(., "clonehd_input/purity.txt", col_names = FALSE)
  
  
  system("bash ../run_clonehd.sh ")
  
  
  setwd("..")
  
}

merge_under <- function(x, e = 0.01) {
  if(nrow(x) == 1) return(x)
  copy <- x
  pi <- x$n_ssms / sum(x$n_ssms)
  idx_to_collapse <- which(pi < e)
  idx_to_retain <- which(pi >= e)
  for (idx in idx_to_collapse){
    nearest_cluster <- which.min(abs(x$proportion[idx] - x$proportion[idx_to_retain]))
    copy$n_ssms[idx_to_retain[nearest_cluster]] <- copy$n_ssms[idx_to_retain[nearest_cluster]] + x$n_ssms[idx]
  }
  
  return(copy[idx_to_retain,])
}

format_for_weme <- function(dir) {
  library(tidyverse)
  library(readr)
  setwd(dir)
  dir.create("example_simulated", showWarnings = F)
  
  # Method 1 - Ccube
  dir.create("example_simulated/method1", showWarnings = F)
  
  load("results_ccube.rda", verbose = T)
  
  results_ccube$ssm$cluster_id <- results_ccube$res$label
  ccube_sc = results_ccube$ssm %>% 
    group_by(cluster_id) %>% 
    summarise(
      n_ssms = n(), 
      proportion = ccube_ccf_mean %>%  unique()
    ) %>%  
    rename(cluster = cluster_id) %>%  
    select(cluster, n_ssms, proportion) %>%  
    mutate(proportion = proportion * PURITY)
  
  readr::write_tsv(
    "example_simulated/method1/example_subclonal_structure.txt",
    x = ccube_sc %>%  
      merge_under()  %>%  
      arrange(-proportion) %>% 
      mutate(cluster = 1:length(cluster))
  )
  
  
  dir.create("example_simulated/method2", showWarnings = F)
  
  cluster_composition_dpclust <- readr::read_tsv("simulation_DPoutput_2000iters_1000burnin_seed123/simulation_2000iters_1000burnin_bestClusterInfo.txt", col_types = readr::cols())
  
  ## dpclust apparently returns the location in terms of total tumour cells so we multiply it by purity
  dpclust_sc <-  cluster_composition_dpclust %>%  
    rename(
      cluster = cluster.no, 
      proportion = location, 
      n_ssms = no.of.mutations) %>%  
    mutate(proportion= proportion * PURITY) %>% 
    select(cluster, n_ssms, proportion)
  
  readr::write_tsv(
    "example_simulated/method2/example_subclonal_structure.txt", 
    x = dpclust_sc %>%  
      merge_under() %>%  
      arrange(-proportion) %>% 
      mutate(cluster = 1:length(cluster))
  )
  
  dir.create("example_simulated/method3", showWarnings = F)
  
  pyclone_out <-
    readr::read_tsv("output_pyclone/simulation.tsv", col_types = cols()) %>%
    group_by(cluster_id, cellular_prevalence) %>%  
    summarize(n_ssms = n())
  
  pyclone_sc <- pyclone_out %>%  
    rename(cluster = cluster_id, proportion = cellular_prevalence) %>%  
    select(cluster, n_ssms, proportion) %>%
    mutate(proportion = proportion * PURITY)
  
  readr::write_tsv(
    "example_simulated/method3/example_subclonal_structure.txt",
    x = pyclone_sc %>%  
      merge_under() %>%  ungroup() %>%  
      arrange(-proportion) %>% 
      mutate(cluster = 1:length(cluster))
  )
  
  dir.create("example_simulated/method4", showWarnings = F)
  
  mclust_out <-
    readr::read_tsv("mclust_out/mclust_subclonal_structure.tsv", col_types = cols())
  
  mclust_out <- mclust_out %>%  
    rename(cluster = id, n_ssms = n_ssm, proportion = loc)
  
  readr::write_tsv(
    "example_simulated/method4/example_subclonal_structure.txt",
    x = mclust_out %>%  
      merge_under %>%  
      arrange(-proportion)
    
  )
  
  dir.create("example_simulated/method5", showWarnings = F)
  
  clip_out <-
    readr::read_tsv("clip_output/clip_subclonal_structure.tsv", col_types = cols())
  
  clip_out <- clip_out %>%  
    rename(cluster = cluster_index, n_ssms = num_SNV, proportion = cellular_prevalence) 
  
  readr::write_tsv(
    "example_simulated/method5/example_subclonal_structure.txt",
    x = clip_out %>%  
      merge_under() %>%  
      arrange(-proportion) %>% 
      mutate(cluster = 1:nrow(clip_out))
  )
  
  dir.create("example_simulated/method6", showWarnings = F)
  
  ctp_out <-
    readr::read_tsv("ctp_output/example_PCAWG_ctp_cluster_assignments.txt", col_types = cols()) %>% 
    group_by(mostLikely) %>%  summarize(n_ssms = n(), proportion = unique(averageFrequency) * PURITY) %>%  ungroup()
  
  ctp_out <- ctp_out %>%  
    rename(cluster = mostLikely)
  
  readr::write_tsv(
    "example_simulated/method6/example_subclonal_structure.txt",
    x = ctp_out %>%  
      merge_under() %>%  
      arrange(-proportion) %>% 
      mutate(cluster = 1:length(cluster))
  )
  
  dir.create("example_simulated/method7", showWarnings = F)
  
  ndt_out <-
    readr::read_tsv("./sample0.cluster_ccfs.txt", col_types = cols()) %>%  
    select(Cluster_ID, postDP_ccf_mean)
  
  ndt_ass <- readr::read_tsv("./sample0.mut_ccfs.txt", col_types = cols()) %>%  group_by(Cluster_Assignment) %>%  summarize(n_ssms = n()) %>%  ungroup()
  
  ndt_ass$proportion <- ndt_out$postDP_ccf_mean * PURITY
  
  ndt_ass <- ndt_ass %>%  
    rename(cluster = Cluster_Assignment) 
  
  readr::write_tsv(
    "example_simulated/method7/example_subclonal_structure.txt",
    x = ndt_ass %>%  
      merge_under() %>%  
      arrange(-proportion) %>% 
      mutate(cluster = 1:length(cluster))
  )
  
  dir.create("example_simulated/method8", showWarnings = F)
  
  wgs_out <-
    readr::read_tsv("wgs_output/final_struct_wgs.tsv", col_types = cols()) 
  
  wgs_out <- wgs_out %>%  
    rename(cluster = id, n_ssms = n_ssm, proportion = loc) 
  
  readr::write_tsv(
    "example_simulated/method8/example_subclonal_structure.txt",
    x = wgs_out %>%  
      merge_under() %>%  
      arrange(-proportion) %>% 
      mutate(cluster = 1:length(cluster))
  )
  
  dir.create("example_simulated/method9", showWarnings = F)
  
  sclust_out <-
    readr::read_tsv("sclust_output/example_mclusters.txt", col_types = cols()) 
  
  sclust_out <- sclust_out %>%  
    rename(cluster = Cluster_ID, n_ssms = Mutations_In_Cluster, proportion = CCF_Cluster ) %>%  mutate(proportion = proportion * PURITY) %>%  select(-Cluster_Peak_Height) %>%  select(cluster, n_ssms, proportion)
  
  readr::write_tsv(
    "example_simulated/method9/example_subclonal_structure.txt",
    x = sclust_out %>%  
      merge_under() %>%  
      arrange(-proportion) %>% 
      mutate(cluster = 1:length(cluster))
  )
  
  dir.create("example_simulated/method10", showWarnings = F)
  
  sv_out <-
    readr::read_tsv("svclone_output/ccube_out/snvs/example_subclonal_structure.txt", col_types = cols()) 
  
  readr::write_tsv(
    "example_simulated/method10/example_subclonal_structure.txt",
    x = sv_out %>%  
      merge_under() %>%  
      arrange(-proportion) %>% 
      mutate(cluster = 1:length(cluster))
  )
  
  dir.create("example_simulated/method11", showWarnings = F)
  
  hd_out <-read.table("clonehd_output/tumorSNV.snv.posterior.txt",skip = 5, header = F, sep = "") 
  clsts <- hd_out[,-c(1,2)] %>% apply(.,1,which.max)
  r_inp <- readRDS("input_data_sim.rds") %>% filter(position %in% hd_out$V2)
  r_inp$hd_c <- clsts
  hd_out <- r_inp %>% group_by(hd_c) %>%  summarize(proportion = mean(NV / DP))
  colnames(hd_out)[1] <- "cluster"
  hd_out$n_ssms <- clsts %>%  table() %>% as.numeric()
  
  readr::write_tsv(
    "example_simulated/method11/example_subclonal_structure.txt",
    x = hd_out %>%  select(cluster, n_ssms, proportion) %>% 
      merge_under() %>%  
      arrange(-proportion) %>% 
      mutate(cluster = 1:length(cluster))
  )
  
  ### make subdirectories 
  
  # 3 methods: pyclone, dpclust, ctpsingle
  dir.create("example_simulated_3", showWarnings = F)
  
  system("cp -r example_simulated/method2/ example_simulated_3")
  system("cp -r example_simulated/method3/ example_simulated_3")
  system("cp -r example_simulated/method6/ example_simulated_3")
  
  # 6 methods: pyclone, dpclust, ctpsingle, sclust, phylogic, ccube
  dir.create("example_simulated_6", showWarnings = F)
  
  system("cp -r example_simulated/method2/ example_simulated_6")
  system("cp -r example_simulated/method3/ example_simulated_6")
  system("cp -r example_simulated/method6/ example_simulated_6")
  system("cp -r example_simulated/method1/ example_simulated_6")
  system("cp -r example_simulated/method7/ example_simulated_6")
  system("cp -r example_simulated/method9/ example_simulated_6")
  
  # 9 methods: pyclone, dpclust, ctpsingle, sclust, phylogic, ccube, phylowgs, clip, mclust
  dir.create("example_simulated_9", showWarnings = F)
  
  system("cp -r example_simulated/method2/ example_simulated_9")
  system("cp -r example_simulated/method3/ example_simulated_9")
  system("cp -r example_simulated/method6/ example_simulated_9")
  system("cp -r example_simulated/method1/ example_simulated_9")
  system("cp -r example_simulated/method7/ example_simulated_9")
  system("cp -r example_simulated/method9/ example_simulated_9")
  system("cp -r example_simulated/method4/ example_simulated_9")
  system("cp -r example_simulated/method5/ example_simulated_9")
  system("cp -r example_simulated/method8/ example_simulated_9")
  
  setwd("..")
  
  return(NULL)
}


calculate_consensus <- function(dir) {
  source("./weme.R")
  
  setwd(dir)
  setwd("example_simulated/")
  sids = find_sids()
  genconsensus(sids,rounddown=FALSE)
  setwd("../")
  setwd("example_simulated_3/")
  sids = find_sids()
  genconsensus(sids,rounddown=FALSE)
  setwd("../")
  setwd("example_simulated_6/")
  sids = find_sids()
  genconsensus(sids,rounddown=FALSE)
  setwd("../")
  setwd("example_simulated_9/")
  sids = find_sids()
  genconsensus(sids,rounddown=FALSE)
  setwd("../")
  setwd("..")
}

harvest_data <- function(dir) {
  library(tidyverse)
  setwd(dir)
  setwd("example_simulated/")
  consensus <- readr::read_tsv("example_subclonal_structure.txt") %>%  
    mutate(tool = "weme_all")
  consensus_3 <- readr::read_tsv("../example_simulated_3/example_subclonal_structure.txt") %>%  
    mutate(tool = "weme_3")
  consensus_6 <- readr::read_tsv("../example_simulated_6/example_subclonal_structure.txt") %>%  
    mutate(tool = "weme_6")
  consensus_9 <- readr::read_tsv("../example_simulated_9/example_subclonal_structure.txt") %>%  
    mutate(tool = "weme_9")
  ccube <- readr::read_tsv("method1/example_subclonal_structure.txt")%>%  
    mutate(tool = "ccube")
  dpclust <- readr::read_tsv("method2/example_subclonal_structure.txt")%>%  
    mutate(tool = "dpclust")
  pyclone <- readr::read_tsv("method3/example_subclonal_structure.txt")%>%  
    mutate(tool = "pyclone")
  mclust <- readr::read_tsv("method4/example_subclonal_structure.txt")%>%  
    mutate(tool = "mclust")
  clip <- readr::read_tsv("method5/example_subclonal_structure.txt")%>%  
    mutate(tool = "clip")
  ctp <- readr::read_tsv("method6/example_subclonal_structure.txt")%>%  
    mutate(tool = "ctpsingle")
  phylogicndt <- readr::read_tsv("method7/example_subclonal_structure.txt")%>%  
    mutate(tool = "phylogicndt")
  phylowgs <- readr::read_tsv("method8/example_subclonal_structure.txt")%>%  
    mutate(tool = "phylowgs")
  sclust <- readr::read_tsv("method9/example_subclonal_structure.txt") %>%  
    mutate(tool = "sclust")
  svclone <- readr::read_tsv("method10/example_subclonal_structure.txt") %>%  
    mutate(tool = "svclone")
  #clonehd <- readr::read_tsv("method11/example_subclonal_structure.txt") %>%  
   # mutate(tool = "clonehd")
  setwd("../..")
  
  return(list(consensus, consensus_3, consensus_6, consensus_9, ccube, dpclust,
              pyclone, mclust, clip, ctp, phylogicndt, phylowgs, sclust, svclone) %>%  do.call(rbind,.) )
}


plot_ass <- function(dir) {
  setwd(dir)
  
  # load_input
  inp = readRDS("input_data_sim.rds") 
  
  # Heterozygous diploid mutations
  inp$CCF <-  inp$VAF * 2
  
  # load ccube fit 
  load("results_ccube.rda", verbose = T) 
  inp$ccube_clusters <- results_ccube$res$label
  
  #load pyclone-vi fit
  pyclone_res <- readr::read_tsv("output_pyclone/simulation.tsv", col_types = readr::cols()) 
  pyclone_res <- pyclone_res[gtools::mixedorder(pyclone_res$mutation_id),]
  inp$pyclone_res <- pyclone_res$cluster_id + 1
  
  inp <- inp %>% arrange(as.numeric(position))
  
  # load dpclust
  dpclust_res <-  readr::read_tsv("simulation_DPoutput_2000iters_1000burnin_seed123/simulation_2000iters_1000burnin_bestConsensusAssignments.bed", col_types = readr::cols()) %>%  arrange(as.numeric(start))
  inp$dpclust_clusters <- dpclust_res$cluster
  

  
  #load mclust
  mclust_clusters <- readr::read_tsv("mclust_out/mclust_clusters.tsv", col_types = readr::cols())
  inp$mclust_res <- mclust_clusters$clusters
  
  #load CliP
  clip_res <- readr::read_tsv("clip_output/clip_clusters.tsv", col_types = readr::cols()) %>%  arrange(as.numeric(position))
  inp$clip_res <- clip_res$cluster_index + 1
  
  #load CTPsingle
  ctp_res <- readr::read_tsv("ctp_output/example_PCAWG_ctp_cluster_assignments.txt", col_types = readr::cols()) %>%  arrange(as.numeric(Position))
  inp$ctp_res <- ctp_res$mostLikely
  
  #load PhylogicNDT
  phylogic_res <- readr::read_tsv("./sample0.mut_ccfs.txt", col_types = readr::cols()) %>%  arrange(Start_position) %>%  select(Chromosome, Start_position, Cluster_Assignment) %>% 
    rename(phylogic_res = Cluster_Assignment, chr = Chromosome, position = Start_position)
  inp <- full_join(phylogic_res, inp)
  
  #load PhyloWGS (simplified binomial assignment)
  wgs_res <- readr::read_tsv("wgs_output/final_struct_wgs.tsv", col_types = readr::cols())
  
  
  mat_lk <- apply(inp,1,  function(x) apply(wgs_res, 1, function(par) dbinom(x = x["NV"],prob = par["loc"]/2,size = x["DP"], log = T) +   log(par["n_ssm"] / nrow(inp)) )) %>%  t()
  
  tot_lk <- apply(mat_lk, 1,function(lk) logSumExp(lk))
  
  ps <- exp(mat_lk - tot_lk)
  
  inp$wgs_res <- apply(ps, 1, which.max)
  
  #load Sclust
  sclust_res <- readr::read_tsv("sclust_output/example_cluster_assignments.txt", col_types = readr::cols()) %>%  arrange(Position)
  inp$sclust_res <- as.numeric(sclust_res$Cluster_Id) + 1
  
  #load SVclone
  svclone_res <- readr::read_tsv("svclone_output/ccube_out/snvs/example_assignment_probability_table.txt", col_types = readr::cols())
  inp$svclone_res <- svclone_res %>% arrange(as.numeric(pos)) %>% select(starts_with("cluster")) %>%  apply(., 1 , which.max)
  #load consensus
  consensus_structure <- readr::read_tsv("example_simulated/example_subclonal_structure.txt", col_types = readr::cols())
  
  RHO = 0.01
  
  # Posterior clustering assignments
  mat_lk <- apply(inp,1,  function(x) apply(consensus_structure, 1, function(par) dbetabinom(x = x["NV"],prob = par["proportion"]/2,size = x["DP"], log = T, rho = RHO)  +  log(par["n_ssms"] / 100) )) %>%  t()
  
  tot_lk <- apply(mat_lk, 1,function(lk) logSumExp(lk))
  
  ps <- exp(mat_lk - tot_lk)
  
  inp$consensus_cluster_all <- apply(ps, 1, which.max)
  
  #load consensus
  consensus_structure <- readr::read_tsv("example_simulated_3/example_subclonal_structure.txt", col_types = readr::cols())
  
  RHO = 0.01
  
  # Posterior clustering assignments
  mat_lk <- apply(inp,1,  function(x) apply(consensus_structure, 1, function(par) dbetabinom(x = x["NV"],prob = par["proportion"]/2,size = x["DP"], log = T, rho = RHO)  +  log(par["n_ssms"] / 100) )) %>%  t()
  
  tot_lk <- apply(mat_lk, 1,function(lk) logSumExp(lk))
  
  ps <- exp(mat_lk - tot_lk)
  
  inp$consensus_cluster_3 <- apply(ps, 1, which.max)
  
  #load consensus
  consensus_structure <- readr::read_tsv("example_simulated_6/example_subclonal_structure.txt", col_types = readr::cols())
  
  RHO = 0.01
  
  # Posterior clustering assignments
  mat_lk <- apply(inp,1,  function(x) apply(consensus_structure, 1, function(par) dbetabinom(x = x["NV"],prob = par["proportion"]/2,size = x["DP"], log = T, rho = RHO)  +  log(par["n_ssms"] / 100) )) %>%  t()
  
  tot_lk <- apply(mat_lk, 1,function(lk) logSumExp(lk))
  
  ps <- exp(mat_lk - tot_lk)
  
  inp$consensus_cluster_6 <- apply(ps, 1, which.max)
  
  #load consensus
  consensus_structure <- readr::read_tsv("example_simulated_9/example_subclonal_structure.txt", col_types = readr::cols())
  
  RHO = 0.01
  
  # Posterior clustering assignments
  mat_lk <- apply(inp,1,  function(x) apply(consensus_structure, 1, function(par) dbetabinom(x = x["NV"],prob = par["proportion"]/2,size = x["DP"], log = T, rho = RHO)  +  log(par["n_ssms"] / 100) )) %>%  t()
  
  tot_lk <- apply(mat_lk, 1,function(lk) logSumExp(lk))
  
  ps <- exp(mat_lk - tot_lk)
  
  inp$consensus_cluster_9 <- apply(ps, 1, which.max)
  
  inp %>% head()
  
  df_b <- inp %>% 
    select(
      CCF, 
      ccube_clusters, 
      dpclust_clusters, 
      pyclone_res, 
      mclust_res,
      clip_res,
      ctp_res,
      phylogic_res,
      wgs_res,
      sclust_res,
      svclone_res,
      consensus_cluster_3,
      consensus_cluster_6,
      consensus_cluster_9,
      consensus_cluster_all) %>%  
    rename(
      pyclone = pyclone_res, 
      ccube = ccube_clusters, 
      dpclust = dpclust_clusters,
      mclust = mclust_res,
      clip = clip_res,
      ctpsingle = ctp_res,
      phylogicndt = phylogic_res,
      phylowgs = wgs_res,
      sclust = sclust_res,
      svclone = svclone_res,
      weme_3 = consensus_cluster_3,
      weme_6 = consensus_cluster_6,
     weme_9 = consensus_cluster_9,
     weme_all = consensus_cluster_all
    ) %>% 
    reshape2::melt(id.vars = "CCF")
  
  setwd("..")
  
  return(df_b)
  
}

data <-  lapply(1:REPS, function(x) generate_data(x))

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_ccube, filter_errors = FALSE, export = "PURITY", cores.ratio = 0.95)

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_pyclone_vi, filter_errors = FALSE, export = "PURITY", cores.ratio = 0.95)

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_dpclust, filter_errors = FALSE, export = "PURITY",  cores.ratio = 0.95)

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_mclust, filter_errors = FALSE, export = "PURITY",  cores.ratio = 0.95)

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_clip, filter_errors = FALSE, export = "PURITY",  cores.ratio = 0.95)

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = best_cluster_clip, filter_errors = FALSE, export = c("calculate_BIC_clip",
                                                                        "PURITY",
                                                                        "collapse_clusters_clip"))

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_ctpsingle, filter_errors = FALSE, export = "PURITY", cores.ratio = 0.5)

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_phylogic, filter_errors = FALSE, export = "PURITY",  cores.ratio = 0.95)

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_phylowgs, filter_errors = FALSE, export = "PURITY",  cores.ratio = 0.95)

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = filter_phylowgs, filter_errors = FALSE, export = "PURITY")

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_sclust, filter_errors = FALSE, export = "PURITY")

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_svclone, filter_errors = FALSE, export = "PURITY")

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list), 
             FUN = run_clonehd, filter_errors = FALSE, export = "PURITY", cores.ratio = 1)

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list),  
             FUN = format_for_weme, filter_errors = FALSE, export = c("PURITY", "merge_under"))

easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list),  
             FUN = calculate_consensus, filter_errors = T)

res <- easypar::run(PARAMS = lapply(paste0("sim_", 1:REPS), list),  
             FUN = harvest_data, filter_errors = T)


res_bind <- do.call(rbind, lapply(1:length(res), FUN = function(i) res[[i]] %>%  mutate(rep = i)))

palette <- list("1" = "forestgreen", "2" = "indianred", "3" = "indianred", 
                "4" = "indianred", "5" = "indianred", "6" = "indianred")

library(patchwork)

p_sim_1 <-  ggplot(res_bind %>% group_by(tool, rep) %>% summarise(nclust = length(unique(cluster))) %>% filter(!grepl(pattern = "weme", tool)) 
                   , aes(nclust %>%  as.character(), fill = nclust %>%  as.character())) + geom_histogram(alpha = 0.85, stat = "count") + 
  facet_wrap(.~tool, ncol = 4) + scale_fill_manual(values = palette) + theme_bw() + xlab("#cluster") +
  ylab("") + ggtitle("PCAWG pipeline methods on simulated data") + theme(legend.position = "None")

p_sim_2 <-  ggplot(res_bind %>% group_by(tool, rep) %>% summarise(nclust = max(cluster)) %>% filter(grepl(pattern = "weme", tool)) %>% mutate(tool = gsub("weme_", "", tool)) 
                   , aes(nclust %>%  as.character(), fill = nclust %>%  as.character())) + geom_histogram(alpha = 0.85, stat = "count") + 
  facet_wrap(.~tool, ncol = 4) + scale_fill_manual(values = palette) + theme_bw() + xlab("#cluster") +
  ylab("") + ggtitle("WEME consensus on simulated data", subtitle = "Consensus considering different number of methods") + theme(legend.position = "None")  

df_b <- plot_ass("sim_15")


new_ids <- df_b %>%  group_by(variable, value) %>%  mutate(pi = n() / sum(df_b$variable == "ccube")) %>%group_by(variable, value, pi) %>% summarize() %>% group_by(variable) %>% mutate(value = order(pi,decreasing = T))

df_b_plot <- df_b %>%  group_by(variable, value) %>%  mutate(pi = n() / sum(df_b$variable == "ccube")) %>% filter(pi > 0.01) %>% ungroup %>%  select(-value)


df_b_plot <- inner_join(df_b_plot, new_ids)

p_sim_3 <- ggplot( df_b_plot, 
               aes(x = CCF, fill = value %>% paste())
) + 
  geom_histogram(binwidth = 0.02) + 
  scale_fill_brewer("Cluster", palette = "Set2") + 
  facet_wrap(.~variable, nrow = 3) + 
  theme_bw() + ggtitle("Clustering assignment simulated data") + theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")


p_final = p_sim_1 / p_sim_2/ p_sim_3 + patchwork::plot_layout(heights = c(4,1,4)) + patchwork::plot_annotation(tag_levels = "a")&
  theme(plot.tag = element_text(face = 'bold'))

p_final %>%  ggsave("figure_2.png", ., device = "png", units = "px", dpi = 300, width = 2400, height = 3600)

p_3_df <- plot_ass("sim_2")

saveRDS(p_sim_1, "p_sim.rds")

