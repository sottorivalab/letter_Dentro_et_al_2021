library(tidyverse)

# load_input
load("input_data.rda", verbose = T)

inp$CCF <-  inp$VAF * 2

# load ccube fit 
load("results_ccube.rda", verbose = T)


inp$ccube_clusters <- results_ccube$res$label

p1 <- ggplot(inp, aes(x = CCF, fill = ccube_clusters %>% paste())) + geom_histogram(binwidth = 0.02) + 
  scale_fill_brewer("Cluster",palette = "Set2") +
  theme_bw() + ggtitle("Ccube (Gaussian-Binomial)")

dpclust_res <-  readr::read_tsv("output_dpclust/sample0_DPoutput_2000iters_1000burnin_seed123/sample0_2000iters_1000burnin_bestConsensusAssignments.bed")

inp$dpclust_clusters <- dpclust_res$cluster

p2 <- ggplot(inp, aes(x = CCF, fill = dpclust_clusters %>% paste())) + geom_histogram(binwidth = 0.02) + 
  scale_fill_brewer("Cluster",palette = "Set2") +
  theme_bw() + ggtitle("Dpclust (Binomial)")

pyclone_res <- readr::read_tsv("output_pyclone/example.tsv")

inp$pyclone_res <- pyclone_res$cluster_id + 1


p3 <- ggplot(inp, aes(x = CCF, fill = pyclone_res %>% paste())) + geom_histogram(binwidth = 0.02) + 
  scale_fill_brewer("Cluster",palette = "Set2") +
  theme_bw() + ggtitle("Pyclone (beta-Binomial)")

consensus_structure <- readr::read_tsv("weme/example_simulated/example_subclonal_structure.txt")

lk_c1 <- dbinom(x = inp$NV, size = inp$DP ,prob = consensus_structure$proportion[1] * 0.65, log = T) + log(consensus_structure$n_ssms[1] / 100) 
lk_c2 <- dbinom(x = inp$NV, size = inp$DP ,prob = consensus_structure$proportion[2] * 0.65, log = T) + log(consensus_structure$n_ssms[2] / 100)


inp$consensus_cluster <- apply(cbind(lk_c1, lk_c2), MARGIN = 1, FUN = function(x) which.max(x))


p4 <- ggplot(inp, aes(x = VAF, fill = consensus_cluster %>% paste())) + geom_histogram(binwidth = 0.02) + 
  scale_fill_brewer("Cluster",palette = "Set2") +
  theme_bw() + ggtitle("Weme consensus clustering")
