library(tidyverse)
library(patchwork)

# load_input
load("input_data.rda", verbose = T)

inp$CCF <-  inp$VAF * 2

# load ccube fit 
load("results_ccube.rda", verbose = T)


inp$ccube_clusters <- results_ccube$res$label

dpclust_res <-  readr::read_tsv("output_dpclust/sample0_DPoutput_2000iters_1000burnin_seed123/sample0_2000iters_1000burnin_bestConsensusAssignments.bed")

inp$dpclust_clusters <- dpclust_res$cluster

pyclone_res <- readr::read_tsv("output_pyclone/example.tsv")

inp$pyclone_res <- pyclone_res$cluster_id + 1

consensus_structure <- readr::read_tsv("example_simulated/example_subclonal_structure.txt")

RHO = 0.01

lk_c1 <- VGAM::dbetabinom(x = inp$NV, size = inp$DP ,prob = consensus_structure$proportion[1] / 2, log = T, rho = RHO) + log(consensus_structure$n_ssms[1] / 100) 
lk_c2 <- VGAM::dbetabinom(x = inp$NV, size = inp$DP ,prob = consensus_structure$proportion[2] / 2, log = T, rho = RHO) + log(consensus_structure$n_ssms[2] / 100)


inp$consensus_cluster <- apply(cbind(lk_c1, lk_c2), MARGIN = 1, FUN = function(x) which.max(x))


p_top <-  readRDS("../BB_vs_B_weme/supp2_ab.rds")

p_middle <- readRDS("p_sim.rds")

df_b <- inp %>% select(CCF, ccube_clusters, dpclust_clusters, pyclone_res, consensus_cluster) %>%  
  rename(weme = consensus_cluster, pyclone = pyclone_res, ccube = ccube_clusters, dpclust = dpclust_clusters) %>% 
  reshape2::melt(id.vars = "CCF")

p_bottom <- ggplot(df_b, aes(x = CCF, fill = value %>% paste())) + geom_histogram(binwidth = 0.02) + 
  scale_fill_brewer("Cluster",palette = "Set2") + facet_wrap(.~variable, nrow = 1) + 
  theme_bw() + ggtitle("Weme consensus clustering with PCAWG tools")

p_S2 <- p_top / p_middle / p_bottom + patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold')) 

p_S2 %>% ggsave(., filename = "figure_S2.png", device = png, height = 2800, width = 2500, units = "px")
