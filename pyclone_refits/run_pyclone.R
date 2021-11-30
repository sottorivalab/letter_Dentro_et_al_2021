run_pyclone <- function(x){
  
  # run the script in each directory and go back to the original wd
  setwd(x)
  if(!file.exists("pyclone_output.tsv"))
    tryCatch( system("bash ../pyclone.sh ")
      ,finally = setwd("..")
    ) 
  else 
    setwd("..")
}

plot_pyclone <- function(x){
  library(patchwork)
  library(tidyverse)
  # rplot ccf histogram in each directory and go back to the original wd
  setwd(x)
  if(!file.exists("pyclone_plot.pdf"))
  {
    
    inp <-  readr::read_tsv("pyclone_output.tsv") %>% separate(col = mutation_id,
                                                               into = c("chr", "from", "to", "ccf", "gene", "is_driver"),
                                                               sep = ":")
    
    # Clusters cellular prevalence plot
    p1 <-  ggplot(inp, aes(x = cellular_prevalence, fill = cluster_id %>% paste)) +
      geom_histogram() + theme_bw() + scale_fill_brewer("Cluster", palette = "Set1") + 
      ggtitle( paste0("Pyclone clustering ", inp$sample_id %>% unique()) )
    
    # As pyclone just returns the CCF per cluster we are using the PCAWG estimated CCF per mutation
    p2 <- ggplot(inp, aes(x = ccf %>% as.numeric, fill = cluster_id %>% paste)) +
      geom_histogram(bins = 100) + theme_bw() + scale_fill_brewer("Cluster", palette = "Set1") + 
      ggtitle( paste0("Pyclone CCF ", inp$sample_id %>% unique()) ) + xlab("CCF") 
      
    # Assamble and save the plots
    p12 <- p1 / p2
    
    p12 %>%  ggplot2::ggsave(filename = "pyclone_plot.pdf", units = "px", device = "pdf", width = 2400, height = 3600)
  }
  setwd("..")
}

# load "names" variable with used sample names
load("samples.rda")

easypar::run(lapply(names , list), FUN = run_pyclone, filter_errors = T, export = NULL,parallel = T,
             cores.ratio = 0.6)
easypar::run(lapply(names , list), FUN = plot_pyclone, filter_errors = F, export = NULL,parallel = T,
             cores.ratio = 0.6)
