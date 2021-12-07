### FUNCTION TO PLOT CCF HISTOGRAMS ###

plot_hist <- function(data) {
  
  library(ggrepel)
  
  clusters <- data %>% group_by(cluster) %>%  summarise(mean_ccf = mean(ccf.y %>% as.numeric)) %>% arrange(-mean_ccf)
  
  clusters$new_clusters <-  seq(1:nrow(clusters)) %>%  paste
  
  data_1 <- inner_join(data,clusters) %>% select(ccf.y, new_clusters, name, is_driver, gene.y, clonality)
  
  data_1$pipeline <-  "PCWAG"
  
  clusters_2 <- data %>% group_by(cluster_id) %>%  summarise(mean_ccf = cellular_prevalence %>% unique()) %>% arrange(-mean_ccf)
  clusters_2$new_clusters <-  seq(1:nrow(clusters_2)) %>%  paste
  
  data_2 <- inner_join(data,clusters_2)%>% select(ccf.y, new_clusters, name, is_driver, gene.y, clonality)
  
  data_2$pipeline <-  "Pyclone"
  
  data_plot <- rbind(data_1,data_2)
  
  ggplot(data = data_plot, aes(x = ccf.y %>% as.numeric, fill = new_clusters %>%  paste)) +
    geom_histogram(bins = 100, position = "stack", alpha = 0.6) +
    geom_text_repel(data = data_plot %>%  filter(is_driver != "NA",  clonality == "subclonal"),
                    aes(label = gene.y, y = Inf, x = ccf.y %>%  as.numeric), colour = "black",
                    label.size = 0.01, alpha = 0.85) +
    geom_vline(data = data_plot %>%  filter(is_driver != "NA", clonality == "subclonal"),
               aes(xintercept = ccf.y %>% as.numeric, colour = new_clusters %>%  paste), lty = 2, show.legend = FALSE) +
    scale_fill_brewer("cluster",palette = "Dark2", aesthetics = c("colour", "fill")) +
    facet_wrap(.~pipeline, scales = "free_y") + ggtitle(data_plot$name %>% unique()) +
    theme_bw() + theme(legend.position = "None", plot.background  = element_rect(colour = "red4", fill=NA, size=2), plot.title = element_text(hjust = 0.5)) +
    xlim(c(0,2)) +
    xlab("CCF") +
    ylab("Counts") 
}


# As this data is closed access we cannot share it, we provide the columns you need from the original PCWAG data release 
# to make this script work, you also need results from pyclone as obtained by our pipeline

#all_pyclone <- readRDS("../overdispersion/letter_pyclone_BB_all.rds")
#all_PCWAG <- data.table::fread("../overdispersion/maf_ccf_withdriver.csv", data.table = F)

# This file is publicly downloadable at this link https://dcc.icgc.org/releases/PCAWG/evolution_and_heterogeneity
annots <- data.table::fread("../overdispersion/icgc_sample_annotations_summary_table.txt", data.table=F) %>% rename(sample_id = tumour_aliquot_id) %>%
  mutate(name = paste0(histology_abbreviation, "\n",icgc_sample_id)) %>%
  select(sample_id, name)

# We have a subset of manually choosen 41 samples
load("samples_letter.rda", verbose = T)

all_pyclone_letter <- all_pyclone %>%  do.call(rbind,.) %>% filter(sample_id %in% nms_letter)

p3_1_df <- all_pyclone_letter %>% separate(col = mutation_id,
                                           into = c("chr", "from", "to", "ccf", "gene", "is_driver"),
                                           sep = ":")
p3_1_df <- p3_1_df %>%  dplyr::inner_join(annots, .)

p3_1_df <- dplyr::inner_join(all_PCWAG %>% mutate(bind_id = paste(paste0("chr",Chromosome), Start_position, Start_position, round(ccf,3), sep = ":")),
                             p3_1_df %>%  mutate(bind_id = paste(chr, from, to, ccf, sep = ":")), by ="bind_id")

p3_1_df <-  p3_1_df %>%  filter(Variant_Type== "SNP")


## Corresponding supplemetnary figure

p3_supp <- ggplot(data = p3_1_df, aes(x = ccf.y %>% as.numeric, fill = cluster_id %>%  paste)) +
  geom_histogram(bins = 100, position = "stack", alpha = 0.6) +
  geom_text_repel(data = p3_1_df %>%  filter(is_driver != "NA",  clonality == "subclonal"),
                  aes(label = gene.y, y = Inf, x = ccf.y %>%  as.numeric), colour = "black",
                  label.size = 0.01, alpha = 0.85) +
  geom_vline(data = p3_1_df %>%  filter(is_driver != "NA", clonality == "subclonal"),
             aes(xintercept = ccf.y %>% as.numeric, colour = cluster_id %>%  paste), lty = 2, show.legend = FALSE) +
  scale_fill_brewer("cluster",palette = "Dark2", aesthetics = c("colour", "fill")) +
  facet_wrap(.~name, scales = "free_y") + 
  theme_bw() + theme(legend.position = "None") +
  xlim(c(0,2)) +
  xlab("CCF") +
  ylab("Counts") + ggtitle("PCAWG consensus clustering")


p3_supp %>% ggsave(.,filename = "figure_S4.png", device = "png", units = "px", width = 3200, height = 2500)



# We assamble the 6 plots of figure a manually

samples_to_plot <- c("Panc-AdenoCA\nSA410399", "Panc-AdenoCA\nSA410978", "Panc-AdenoCA\nSA520252", "Panc-AdenoCA\nSA520283",
                     "Biliary-AdenoCA\nSA543550", "Liver-HCC\nSA501649")

p_3_1_1 <- plot_hist(p3_1_df %>% filter(name == samples_to_plot[1]))

p_3_1_2 <- plot_hist(p3_1_df %>% filter(name == samples_to_plot[2]))

p_3_1_3 <- plot_hist(p3_1_df %>% filter(name == samples_to_plot[3]))

p_3_1_4 <- plot_hist(p3_1_df %>% filter(name == samples_to_plot[4]))

p_3_1_5 <- plot_hist(p3_1_df %>% filter(name == samples_to_plot[5]))

p_3_1_6 <- plot_hist(p3_1_df %>% filter(name == samples_to_plot[6]))

library(patchwork)

p3_1 <- (p_3_1_1 | p_3_1_2 |p_3_1_3 ) / (p_3_1_4 | p_3_1_5 | p_3_1_6)

#p5_1 %>%  ggsave(filename = "p5_1.pdf", device = "pdf")

p3_2_df1 <-  all_PCWAG %>% filter(isdriver == T,mut_type == "SNV" ) %>%  group_by(clonality) %>%  summarize(N = n()) %>% filter(clonality %in% c("clonal", "subclonal"))

p3_2_df2 <- all_pyclone %>%  do.call(rbind,.) %>% separate(col = mutation_id,
                                                           into = c("chr", "from", "to", "ccf", "gene", "is_driver"),  sep = ":")

p3_2_df2 <- dplyr::inner_join(all_PCWAG %>% select(-clonality, - gene )%>%  mutate(bind_id = paste(paste0("chr",Chromosome), Start_position, Start_position, round(ccf,3), sep = ":")),
                              p3_2_df2 %>%  mutate(bind_id = paste(chr, from, to, ccf, sep = ":")), by ="bind_id")


# Here we somehow copy the assignment criteria of the Dentro paper, a pyclone cluster is clonal if it is the max(cellular prevalence) clsuter
# or it has cellular > 90%  

p3_2_df2 <- p3_2_df2 %>% group_by(sample_id) %>% mutate(clonal_ccfs = max(cellular_prevalence)) %>%  filter(Variant_Type == "SNP")                                                                           
p3_2_df2 <-  p3_2_df2 %>% mutate(clonality = if_else(cellular_prevalence > 0.9, "clonal", "subclonal")) %>% mutate(clonality = if_else(cellular_prevalence == clonal_ccfs,
                                                                                                                                       "clonal", clonality))
p3_2_df2_aggr <- p3_2_df2 %>% filter(is_driver != "NA") %>%  group_by(clonality) %>%  summarize(N = n()) %>% filter(clonality %in% c("clonal", "subclonal"))

p3_2_df2_aggr$pipeline <- "PyClone"
p3_2_df1$pipeline <-  "PCWAG"

p3_2 <- ggplot(data =  rbind(p3_2_df2_aggr, p3_2_df1)%>% filter(clonality == "subclonal"),
               mapping = aes(x = pipeline , y = N)) +
  geom_col(colour = "black", size = 0.3, fill = "gainsboro", alpha = 0.8) +
  theme_bw() + theme(legend.position = "None") +
  xlab("") + ylab("#drivers") + 
  ggtitle("Subclonal SNV drivers in 1340 samples")


#p3_2 %>%  ggsave(filename = "p3_2.pdf", device = "pdf")

p3_3_df1 <- p3_2_df2 %>%  filter(sample_id %in% nms_letter, is_driver != "NA", mut_type == "SNV" ) %>%  group_by(clonality) %>%  summarize(N = n()) %>% filter(clonality %in% c("clonal", "subclonal"))
p3_3_df2 <- all_PCWAG %>% filter(isdriver == T,mut_type == "SNV" ) %>%  filter(Tumor_Sample_Barcode %in% nms_letter) %>%  group_by(clonality) %>%  summarize(N = n()) %>% filter(clonality %in% c("clonal", "subclonal"))

p3_3_df1$pipeline <- "PyClone"
p3_3_df2$pipeline <-  "PCWAG"

p3_3 <- ggplot(data =  rbind(p3_3_df1, p3_3_df2) %>% filter(clonality == "subclonal"),
               mapping = aes(x = pipeline , y = N )) +
  geom_col(colour = "black", size = 0.3, alpha = 0.8, fill = "gainsboro") +
  theme_bw() + theme(legend.position = "None") +
  xlab("") + ylab("#drivers") + 
  ggtitle("Subclonal SNV drivers in 41 samples with clonal split")

#p3_3 %>%  ggsave(filename = "p3_3.pdf", device = "pdf")

p3_4_1 <-  all_pyclone %>%  do.call(rbind,.) 

p3_4_1 <-  p3_4_1 %>% group_by(sample_id) %>% mutate(mut_tot = n()) 


### HERE WE REMOVE THE PYCLONE CLUSTERS WITH MIXING PROPORTION LESS THAN 0.01 AS IN THE DENTRO ANALYSIS ###
p3_4_1 <-  p3_4_1 %>% group_by(sample_id, cluster_id) %>% mutate(mixing_prop = n() / mut_tot) %>%  mutate(to_remove = if_else(mixing_prop < 0.01 , T, F)) 

p3_4_1_aggr <-  p3_4_1 %>% filter(!to_remove) %>% group_by(sample_id) %>% summarize(N = length(unique(cluster_id)))

p3_4_2_aggr <-  all_PCWAG %>% group_by(Tumor_Sample_Barcode) %>% summarize(N = length(unique(cluster)))

p3_4_2_aggr <-  p3_4_2_aggr %>% rename(sample_id = Tumor_Sample_Barcode)

p3_4_2_aggr$pipeline <- "PCWAG"
p3_4_1_aggr$pipeline <- "PyClone"


plot_df <-  rbind(p3_4_2_aggr, p3_4_1_aggr) %>%  mutate(N = N -1) %>% 
  mutate(N = if_else(N > 3, ">3", N %>%  paste))

plot_df$N <- factor(plot_df$N, levels = c("0","1", "2", "3", ">3"))
p3_4 <- ggplot(data =  plot_df, mapping = aes(x = N, fill = N %>%  paste())) +
  geom_bar(colour = "white") +
  scale_fill_brewer("Clonality", palette = "Blues") +
  theme_bw() + theme(legend.position = "None") +
  facet_wrap(.~pipeline) +
  xlab("") + ylab("#subclones") + 
  ggtitle("Number of subclones in 1340 samples")


snvs_evo_bind_drivers <- readRDS("../overdispersion/drivers_as_prob.rds") %>%  mutate(prob = if_else(clonality == "clonal", p_clonal, 1 - p_clonal))

p3_5 <- ggplot(snvs_evo_bind_drivers , aes(x = clonality, y = prob, color = clonality)) + 
  geom_quasirandom(alpha = 0.8) + theme_minimal() + 
  scale_color_brewer(palette = "Set2") + xlab("") + 
  stat_compare_means(label.x = 1.3,label.y = 1.08) + ylim(0.35,1.15) + ggtitle("Cluster probability distribution")

df3_6 <- snvs_evo_bind_drivers %>%  filter(clonality == "subclonal") %>%  
  mutate(">50" = if_else(prob > 0.5 , T, F),
         ">60" = if_else(prob > 0.6 , T, F),
         ">70" = if_else(prob > 0.7 , T, F),
         ">80" = if_else(prob > 0.8 , T, F),
         ">90" = if_else(prob > 0.9 , T, F)) %>% select(starts_with(">"), gene, clonality) %>%  reshape2::melt(id.vars	
                                                                                   = c("gene", "clonality"))

p3_6 <- ggplot(df3_6 %>%  group_by(clonality, variable) %>%  summarize(N = sum(value)),
               aes(x = variable, y = N)) + 
  theme_bw() + geom_point(color = "black", alpha = 0.85) + xlab("Probability threshold") + ylab("#drivers") + ggtitle("Number of drivers with different probability thresholds")



p3 = p3_1 / (p3_2 | p3_3 | p3_4) / (p3_5 | p3_6) + patchwork::plot_layout(nrow =  4) + patchwork::plot_annotation(tag_levels = "a") &
                                                                                theme(plot.tag = element_text(face = 'bold')) 

p3 %>%  ggsave(.,filename = "figure_3.png", device = "png", units = "px", width = 4800, height = 4800)

