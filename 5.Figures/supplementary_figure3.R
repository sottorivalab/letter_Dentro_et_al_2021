### Figure S2: Heathmap driver clustering probabilities
library(tidyverse)

cluster_probs_41_cases <- readr::read_csv("../overdispersion/cluster_probs_41cases.csv")

snvs_evo <-  readRDS("../overdispersion/letter_snvs.rds")

snvs_evo_bind <- do.call(rbind, snvs_evo)

snvs_evo_bind_drivers <- snvs_evo_bind %>% filter(isdriver) %>% 
  mutate(gene = paste0(gene,"_",substr(Tumor_Sample_Barcode,1,4) )) %>% 
  select(gene, Chromosome, Start_position, clonality, Tumor_Sample_Barcode) %>%  rename(chromosome = Chromosome, position = Start_position, sample = Tumor_Sample_Barcode)


snvs_evo_bind_drivers <- inner_join(cluster_probs_41_cases, snvs_evo_bind_drivers)

plot_inp <- snvs_evo_bind_drivers %>%  filter(mut_type == "SNV") %>%  select(starts_with("cluster"), gene, clonality) %>%  reshape2::melt(id.vars = c("gene", "clonality")) %>% 
  filter(variable != "cluster_4")
p4_1_1 <- ggplot(plot_inp %>%  filter(clonality == "clonal"), aes(fill = value, x = gene, y = variable)) + 
  geom_tile() + theme_minimal() + 
  scale_fill_distiller("prob",palette = "RdBu") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Clonal drivers cluster probability")

p4_1_2 <- ggplot(plot_inp %>%  filter(clonality == "subclonal"),  aes(fill = value, x = gene, y = variable)) + 
  geom_tile() + theme_minimal() + 
  scale_fill_distiller("prob",palette = "RdBu") + xlab("") + ylab("") +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Subclonal drivers cluster probability")


p4 <- p4_1_1 / p4_1_2  + plot_annotation(tag_levels = 'a')

p4  %>%  ggsave(filename = "figure_S2.png", device = "png", plot = ., width = 12, height = 6.9)



