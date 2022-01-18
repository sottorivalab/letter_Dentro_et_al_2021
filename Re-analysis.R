require(tidyverse)

# Results from then analysis
x =  readRDS("~/Downloads/all_1340_PCWAG_pycloneBB_clustered.rds")
x = x %>% as_tibble
x$ccf = as.numeric(x$ccf)
x$from = as.numeric(x$from)
x$to = as.numeric(x$to)
x$cluster_id = as.character(x$cluster_id)

# Plotting function for a single sample
single_plot = function(sample_id = "0bfe2ac9-0af5-c248-e050-11ac0d487e1c")
{
  x %>%
    filter(sample_id == !!sample_id) %>%
    ggplot(aes(ccf, fill = cluster_id)) +
    geom_histogram(binwidth = 0.01) +
    facet_wrap(~icgc_sample_id, scales = 'free')
}

# All samples analysed
samples = x$sample_id %>% unique

# Same random samples plot
# original_plot = x %>%
#   filter(sample_id %in% samples[1:100]) %>%
#   ggplot(aes(ccf, fill = cluster_id)) +
#   geom_histogram(binwidth = 0.01) +
#   facet_wrap(~icgc_sample_id, scales = 'free')

samples_tomerge = x %>%
  distinct(sample_id, cellular_prevalence, cluster_id) %>%
  filter(cellular_prevalence >.95) %>%
  group_by(sample_id) %>%
  filter(n() > 1) %>%
  pull(sample_id)

# Merging clones with CCF >.95
for(s in samples_tomerge)
{
  clusters_tomerge = x %>%
    filter(sample_id == s) %>%
    distinct(cellular_prevalence, cluster_id) %>%
    filter(cellular_prevalence >.95) %>%
    arrange(cellular_prevalence %>% desc) %>%
    pull(cluster_id)

  if(length(clusters_tomerge) > 1) {
    print(s)

    to_be = clusters_tomerge[1]
    to_swap = clusters_tomerge[2:length(clusters_tomerge)]

    new_cellularity = which(x$sample_id == s & x$cluster_id == to_be)
    new_cellularity = x$cellular_prevalence[new_cellularity[1]]
    where_about = which(x$sample_id == s & x$cluster_id %in% to_swap)

    x$cluster_id[where_about] = to_be
    x$cellular_prevalence[where_about] = new_cellularity
  }
}

# Size of clusters from mutation assignments, and CCF peaks from PyClone-VI
clusters_stat = x %>%
  group_by(sample_id, cluster_id, icgc_sample_id) %>%
  mutate(ccf = cellular_prevalence) %>%
  summarise(n = n(), ccf = ccf[1]) %>%
  group_by(sample_id) %>%
  mutate(prop = n/sum(n)) %>%
  arrange(sample_id, ccf)

# genuine small clonal clusters:  cases that have a clonal CCF peak (>.9)
# which is small (<5% of mutations), but that is unique meaning that the closest
# subclonal peak is <.8.
small_clonal = NULL

for(sample in clusters_stat$sample_id)
{
  this_sample = clusters_stat %>% filter(sample_id == sample)

  if(nrow(this_sample) == 1) next

  has_small_superclonal = (this_sample[nrow(this_sample), 'prop'] < 0.05) %>% as.logical()
  subclonal_peak_is_far = (
    abs(
      this_sample[nrow(this_sample) - 1, 'ccf'] - this_sample[nrow(this_sample), 'ccf']
      ) > 0.2
    ) %>% as.logical()

  if(has_small_superclonal & subclonal_peak_is_far)
    small_clonal = c(small_clonal, sample)
}

pdf("1.small_clonal_cluster.pdf", width = 5, height = 3)
sapply(small_clonal, function(x) single_plot(x) %>% print)
dev.off()

# Filter: remove clusters < 5% of overall mutation burden, unless they are small
# clonal clusters identified above
clusters_stat = clusters_stat %>%
  filter(prop > 0.05 | sample_id %in% small_clonal)

# CCF peaks distribution
clusters_stat %>%
  group_by(sample_id) %>%
  ggplot(aes(ccf)) +
  geom_histogram(binwidth = 0.01)
ggsave("2.1.distribution_CCF_peaks.pdf", scale = .5)

# Filter: annotate low-frequency clusters (tail), clonal clusters and all the
# rest as subclonal clusters
clusters_stat = clusters_stat %>%
  group_by(sample_id) %>%
  mutate(cluster = case_when(
    row_number() == n() ~ "clonal",
    row_number() == 1 & ccf <= 0.5 ~ "tail",
    TRUE ~ "subclonal"
  ))

# Statistics - get tumour types, remove histologies with <5 cases
t_types = readxl::read_xlsx('~/Downloads/pcawg_specimen_histology_August2016_v9.xlsx') %>%
  select(icgc_sample_id, histology_tier2, histology_tier3, project_code) %>%
  mutate(histology = paste(histology_tier2, histology_tier3)) %>%
  select(-histology_tier2, -histology_tier3)

renaming = c(
  `Kidney Renal cell carcinoma (proximal tubules)` = "Kidney renal cell carcinoma",
  `Bone/SoftTissue Sarcoma, bone` = "Sarcoma",
  `CNS Non-diffuse glioma` = "Non-diffuse glioma",
  `CNS Medulloblastoma` = "Medulloblastoma"
)

for(t in 1:nrow(t_types))
{
  if(t_types$histology[t] %in% names(renaming))
    t_types$histology[t] = renaming[t_types$histology[t]]
}

# t_types = x %>% distinct(sample_id, histology_abbreviation)

(clusters_stat$icgc_sample_id %in% t_types$icgc_sample_id) %>% all # must be TRUE

sample_stat = clusters_stat %>%
  group_by(sample_id, icgc_sample_id) %>%
  summarise(type = ifelse('subclonal' %in% cluster, "polyclonal", "monoclonal"), n_subclones = sum(cluster == 'subclonal')) %>%
  left_join(t_types, by = 'icgc_sample_id')

retain_samples = sample_stat %>%
  distinct(sample_id, histology) %>%
  group_by(histology) %>%
  filter(n() > 10) %>%
  pull(sample_id) %>%
  unique()

clusters_stat = clusters_stat %>% filter(sample_id %in% retain_samples)
sample_stat = sample_stat %>% filter(sample_id %in% retain_samples)
x = x %>% filter(sample_id %in% retain_samples)

sample_stat$type %>% table()
sample_stat$n_subclones %>% table()

x = x %>%
  left_join(
    clusters_stat %>% select(sample_id, cluster, cluster_id),
    by = c('sample_id', 'cluster_id')
  )

x = x %>%
  left_join(
    sample_stat %>% select(sample_id, histology)
  )

# Save drive assignment list
x %>%
  mutate(is_driver = as.logical(is_driver)) %>%
  filter(is_driver) %>%
  select(chr, from, to, cellular_prevalence, gene, sample_id, cluster_id, cluster_assignment_prob, histology, cluster) %>%
  saveRDS("Drivers_PyCloneVI.rds")

# CCF peaks distribution 0 by cluster type
clusters_stat %>%
  group_by(sample_id) %>%
  ggplot(aes(ccf, fill = cluster)) +
  geom_histogram(binwidth = 0.01)
ggsave("2.2.distribution_CCF_peaks_by_cluster.pdf", scale = .5)

# Plot all fits by cluster type
pdf("3.Fits_CCF_PCAWG_by_cluster.pdf", width = 5, height = 3)
for(s in samples)
  print(
    x %>%
      filter(sample_id == !!s) %>%
      ggplot(aes(ccf, fill = cluster)) +
      geom_histogram(binwidth = 0.01) +
      facet_wrap(~icgc_sample_id, scales = 'free') +
      scale_fill_manual(
        values = c(`tail` = 'red', `clonal` = 'black', `subclonal` = 'yellow')
      ) +
      CNAqc:::my_ggplot_theme()
  )
dev.off()

# Plot all fits with subclone
pdf("4.1.2_subclones.pdf", width = 5, height = 3)
sapply(
  clusters_stat %>%
    group_by(sample_id) %>%
    summarise(clones = 1 + sum(cluster == 'subclonal')) %>%
    filter(clones == 3) %>%
    pull(sample_id),
  function(x) single_plot(x) %>% print)
dev.off()

pdf("4.2.1_subclones.pdf", width = 5, height = 3)
sapply(
  clusters_stat %>%
    group_by(sample_id) %>%
    summarise(clones = 1 + sum(cluster == 'subclonal')) %>%
    filter(clones == 2) %>%
    pull(sample_id),
  function(x) single_plot(x) %>% print)
dev.off()

# Summary statistics - cluster types
# colors_clusters
donut_style = CNAqc:::my_ggplot_theme() +
  theme(
    legend.position = 'bottom',
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    plot.background = element_blank()
  )

cluster_donut_df = table(clusters_stat$cluster) %>%
  as.data.frame() %>%
  rename(class = Var1, n = Freq) %>%
  mutate(prop = n / sum(n), class = as.character(class)) %>%
  mutate(class = ifelse(class == 'tail', 'neutral tail or \nlow-frequency subclone', class)) %>%
  # mutate(class = ifelse(class == 'subclonal', 'subclonal', class)) %>%
  mutate(class = ifelse(class == 'subclonal', 'subclone', class)) %>%
  mutate(class) %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5 * prop) %>%
  mutate(class = paste0(class, ' (n = ', n, ')')) %>%
  arrange(class)

clusters_donut = cluster_donut_df %>%
  ggplot(aes(x = 2, y = prop, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(y = lab.ypos, label = paste0(round(prop, 2) * 100, '%')), color = "black") +
  scale_fill_manual(
    values = pio:::nmfy(
      cluster_donut_df$class,
      ggsci::pal_uchicago()(3)
    )
  ) +
  donut_style +
  labs(title = paste0(
    'Type of cluster'
  )) +
  guides(fill = guide_legend(paste0('n = ', nrow(clusters_stat)), ncol = 1))


# Summary statistics - tumour types
architecture_donut_df = table(sample_stat$type) %>%
  as.data.frame() %>%
  rename(class = Var1, n = Freq) %>%
  mutate(prop = n/sum(n)) %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) %>%
  mutate(class = paste0(class, ' (n = ', n, ')')) %>%
  arrange(class)

architecture_donut = architecture_donut_df %>%
  ggplot(aes(x = 2, y = prop, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(y = lab.ypos, label = paste0(round(prop, 2) * 100, '%')), color = "black")+
  scale_fill_manual(
    values = pio:::nmfy(
      architecture_donut_df$class,
      c( "#0073C2FF",  "#EFC000FF")
    )
  ) +
  donut_style +
  # xlim(0.5, 2.5) +
  labs(title = paste0('Tumour classification')) +
  guides(fill = guide_legend(paste0(
    'n = ', architecture_donut_df$n %>% sum
  ), ncol = 1))

# Summary statistics - n subclones
tumour_donut_df = table(sample_stat$n_subclones) %>%
  as.data.frame() %>%
  rename(class = Var1, n = Freq) %>%
  mutate(prop = n/sum(n)) %>%
  arrange(desc(class)) %>%
  mutate(class = paste0(class, ' (n = ', n, ')')) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop) %>%
  arrange(class)

tumour_donut =  tumour_donut_df %>%
  ggplot(aes(x = 2, y = prop, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(y = lab.ypos, label = paste0(round(prop, 2) * 100, '%')), color = "black", size = 3)+
  scale_fill_manual(values = ggsci::pal_locuszoom()(3) %>% rev) +
  donut_style +
  # xlim(0.5, 2.5) +
  labs(title = paste0('Number of subclones')) +
  guides(fill = guide_legend(paste0(
    'n = ', tumour_donut_df$n %>% sum
  ), ncol = 1))

# Histology-based classification
histology_ordering = sample_stat %>%
  group_by(histology) %>%
  summarise(n = n()) %>%
  arrange((n)) %>%
  pull(histology)

left_panel = sample_stat %>%
  group_by(histology, type) %>%
  summarise(n = n()) %>%
  group_by(histology) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x = histology, prop, fill = type)) +
  geom_bar(stat = 'identity') +
  # facet_wrap(~histology_abbreviation) +
  CNAqc:::my_ggplot_theme() +
  scale_fill_manual(values = c(`monoclonal` = "#0073C2FF", `polyclonal` = "#EFC000FF"))  +
  coord_flip() +
  scale_x_discrete(limits = histology_ordering) +
  labs(x = 'Histology', y = "Proportion", title =  "Clonality by histology") +
  guides(fill = guide_legend('Clonality', ncol = 1))

right_panel = sample_stat %>%
  ggplot(aes(x = histology, fill = n_subclones %>% paste())) +
  geom_bar() +
  # facet_wrap(~histology_abbreviation) +
  CNAqc:::my_ggplot_theme() +
  scale_fill_manual(values = mycols)  +
  coord_flip() +
  scale_x_discrete(limits = histology_ordering) +
  guides(fill = guide_legend('Number of \nsubclones', ncol = 1)) +
  labs(y = "Cases", x = NULL)  +
  theme(axis.text.y = element_blank())

heatmap_histology = cowplot::plot_grid(left_panel, right_panel, ncol = 2, align = 'h', axis = 'tb', rel_widths = c(1,.3))

figure = cowplot::plot_grid(
    cowplot::plot_grid(
      clusters_donut,
      architecture_donut,
      tumour_donut,
      ncol = 3, align = 'h', axis = 'tb', labels = c('a', 'b', 'c')),
    heatmap_histology,
    ncol = 1,
    rel_heights = c(1, 3),
    labels = c('', 'd')
)

ggsave(figure, filename = 'Figure_clustering.png', width = 12, height = 10)


# Get drivers
drivers = readRDS("Drivers_PyCloneVI.rds") %>%
  filter(!is.na(cluster)) %>%
  mutate(cluster = ifelse(cluster == 'tail', 'neutral tail or \nlow-frequency subclone', cluster))

# Split by group
counts_hc_drivers = drivers %>%
  group_by(gene, cluster, histology) %>%
  summarise(n = n())

high_freq = counts_hc_drivers %>% filter(n > 100) %>% pull(gene)
mid_freq = counts_hc_drivers %>% filter(n > 20, n <= 100, !(gene %in% c(high_freq))) %>% pull(gene)
low_freq = counts_hc_drivers %>% filter(n <= 20, n > 5, !(gene %in% c(mid_freq, high_freq))) %>% pull(gene)

# Lowfreq
p1 = ggplot(counts_hc_drivers %>% filter(gene %in% low_freq)) +
  geom_tile(
    aes(x = gene, y = histology, fill = n)
  ) +
  coord_flip() +
  CNAqc:::my_ggplot_theme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~cluster, nrow = 1) +
  scale_fill_viridis_c(direction = -1, option = "B") +
  # scale_fill_gradientn(colours = c("darkgoldenrod2", "orange", 'indianred3')) +
  # guides(fill = guide_colorbar(barwidth = unit(6, 'cm')))+
  # theme(axis.text.x = element_text(angle = 45))
  guides(fill = guide_colorbar("Cases", barheight = unit(3, 'cm'))) +
  theme(axis.text.x = element_text(angle = 45), legend.position = 'right') +
  labs(x = "Low-frequency\n(>5 and <=20)", y = "Tumour histology")

# Mid
p2 = ggplot(counts_hc_drivers %>% filter(gene %in% mid_freq)) +
  geom_tile(
    aes(x = gene, y = histology, fill = n)
  ) +
  coord_flip() +
  CNAqc:::my_ggplot_theme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~cluster, nrow = 1) +
  scale_fill_viridis_c(direction = -1, option = "B") +
  # scale_fill_gradientn(colours = c("darkgoldenrod2", "orange", 'indianred3')) +
  guides(fill = guide_colorbar("Cases", barheight = unit(3, 'cm'))) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = 'right') +
  labs(x = "Mid-frequency\n(>20 and <=100)")

# High
p3 = ggplot(counts_hc_drivers %>% filter(gene %in% high_freq)) +
  geom_tile(
    aes(x = gene, y = histology, fill = n)
  ) +
  coord_flip() +
  CNAqc:::my_ggplot_theme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~cluster, nrow = 1) +
  scale_fill_viridis_c(direction = -1, option = "B") +
  # scale_fill_viridis_c(direction = -1, option = "B") +
  # scale_fill_gradientn(colours = c("forestgreen", "orange", "indianred3", 'purple4')) +
  # scale_fill_gradientn(colours = c("darkgoldenrod2", "orange", 'indianred3')) +
  guides(fill = guide_colorbar("Cases", barheight = unit(3, 'cm'))) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = 'right') +
  labs(x = "High-frequency\n(>100)")

figure = cowplot::plot_grid(p3,p2,p1, ncol = 1, rel_heights = c(1.5,2.5,10), axis = 'lr', align = 'v')
ggsave(figure, filename = "Figure_drivers_heatmap.png", height = 14, width = 12)

# drivers_clus_freq_df %>%
#   ggplot() +
#   geom_histogram(aes(n, fill = cluster), binwidth =  5)+
#   CNAqc:::my_ggplot_theme() +
#   labs(title = "Driver occurrence (pan-cancer)", x = 'Occurrences', y = "Observations") +
#   scale_fill_manual(
#     values = pio:::nmfy(
#       drivers$cluster %>% unique(na.rm = T) %>% sort,
#       ggsci::pal_uchicago()(3)
#     )
#   ) +
#   ggrepel::geom_label_repel(
#     data = drivers_clus_freq_df %>% filter(n>50),
#     aes(x = n, y = 50, label = gene, fill = cluster),
#     size = 3, nudge_y = 0, nudge_x = 0, segment.color = 'black'
#   )
# drivers_mco_his = drivers %>%
#   filter(
#     gene %in% mobster::cancer_genes_dnds$Martincorena_drivers
#   ) %>%
#   group_by(gene, cluster, histology) %>%
#   summarise(n = n())
#
# selected = drivers_mco_his %>% filter(n > 5) %>% pull(gene)
#
#
# drivers_clus_freq_df
#
# ggplot(
#   drivers_mco_his %>% filter(gene %in% selected)
#   ) +
#   geom_tile(
#     aes(x = gene, y = histology, fill = n)
#   ) +
#   coord_flip() +
#   CNAqc:::my_ggplot_theme() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   facet_wrap(~cluster, nrow = 1) +
#   scale_fill_viridis_c(direction = -1) +
#   # scale_fill_gradientn(colours = c("darkgoldenrod2", "orange", 'indianred3')) +
#   guides(fill = guide_colorbar(barwidth = unit(6, 'cm')))
#
#
#
# ggplot(
#   drivers_mco_his %>% filter(gene %in% selected)
# ) +
#   geom_bar(
#     aes(x = gene, n , fill = cluster),
#     stat = 'identity'
#   ) +
#   coord_flip() +
#   CNAqc:::my_ggplot_theme() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   facet_wrap(~histology_abbreviation, scales = 'free') +
#   # scale_fill_viridis_c(direction = -1) +
#   # scale_fill_gradientn(colours = c("darkgoldenrod2", "orange", 'indianred3')) +
#   guides(fill = guide_colorbar(barwidth = unit(6, 'cm')))
#
# drivers_mco_his %>%
#   filter(gene %in% selected) %>%
#   ggplot() +
#     geom_bar(aes(x = histology_abbreviation, n, fill = cluster), stat = 'identity') +
#     facet_wrap(~gene) +
#     CNAqc:::my_ggplot_theme() +
#     coord_flip()

# CCF by subclone
p1 = drivers %>%
  # filter(cluster == 'subclonal') %>%
  ggplot(aes(x = cellular_prevalence, y = cluster_assignment_prob, color = cluster)) +
  geom_hline(yintercept = 0.8, color = 'red', linetype = 'dashed')  +
  geom_point() +
  CNAqc:::my_ggplot_theme() +
  scale_color_manual(
    values = pio:::nmfy(
      drivers$cluster %>% unique %>% sort,
      ggsci::pal_uchicago()(3)
    )
  ) +
  ylim(0,1) +
  xlim(0,1) +
  labs(x = "Clone peak", y = "Assignment probability") +
  guides(color = guide_legend('Cluster'))

p2 = drivers %>%
  ggplot(aes(cluster_assignment_prob, fill = cluster)) +
  # geom_hline(yintercept = 0.8, color = 'red', linetype = 'dashed')  +
  geom_histogram(binwidth = 0.01) +
  CNAqc:::my_ggplot_theme() +
  scale_fill_manual(
    values = pio:::nmfy(
      drivers$cluster %>% unique %>% sort,
      ggsci::pal_uchicago()(3)
    )
  ) +
  # ylim(0,1) +
  xlim(-0.1,1.01) +
  labs(x = "Clone peak", y = "Assignment probability") +
  guides(color = guide_legend('Cluster'))  +
  labs(y = "Counts", x = "Assignment probability") +
  guides(fill = guide_legend('Cluster'))

cutf = function(m){
  nc = table(drivers$cluster)
  drivers %>%
    filter(cluster_assignment_prob > m) %>%
    group_by(cluster) %>%
    summarise(n = n()) %>%
    mutate(cutoff = m, prop = n/nc[cluster])
}

p3 = Reduce(bind_rows, lapply(seq(0, .99, 0.1), cutf)) %>%
  ggplot(aes(x = cutoff, y = prop, color = cluster)) +
  geom_vline(xintercept = 0.8, color = 'red', linetype = 'dashed')  +
  geom_line() +
  geom_point(size = 4) +
  CNAqc:::my_ggplot_theme() +
  scale_color_manual(
    values = pio:::nmfy(
      drivers$cluster %>% unique %>% sort,
      ggsci::pal_uchicago()(3)
    )
  ) +
  # ylim(0,1) +
  xlim(-0.1,1.01) +
  labs(x = "Clone peak", y = "Assignment probability") +
  guides(color = guide_legend('Cluster'))  +
  labs(y = "Retained drivers (%)", x = "Assignment probability cutoff") +
  guides(fill = guide_legend('Cluster'))

# p3 = Reduce(bind_rows, lapply(seq(0, 1, 0.1), cutf)) %>%
#   ggplot(aes(x = cutoff, y = n, color = cluster)) +
#   geom_vline(xintercept = 0.8, color = 'red', linetype = 'dashed')  +
#   geom_point() +
#   CNAqc:::my_ggplot_theme() +
#   scale_color_manual(
#     values = pio:::nmfy(
#       drivers$cluster %>% unique %>% sort,
#       ggsci::pal_uchicago()(3)
#     )
#   ) +
#   # ylim(0,1) +
#   xlim(-0.1,1.01) +
#   facet_wrap(~cluster, scales = 'free', ncol = 1) +
#   labs(x = "Clone peak", y = "Assignment probability") +
#   guides(color = guide_legend('Cluster'))  +
#   labs(y = "Counts", x = "Assignment probability") +
#   guides(fill = guide_legend('Cluster'))

# x %>% filter(is_driver) %>% group_by(sample_id, cluster)

figure = cowplot::plot_grid(p1,p2,p3, ncol = 3, labels = c('a', 'b', 'c'))
ggsave(figure, filename = "Figure_drivers_uncertainty.png", height = 4, width = 12)


"b37d6283-6f95-4975-a794-f3d5c4bbc7b3" %>% single_plot()
"2d0e4b82-c623-11e3-bf01-24c6515278c0" %>% single_plot()
"72f82fbd-9838-4082-b605-bc3d80226f16" %>% single_plot()


# Unexplained cases
subclonal_cases = x %>%
  filter(cluster =='subclonal') %>%
  pull(sample_id) %>%
  unique()

subclonal_cases_dr = drivers %>%
  filter(cluster =='subclonal') %>%
  pull(sample_id) %>%
  unique()

all(subclonal_cases_dr %in% subclonal_cases) # TRUE

length(subclonal_cases_dr)/length(subclonal_cases)

non_explained = setdiff(subclonal_cases, subclonal_cases_dr)


pdf("5.no_driver.pdf", width = 5, height = 3)
sapply(non_explained, function(x) single_plot(x) %>% print)
dev.off()
