set.seed(3)

source('sampler.R')

# Simulate data
simulated_data = simulate_data(overdispersion = 0.015, N = 5000, coverage = rpois(5000, 45))

cov_plot = ggplot(simulated_data %>% as_tibble) + 
  geom_histogram(aes(DP), bins = 30) +
  labs(x = 'Sequencing depth', y = 'Counts', title = "Coverage distribtion (45x)") +
  CNAqc:::my_ggplot_theme()

vaf_plot = ggplot(simulated_data %>% as_tibble) + 
  geom_histogram(aes(VAF), bins = 80) +
  labs(x = 'VAF', y = 'Counts', title = bquote("Variant Allele Frequencies ("*rho~'= 0.015)')) +
  CNAqc:::my_ggplot_theme() +
  xlim(0, 1)

# Fit Binomial
binomial_fit = fit_Binomial(simulated_data %>% select(NV, DP), return_pars = T)

bin_fit_plot = simulated_data %>% mutate(cluster = if_else(binomial_fit$Z >= 0.5, "B1" , "B2")) %>% 
  ggplot(aes(VAF, fill = cluster)) +
  geom_histogram(binwidth = 0.01)+
  scale_fill_brewer(palette = "Set1") +
  labs(subtitle = NULL, 
       title = bquote("Beta-Binomial clusters (K"["B"]~'= 2)'))


dir.create("method1")

tribble(
  ~"cluster", ~"n_ssms", ~"proportion",
  1,   round(binomial_fit$phi[1] * 5000), binomial_fit$p[1],
  2,   round(binomial_fit$phi[2] * 5000), binomial_fit$p[2]
) %>% 
  write_tsv("method1/test_subclonal_structure.txt")

# Fit BetaBinomial
bbinomial_fit = fit_BetaBinomial(simulated_data %>% select(NV, DP), return_pars = TRUE)

bbin_fit_plot = simulated_data %>% mutate(cluster = "BB1") %>% 
  ggplot(aes(VAF, fill = cluster)) +
  geom_histogram(binwidth = 0.01)+
  scale_fill_brewer(palette = "Set1") +
  labs(subtitle = NULL, 
       title = bquote("Beta-Binomial clusters (K"["BB"]~'= 1)'))


dir.create("method2")


tribble(
  ~"cluster", ~"n_ssms", ~"proportion",
  1,   5000, bbinomial_fit$p,
) %>% 
  write_tsv("method2/test_subclonal_structure.txt")

# Source WeMe
source("weme.R")

# Run WeMe
sids = find_sids()
sids %>% print()

# Modified WeMe version that returns the plot object and dumps and extra rda
wemw_consensus = genconsensus(sids, rounddown = FALSE)
wemw_consensus = wemw_consensus[[1]] + CNAqc:::my_ggplot_theme() +
  labs(
    title = "WeMe consensus K = 2"
  )
test_plot = readRDS('test_plot.rds')

cowplot::plot_grid(
  cov_plot,
  vaf_plot,
  bin_fit_plot,
  bbin_fit_plot,
  wemw_consensus,
  nrow = 2,
  ncol = 3,
  align = 'h',
  axis = 'tb'
)

panel_one = cowplot::plot_grid(
  cowplot::plot_grid(bin_fit_plot, bbin_fit_plot, ncol = 1, labels = c('d', 'e')),
  test_plot,
  nrow = 1,
  ncol = 2, 
  rel_widths = c(1, 2),
  labels = c('', 'f')
)

panel_two = cowplot::plot_grid(
  cov_plot,
  vaf_plot,
  wemw_consensus,
  nrow = 1,
  ncol = 3,
  align = 'h',
  axis = 'tb',
  labels = c('a', 'b', 'c')
)

image_plot = cowplot::plot_grid(panel_two, panel_one, nrow = 2, ncol = 1,
                   rel_heights = c(1, 2))

ggsave(plot = image_plot, 'Main_Text_figure.png', width = 11, height = 8)
