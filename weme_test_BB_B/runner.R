set.seed(3)

library(fitdistrplus)

# Real values from the cohort SNPs in tumour BAMs
dispersion_PCAWG = readr::read_tsv("../rho_selected_cases.tsv")
dispersion_PCAWG = dispersion_PCAWG %>% pull(rho)
dispersion_PCAWG = dispersion_PCAWG[dispersion_PCAWG > 1e-8]

hist(dispersion_PCAWG, breaks = 100)

saveRDS(dispersion_PCAWG, file = "dispersion_PCAWG.rds")

# Parameters 
grid_samples = expand.grid(
  n_binomial_methods = 1:7,
  n_betabinomial_methods = 1:7,
  sample = 1:30
)

N_samples = grid_samples %>% nrow

### Fit gamma distribution to the rho to sample from

gamma <-  fitdistrplus::fitdist(dispersion_PCAWG, distr = "gamma")

shape <-  gamma$estimate[1]

rate <-  gamma$estimate[2]


hist(dispersion_PCAWG, breaks = 100, probability = T)
curve(dgamma(x, shape = shape, rate = rate),from = 0.005,to =  0.035, add = T, col = "red")

table_runs = easypar::run(FUN = 
   function(i)
  {
    set.seed(i)
    require(dplyr)

    source("sampler.R")

    overdispersion = runif(1, 0.001, qgamma(0.85, shape = shape, rate = rate))
    #overdispersion = rgamma(1, shape = shape, rate = rate)
    N = runif(1, 500, 5000) %>% round()
    coverage = rpois(1,lambda = 45) %>% round()

    smpl <- sampler(
      overdispersion = overdispersion,
      N = N,
      coverage = rpois(N, coverage),
      n_binomial_methods = grid_samples$n_binomial_methods[i],
      n_betabinomial_methods = grid_samples$n_betabinomial_methods[i],
      seed = i
    ) %>%
      mutate(
        overdispersion = overdispersion,
        N = N,
        coverage = coverage
      )
    smpl <- smpl %>% as_tibble() %>%  mutate(i = i)
    
    smpl
  },PARAMS = lapply(1:nrow(grid_samples), list),parallel = TRUE,cache = FALSE,
  export = c("grid_samples", "shape", "rate"), filter_errors = TRUE, 
  cores.ratio = 0.6)

saveRDS(table_runs, 'Table_runs.rds')

# Assembly results, remove errors?
# table_runs = readRDS('Table_runs.rds')

w_ok = sapply(table_runs, function(x) !inherits(x, 'error'))
table_runs = table_runs[w_ok]
table_runs = do.call(rbind, table_runs)

consensus_table = table_runs %>% 
  filter(method == "consensus") %>% 
  as_tibble() %>% 
  group_by(i) %>% 
  summarise(K = n())

grid_samples = grid_samples[w_ok, ]

grid_samples$K = consensus_table$K

grid_samples = grid_samples %>% 
  group_by(n_binomial_methods, n_betabinomial_methods) %>% 
  summarise(p = sum(K == 1)/n(), sd = sd(K))

test_plot = ggplot(grid_samples,
       aes(x = n_binomial_methods, y = n_betabinomial_methods, fill = p)) +
  geom_tile(width = .9, height = .9) +
  geom_text(aes(label = round(p * 100, 1)), size = 3) +
  # geom_text(aes(x = 7, y = 2), label = "0", color = 'orange') +
  scale_fill_distiller(palette = 'Spectral')+
  labs(
    title = paste0("WeMe consensus clustering (n = ", N_samples, ")"),
    subtitle = "30 simulations per square",
    x = "Binomial methods",
    y = "Beta-Binomial methods"
  ) +
  CNAqc:::my_ggplot_theme() +
  guides(fill = guide_colorbar("Monoclonal (%)   ", barwidth = unit(3,'cm')))  +
  scale_x_discrete(limits = 1:7) +
  scale_y_discrete(limits = 1:7) 

ggsave(filename = "weme_test.png", width = 5, height = 5)
  
saveRDS(test_plot, "test_plot.rds")

