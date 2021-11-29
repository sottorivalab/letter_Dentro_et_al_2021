simulate_data = function(overdispersion = 0.01,
                         N = 1000,
                         coverage = rpois(N, 45))
{
  data.frame(
    DP = coverage,
    NV = VGAM::rbetabinom(
      n = N,
      size = coverage,
      prob = 0.5,
      rho = overdispersion
    )
  ) %>%
    mutate(VAF = NV / DP)
}

fit_Binomial = function(x)
{
  toString_Binomial_mix = BMix::bmixfit(
    x %>% dplyr::select(NV, DP) %>% as.data.frame(), 
    K.Binomials = 1:2, 
    K.BetaBinomials = 0, 
    score = "BIC"
  ) %>% 
    BMix::to_string() 
  
  if(!is.null(toString_Binomial_mix$Pi_Bin_2))
  tribble(
    ~"cluster", ~"n_ssms", ~"proportion",
    1,   toString_Binomial_mix$N_Bin_1, toString_Binomial_mix$Mean_Bin_1,
    2,   toString_Binomial_mix$N_Bin_2, toString_Binomial_mix$Mean_Bin_2
  )
  
else
  tribble(
    ~"cluster", ~"n_ssms", ~"proportion",
    1,   toString_Binomial_mix$N_Bin_1, toString_Binomial_mix$Mean_Bin_1,
  )
}

fit_BetaBinomial = function(x)
{
  toString_BetaBinomial_mix = BMix::bmixfit(
    x %>% dplyr::select(NV, DP) %>% as.data.frame(), 
    K.Binomials = 0, 
    K.BetaBinomials = 1, 
    score = "BIC"
  ) %>% 
    BMix::to_string() 
  
  tribble(
    ~"cluster", ~"n_ssms", ~"proportion",
    1,   toString_BetaBinomial_mix$N_BBin_1, toString_BetaBinomial_mix$Mean_BBin_1,
  ) 
}

sampler = function(
  overdispersion = 0.01,
  N = 1000,
  coverage = rpois(N, 45),
  n_binomial_methods = 7,
  n_betabinomial_methods = 2
)
{
  rndm_folder = sample(LETTERS, 16) %>% paste(collapse = '')
  dir.create(rndm_folder)
  
  # Simulate data
  dataset = simulate_data(overdispersion, N, coverage)
  
  # Fit Binomial
  binomial_methods = lapply(
    1:n_binomial_methods,
    function(x) fit_Binomial(dataset)
  )
  
  # Fit BetaBinomial
  betabinomial_methods = lapply(
    1:n_betabinomial_methods,
    function(x) fit_BetaBinomial(dataset)
  )
  
  n_methods = n_betabinomial_methods + n_binomial_methods
  
  # Dump data
  sapply(1:n_methods, function(i) {
    output_folder = paste0('./', rndm_folder, "/method", i)
    output_file = paste0(output_folder, '/test_subclonal_structure.txt')
    
    dir.create(output_folder, recursive = TRUE)
    
    if(i <= n_binomial_methods)
      readr::write_tsv(binomial_methods[[i]], output_file)
    else
      readr::write_tsv(betabinomial_methods[[i- n_binomial_methods]], output_file)
  })
  
  # Source WeMe
  source("/Users/salvatore.milite/work/Trieste/analysis/overdispersion/weme_test_BB_B/weme.R")
  
  # Run WeMe
  setwd(rndm_folder)
  
  sids = find_sids()

  genconsensus(sids, rounddown = FALSE, ncores = 1)
  
  load("./results.rda")
  
  setwd('..')
  
  result = result  %>% distinct(phi, method)
  
  # result$model = NA
  # result$model[1:n_binomial_methods] = 'Binomial'
  # result$model[n_binomial_methods + 1:n_betabinomial_methods] = 'BetaBinomial'
  
  result
}
