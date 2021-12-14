simulate_data = function(overdispersion = 0.01,
                         N = 1000,
                         coverage = rpois(N, 45), seed = 3, purity = 1)
{
  data.frame(
    DP = coverage,
    NV = VGAM::rbetabinom(
      n = N,
      size = coverage,
      prob = 0.5 * purity,
      rho = overdispersion
    )
  ) %>%
    mutate(VAF = NV / DP)
}

fit_Binomial = function(x, seed = 3, return_pars = FALSE)
{
  
  library(matrixStats)
  library(VGAM)
  
  set.seed(seed)
  
  input <- x %>% dplyr::select(NV, DP) %>% as.data.frame()
  
  LL <- function(prob) {
    R = dbinom(x = input$NV, prob = prob, size = input$DP)
    
    -sum(log(R))
  }
  
  b1 <- stats4::mle(LL, start = list(prob = 0.5), lower = 0.001, upper = 0.999)
  
  p_b1 <-  b1@coef[1]
  
  phi1 = runif(1)
  
  phi = c(phi1, 1 - phi1)
  
  p = runif(2)
  
  E_tol = 1e-2
  
  S = x$NV
  N = x$DP
  
  
  lk1 = dbinom(S, p[1], size = N, log = T) + log(phi[1])  
  lk2 = dbinom(S, p[2], size = N, log = T) + log(phi[2])
  lk = apply(cbind(lk1, lk2), 1,  logSumExp, simplify = TRUE)
  tol = Inf
  
  idx = 1
  
  while(tol > E_tol | idx < 3) {
    

    # E step
    z = lk1 - lk
    z = exp(z)
    
    #M step
    
    phi1 = sum(z) / length(N)  
    phi = c(phi1, 1 - phi1)
    
    p1 = sum(S * z) / sum(N * z)
    p2 = sum(S * (1 - z)) / sum(N * (1 - z))
    
    p = c(p1, p2)
    
    lk1 = dbinom(S, p[1], size = N, log = T) + log(phi[1])  
    lk2 = dbinom(S, p[2], size = N, log = T) + log(phi[2]) 
    
    lk_old = lk
    
    lk = apply(cbind(lk1, lk2), 1,  logSumExp, simplify = TRUE)
    
    tol = abs(sum(lk_old) - sum(lk))
    
    idx = idx + 1
  }
  
  BIC_b1 <- log(length(N)) + 2 * b1@min 
  BIC_b2 <- 3 * log(length(N)) - 2 * sum(lk)
  
  if(return_pars) {
    return(list(p = p, phi = phi, Z = z))
  }
  
  if(BIC_b2 < BIC_b1){
    tribble(
      ~"cluster", ~"n_ssms", ~"proportion",
      1,   round(length(N) * phi[1]) , p[1],
      2,   round(length(N) * phi[2]), p[2]
    )
  } else {
    tribble(
      ~"cluster", ~"n_ssms", ~"proportion",
      1,   length(N), p_b1,
    )
  }

}

fit_BetaBinomial = function(df, seed = 3, return_pars = FALSE)
{
  set.seed(seed)
  library(VGAM)
  input <- df %>% dplyr::select(NV, DP) %>% as.data.frame()
  
  LL <- function(a, b) {
         R = dbetabinom.ab(x = input$NV, shape1 = a, shape2 = b, size = input$DP, log = T)
         
           -sum(R)
     }
  
  bb <- stats4::mle(LL, start = c( 50, 50), lower = c(0,0), upper = c(Inf, Inf))
  
  if(return_pars) {
    return(list(p = bb@coef[1], rho = bb@coef[2]))
  }
  
  tribble(
    ~"cluster", ~"n_ssms", ~"proportion",
    1,   length(input$NV) , bb@coef[1]  / sum(bb@coef[1] + bb@coef[2]) ,
  ) 
}

sampler = function(
  overdispersion = 0.01,
  N = 1000,
  coverage = rpois(N, 45),
  n_binomial_methods = 7,
  n_betabinomial_methods = 2,
  seed = 3
)
{
  set.seed(seed)
  rndm_folder = sample(LETTERS, 16) %>% paste(collapse = '')
  dir.create(rndm_folder)
  
  # Simulate data
  dataset = simulate_data(overdispersion, N, coverage, seed)
  
  # Fit Binomial
  binomial_methods = lapply(
    1:n_binomial_methods,
    function(x) fit_Binomial(dataset,seed = x)
  )
  
  # Fit BetaBinomial
  betabinomial_methods = lapply(
    1:n_betabinomial_methods,
    function(x) fit_BetaBinomial(dataset,seed = x)
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
