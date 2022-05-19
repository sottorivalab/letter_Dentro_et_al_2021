require(dplyr)
require(patchwork)

set.seed(3)

source("sampler.R")


fit_Binomial_fixed_p = function(x,p1, p2, seed = 3)
{
  
  library(matrixStats)
  library(VGAM)
  
  set.seed(seed)
  
  phi1 = runif(1)
  
  p = c(p1,p2)
  
  phi = c(phi1, 1 - phi1)
  
  E_tol = 1e-3
  
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
    

    lk1 = dbinom(S, p[1], size = N, log = T) + log(phi[1])  
    lk2 = dbinom(S, p[2], size = N, log = T) + log(phi[2]) 
    
    lk_old = lk
    
    lk = apply(cbind(lk1, lk2), 1,  logSumExp, simplify = TRUE)
    
    tol = abs(sum(lk_old) - sum(lk))
    
    idx = idx + 1
  }
  
  return(c(-sum(lk),  phi[1]))
  
}

overdispersion = 0.015
#overdispersion = rgamma(1, shape = shape, rate = rate)
N = 5000
coverage = rpois(1,lambda = 45) %>% round()

smpl <- simulate_data(
  overdispersion = overdispersion,
  N = N,
  coverage = rpois(N, coverage),
  seed = i
)


grid = expand.grid(p1 = seq(0.1,0.9,by = 0.01),p2 = seq(0.1,0.9,by = 0.01))


res = easypar::run(FUN = function(x) {
  fit_Binomial_fixed_p(smpl, grid[x,1], grid[x,2])
}
    
  ,PARAMS = lapply(1:nrow(grid), list),parallel = FALSE,cache = FALSE,
  export = c("grid", "smpl", "fit_Binomial_fixed_p"), filter_errors = FALSE, 
  cores.ratio = 0.6
)

lks <- res %>% do.call(cbind, .)



grid$lk <-  lks[1,]
grid$phi <- lks[2,]

grid %>%  saveRDS("grid_ps.rds")

grid$phi <-  pmax(grid$phi , 1 - grid$phi)

p1 <- ggplot(grid, aes(p1, p2, fill= lk / nrow(grid))) + 
  geom_tile() + scale_fill_distiller("NLL/N", palette = "Spectral", direction = -1) +
  theme_minimal() + ggtitle("NLL for fixed values of prob") + geom_contour( aes(z = ifelse(lk < 15000, lk, NA) ), size = 0.5,alpha = 0.45, bins = 10, color = "black", lty = "dashed")

p2 <- ggplot(grid, aes(p1, p2, fill= phi)) + 
  geom_tile() + scale_fill_distiller("maxprob", palette = "Spectral", direction = -1) +
  theme_minimal() + ggtitle("Max phi for fixed values of prob") 

(p1 | p2) %>%  saveRDS(file = "supp2_ab.rds")

(p1 | p2) %>%  ggsave(plot = .,filename = "supp2_ab.png", device = "png", scale = "px", width = 1600, height = 800) 

