library(tidyverse)
library(patchwork)

### Figure S1 panel A: example of fit

data <- readRDS("../overdispersion/germline_analysis.rds")

germline_SNPs <- readRDS("../overdispersion/letter_snps.rds")

germline_SNPs <- germline_SNPs[names(data)[1]]

df2 <-  data[[1]]

x <-  seq(0,100,by = 1)

COV <- 45

bin1 <-  dbinom(x, size = COV, prob = df2$BMIX$B.params[1]) %>%  data.frame(density = .) %>%
  mutate(weight = 1, type = "bin1", x = x)
bbin <- VGAM::dbetabinom(x, size = COV, prob = df2$BBMIX$BB.params[1,1], rho = df2$BBMIX$BB.params[2,1]) %>%  data.frame(density = .)%>%
  mutate(weight = 1, type = "bbin", x = x)

density <- rbind(bin1, bbin)

hist <-  germline_SNPs[[1]] %>% separate(X4, into = c("NR", "NV", "NO")) %>%
  select(NR, NV) %>%  mutate(NR = as.numeric(NR), NV = as.numeric(NV))

#hist <- hist %>%  mutate(NV = (NV / (NR + NV)) * 45 %>% round())

p2_1 <-
  ggplot(data = density , aes(x = x, y = density , color = type)) + geom_line(size = 1.5) + 
  theme_bw() + xlim(0, 45) + xlab("") +
  scale_color_brewer("model", palette = "Set1") + ylab("density") + ggtitle(paste0("bbin vs bin fits (rho= ",round(df2$BBMIX$BB.params[2],3),") for sample 02c97e2b−914e−4afc−bf50−78f0cfbfa67b"))

p2_2 <- 
  ggplot(data = hist %>% filter((NV / (NR + NV)) <  0.88, (NV / (NR + NV)) >  0.12), aes(x = NV)) + 
  geom_histogram( alpha = .85, bins = 46) + theme_bw() + xlim(0, 45) + xlab("") +
  scale_color_brewer("model", palette = "Set1") + ylab("density") + ggtitle(paste0("Number of variant reads for sample 02c97e2b−914e−4afc−bf50−78f0cfbfa67b"))


BB_S <- data.frame(model = "bbin", NLL2 = df2$BBMIX$NLL * 2, AIC =2 * 2 + df2$BBMIX$NLL * 2,  BIC = df2$BBMIX$BIC)

B1_S <- data.frame(model = "bin1", NLL2 = df2$BMIX$NLL * 2, AIC =2 * 1  +df2$BMIX$NLL * 2,  BIC = df2$BMIX$BIC)

# B2_S <- data.frame(model = "bin2", NLL2 = BMIX_fit_B2$NLL * 2, ICL = BMIX_fit_B2$ICL,  BIC = BMIX_fit_B2$BIC)

score_df <-  rbind(BB_S, B1_S)

score_df_l <-  reshape2::melt(score_df, id.vars = c("model"))

p2_3 <- ggplot(score_df_l , aes(x = variable %>%  as.factor(), y = value %>%  as.numeric(), fill = model %>%  as.factor())) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw() + xlab("") +
  scale_fill_brewer("model", palette = "Set1") + ylab("score") + ggtitle("Fits score comparison for sample 02c97e2b−914e−4afc−bf50−78f0cfbfa67b")




### Figure 2 panel B: Rho distribution (histogram)

rho <-  readr::read_tsv("../overdispersion/rho_selected_cases.tsv")

colnames(rho) <-  c("sample", "prob", "rho")

p2_4 <- ggplot(rho, aes(x = rho)) + geom_histogram(bins = 60) + theme_bw() + ggtitle("Rho distribution in 41 selected samples")

### Figure 2 panel C: Scatter of NLL and BIC values of B vs BB

all_scores <-  lapply(data, function(x) x$scores) %>%  do.call(rbind,.)

# p2_4 <- ggplot(all_scores %>%  select(model, sample, NLL2) %>% spread(data = ., key = model, value = c(NLL2))  ,
#                aes(x = bbin / 2, y = bin1 / 2) ) + geom_point() + theme_bw() + xlim(1.25e5, 1.9e5)+ ylim(1.25e5, 1.9e5) + 
#   geom_abline(color = "darkred", linetype = 2) + xlab("beta-binomial") + ylab("binomial") + ggtitle("NLL values for binomial and beta-binomial fits on het-SNPs")
# 
# 

p2_5 <- ggplot(all_scores %>%  select(model, sample, AIC) %>% spread(data = ., key = model, value = c(AIC))  ,
               aes(x = bbin , y = bin1 ) ) + geom_point() + theme_bw() +  xlim(2.5e5, 3.7e5)+ ylim(2.5e5, 3.7e5) +
  geom_abline(color = "darkred", linetype = 2) + xlab("beta-binomial") + ylab("binomial") + ggtitle("AIC values for binomial and beta-binomial fits on het-SNPs")


p2_6 <- ggplot(all_scores %>%  select(model, sample, BIC) %>% spread(data = ., key = model, value = c(BIC))  ,
               aes(x = bbin, y = bin1 ) ) + geom_point() + theme_bw() + xlim(2.5e5, 3.7e5)+ ylim(2.5e5, 3.7e5) + 
  geom_abline(color = "darkred", linetype = 2) + xlab("beta-binomial") + ylab("binomial") + ggtitle("BIC values for binomial and beta-binomial fits on het-SNPs")


figS1 <- (p2_1 | p2_2) / (p2_3 | p2_4) / (p2_5 | p2_6)  + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold'))

figS1 %>% ggsave(filename = "figure_S1.png", device = "png", units = "px", width = 4800, height = 6000)

