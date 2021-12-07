generate_report_overdispersion <- function(snps, snvs, cna, name, meta){

  library(ggplot2)
  library(tidyverse)
  library(ggrepel)
  
  set.seed(3)


  print(paste0("Processing sample ", name))

  colnames(snps) <-  c("chr", "ref", "alt", "AD", "DP")

  # Separate the AD filed to obtain DR and DV
  snps <-
    snps %>% tidyr::separate(col = "AD",
                             into = c("DR", "DV", "DO"),
                             sep = ",") %>%
    separate(col = "alt",
             into = c("alt", "other"),
             sep = ",")


  # Filter multiallelic SNPs
  snps_filt <- snps %>%  filter(DO == 0)

  # Filter hetero SNPs
  L_FILT <- 0.12
  H_FILT <-  0.88


  snps_het <-
    snps_filt %>% mutate(VAF = as.numeric(DV) / as.numeric(DP)) %>%  filter(VAF > L_FILT, VAF < H_FILT)


  T_test <-  Tarone.test(snps_het$DP %>% as.numeric()
                         ,snps_het$DV  %>% as.numeric())

  bmix_data <- data.frame(snps_het$DV %>% as.numeric(), snps_het$DP  %>% as.numeric())

  colnames(bmix_data) <- c("DV", "DP")

  # For efficiency reasons we subsample the data
  if(nrow(bmix_data) > 50000) {

    bmix_data_s <- bmix_data %>% sample_n(50000)

  } else {
    bmix_data_s <-  bmix_data
  }

  # Standard BMix fit, with EM
  BMIX_fit_BB <- BMix::bmixfit(K.Binomials = 0, K.BetaBinomials = 1, 
                               samples = 1, score = "NLL",data = bmix_data_s)
  
  BMIX_fit_B1 <- BMix::bmixfit(data = bmix_data_s, K.Binomials = 1, K.BetaBinomials = 0,
                               samples = 1, score = "NLL")

  
  x <-  seq(0,100,by = 1)

  COV <- median(snps_filt$DP) %>%  round

  bin1 <-  dbinom(x, size = COV, prob = BMIX_fit_B1$B.params[1]) %>%  data.frame(density = .) %>%
    mutate(weight = 1, type = "bin1", x = x)


  bbin <- VGAM::dbetabinom(x, size = COV, prob = BMIX_fit_BB$BB.params[1,1], rho = BMIX_fit_BB$BB.params[2,1]) %>%  data.frame(density = .)%>%
    mutate(weight = 1, type = "bbin", x = x)

  density <- rbind(bin1, bbin)


  dens_plot <-
    ggplot(density , aes(x = x, y = density * weight, color = type)) + geom_line() + theme_bw() + xlim(0, 100) + xlab("") +
    scale_color_brewer("model", palette = "Set1") + ylab("density") + ggtitle(paste0("bbin vs bin fits (rho= ",round(BMIX_fit_BB$BB.params[2],3),")"))


  BB_S <- data.frame(model = "bbin", NLL2 = BMIX_fit_BB$NLL * 2, AIC =2 * 2 + BMIX_fit_BB$NLL * 2,  BIC = BMIX_fit_BB$BIC)

  B1_S <- data.frame(model = "bin1", NLL2 = BMIX_fit_B1$NLL * 2, AIC =2 * 1  +BMIX_fit_B1$NLL * 2,  BIC = BMIX_fit_B1$BIC)
  
  score_df <-  rbind(BB_S, B1_S)

  score_df_l <-  reshape2::melt(score_df, id.vars = c("model"))

  score_plot <- ggplot(score_df_l, aes(x = variable %>%  as.factor(), y = value %>%  as.numeric(), fill = model %>%  as.factor())) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw() + xlab("") +
    scale_fill_brewer("model", palette = "Set1") + ylab("score") + ggtitle("Fits score comparison")

  vaf_plot <-
    ggplot(data = snps_het, aes(x = as.numeric(DV) / as.numeric(DP))) +
    geom_histogram(binwidth =  0.01, alpha = 0.8) + theme_bw() + xlab("VAF") + ggtitle("Germline SNPs VAF plot")


  somatic_plot <-
  ggplot(data = snvs, aes(x = ccf, fill = cluster %>%  paste)) +
    geom_histogram(bins = 100, position = "stack", alpha = 0.6) +
    geom_text_repel(data = snvs %>%  filter(!is.na(isdriver)),
                    aes(label = Hugo_Symbol, y = Inf, x = ccf), colour = "black",
                    label.size = 0.01, alpha = 0.85) +
    geom_vline(data = snvs %>%  filter(!is.na(isdriver)),
               aes(xintercept = ccf, colour = cluster), lty = 2, show.legend = FALSE) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Dark2"), breaks = paste0("cluster", 1:6)) +
    scale_fill_manual("cluster",values = RColorBrewer::brewer.pal(6, "Dark2"), breaks = paste0("cluster", 1:6)) +
    theme_bw() +
    xlim(c(0,2)) +
    xlab("CCF") +
    ylab("Counts") + ggtitle("Somatic SNVs PCAWG CCF and clustering")

  # Tarone statistics, beware of how this statistics needs to be interpreted together with the 
  # actal effect size (rho value)
  tarone_plot <- Tarone_Boostrap(bmix_data_s %>% mutate(NV = DV, VAF = NV / DP))

  cna_p <-  cna %>% filter(!is.na(Major)) %>% mutate(karyotype = paste0(Major, ":", minor), length = to - from) %>%
    mutate(is_diploid = ifelse(karyotype == "1:1", "diploid", "aneuploid"), length = to - from ) %>%
    dplyr::select(length, is_diploid) %>% group_by(is_diploid) %>%  summarize(t_length = sum(length)) %>%
    mutate(t_length =(t_length / sum(t_length)) %>% round(2))


  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )

  cna_p$position <-  c(0,cumsum(cna_p$t_length)[1]) + cna_p$t_length/2

  cna_p$is_diploid <-  factor(cna_p$is_diploid, levels = c("diploid", "aneuploid"))

  cna_plot <-  ggplot(cna_p, aes(x="", y=t_length, fill=is_diploid) ) +
    geom_bar(width = 1, stat = "identity") + coord_polar("y", direction = 1) +
     blank_theme +
    theme(axis.text.x=element_blank()) +
    geom_text(aes(y = position,
                  label = scales::percent(t_length)), size=4) + ggtitle("Percentage of diploid genome") +
    scale_fill_manual(values = c("orange", "gainsboro"))

  library(cowplot)

  pfin <- cowplot::plot_grid(dens_plot, score_plot, vaf_plot, somatic_plot,
                             tarone_plot, cna_plot,
                             axis = "rlbt", align = "hv", ncol = 2, labels = "auto")
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(name, "  ttype=", snvs$ctype %>%  unique, "  purity=", meta$purity, "  ploidy=", meta$ploidy),
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  psave = cowplot::plot_grid(
    title, pfin,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.09, 1)
  )

  psave %>%  ggsave(paste0("./plots/",name, "_overdispersion_plot.pdf"), plot = ., width = 10, height = 12)

  return(list(scores = score_df %>%  mutate(sample = name),plts = psave, BBMIX = BMIX_fit_BB, BMIX = BMIX_fit_B1))

}



Tarone.test <- function(N, M) {
  #Check validity of inputs
  if (!(all(N == as.integer(N)))) {
    stop("Error: Number of trials should be integers")

  }
  if (min(N) < 1) {
    stop("Error: Number of trials should be positive")

  }
  if (!(all(M == as.integer(M)))) {
    stop("Error: Count values should be integers")

  }
  if (min(M) < 0) {
    stop("Error: Count values cannot be negative")

  }
  if (any(M > N)) {
    stop("Error: Observed count value exceeds number of trials")

  }
  #Set description of test and data
  method      <- "Tarone's Z test"

  data.name   <- paste0(deparse(substitute(M)),
                        " successes from ",
                        deparse(substitute(N)),
                        " trials")

  #Set null and alternative hypotheses
  null.value  <- 0

  attr(null.value, "names") <- "dispersion parameter"

  alternative <- "greater"

  #Calculate test statistics
  estimate    <- sum(M) / sum(N)

  attr(estimate, "names") <- "proportion parameter"

  S           <- ifelse(estimate == 1, sum(N),
                        sum((M - N * estimate) ^ 2 / (estimate * (1 - estimate))))

  statistic   <- (S - sum(N)) / sqrt(2 * sum(N * (N - 1)))

  attr(statistic, "names") <- "z"

  #Calculate p-value
  p.value     <- 2 * pnorm(-abs(statistic), 0, 1)

  attr(p.value, "names") <- NULL

  #Create htest object
  TEST        <- list(
    method = method,
    data.name = data.name,
    null.value = null.value,
    alternative = alternative,
    estimate = estimate,
    statistic = statistic,
    p.value = p.value
  )

  class(TEST) <- "htest"

  TEST

}


# TARONE bootstrap test

Tarone_Boostrap <-
  function(x,
           num_samples = 1000)  {


    N = x$DP %>% as.numeric()
    M = x$NV %>% as.numeric()
    z_scores = c()

    for (i in 1:num_samples) {
      index = sample(1:length(N), size = length(N), replace = TRUE)
      n = N[index]
      m = M[index]

      p = sum(m) / sum(n)
      S = ((m - n * p) ^ 2)/(p * (1 - p))
      S = sum(S)
      Z = (S - sum(n)) / sqrt(2 * sum(n * (n - 1)))
      z_scores = c(z_scores, Z)
    }

    log_p_value = dnorm(z_scores, log = TRUE) %>% sum()

    null_hyp = rnorm(num_samples)
    y = data.frame(value = z_scores, data = "Z-score")
    w = data.frame(value = null_hyp, data = "Binomial (sim.)")
    data = rbind(y, w)

    # Create histogram plot
    p_label = format(log_p_value %>% exp(), digits = 3)
    p_label = ifelse(p_label == "0", "p<1e-12", paste("p = ", p_label))

    plot_tarone <-
      ggplot(data, aes(x = value, fill = data)) +
      geom_histogram(
        aes(y = ..density..),
        alpha = 0.6,
        position = "dodge",
        bins = 50
      ) +
      geom_vline(xintercept = 0,
                 linetype = "dashed",
                 size = 0.7) +
      stat_function(
        fun = dnorm,
        colour = "cornflowerblue",
        size = 1,
        args = list(mean = 0, sd = 1)
      ) +
      labs(
        title = paste0("Tarone's Z score; ", p_label),
        x = "log10(Z) (overdispersion)",
        y = "Density"
      ) +
      scale_fill_manual(values = c('Binomial (sim.)' = "#377EB8", "Z-score" = "#E41A1C")) +
      theme_linedraw(base_size = 9) +
      guides(fill = guide_legend("")) +
      theme(legend.position = "bottom")

    plot_tarone

  }

# We assume to have an rds with the following information for each sample: CNAs, SNVs, SNPs, and metadata (purity)
# The expected input type is a list
# CNAs can be downloaded from https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_cnv/consensus.20170119.somatic.cna.annotated.tar.gz
# SNVs can be obtained from the mutation assignment files from the official PCWAG release
# Meta can be downloaded from consensus.20170217.purity.ploidy.txt.gz
# SNPs need to be extracted from closed access BAM files, they need to have the following columns: ["chr", "ref", "alt", "AD", "DP"].

<<<<<<< HEAD
# cnas <-  readRDS("letter_cnas.rds")
# 
# snvs <-  readRDS("letter_snvs.rds")
# 
# snps <-  readRDS("letter_snps.rds")
# 
# meta <- readRDS("letter_metadata.rds")



#res <-  mapply(snps, snvs, cnas,names(cnas), meta, 
#                           FUN = function(x,y,z,n,m) 
#                             try(generate_report_overdispersion(x,y,z,n, m)),
#                           SIMPLIFY = F)

#res %>% saveRDS("germline_analysis.rds")


