library(gtools)
library(DESeq2)
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(dirmult)

setwd('~/Desktop/research/microbial/')

permutation_t_test <- function(x, group, num_perm=999) {
  # group: a boolean
  n1 <- sum(group)
  n2 <- sum(!group)
  t_stat <- abs((mean(x[group]) - mean(x[!group])) / 
    sqrt(var(x[group]) / n1 + var(x[!group]) / n2))
  perm_t_stat <- rep(NA, num_perm)
  for (i in 1:num_perm) {
    group_perm <- sample(group)
    perm_t_stat[i] <- abs((mean(x[group_perm]) - 
                      mean(x[!group_perm])) / 
      sqrt(var(x[group_perm]) / n1 + 
             var(x[!group_perm]) / n2))
  }
  (1 + sum(perm_t_stat >= t_stat)) / (1 + num_perm)
}

sim_richness_null <- function(num_sim=100,
                         lib1, lib2,
                         alpha, prob_vector,
                         n1=30, n2=30,
                         num_perm=99,
                         mc.cores=12) {
  perm_p_values <- rbindlist(mclapply(1:num_sim, function(i) {
    L1 <- sample(lib1, n1)
    L2 <- sample(lib2, n2)
    x1 <- sapply(L1, function(L) {
      rmultinom(n=1, size=L, prob=rdirichlet(1, alpha * prob_vector))
    })
    x2 <- sapply(L2, function(L) {
      rmultinom(n=1, size=L, prob=rdirichlet(1, alpha * prob_vector))
    })
    richness1 <- colSums(x1 > 0)
    richness2 <- colSums(x2 > 0)
    group <- c(rep(T, n1), rep(F, n2))
    original <- permutation_t_test(c(richness1, richness2), 
                       group,
                       num_perm=num_perm)
    L_star <- min(c(L1, L2))
    x1_rrf <- t(rrarefy(t(x1), L_star))
    x2_rrf <- t(rrarefy(t(x2), L_star))
    richness1_rrf <- colSums(x1_rrf > 0)
    richness2_rrf <- colSums(x2_rrf > 0)
    rarefied <- permutation_t_test(c(richness1_rrf, richness2_rrf), 
                                   group=group,
                                   num_perm=num_perm)
    data.frame(p_val=c(original, rarefied),
               data=rep(c('original', 'rarefied')))
  }, mc.cores=12))
  perm_p_values
}

# Rob Knight's data ----
library(MicrobeDS)
data('qa10394')
sample_data_10394 <- read.table('sample_data_10394.txt', 
                                header=TRUE, sep='\t')
sample_data_10394$name <- paste0('10394.', sample_data_10394$description)
sample_df <- data.frame(name=colnames(otu_table(qa10394)), 
                        stringsAsFactors=F)
sample_df <- sample_df %>% left_join(sample_data_10394, by='name')
well_defined <- !grepl('BLANK', sample_df$host_subject_id) & 
  !grepl('FTA', sample_df$host_subject_id) &
  !grepl('mistake', sample_df$host_subject_id) &
  colSums(otu_table(qa10394)) >= 10000

sample_df_subset <- sample_df[well_defined, ]
otu_table_subset <- otu_table(qa10394)[, well_defined]
library_sizes <- colSums(otu_table_subset)
freezethaw <- library_sizes[
  sample_df_subset$sample_storage_temp_treatment == 'freezethaw'
  ]
fourC <- library_sizes[
  sample_df_subset$sample_storage_temp_treatment == '4C'
  ]

total_vector <- rowSums(otu_table_subset)
prob_vector <- total_vector / sum(total_vector)
theta <- weirMoM(t(otu_table_subset))
alpha <- (1 - theta) / theta

set.seed(2000)
num_sim <- 400
ptm <- proc.time()
perm_p_values_all <- sim_richness_null(
  num_sim=num_sim,
  lib1=freezethaw, lib2=fourC,
  alpha=alpha, prob_vector=prob_vector,
  n1=30, n2=30,
  num_perm=199,
  mc.cores=12)
proc.time() - ptm

richness_results <- perm_p_values_all %>% group_by(data) %>%
  summarize(typeIerror=mean(p_val <= 0.05),
            U=mean(p_val <= 0.05) + 
              1.96 * sd(p_val <= 0.05) / sqrt(num_sim),
            L=mean(p_val <= 0.05) - 
              1.96 * sd(p_val <= 0.05) / sqrt(num_sim))

pdf('richness.pdf', width=8, height=8)
ggplot(richness_results, aes(x=data, y=typeIerror)) +
  geom_point() +
  geom_errorbar(aes(ymax = U, ymin = L), width = 0.5) +
  ylab('95% CI of Type I error rate') +
  geom_hline(yintercept=0.05, lty=2) +
  ggtitle('Type I error rates for species richness\n(Permutation t-test)') +
  theme(text=element_text(size=18)) +
  ylim(c(0, 0.20))
dev.off()
