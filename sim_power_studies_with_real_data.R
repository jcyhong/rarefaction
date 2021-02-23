# Simulation (realistic)
# Correctly specified

library(gtools)
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(MASS)
library(MicrobeDS)
data('HMPv35')
setwd('~/Desktop/research/microbial/')

REI <- function(otu_matrix, rarefied_depth, group, taxa_are_rows=FALSE) {
  # Compute rarefaction efficiency indices based on a given OTU matrix
  # and the rarefied depth.
  #
  # Input(s):
  #  - otu_matrix: a matrix of counts (taxa are columns)
  #  - rarefied_depth: a scalar indicating the rarefied depth
  #  - group: a vector indicating the group memberships
  #  - taxa_are_rows: a boolean, TRUE if each row represents a taxon
  #
  # Returns:
  #  - a vector of REIs
  
  if (taxa_are_rows) {
    otu_matrix <- t(otu_matrix)
  }
  library_size <- rowSums(otu_matrix)
  
  # Drop all observations with library sizes less than rarefied depth.
  otu_matrix <- otu_matrix[library_size >= rarefied_depth, ]
  group <- group[library_size >= rarefied_depth]
  library_size <- library_size[library_size >= rarefied_depth]
  group_lvl <- unique(group)
  n1 <- sum(group == group_lvl[1])
  n2 <- sum(group == group_lvl[2])
  
  var_prop1 <- apply(otu_matrix[group == group_lvl[1], ] / 
                       library_size[group == group_lvl[1]], 
                     2, function(col) {
                       var(col)
                     })
  var_prop_rrf1 <- apply(otu_matrix[group == group_lvl[1], ] / 
                           library_size[group == group_lvl[1]], 
                         2, function(col) {
                           var(col) + 
                             1 / rarefied_depth * mean(col * (1 - col) * (library_size[group == group_lvl[1]] - rarefied_depth) / 
                                                         (library_size[group == group_lvl[1]] - 1))
                         })
  
  var_prop2 <- apply(otu_matrix[group == group_lvl[2], ] / 
                       library_size[group == group_lvl[2]], 
                     2, function(col) {
                       var(col)
                     })
  var_prop_rrf2 <- apply(otu_matrix[group == group_lvl[2], ]  / 
                           library_size[group == group_lvl[2]], 
                         2, function(col) {
                           var(col) + 
                             1 / rarefied_depth * mean(col * (1 - col) * (library_size[group == group_lvl[2]] - rarefied_depth) / 
                                                         (library_size[group == group_lvl[2]] - 1))
                         })
  
  (var_prop1 / n1 + var_prop2 / n2) / 
    (var_prop_rrf1 / n1 + var_prop_rrf2 / n2)
}

gen_sim2 <- function(FC=2,
                    num_obs1=100, 
                    num_obs2=100,
                    libs1, libs2, 
                    prob1,
                    alpha=40) {
  num_spp <- length(prob1)
  # alpha: concentration parameter
  L1 <- sample(libs1, size=num_obs1, replace=TRUE)
  L2 <- sample(libs2, size=num_obs2, replace=TRUE)
  prob2 <- c(prob1[1] * FC, (1 - prob1[1] * FC) * prob1[2:num_spp] / sum(prob1[2:num_spp]))
  
  x1 <- sapply(L1, function(L) {
    rmultinom(n=1, size=L, prob=rdirichlet(1, alpha * prob1))
  })
  x2 <- sapply(L2, function(L) {
    rmultinom(n=1, size=L, prob=rdirichlet(1, alpha * prob2))
  })
  t_no_rrf <- t.test(x1[1, ] / colSums(x1), x2[1, ] / colSums(x2))$p.value < 0.05
  
  df <- data.frame(y=c(x1[1, ], x2[1, ]), 
                   x=c(rep(0, num_obs1), rep(1, num_obs2)),
                   L=c(L1, L2))
  nb_no_rrf <- summary(glm.nb(y ~ offset(log(L)) + x, data=df))$coefficients[2, 4] < 0.05
  
  L_star <- min(c(L1, L2))
  x1_rrf <- t(rrarefy(t(x1), L_star))
  x2_rrf <- t(rrarefy(t(x2), L_star))
  t_rrf <- t.test(x1_rrf[1, ] / colSums(x1_rrf), x2_rrf[1, ] / colSums(x2_rrf))$p.value < 0.05
  
  df_rrf <- data.frame(y=c(x1_rrf[1, ], x2_rrf[1, ]),
                       x=c(rep(0, ncol(x1_rrf)), rep(1, ncol(x2_rrf))),
                       L=L_star)
  nb_rrf <- summary(glm.nb(y ~ offset(log(L)) + x, data=df_rrf))$coefficients[2, 4] < 0.05
  
  otu_matrix <- cbind(x1, x2)
  group <- c(rep(1, ncol(x1)), rep(2, ncol(x2)))
  
  data.frame(test=c('t test', 't test', 'NB test', 'NB test'),
             data=c('original', 'rarefied', 'original', 'rarefied'),
             reject=c(t_no_rrf, t_rrf, nb_no_rrf, nb_rrf),
             REI1=REI(otu_matrix, rarefied_depth=L_star,
                      group=group, taxa_are_rows=TRUE)[1])
}

otu_table_common <- otu_table(HMPv35)[rowSums(otu_table(HMPv35) > 0) > 600, ]
otu_table_common <- otu_table_common[, colSums(otu_table_common) > 1000]
base_libs <- colSums(otu_table_common)
base_prob <- rowSums(otu_table_common)
base_prob <- base_prob / sum(base_prob)

num_sim <- 500
# num_sim <- 12
num_obs_list <- c(20, 100)
factor1 <- 1
factor2 <- c(2, 10)

small_alpha_small_sample <- expand.grid(20, 
                           factor1,
                           factor2,
                           100,
                           seq(1, 9, 1))
medium_alpha_small_sample <- expand.grid(20, 
                            factor1,
                            factor2,
                            1000,
                            seq(1, 4, 0.3))
large_alpha_small_sample <- expand.grid(20, 
                           factor1,
                           factor2,
                           10000,
                           seq(1, 3, 0.25))
small_alpha_large_sample <- expand.grid(100, 
                                        factor1,
                                        factor2,
                                        100,
                                        seq(1, 5, 0.5))
medium_alpha_large_sample <- expand.grid(100, 
                                         factor1,
                                         factor2,
                                         1000,
                                         seq(1, 2.4, 0.2))
large_alpha_large_sample <- expand.grid(100, 
                                        factor1,
                                        factor2,
                                        10000,
                                        seq(1, 1.8, 0.08))
config_list <- rbind(small_alpha_small_sample, 
                     medium_alpha_small_sample, 
                     large_alpha_small_sample,
                     small_alpha_large_sample, 
                     medium_alpha_large_sample, 
                     large_alpha_large_sample)

colnames(config_list) <- c('num_obs', 'factor1', 'factor2', 'alpha', 'FC')
set.seed(3000) # fix seed for reproducibility
ptm <- proc.time()
results <- rbindlist(mclapply(1:nrow(config_list), function(i) {
  num_obs <- config_list[i, ]$num_obs
  alpha <- config_list[i, ]$alpha
  FC <- config_list[i, ]$FC
  factor1 <- config_list[i, ]$factor1
  factor2 <- config_list[i, ]$factor2
  results <- rbindlist(
    replicate(num_sim, 
              gen_sim2(num_obs1=num_obs,
                      num_obs2=num_obs,
                      libs1=factor1 * base_libs,
                      libs2=factor2 * base_libs,
                      prob1=base_prob,
                      FC=FC, 
                      alpha=alpha), 
              simplify=FALSE)
  )
  results %>% group_by(test, data) %>% summarize(power=mean(reject), 
                                                 REI=mean(REI1)) %>%
    mutate(num_obs=num_obs, factor1=factor1, factor2=factor2, alpha=alpha, FC=FC)
}, mc.cores=12))
ptm <- proc.time() - ptm

write.csv(results, file='sim_results/realistic/sim_realistic.csv')

plots <- list()

for (num_obs_value in num_obs_list) {
  for (factor2_value in factor2) {
    df_sub <- results %>% filter(num_obs==num_obs_value & factor2==factor2_value)
    rei_text <- df_sub %>% group_by(alpha, data) %>% 
      summarize(REI=mean(REI, na.rm=TRUE))
    rei_text$REI <- round(rei_text$REI, 2)
    colnames(rei_text)[3] <- 'Average REI'
    df_sub <- df_sub %>% left_join(rei_text, by=c('alpha', 'data'))
    plots <- c(plots, list(ggplot(df_sub, 
                 aes(x=FC, y=power, color=data, lty=data)) +
            geom_line() +
            facet_grid(test~`Average REI`, labeller = label_both,
                       scales='free') +
            #ggtitle(sprintf('n=%s; %sx library size difference; p=%s',
            #                num_obs_value, factor2_value, round(base_prob[1], 4))) +
            ggtitle(sprintf('n=%s; %sx library size difference',
                            num_obs_value, factor2_value)) +
            xlab('fold change') +
            theme(text = element_text(size=15))))
  }
}

pdf('sim_results/realistic/all.pdf', width=14, height=14)
gridExtra::grid.arrange(plots[[1]], plots[[2]],
                        plots[[3]], plots[[4]], nrow=2)
dev.off()


