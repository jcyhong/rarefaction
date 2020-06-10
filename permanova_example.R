# PERMANOVA example

setwd('~/Desktop/research/microbial/')

library(ggplot2)
library(vegan)
library(dplyr)
library(parallel)
library(data.table)
library(MCMCpack)
library(gridExtra)
library(grid)
library(dirmult)

lambda1 <- 200
lambda2 <- 800
num_obs1 <- 50
num_obs2 <- 50
num_sim <- 1000
num_spp <- 100

normalize_by_library_size <- function(mat) {
  # Each row is a site.
  t(apply(mat, 1, function(row) {row/sum(row)}))
}

p_val <- mclapply(1:num_sim, function(i) {
  L1 <- rpois(n=num_obs1, lambda=lambda1)
  L2 <- rpois(n=num_obs2, lambda=lambda2)
  prob <- rep(1, num_spp)
  x1 <- sapply(L1, function(L) {
    rmultinom(n=1, size=L, prob=rdirichlet(1, prob))
  })
  x2 <- sapply(L2, function(L) {
    rmultinom(n=1, size=L, prob=rdirichlet(1, prob))
  })
  df <- data.frame(normalize_by_library_size(
    rbind(t(x1), 
          t(x2))
    ))
  ind <- c(rep(1, num_obs1), rep(2, num_obs2))

  L_star <- min(c(L1, L2))
  df_rrf <- data.frame(normalize_by_library_size(
    rbind(rrarefy(t(x1), L_star), 
          rrarefy(t(x2), L_star))
    ))
  data.frame(original=adonis(df ~ ind)$aov.tab[1, 6], 
             rarefied=adonis(df_rrf ~ ind)$aov.tab[1, 6])
}, mc.cores=12)
p_val <- rbindlist(p_val)



p_val_df <- data.frame(
  p_value=c(p_val$original, p_val$rarefied), 
  data=c(rep('original', nrow(p_val)),
            rep('rarefied', nrow(p_val))))

plot_title_size <- 20
plot_text_size <- 16

p1 <- ggplot(p_val_df %>% filter(data=='original'), aes(x=p_value)) +
  geom_histogram(bins=20, colour='red') + 
  ggtitle('PERMANOVA with original data') + 
  theme(text=element_text(size=plot_text_size),
        plot.title = element_text(size=plot_title_size)) +
  xlab('permutation p-value')

p2 <- ggplot(p_val_df %>% filter(data=='rarefied'), aes(x=p_value)) +
  geom_histogram(colour='red', breaks=seq(0, 1, 0.05)) + 
  ggtitle('PERMANOVA with rarefied data') + 
  theme(text=element_text(size=plot_text_size),
        plot.title = element_text(size=plot_title_size)) +
  xlab('permutation p-value')

grid.arrange(p1, p2)

pdf('PERMANOVA_p_values.pdf')
grid.arrange(p1, p2)
dev.off()

# Rob Knight's data ----
library(MicrobeDS)
data('qa10394')
sample_data_10394 <- read.table('sample_data_10394.txt', header=TRUE, sep='\t')
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
boxplot(library_sizes ~ sample_df_subset$sample_storage_temp_treatment)
boxplot(library_sizes ~ sample_df_subset$sample_preservation_method)
boxplot(library_sizes ~ sample_df_subset$duration_of_storage)

freezethaw <- library_sizes[
  sample_df_subset$sample_storage_temp_treatment == 'freezethaw'
    ]
fourC <- library_sizes[
  sample_df_subset$sample_storage_temp_treatment == '4C'
  ]
boxplot(freezethaw, fourC)
summary(freezethaw)
summary(fourC)

# Estimate the alpha (concentration) parameter via the method of moments.
theta <- weirMoM(t(otu_table_subset))
conc <- (1 - theta) / theta

# Simulation ----
prob_vector <- rowSums(otu_table_subset) / 
  sum(rowSums(otu_table_subset))
ptm <- proc.time()
p_val <- mclapply(1:num_sim, function(i) {
  L1 <- sample(freezethaw, num_obs1, replace=T)
  L2 <- sample(fourC, num_obs2, replace=T)
  prob <- prob_vector
  x1 <- sapply(L1, function(L) {
    rmultinom(n=1, size=L, prob=rdirichlet(1, conc * prob))
  })
  x2 <- sapply(L2, function(L) {
    rmultinom(n=1, size=L, prob=rdirichlet(1, conc * prob))
  })
  df <- data.frame(normalize_by_library_size(
    rbind(t(x1), 
          t(x2))
  ))
  ind <- c(rep(1, num_obs1), rep(2, num_obs2))
  
  L_star <- min(c(L1, L2))
  df_rrf <- data.frame(normalize_by_library_size(
    rbind(rrarefy(t(x1), L_star), 
          rrarefy(t(x2), L_star))
  ))
  data.frame(original=adonis(df ~ ind, permutations=199)$aov.tab[1, 6], 
             rarefied=adonis(df_rrf ~ ind, permutations=199)$aov.tab[1, 6])
}, mc.cores=12)
proc.time() - ptm

p_val <- rbindlist(p_val)

p_val_df <- data.frame(
  p_value=c(p_val$original, p_val$rarefied), 
  data=c(rep('original', nrow(p_val)),
         rep('rarefied', nrow(p_val))))

plot_title_size <- 20
plot_text_size <- 16

p1 <- ggplot(p_val_df %>% filter(data=='original'), aes(x=p_value)) +
  geom_histogram(bins=20, colour='red') + 
  ggtitle('PERMANOVA with original data') + 
  theme(text=element_text(size=plot_text_size),
        plot.title = element_text(size=plot_title_size)) +
  xlab('permutation p-value')

p2 <- ggplot(p_val_df %>% filter(data=='rarefied'), aes(x=p_value)) +
  geom_histogram(colour='red', breaks=seq(0, 1, 0.05)) + 
  ggtitle('PERMANOVA with rarefied data') + 
  theme(text=element_text(size=plot_text_size),
        plot.title = element_text(size=plot_title_size)) +
  xlab('permutation p-value')

grid.arrange(p1, p2)

pdf('PERMANOVA_p_values_realistic.pdf')
grid.arrange(p1, p2)
dev.off()

lib_size_df <- rbind(
  data.frame(lib_size=freezethaw, treatment='freezethaw'),
  data.frame(lib_size=fourC, treatment='4C')
)
plot_lib_size <- ggplot(lib_size_df,
                        aes(x=lib_size, 
                            col=treatment, 
                            lty=treatment)) +
  geom_density() +
  xlab('library size') +
  ggtitle('Library size distributions') + 
  theme(text=element_text(size=plot_text_size),
        plot.title = element_text(size=plot_title_size))

sim_lib_size_df <- rbind(
  data.frame(lib_size=rpois(10000, 200), distribution='Poisson(200)'),
  data.frame(lib_size=rpois(10000, 800), distribution='Poisson(800)')
)
plot_sim_lib_size <- ggplot(sim_lib_size_df, 
       aes(x=lib_size, col=distribution, lty=distribution)) +
  geom_density() +
  xlab('library size') +
  ggtitle('Library size distributions') + 
  theme(text=element_text(size=plot_text_size),
        plot.title = element_text(size=plot_title_size))

pdf('lib_size_distributions.pdf', width=9, height=7)
grid.arrange(plot_sim_lib_size, plot_lib_size, nrow=2)
dev.off()








