# Simulation studies for DESeq2 motivating examples
# 1. Various hypothetical data scenarios
# 2. Simulation based on realistic data (Rob Knight's data)

library(gtools)
library(DESeq2)
library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)

setwd('~/Desktop/research/microbial/')
results_dir <- 'deseq_results/'

nspp <- 2000

nsample <- 10
lambda1 <- 10000
lambda2 <- 50000
p1 <- rep(1/nspp, nspp)
p2 <- c(rep(1.5/nspp, nspp / 2), rep(0.5/nspp, nspp / 2))
pi0 <- 0.2
alpha1 <- 30
alpha2 <- 50
num_sim <- 100

# Asymptotic distribution ----

run_no_perm_sim <- function(num_sim,
                            nsample,
                            lambda1, lambda2,
                            p1, p2, pi0, alpha1, alpha2) {
  rbindlist(mclapply(1:num_sim, function(i) {
    L1 <- rpois(n=nsample, lambda=lambda1)
    L2 <- rpois(n=nsample, lambda=lambda2)
    x1 <- sapply(L1, function(L) {
      p <- rbinom(n=1, size=1, prob=pi0)
      rmultinom(n=1, size=L, prob=p * rdirichlet(1, alpha1 * p1) + 
                  (1 - p) * rdirichlet(1, alpha2 * p2))
    })
    x2 <- sapply(L2, function(L) {
      p <- rbinom(n=1, size=1, prob=pi0)
      rmultinom(n=1, size=L, prob=p * rdirichlet(1, alpha1 * p1) + 
                  (1 - p) * rdirichlet(1, alpha2 * p2))
    })
    
    cnts <- cbind(x1, x2)
    cond <- factor(rep(1:2, each=nsample))
    dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
    dds <- estimateSizeFactors(dds, type="poscounts")
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res1 <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
    
    cnts_rrf <- t(rrarefy(t(cnts), min(L1, L2)))
    cond <- factor(rep(1:2, each=nsample))
    dds <- DESeqDataSetFromMatrix(cnts_rrf, DataFrame(cond), ~ cond)
    dds <- estimateSizeFactors(dds, type="poscounts")
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res2 <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
    
    p_val <- c(res1$pvalue, res2$pvalue)
    species <- c(1:nspp, 1:nspp)
    data <- c(rep('original', nspp), rep('rarefied', nspp))
    df <- data.frame(p_val=p_val,
                     species=species,
                     data=data)
    df
  }, mc.cores=12))
}

compute_typeI_error <- function(result, data_type='original') {
  result$sim_num <- rep(1:num_sim, each=nspp * 2)
  (result %>% 
     filter(data==data_type) %>%
     group_by(sim_num) %>% 
     summarize(typeIerror=mean(p_val < 0.05, na.rm=TRUE)))$typeIerror
}

# 0. large sample sizes + equal lib size distr. + dirichlet latent comp.
ptm <- proc.time()
results_no_perm <- 
  run_no_perm_sim(num_sim=100, 
                  nsample=30,
                  lambda1=10000, lambda2=10000,
                  p1=rep(1/nspp, nspp), 
                  p2=rep(1/nspp, nspp), 
                  pi0=0, 
                  alpha1=30, alpha2=30)
ptm <- proc.time() - ptm
fwrite(results_no_perm,
       paste0(results_dir, 'no_perm.csv'))
results_no_perm <- read.csv(paste0(results_dir, 'no_perm.csv'))
typeIerror_no_perm <- compute_typeI_error(result=results_no_perm,
                       data_type='original')
typeIerror_no_perm_rrf <- compute_typeI_error(result=results_no_perm,
                                          data_type='rarefied')
summary(typeIerror_no_perm)
summary(typeIerror_no_perm_rrf)

# 1. small sample sizes (equal lib size distr. + dirichlet latent comp.)
results_small_sample_sizes_no_perm <- 
  run_no_perm_sim(num_sim=100, 
                  nsample=10,
                  lambda1=10000, lambda2=10000,
                  p1=rep(1/nspp, nspp), 
                  p2=rep(1/nspp, nspp), 
                  pi0=0, 
                  alpha1=30, alpha2=30)
fwrite(results_small_sample_sizes_no_perm,
       paste0(results_dir, 'small_sample_sizes_no_perm.csv'))
results_small_sample_sizes_no_perm <- read.csv(
  paste0(results_dir, 'small_sample_sizes_no_perm.csv')
)
typeIerror_small_sample_sizes_no_perm <- 
  compute_typeI_error(result=results_small_sample_sizes_no_perm,
                      data_type='original')
typeIerror_small_sample_sizes_no_perm_rrf <- 
  compute_typeI_error(result=results_small_sample_sizes_no_perm,
                      data_type='rarefied')
summary(typeIerror_small_sample_sizes_no_perm)
summary(typeIerror_small_sample_sizes_no_perm_rrf)

# 2. diff lib size distr. (large sample sizes + dirichlet latent comp.)
results_diff_lib_sizes_no_perm <- 
  run_no_perm_sim(num_sim=100, 
                  nsample=30,
                  lambda1=10000, lambda2=50000,
                  p1=rep(1/nspp, nspp), 
                  p2=rep(1/nspp, nspp), 
                  pi0=0, 
                  alpha1=30, alpha2=30)
fwrite(results_diff_lib_sizes_no_perm,
       paste0(results_dir, 'diff_lib_sizes_no_perm.csv'))
results_diff_lib_sizes_no_perm <- read.csv(
  paste0(results_dir, 'diff_lib_sizes_no_perm.csv')
)
typeIerror_diff_lib_sizes_no_perm <- 
  compute_typeI_error(result=results_diff_lib_sizes_no_perm,
                      data_type='original')
typeIerror_diff_lib_sizes_no_perm_rrf <- 
  compute_typeI_error(result=results_diff_lib_sizes_no_perm,
                      data_type='rarefied')
summary(typeIerror_diff_lib_sizes_no_perm)
summary(typeIerror_diff_lib_sizes_no_perm_rrf)

# 3. mixture latent comp. (large sample sizes + equal lib size distr.)
results_mixture_no_perm <- 
  run_no_perm_sim(num_sim=100, 
                  nsample=30,
                  lambda1=10000, lambda2=10000,
                  p1=rep(1/nspp, nspp), 
                  p2=c(rep(1.5/nspp, nspp / 2), rep(0.5/nspp, nspp / 2)), 
                  pi0=0.3, 
                  alpha1=30, alpha2=3)
fwrite(results_mixture_no_perm,
       paste0(results_dir, 'mixture_no_perm.csv'))
results_mixture_no_perm <- read.csv(
  paste0(results_dir, 'mixture_no_perm.csv')
)
typeIerror_mixture_no_perm <- 
  compute_typeI_error(result=results_mixture_no_perm,
                      data_type='original')
typeIerror_mixture_no_perm_rrf <- 
  compute_typeI_error(result=results_mixture_no_perm,
                      data_type='rarefied')
summary(typeIerror_mixture_no_perm)
summary(typeIerror_mixture_no_perm_rrf)

get_stat_DESeq2 <- function(cnts, cond) {
  dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
  dds <- estimateSizeFactors(dds, type="poscounts")
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
  res$stat
}

# Permutation distribution ----

run_perm_sim <- function(num_sim,
                         nsample,
                         lambda1, lambda2,
                         p1, p2, pi0, alpha1, alpha2,
                         num_perm) {
  rbindlist(mclapply(1:num_sim, function(i) {
    L1 <- rpois(n=nsample, lambda=lambda1)
    L2 <- rpois(n=nsample, lambda=lambda2)
    x1 <- sapply(L1, function(L) {
      p <- rbinom(n=1, size=1, prob=pi0)
      rmultinom(n=1, size=L, prob=p * rdirichlet(1, alpha1 * p1) + 
                  (1 - p) * rdirichlet(1, alpha2 * p2))
    })
    x2 <- sapply(L2, function(L) {
      p <- rbinom(n=1, size=1, prob=pi0)
      rmultinom(n=1, size=L, prob=p * rdirichlet(1, alpha1 * p1) + 
                  (1 - p) * rdirichlet(1, alpha2 * p2))
    })
    cnts <- cbind(x1, x2)
    L_star <- min(c(L1, L2))
    cnts_rrf <- t(data.frame(rbind(rrarefy(t(x1), L_star), 
                                   rrarefy(t(x2), L_star))))
    cond <- factor(rep(1:2, each=nsample))
    stat <- abs(get_stat_DESeq2(cnts, cond))
    stat_rrf <- abs(get_stat_DESeq2(cnts_rrf, cond))
    perm_stat <- rbindlist(lapply(1:num_perm, function(b) {
      cond_shuffled <- sample(cond)
      data.frame(original=abs(get_stat_DESeq2(cnts, cond_shuffled)),
                 rarefied=abs(get_stat_DESeq2(cnts_rrf, cond_shuffled)),
                 species=1:nspp)
    }))
    perm_pval <- sapply(1:nspp, function(i) {
      as.numeric(perm_stat %>% filter(species == i) %>% 
        summarize(p_val=mean(c(stat[i], original) >= stat[i], na.rm=T)))
      })
    perm_pval_rrf <- sapply(1:nspp, function(i) {
      as.numeric(perm_stat %>% filter(species == i) %>% 
        summarize(p_val=mean(c(stat_rrf[i], rarefied) >= stat_rrf[i], na.rm=T)))
      })
    
    p_val <- c(perm_pval, perm_pval_rrf)
    species <- c(1:nspp, 1:nspp)
    data <- c(rep('original', nspp), rep('rarefied', nspp))
    df <- data.frame(p_val=p_val,
                     species=species,
                     data=data)
    df
  }, mc.cores=12))
}

B <- 99
# 0. large sample sizes + equal lib size distr. + dirichlet latent comp.
ptm <- proc.time()
results_perm <- 
  run_perm_sim(num_sim=100, 
               nsample=30,
               lambda1=10000, lambda2=10000,
               p1=rep(1/nspp, nspp), 
               p2=rep(1/nspp, nspp), 
               pi0=0, 
               alpha1=30, alpha2=30,
               num_perm=B)
proc.time() - ptm
fwrite(results_perm, paste0(results_dir, 'perm.csv'))

results_perm <- read.csv(
  paste0(results_dir, 'perm.csv')
)
typeIerror_perm <- compute_typeI_error(result=results_perm,
                                          data_type='original')
typeIerror_perm_rrf <- compute_typeI_error(result=results_perm,
                                              data_type='rarefied')
summary(typeIerror_perm)
summary(typeIerror_perm_rrf)

# 1. small sample sizes (equal lib size distr. + dirichlet latent comp.)
results_small_sample_sizes_perm <- 
  run_perm_sim(num_sim=100, 
               nsample=10,
               lambda1=10000, lambda2=10000,
               p1=rep(1/nspp, nspp), 
               p2=rep(1/nspp, nspp), 
               pi0=0, 
               alpha1=30, alpha2=30,
               num_perm=B)
fwrite(results_small_sample_sizes_perm, 
       paste0(results_dir, 'small_sample_sizes_perm.csv'))
results_small_sample_sizes_perm <- read.csv(
  paste0(results_dir, 'small_sample_sizes_perm.csv')
)
typeIerror_small_sample_sizes_perm <- 
  compute_typeI_error(result=results_small_sample_sizes_perm,
                      data_type='original')
typeIerror_small_sample_sizes_perm_rrf <- 
  compute_typeI_error(result=results_small_sample_sizes_perm,
                      data_type='rarefied')
summary(typeIerror_small_sample_sizes_perm)
summary(typeIerror_small_sample_sizes_perm_rrf)

# 2. diff lib size distr. (large sample sizes + dirichlet latent comp.)
results_diff_lib_sizes_perm <- 
  run_perm_sim(num_sim=100, 
               nsample=30,
               lambda1=10000, lambda2=50000,
               p1=rep(1/nspp, nspp), 
               p2=rep(1/nspp, nspp), 
               pi0=0, 
               alpha1=30, alpha2=30,
               num_perm=B)
fwrite(results_diff_lib_sizes_perm, 
       paste0(results_dir, 'diff_lib_sizes_perm.csv'))
results_diff_lib_sizes_perm <- read.csv(
  paste0(results_dir, 'diff_lib_sizes_perm.csv')
)
typeIerror_diff_lib_sizes_perm <- 
  compute_typeI_error(result=results_diff_lib_sizes_perm,
                      data_type='original')
typeIerror_diff_lib_sizes_perm_rrf <- 
  compute_typeI_error(result=results_diff_lib_sizes_perm,
                      data_type='rarefied')
summary(typeIerror_diff_lib_sizes_perm)
summary(typeIerror_diff_lib_sizes_perm_rrf)

# 3. mixture latent comp. (large sample sizes + equal lib size distr.)
results_mixture_perm <- 
  run_perm_sim(num_sim=100, 
               nsample=30,
               lambda1=10000, lambda2=10000,
               p1=rep(1/nspp, nspp), 
               p2=c(rep(1.5/nspp, nspp / 2), rep(0.5/nspp, nspp / 2)), 
               pi0=0.3, 
               alpha1=30, alpha2=3,
               num_perm=B)
fwrite(results_mixture_perm, 
       paste0(results_dir, 'mixture_perm.csv'))
results_mixture_perm <- read.csv(
  paste0(results_dir, 'mixture_perm.csv')
)
typeIerror_mixture_perm <- 
  compute_typeI_error(result=results_mixture_perm,
                      data_type='original')
typeIerror_mixture_perm_rrf <- 
  compute_typeI_error(result=results_mixture_perm,
                      data_type='rarefied')
summary(typeIerror_mixture_perm)
summary(typeIerror_mixture_perm_rrf)


asymp_df <- data.frame(
  typeIerror=c(typeIerror_no_perm, typeIerror_no_perm_rrf,
               typeIerror_small_sample_sizes_no_perm,
               typeIerror_small_sample_sizes_no_perm_rrf,
               typeIerror_diff_lib_sizes_no_perm,
               typeIerror_diff_lib_sizes_no_perm_rrf,
               typeIerror_mixture_no_perm,
               typeIerror_mixture_no_perm_rrf),
  data=c(rep('original', length(typeIerror_no_perm)),
         rep('rarefied', length(typeIerror_no_perm_rrf)),
         rep('original', length(typeIerror_small_sample_sizes_no_perm)),
         rep('rarefied', length(typeIerror_small_sample_sizes_no_perm_rrf)),
         rep('original', length(typeIerror_diff_lib_sizes_no_perm)),
         rep('rarefied', length(typeIerror_diff_lib_sizes_no_perm_rrf)),
         rep('original', length(typeIerror_mixture_no_perm)),
         rep('rarefied', length(typeIerror_mixture_no_perm_rrf))),
  condition=factor(c(rep('baseline', length(typeIerror_no_perm) + 
                    length(typeIerror_no_perm_rrf)), 
              rep('small sample sizes', 
                  length(typeIerror_small_sample_sizes_no_perm) + 
                    length(typeIerror_small_sample_sizes_no_perm_rrf)),
              rep('diff. mean library sizes', 
                  length(typeIerror_small_sample_sizes_no_perm) + 
                    length(typeIerror_small_sample_sizes_no_perm_rrf)),
              rep('mixture', 
                  length(typeIerror_mixture_no_perm) + 
                    length(typeIerror_mixture_no_perm_rrf))),
              levels=c('baseline', 'small sample sizes',
                       'diff. mean library sizes',
                       'mixture')),
  inference='DESeq2 (likelihood asymptotics)'
  )
perm_df <- data.frame(
  typeIerror=c(typeIerror_perm, 
               typeIerror_perm_rrf,
               typeIerror_small_sample_sizes_perm,
               typeIerror_small_sample_sizes_perm_rrf,
               typeIerror_diff_lib_sizes_perm,
               typeIerror_diff_lib_sizes_perm_rrf,
               typeIerror_mixture_perm,
               typeIerror_mixture_perm_rrf),
  data=c(rep('original', length(typeIerror_perm)),
         rep('rarefied', length(typeIerror_perm_rrf)),
         rep('original', length(typeIerror_small_sample_sizes_perm)),
         rep('rarefied', length(typeIerror_small_sample_sizes_perm_rrf)),
         rep('original', length(typeIerror_diff_lib_sizes_perm)),
         rep('rarefied', length(typeIerror_diff_lib_sizes_perm_rrf)),
         rep('original', length(typeIerror_mixture_perm)),
         rep('rarefied', length(typeIerror_mixture_perm_rrf))),
  condition=factor(c(rep('baseline', length(typeIerror_perm) + 
                           length(typeIerror_perm_rrf)), 
                     rep('small sample sizes', 
                         length(typeIerror_small_sample_sizes_perm) + 
                           length(typeIerror_small_sample_sizes_perm_rrf)),
                     rep('diff. mean library sizes', 
                         length(typeIerror_small_sample_sizes_perm) + 
                           length(typeIerror_small_sample_sizes_perm_rrf)),
                     rep('mixture', 
                         length(typeIerror_mixture_perm) + 
                           length(typeIerror_mixture_perm_rrf))),
                   levels=c('baseline', 'small sample sizes',
                            'diff. mean library sizes',
                            'mixture')),
  inference='DESeq2 + permutation inference'
)
typeIerror_df <- rbind(asymp_df, perm_df)

typeIerror_plot <- ggplot(typeIerror_df %>% group_by(data, inference, condition) %>%
                            summarize(mean_typeIerror=mean(typeIerror),
                                      U=mean(typeIerror) + 1.96 * sd(typeIerror) / sqrt(length(typeIerror)),
                                      L=mean(typeIerror) - 1.96 * sd(typeIerror) / sqrt(length(typeIerror))), 
                          aes(x=condition, y=mean_typeIerror, 
                              color=condition, shape=condition)) + 
  geom_point() +
  geom_errorbar(aes(ymax = U, ymin = L), width = 0.5) +
  facet_grid(data ~ inference) +
  ylab('95% CI of Type I error rate') +
  geom_hline(yintercept=0.05, lty=2) +
  ggtitle('Type I error rates for differential abundance testing') +
  theme(text=element_text(size=18)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  ylim(c(0, 0.25))

pdf('deseq_results/deseq_type_I_error.pdf',
    width=12)
typeIerror_plot
dev.off()

# Realistic ----

run_perm_sim_realistic <- function(num_sim=100,
                                   lib1, lib2,
                                   latent1, latent2,
                                   n1=30, n2=30,
                                   num_perm=99) {
  k1 <- nrow(latent1)
  k2 <- nrow(latent2)
  nspp <- ncol(latent1)
  rbindlist(mclapply(1:num_sim, function(i) {
    L1 <- sample(lib1, n1)
    L2 <- sample(lib2, n2)
    x1 <- sapply(L1, function(L) {
      rmultinom(n=1, size=L, prob=latent1[sample(1:k1, 1), ])
    })
    x2 <- sapply(L2, function(L) {
      rmultinom(n=1, size=L, prob=latent2[sample(1:k2, 1), ])
    })
    cnts <- cbind(x1, x2)
    L_star <- min(c(L1, L2))
    cnts_rrf <- t(data.frame(rbind(rrarefy(t(x1), L_star), 
                                   rrarefy(t(x2), L_star))))
    cond <- factor(c(rep(1, n1), rep(2, n2)))
    stat <- abs(get_stat_DESeq2(cnts, cond))
    stat_rrf <- abs(get_stat_DESeq2(cnts_rrf, cond))
    perm_stat <- rbindlist(lapply(1:num_perm, function(b) {
      cond_shuffled <- sample(cond)
      data.frame(original=abs(get_stat_DESeq2(cnts, cond_shuffled)),
                 rarefied=abs(get_stat_DESeq2(cnts_rrf, cond_shuffled)),
                 species=1:nspp)
    }))
    perm_pval <- sapply(1:nspp, function(i) {
      as.numeric(perm_stat %>% filter(species == i) %>% 
                   summarize(p_val=mean(c(stat[i], original) >= stat[i], na.rm=T)))
    })
    perm_pval_rrf <- sapply(1:nspp, function(i) {
      as.numeric(perm_stat %>% filter(species == i) %>% 
                   summarize(p_val=mean(c(stat_rrf[i], rarefied) >= stat_rrf[i], na.rm=T)))
    })
    
    p_val <- c(perm_pval, perm_pval_rrf)
    species <- c(1:nspp, 1:nspp)
    data <- c(rep('original', nspp), rep('rarefied', nspp))
    df <- data.frame(p_val=p_val,
                     species=species,
                     data=data)
    df
  }, mc.cores=12))
}


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
freezethaw <- library_sizes[
  sample_df_subset$sample_storage_temp_treatment == 'freezethaw'
  ]
fourC <- library_sizes[
  sample_df_subset$sample_storage_temp_treatment == '4C'
  ]
otu_freezethaw <- t(otu_table_subset[
  , sample_df_subset$sample_storage_temp_treatment == 'freezethaw'
])
otu_fourC <- t(otu_table_subset[
  , sample_df_subset$sample_storage_temp_treatment == '4C'
  ])

# Remove OTU with zero counts and normalize.
otu_latent <- otu_freezethaw[, colSums(otu_freezethaw) > 0]
otu_latent <- otu_latent / rowSums(otu_latent)
num_sim <- 100
nspp <- ncol(otu_latent)
set.seed(2000)
ptm <- proc.time()
results_realistic <- run_perm_sim_realistic(
  n1=30, n2=30,
  num_sim=num_sim, num_perm=99,
  lib1=freezethaw, lib2=fourC,
  latent1=otu_latent, latent2=otu_latent)
proc.time() - ptm

fwrite(results_realistic, file=paste0(results_dir, 'no_perm.csv'))

typeIerror_realistic_perm <- 
  compute_typeI_error(result=results_realistic,
                      data_type='original')
typeIerror_realistic_perm_rrf <- 
  compute_typeI_error(result=results_realistic,
                      data_type='rarefied')

df_realistic <- data.frame(data=rep(c('original', 'rarefied'), each=100),
                           typeIerror=c(typeIerror_realistic_perm,
                                          typeIerror_realistic_perm_rrf))
df_realistic_summary <- df_realistic %>%
  group_by(data) %>%
  summarize(mean_typeIerror=mean(typeIerror),
            U=mean(typeIerror) + 1.96 * sd(typeIerror) / sqrt(length(typeIerror)),
            L=mean(typeIerror) - 1.96 * sd(typeIerror) / sqrt(length(typeIerror)))
pdf('deseq_results/deseq_realistic_type_I_error.pdf',
    width=12)
ggplot(df_realistic_summary, aes(x=data, y=mean_typeIerror)) +
  geom_point() +
  geom_errorbar(aes(ymax = U, ymin = L), width = 0.5) +
  ylab('95% CI of Type I error rate') +
  geom_hline(yintercept=0.05, lty=2) +
  ggtitle('Type I error rates for differential abundance testing\nDESeq2 + permutation inference w/ simluations based on real data') +
  theme(text=element_text(size=18)) +
  ylim(c(0, 0.25))
dev.off()

