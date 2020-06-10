# Check if results are different if number of observations is increased to 100.
# Conclusions are the same as using n = 30.

ptm <- proc.time()
results_no_perm_100 <- 
  run_no_perm_sim(num_sim=100, 
                  nsample=100,
                  lambda1=10000, lambda2=10000,
                  p1=rep(1/nspp, nspp), 
                  p2=rep(1/nspp, nspp), 
                  pi0=0, 
                  alpha1=30, alpha2=30)
runtime_100 <- proc.time() - ptm
fwrite(results_no_perm_100,
       paste0(results_dir, 'no_perm_100.csv'))
typeIerror_no_perm_100 <- compute_typeI_error(result=results_no_perm_100,
                                          data_type='original')
typeIerror_no_perm_rrf_100 <- compute_typeI_error(result=results_no_perm_100,
                                              data_type='rarefied')
tmp <- data.frame(typeIerror=c(typeIerror_no_perm_100,
                               typeIerror_no_perm_rrf_100),
                  data=c(rep('original', length(typeIerror_no_perm_100)),
                         rep('rarefied', length(typeIerror_no_perm_100))))
pdf('deseq_results/deseq_100.pdf')
ggplot(tmp, aes(x=data, y=typeIerror, color=data)) +
  geom_boxplot() + ggtitle('DESEq2, Wald test, n=100, baseline')
dev.off()



ptm <- proc.time()
results_diff_lib_sizes_no_perm_100 <- 
  run_no_perm_sim(num_sim=100, 
                  nsample=100,
                  lambda1=10000, lambda2=50000,
                  p1=rep(1/nspp, nspp), 
                  p2=rep(1/nspp, nspp), 
                  pi0=0, 
                  alpha1=30, alpha2=30)
runtime_diff_lib_sizes_100 <- proc.time() - ptm

fwrite(results_diff_lib_sizes_no_perm_100,
       paste0(results_dir, 'diff_lib_sizes_no_perm_100.csv'))
typeIerror_diff_lib_sizes_no_perm_100 <- 
  compute_typeI_error(result=results_diff_lib_sizes_no_perm_100,
                      data_type='original')
typeIerror_diff_lib_sizes_no_perm_rrf_100 <- 
  compute_typeI_error(result=results_diff_lib_sizes_no_perm_100,
                      data_type='rarefied')
tmp_diff <- data.frame(typeIerror=c(typeIerror_diff_lib_sizes_no_perm_100,
                               typeIerror_diff_lib_sizes_no_perm_rrf_100),
                  data=c(rep('original', length(typeIerror_diff_lib_sizes_no_perm_100)),
                         rep('rarefied', length(typeIerror_diff_lib_sizes_no_perm_100))))
pdf('deseq_results/deseq_100_diff_lib_sizes.pdf')
ggplot(tmp_diff, aes(x=data, y=typeIerror, color=data)) +
  geom_boxplot() + ggtitle('DESEq2, Wald test, n=100, different mean library sizes')
dev.off()