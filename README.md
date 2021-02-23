# Code for "To rarefy or not to rarefy: robustness and efficiency trade-offs of rarefying microbiome data"

This repository contains scripts for simluation studies and data analysis for the paper "To rarefy or not to rarefy: robustness and efficiency trade-offs of rarefying microbiome data".


- `permanova_example.R`: an example of how varying library sizes can lead to inflated Type I error rate for PERMANOVA (permutational analysis of variance)
- `rarefaction_efficiency_index.R`: contains a function computing the (sample) rarefaction efficiency index (REI) and a usage example of REI
- `sim_deseq.R`: Type I error simulation studies for DESeq2 with hypothetical scenarios and realistic data; remark: this script takes almost a day to run
- `sim_deseq_100.R`: Type I error simulation studies with sample size n = 100; check if the results are similar to n = 30
- `sim_richness.R`: Type I error simulation studies for species richness comparisons
- `sim_power_studies_with_real_data.R`: power simulation studies using real world data
