# Obtain datasets from the R package MicrobeDS. 
# devtools::install_github("twbattaglia/MicrobeDS")

library(MicrobeDS)
library(gridExtra)
library(grid)

plot_text_size <- 25
annotate_text_size <- 10
title_size <- 30

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

# Effects of preservation and storage conditions on the fecal microbiome ----
data('qa10394')
sample_data_10394 <- read.table('sample_data_10394.txt', header=TRUE, sep='\t')
sample_data_10394$name <- paste0('10394.', sample_data_10394$description)
sample_df <- data.frame(name=colnames(otu_table(qa10394)), 
                        stringsAsFactors=F)
sample_df <- sample_df %>% left_join(sample_data_10394, by='name')
well_defined <- !grepl('BLANK', sample_df$host_subject_id) & 
  !grepl('FTA', sample_df$host_subject_id) &
  !grepl('mistake', sample_df$host_subject_id)

L_star <- 30000

# treatments as groups ----
# Use only 95etoh+20C+4week and 95etoh+4C+4week
bool <- well_defined &
  sample_df$sample_preservation_method == '95etoh' & 
  (sample_df$sample_storage_temp_treatment == '20C' | 
     sample_df$sample_storage_temp_treatment == '4C') &
  sample_df$duration_of_storage == '4weeks'
otu_subset <- otu_table(qa10394)[, bool]
group <- sample_df$sample_storage_temp_treatment[bool]
REIs <- REI(otu_matrix=otu_subset, rarefied_depth=L_star, group=group,
    taxa_are_rows=TRUE)
mean(REIs, na.rm=TRUE)

# subjects as groups ----
bool <- well_defined & 
  (sample_df$host_subject_id == 'H4' |
     sample_df$host_subject_id == 'H5')
otu_subset <- otu_table(qa10394)[, bool]
group <- sample_df$host_subject_id[bool]
REIs <- REI(otu_matrix=otu_subset, rarefied_depth=L_star, group=group,
            taxa_are_rows=TRUE)
mean(REIs, na.rm=TRUE)

