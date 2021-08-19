options(stringsAsFactors=FALSE) # default in R 4.0
source("optimise.R")


par0_path <- "results_par0" # see results_par0/README.md (don't change this)


# https://github.com/gagolews/clustering_benchmarks_v1
benchmarks_path <- "~/clustering/clustering_benchmarks_v1"



# during development/testing, you can omit the processing of large datasets
max_n <- 10000 # Inf for no limit  (number of points in the dataset)
max_K <- Inf # Inf for no limit  (number of reference clusters)

# regexes to match datasets' paths:
datasets_subset_regex <- list(
  ".*", # all
  "fcps|graves|other|sipu|uci|wut",
  "g2mg",  # a lower-priority one
  "h2mg",  # a lower-priority one
  "mnist"  # too big
)[[1]]


# logical/integer indices of datasets to select for processing
datasets_subset_indices <- list(
  TRUE, # all
  c(TRUE,  FALSE), # every 2nd - batch 1
  c(FALSE, TRUE) # every 2nd - batch 2
)[[3]]


# your method's name goes here:
output_path <- "results_pso_many_centroids"

# cluster validity index
cluster_validity_index <- 2
CVI_name <- c("CalinskiHarabasz", "DaviesBouldin", "Silhouette")[cluster_validity_index]
fn_scale <- c(1, 1, NA)[cluster_validity_index]
result_scale <- c(-1, -1, NA)[cluster_validity_index]




# @title Your custom label vector optimiser
#
# @param par initial label vector
# @param CVI_ptr will be created by do_optimise_all, actual cluster validity index
# @param X data matrix
# @param K number of clusters
#
# @return A list with at least the following named components, see ?optim:
#     * par    output label vector
#     * value  value of CVI at the result
optimiser <- function(par, CVI_ptr, X, K)
{
  res <- pso_many_centroids(X, K, 5L, labels=par, control=list(maxit=1000, reltol=0.2, fnscale=fn_scale), CVI_ptr=CVI_ptr)

  list(
    par=res$cluster_indices,
    value= result_scale * res$res$value
  )
}


do_optimise_all(
  optimiser,
  CVI_name,
  output_path,
  benchmarks_path,
  par0_path,
  max_n,
  max_K,
  datasets_subset_regex=datasets_subset_regex,
  datasets_subset_indices=datasets_subset_indices,
  labels_subset_indices = 1:5
)
