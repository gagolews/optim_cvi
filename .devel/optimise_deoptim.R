options(stringsAsFactors=FALSE) # default in R 4.0
source("optimise.R")



par0_path <- "results_par0" # see results_par0/README.md (don't change this)


# https://github.com/gagolews/clustering_benchmarks_v1
benchmarks_path <- "~/clustering/clustering_benchmarks_v1"



# during development/testing, you can omit the processing of large datasets
max_n <- 1000 # Inf for no limit  (number of points in the dataset)
max_K <- 5 # Inf for no limit  (number of reference clusters)


# your method's name goes here:
output_path <- "results_deoptim"

# cluster validity index
CVI_name <- "CalinskiHarabasz"





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
  res <- deoptim(X, K, labels=par, control=list(itermax=100, trace=FALSE), CVI_ptr=CVI_ptr)

  list(
    par=as.integer(res$optim$bestmem),
    value=-res$optim$bestval
  )
}


do_optimise_all(
  optimiser,
  CVI_name,
  output_path,
  benchmarks_path,
  par0_path,
  max_n,
  max_K
)
