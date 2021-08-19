#!/usr/bin/env Rscript

# Copyleft 2020, Marek Gagolewski


options(stringsAsFactors=FALSE) # default in R 4.0
source("optimise.R")

par0_path <- "results_par0" # see results_par0/README.md (don't change this)

# your method's name goes here:
output_path <- "results_dummy"

# https://github.com/gagolews/clustering_benchmarks_v1
benchmarks_path <- "~/Projects/clustering_benchmarks_v1"


args <- optimise_read_args()
cat(sprintf("args=(%s); output_path=%s\n", stri_flatten(unlist(args), ", "), output_path))
attach(args)  # sets CVI_name, datasets_subset_regex, datasets_subset_indices, max_n, max_K





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
    par <- identity(par) # LOL, we're not doing anything, this is just an example :P

    .CVI_set_labels(CVI_ptr, par)
    value <- .CVI_compute(CVI_ptr)

    stopifnot(!is.na(value))

    list(
        par=par,
        value=value
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
    datasets_subset_indices=datasets_subset_indices
)
