#!/usr/bin/env Rscript

# Copyleft 2020-2021, Marek Gagolewski


# Tabu Turbo --
# For each initial partition in a given set:
#     if the partition is on the Tabu list, ignore it.
#
# Even for multiple starting points, it returns a single result.
#
# Explores the search space more extensively than independent restarts
# from different points, as the Tabu list is shared between iterations.
#
# Pros: will identify the (local) maximum faster
# Cons: we won't know where we converge from the individual starting points,
#       but we don't need this info in our current setting.



options(stringsAsFactors=FALSE) # default in R 4.0
source("optimise.R")

par0_path <- "results_par0" # see results_par0/README.md (don't change this)


# your method's name goes here:
output_path <- "results_tabu_turbo"


# https://github.com/gagolews/clustering_benchmarks_v1
benchmarks_path <- "~/Projects/clustering_benchmarks_v1"


args <- optimise_read_args()
cat(sprintf("args=(%s); output_path=%s\n", stri_flatten(unlist(args), ", "), output_path))
attach(args)  # sets CVI_name, datasets_subset_regex, datasets_subset_indices, max_n, max_K






do_optimise_all_tabu_turbo <- function(
        optimiser,
        CVI_name,
        output_path,
        benchmarks_path,
        par0_path="results_par0",
        max_n=Inf,
        max_K=Inf,
        datasets_subset_regex=".*",
        datasets_subset_indices=TRUE
)
{
    datasets <- stri_sub(list.files(benchmarks_path, glob2rx("*.data.gz"),
        full.names=FALSE, recursive=TRUE), to=-9)
    datasets <- stri_sort(datasets, numeric=TRUE) # natural sorting of numbers
    datasets <- stri_subset_regex(datasets, datasets_subset_regex)
    datasets <- datasets[datasets_subset_indices]
    datasets <- na.omit(datasets)

#     cat(CVI_name, "\n\n")

    for (dataset in datasets)
    {
        set.seed(123)
        X <- load_dataset(dataset, benchmarks_path)
        n <- nrow(X)
        d <- ncol(X)
        if (n > max_n) {
#             cat(sprintf("[%24s] %40s: n=%6d, d=%3d        <SKIPPING (n>max_n)>\n", CVI_name, dataset, n, d))
            next
        }

        path <- sprintf("%s.result*.gz", file.path(par0_path, dataset))
        files <- list.files(dirname(path), glob2rx(basename(path)),
            recursive=TRUE,
            full.names=TRUE)

        Ks <- as.numeric(stri_match_first_regex(files, "result([0-9]+)\\.gz$")[,2])
        stopifnot(length(Ks) > 0, !is.na(Ks), Ks > 1)

        for (K in Ks) {
            if (K > max_K) {
#                 cat(sprintf("[%24s] %40s: n=%6d, d=%3d, K=%3d <SKIPPING (K>max_K)>\n", CVI_name, dataset, n, d, K))
                next
            }

            Y <- load_pred_labels(dataset, ".", par0_path, K)
            stopifnot(n == nrow(Y), ncol(Y) >= 1, min(Y) == 1, max(Y) == K)

            Y <- Y[, !duplicated(Y, MARGIN=2), drop=FALSE]  # remove duplicates!
            # TODO: in fact, we could try removing more labels,
            #       as (1, 2, 1) === (2, 1, 2)

            out_fname <- file.path(output_path, CVI_name,
                sprintf("%s.result%d.gz", dataset, K))
            out_fname2 <- file.path(output_path, CVI_name,
                sprintf("%s.result%d.csv", dataset, K))
            if (!dir.exists(dirname(out_fname)))
                dir.create(dirname(out_fname), recursive=TRUE)

            if (file.exists(out_fname)) {
#                 cat(sprintf("[%24s] %40s: n=%6d, d=%3d, K=%3d     <ALREADY EXISTS>\n", CVI_name, dataset, n, d, K))
                next
            }

            cat(sprintf("[%24s] %40s: n=%6d, d=%3d, K=%3d                     \n", CVI_name, dataset, n, d, K))

            tryCatch({
                CVI_ptr <- .CVI_create(CVI_name, X, K)
                Y <- Y[, order(sapply(seq_len(ncol(Y)), function(i) {
                    .CVI_set_labels(CVI_ptr, Y[,i])
                    .CVI_compute(CVI_ptr)
                }), decreasing=TRUE), drop=FALSE]

                out_pars <- matrix(NA_real_, nrow=nrow(Y), ncol=1,
                    dimnames=list(NULL, "tabu_turbo"))
                out_vals <- rep(NA_real_, 1)

                res <- .CVI_improve_turbo(CVI_ptr, Y, verbose=TRUE)
                stopifnot(min(res$par)==1)
                stopifnot(max(res$par)==K)
                stopifnot(length(res$par)==n)
                #stopifnot(is.finite(res$value))
                out_vals <- res$value
                out_pars <- matrix(res$par, ncol=1, dimnames=list(NULL, "tabu_turbo"))

                f <- gzfile(out_fname, "w")
                write.csv(out_pars, f, row.names=FALSE)
                close(f)

                write.csv(data.frame(
                    par=dimnames(out_pars)[[2]],
                    value=out_vals
                ), out_fname2, row.names=FALSE)
            },
            error=function(e) {
                print(e)
            })
        }
    }
}


do_optimise_all_tabu_turbo(
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
