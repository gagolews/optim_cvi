#!/usr/bin/env Rscript

# Copyleft 2020-2021, Marek Gagolewski

#pdf("/tmp/Rplots.pdf")
library("CVI")
library("stringi")
set.seed(123)
source("load_data.R")


# https://github.com/gagolews/clustering_benchmarks_v1
benchmarks_path <- "~/Projects/clustering_benchmarks_v1"

# https://github.com/gagolews/clustering_benchmarks_v1_results
results_path <- "~/Projects/clustering_benchmarks_v1_results/original"

output_path <- "results_par0"



sub_paths <- list.dirs(results_path, full.names=FALSE, recursive=FALSE)

datasets <- stri_sub(list.files(benchmarks_path, glob2rx("*.data.gz"),
    full.names=FALSE, recursive=TRUE), to=-9)

for (dataset in datasets) {

    X <- load_dataset(dataset, benchmarks_path)
    n <- nrow(X)
    d <- ncol(X)

    Y_true <- load_true_labels(dataset, benchmarks_path)
    Y_true <- apply(Y_true, 2, merge_noise_points_with_nearest_clusters, X)
    Ks <- apply(Y_true, 2, max)

    for (K in unique(Ks)) {
        cat(sprintf("%20s: n=%6d, d=%3d, K=%3d\n", dataset, n, d, K))

        out_fname <- file.path(output_path,
            sprintf("%s.result%d.gz", dataset, K))
        if (!dir.exists(dirname(out_fname)))
            dir.create(dirname(out_fname))

        if (file.exists(out_fname)) next

        Y_start <- do.call(cbind, lapply(sub_paths, function(sub_path)
            load_pred_labels(dataset, results_path, sub_path, K)))

        Y_start <- cbind(Y_true[, K==Ks, drop=FALSE], Y_start)

        stopifnot(apply(Y_start, 2, min)==1)
        stopifnot(apply(Y_start, 2, max)==K)

        f <- gzfile(out_fname, "w")
        write.csv(Y_start, f, row.names=FALSE)
        close(f)
        rm(f)
    }
}
