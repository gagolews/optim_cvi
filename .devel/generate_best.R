#!/usr/bin/env Rscript

# Copyleft 2020-2021, Marek Gagolewski
#
# For each dataset, K, and CVI:
# a) fetch all the label vectors and the corresponding CVI values,
# b) identify the(a) label vector that yields the greatest CVI,
# c) store it as the output of the "Maximise the CVI" algorithm.

# This script is meant to be executed frequently, preferably
# after each batch of of some optimiser's run,
# e.g., optimise_tabu_turbo on a few datasets.



suppressMessages(library("CVI"))
suppressMessages(library("stringi"))
set.seed(123)

warn_missing <- TRUE

ignore_results <- c(
    # datasets which are accompanied by labels that yield n*(k-1) > 50000:
    "^mnist/digits",
    "^mnist/fashion",
    "^other/chameleon_t7_10k",
    "^other/chameleon_t8_8k",
    "^sipu/a1",
    "^sipu/a2",
    "^sipu/a3",
    "^sipu/birch1",
    "^sipu/birch2",
    "^sipu/d31",
    "^sipu/s1",
    "^sipu/s2",
    "^sipu/s3",
    "^sipu/s4",
    "^sipu/worms_2",
    "^sipu/worms_64",
    # the 1-dimensional datasets:
    "^g2mg/g2mg_1_[0-9]+",
    "^h2mg/h2mg_1_[0-9]+",
    # in fact, g2mg and h2mg are left for further research:
    "^g2mg/g2mg",
    "^h2mg/h2mg",
    # datasets with clusters of less than 25 elements:
    "^uci/glass"
)

ignore_CVIs <- c(
    "^DELME.*",
    "^WCSS$",       # == CalinskiHarabasz (when maximised)
    "^WCNN_(1|10)$",
    "^Dunn$",       # == GDunn_d1_D1
    "^Gamma$",      # too slow
    "^GDunn_d6_D.$" # too slow
)


input_paths <- c(
    "results_dummy",
    "results_tabu_turbo",
    "results_tabu_push",
    "results_tabu_random",
    "results_tabu_random2",
    "results_deoptim_many_centroids",
    "results_pso_many_centroids"
    #"results_deoptim",  # abandoned at some point
    #"results_pso", # abandoned at some point
    #"results_tabu" # abandoned at some point
)


output_path <- "results_best"

if (!dir.exists(output_path)) {
    dir.create(output_path)
} else {
    stop("output directory already exists.")
}

writeLines("
README
======

This directory and all its subdirectories were generated automatically as
follows. For each dataset, number of clusters K, and cluster validity
measure CVI:

1. fetch all the candidate label vectors and the corresponding CVI values,
2. identify the(a) label vector that yields the greatest CVI,
    using the procedure described in the paper by
    Gagolewski, Bartoszuk, and Cena: 'Are Cluster Validity Measures (In)valid?',
3. store it as a new column in the output CSV file.

Each file is named like 'test_battery/dataset.resultsK.gz', where
K is the number of clusters to detect.

Each file is in a CSV format -- (named) columns give the partitions
with elements in 1..K.

", file.path(output_path, "README.md"))





all_files <- unlist(list.files(input_paths, glob2rx("*.result*.csv"),
    full.names=TRUE, recursive=TRUE))
all_files <- stri_match_first_regex(all_files, "^([^/]+)/([^/]+)/([^/]+/[^/]+)\\.csv$")
stopifnot(all(!is.na(all_files)))
dimnames(all_files) <- list(NULL, c("filename", "optim", "CVI", "result"))

regex_omitter <- function(x, y) stri_subset_regex(x, y, negate=TRUE)

all_results <- stri_sort(Reduce(regex_omitter, ignore_results, stri_unique(all_files[, "result"])), numeric=TRUE)
all_CVIs <- stri_sort(Reduce(regex_omitter, ignore_CVIs, stri_unique(all_files[, "CVI"])), numeric=TRUE)
all_optims <- stri_sort(stri_unique(all_files[, "optim"]), numeric=TRUE)



all_datasets <- stri_match_first_regex(all_results, "^((.*)/(.*))\\.result(.*)$")[, -1]
all_datasets <- unlist(Map(stri_flatten, split(all_datasets[, 4], all_datasets[, 1]), collapse=", "))
all_datasets <- all_datasets[stri_order(names(all_datasets), numeric=TRUE)]
all_datasets <- all_datasets[order(match(stri_sub(names(all_datasets), 1, 4), c("g2mg", "h2mg"), 0))]
cat("      dataset                          Ks\n")
cat(sprintf("%5d %-32s (%s)", seq_along(all_datasets), names(all_datasets), all_datasets), sep="\n")
cat("\n")
cat("      CVI\n")
cat(sprintf("%5d %s", seq_along(all_CVIs), all_CVIs), sep="\n")


for (result_i in seq_along(all_results)) {
    result <- all_results[result_i]
    num_missing <- 0
    all_best_labels <- structure(
        lapply(all_CVIs, function(CVI) {
            vals <- lapply(all_optims, function(optim) {
                fname <- file.path(optim, CVI, result %s+% ".csv")
                if (!file.exists(fname)) {
                    num_missing <<- num_missing+1
                    if (warn_missing) cat(sprintf("   missing file: %s\n", fname))
                    return(NULL)
                }
                cbind(optim=optim, CVI=CVI, read.csv(fname))
            })
            vals <- do.call(rbind.data.frame, vals)
            if (length(vals) == 0) return(NULL)
            best_val    <- vals[which.max(vals$value), ]
            best_labels <- read.csv(file.path(best_val$optim, best_val$CVI, result %s+% ".gz"))[[best_val$par]]
        }), names=all_CVIs)
    all_best_labels <- do.call(cbind, all_best_labels)

    out_fname <- file.path(output_path, result %s+% ".gz")
    if (!dir.exists(dirname(out_fname)))
        dir.create(dirname(out_fname))

    f <- gzfile(out_fname, "w")
    write.csv(all_best_labels, f, row.names=FALSE)
    close(f)
    rm(f)

    cat(sprintf("(%3d/%3d) %50s %s\n",
        result_i, length(all_results), result,
        if (num_missing==0) "" else sprintf("(%d missing)", num_missing)))
}
