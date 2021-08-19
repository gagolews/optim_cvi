options(stringsAsFactors=FALSE) # default in R 4.0
suppressMessages(library("CVI"))
suppressMessages(library("stringi"))
source("load_data.R")




# processes all benchmark sets for all benchmark Ks
#
# @param optimiser function
# @param CVI_name string
# @param output_path
# @param benchmarks_path
# @param par0_path
# @param max_n
# @param max_K
# @param datasets_subset_regex regex to match benchmark paths
# @param datasets_subset_indices a vector of indices or a logical vector
#
# @return nothing
do_optimise_all <- function(
        optimiser,
        CVI_name,
        output_path,
        benchmarks_path,
        par0_path="results_par0",
        max_n=Inf,
        max_K=Inf,
        datasets_subset_regex=".*",
        datasets_subset_indices=TRUE,
        labels_subset_indices = NULL
)
{
    datasets <- stri_sub(list.files(benchmarks_path, glob2rx("*.data.gz"),
        full.names=FALSE, recursive=TRUE), to=-9)
    datasets <- stri_sort(datasets, numeric=TRUE) # natural sorting of numbers
#     print(stri_subset_regex(datasets, datasets_subset_regex, negate=TRUE))
    datasets <- stri_subset_regex(datasets, datasets_subset_regex)
    datasets <- datasets[datasets_subset_indices]
    datasets <- na.omit(datasets)



    for (dataset in datasets)
    {
        set.seed(123)
        X <- load_dataset(dataset, benchmarks_path)
        n <- nrow(X)
        d <- ncol(X)
        if (n > max_n) {
#             cat(sprintf("[%24s] %24s: n=%6d, d=%3d <SKIPPING>\n", CVI_name, dataset, n, d))
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
#                 cat(sprintf("[%24s] %24s: n=%6d, d=%3d, K=%3d <SKIPPING>\n", CVI_name, dataset, n, d, K))
                next
            }

            out_fname <- file.path(output_path, CVI_name,
                sprintf("%s.result%d.gz", dataset, K))
            out_fname2 <- file.path(output_path, CVI_name,
                sprintf("%s.result%d.csv", dataset, K))
            if (!dir.exists(dirname(out_fname)))
                dir.create(dirname(out_fname), recursive=TRUE)

            if (file.exists(out_fname)) {
#                 cat(sprintf("[%24s] %24s: n=%6d, d=%3d, K=%3d <ALREADY EXISTS>\n", CVI_name, dataset, n, d, K))
                next
            }

            cat(sprintf("[%24s] %24s: n=%6d, d=%3d, K=%3d\n", CVI_name, dataset, n, d, K))

            CVI_ptr <- .CVI_create(CVI_name, X, K)

            Y <- load_pred_labels(dataset, ".", par0_path, K)
            stopifnot(n == nrow(Y), ncol(Y) >= 1, min(Y) == 1, max(Y) == K)

            if (!is.null(labels_subset_indices) && length(labels_subset_indices) <= ncol(Y)) {
              if (is.numeric(labels_subset_indices)) {
                labels_subset_indices_ordered <- order(sapply(seq_len(ncol(Y)), function(i) {
                  .CVI_set_labels(CVI_ptr, Y[,i])
                  .CVI_compute(CVI_ptr)
                }), decreasing=TRUE)[labels_subset_indices]
              }
              Y <- Y[, labels_subset_indices_ordered]
            }

            Y_dups <- duplicated(Y, MARGIN=2)

            out_pars <- matrix(NA_real_, nrow=nrow(Y), ncol=ncol(Y),
                dimnames=dimnames(Y))
            out_vals <- rep(NA_real_, ncol(Y))
            for (i in which(!Y_dups)) {
                res <- optimiser(Y[,i], CVI_ptr, X, K)
                #print(res)
                stopifnot(min(res$par)==1)
                stopifnot(max(res$par)==K)
                stopifnot(length(res$par)==n)
                out_vals[i]  <- res$value
                out_pars[,i] <- res$par
            }
            for (i in which(Y_dups)) {
                # duplicated label vector -> fetch from "cache"
                w <- 1; while (any(Y[,i] != Y[,w])) { w <- w+1; stopifnot(w<=ncol(Y)) }
                out_vals[i]  <- out_vals[w]
                out_pars[,i] <- out_pars[,w]
            }
            stopifnot(!is.na(out_vals))

            f <- gzfile(out_fname, "w")
            write.csv(out_pars, f, row.names=FALSE)
            close(f)

            write.csv(data.frame(
                par=dimnames(Y)[[2]],
                value=out_vals
            ), out_fname2, row.names=FALSE)
        }
    }
}



optimise_read_args <- function() {

    # regexes to match datasets' paths:
    .datasets_subset_regex <- c(
        "fcps|graves|other/(?!chameleon_t[78]_)|sipu/(?!(a[123]|birch[12]|d31|s[1234]|worms_.*)$)|uci|wut",
        "g2mg/g2mg_(?!1_)",
        "h2mg/h2mg_(?!1_)",
        ".*" # all
    )

    # logical/integer indices of datasets to select for processing
    .datasets_subset_indices <- list(
        TRUE, # all
        c(TRUE,  FALSE, FALSE), # every 3rd - batch 1
        c(FALSE, TRUE,  FALSE), # every 3rd - batch 2
        c(FALSE, FALSE, TRUE)   # every 3rd - batch 3
    )


    # cluster validity index
    .CVI_name <- c(
        "CalinskiHarabasz", "DaviesBouldin", "Silhouette", "SilhouetteW",
        "Dunn", "WCSS", "BallHall", "Gamma",
        crossprintf("GDunn_d%d_D%d", 1:6, 1:3),
        crossprintf("WCNN_%d", c(1, 5, 10, 25)),
        crossprintf("DuNN_%d_%s_%s", c(5, 25), c("Min", "Max", "Mean"), c("Min", "Max", "Mean", "Const")),
        crossprintf("DuNN_%d_%s_%s", c(25), c("SMax:5"), c("SMin:5", "Min", "Const")),
        crossprintf("DuNN_%d_%s_%s", c(25), c("SMin:5"), c("SMax:5", "Max", "Const"))
    )

    .max_n <- Inf
    .max_K <- Inf


    args <- commandArgs(TRUE)
    if (length(args) < 1) {
        cat("Usage: <this_script> CVI_name [datasets_subset_regex [max_n [max_K]]]\n")
        cat("where CVI_name and datasets_subset_regex are integer indexes:\n")
        print(list(
            CVI_name=as.matrix(.CVI_name),
            datasets_subset_regex=as.matrix(.datasets_subset_regex),
            max_n=.max_n,
            max_K=.max_K
    #         datasets_subset_indices=as.matrix(.datasets_subset_indices)
        ))
        cat("For example: <this_script> 42 1 1000 5\n")
        stop("Incorrect number of command args.")
    }


    CVI_name <- .CVI_name[[as.numeric(args[1])]]
    datasets_subset_regex <- .datasets_subset_regex[[1]]
    max_n <- .max_n
    max_K <- .max_K
    datasets_subset_indices <- .datasets_subset_indices[[1]]

    if (length(args) >= 2)
        datasets_subset_regex <- .datasets_subset_regex[[as.numeric(args[2])]]

    if (length(args) >= 3)
        max_n <- as.numeric(args[3]) # Inf for no limit  (number of points in the dataset)

    if (length(args) >= 4)
        max_K <- as.numeric(args[4]) # Inf for no limit  (number of reference clusters)


    list(
        CVI_name                 = CVI_name,
        datasets_subset_regex    = datasets_subset_regex,
        datasets_subset_indices  = datasets_subset_indices,
        max_n                    = max_n,
        max_K                    = max_K
    )
}
