# Copyleft 2020-2021, Marek Gagolewski


load_dataset <- function(dataset, benchmarks_path) {
    data_file <- file.path(benchmarks_path, paste0(dataset, ".data.gz"))
    X    <- read.table(data_file)
    X <- X[, apply(X, 2, var) > 0, drop=FALSE] # remove all columns of 0 variance
    n <- nrow(X)
    d <- ncol(X)

    # add a tiny bit of random noise, center and scale
    X <- t(X)
    X[,] <- X + rnorm(n*d, 0, apply(X, 1, sd)*1e-6)
    X[,] <- X - apply(X, 1, mean)
    X[,] <- X / sd(X) # this is not standardisation
    X <- t(X)

    X
}


load_true_labels <- function(dataset, benchmarks_path) {
    path <- file.path(benchmarks_path, sprintf("%s.labels*.gz", dataset))
    files <- list.files(dirname(path),
        glob2rx(basename(path)),
        recursive=TRUE,
        full.names=TRUE)
    all_labels <- lapply(files, read.table)
    Y <- as.matrix(do.call(cbind, all_labels))
    dimnames(Y) <- list(NULL, paste0("labels", 0:(ncol(Y)-1)))
    Y
}


load_pred_labels <- function(dataset, results_path, sub_path, K) {
    path <- file.path(results_path, sub_path,
        sprintf("%s.result%d.gz", dataset, K))
    if (file.exists(path)) {
        as.matrix(read.csv(path))
    } else {
        cat(sprintf("file does not exist: %s\n", path))
        NULL
    }
}


merge_noise_points_with_nearest_clusters <- function(y, X) {
    stopifnot(length(y) == nrow(X))

    wh_noise <- which(y == 0)
    if (length(wh_noise) == 0) return(y)

    nns <- FNN::get.knnx(X[-wh_noise, , drop=FALSE],
                         X[ wh_noise, , drop=FALSE], 1)$nn.index
    y[wh_noise] <- y[-wh_noise][nns]

    y
}


crossprintf <- function(fmt, ...) {
    args <- c(rev(list(...)), stringsAsFactors=FALSE)
    x <- rev(do.call(expand.grid, args))
    sapply(seq_len(nrow(x)), function(i)
        do.call(sprintf, c(fmt, x[i, , drop=FALSE])))
}
