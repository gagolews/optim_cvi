library("testthat")
library("CVI")
context("internal")

test_that("iris", {
    library("datasets")
    data("iris")
    set.seed(123)

    X <- as.matrix(iris[,1:4])
    X[,] <- jitter(X) # otherwise we get a non-unique solution

    Ys <- list(
        as.integer(iris[[5]]),
        sample(4, nrow(X), replace=TRUE, prob=c(0.5, 0.3, 0.15, 0.05))
    )

    for (y in Ys) {

        K <- max(y)

        expect_error(.CVI_create("UnknownIndex", X, K))

        funs <- list(CVI_CalinskiHarabasz, CVI_DaviesBouldin, CVI_Silhouette,
            CVI_SilhouetteW, CVI_Dunn, CVI_WCSS, CVI_BallHall, CVI_Gamma)
        nams <- c("CalinskiHarabasz", "DaviesBouldin", "Silhouette",
            "SilhouetteW", "Dunn", "WCSS", "BallHall", "Gamma")

        for (u in seq_along(funs)) {
            i1 <- funs[[u]](X, y, K)
            cvi_ptr <- .CVI_create(nams[u], X, K)
            .CVI_set_labels(cvi_ptr, y)
            i2 <- .CVI_compute(cvi_ptr)

            expect_true(i1 == i2)

            for (i in 1:nrow(X)) {
                for (j in 1:K) {
                    if (y[i] == j) next

                    .CVI_modify(cvi_ptr, i, j)
                    i3 <- .CVI_compute(cvi_ptr)
                    y2 <- y
                    y2[i] <- j
                    i5 <- funs[[u]](X, y2, K)
                    expect_true(abs(i3 - i5)<1e-9)

                    .CVI_undo(cvi_ptr)
                    i4 <- .CVI_compute(cvi_ptr)
                    expect_true(abs(i2 - i4)<1e-9)
                }
            }
        }
    }
})

