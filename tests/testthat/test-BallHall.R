source("CVI_test_proc1.R")

CVI_fun <- CVI_BallHall
CVI_name <- "BallHall"

reference_fun <- function(X, y) {
    wcss <- 0.0
    K <- max(y)
    for (i in 1:K) {
        wi <- which(y==i)
        if (length(wi) <= 1) next
        d <- dist(X[wi, ])^2
        # d gives only the unique pairs - object of class dist
        wcss <- wcss + sum(d)/length(wi)**2
    }
    -wcss
}

CVI_test_proc1(CVI_name, CVI_fun, reference_fun)
