source("CVI_test_proc1.R")

CVI_fun <- CVI_Silhouette
CVI_name <- "Silhouette"

reference_fun <- function(X, y) {
    clusterSim::index.S(dist(X), y)
}

CVI_test_proc1(CVI_name, CVI_fun, reference_fun)
