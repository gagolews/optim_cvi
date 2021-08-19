source("CVI_test_proc1.R")

CVI_fun <- CVI_SilhouetteW
CVI_name <- "SilhouetteW"

reference_fun <- function(X, y) {
    clusterCrit::intCriteria(X, y, "Silhouette")[[1]] # returns NA/NaN on singletons
}

CVI_test_proc1(CVI_name, CVI_fun, reference_fun)
