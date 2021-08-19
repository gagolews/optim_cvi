source("CVI_test_proc1.R")

CVI_fun <- CVI_Gamma
CVI_name <- "Gamma"

reference_fun <- function(X, y) {
    clusterCrit::intCriteria(X, y, "Gamma")[[1]]
}

CVI_test_proc1(CVI_name, CVI_fun, reference_fun)
