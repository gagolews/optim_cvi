source("CVI_test_proc1.R")

CVI_fun <- CVI_Dunn
CVI_name <- "Dunn"

reference_fun <- function(X, y) {
    clusterCrit::intCriteria(X, y, "Dunn")[[1]]
}

# reference_fun <- if (require("clValid")) {
#     function(X, y) dunn(dist(X), y)
# } else NULL

CVI_test_proc1(CVI_name, CVI_fun, reference_fun)
