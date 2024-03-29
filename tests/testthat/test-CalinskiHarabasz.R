source("CVI_test_proc1.R")

CVI_fun <- CVI_CalinskiHarabasz
CVI_name <- "CalinskiHarabasz"

reference_fun <- function(X, y) {
    clusterCrit::intCriteria(X, y, "Calinski_Harabasz")[[1]]
}

# reference_fun <- if (require("clusterSim")) {
#     function(X, y) {
#         clusterSim::index.G1(X, y)
#     }
# } else NULL

CVI_test_proc1(CVI_name, CVI_fun, reference_fun)
