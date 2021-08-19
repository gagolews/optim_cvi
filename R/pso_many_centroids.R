## Many centroids PSO algorithm
##
## Copyright (C) 2020,
##     Maciej Bartoszuk (http://bartoszuk.rexamine.com),
##     Anna Cena (http://cena.rexamine.com),
##     Marek Gagolewski (https://www.gagolewski.com).
##
## This program is free software licensed under the European Union Public
## License (EUPL) v. 1.2 or – as soon they will be approved by the European
## Commission - subsequent versions of the EUPL (the "Licence"). You may not
## use this work except in compliance with the Licence. You may obtain a copy
## of the Licence at: https://joinup.ec.europa.eu/software/page/eupl.
##
## This project is a work in progress, which might be subject to continuous
## improvement by numerous contributors. It is not a finished work and may
## therefore contain defects or "bugs" inherent to this type of development.
## It is distributed in the hope that it will be useful. Unless required
## by applicable law or agreed to in writing, software distributed under
## the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR
## CONDITIONS OF ANY KIND, either express or implied. See the Licence for the
## specific language governing permissions and limitations under the Licence.

#' @title
#' The Many Centroids PSO algorithm
#' 
#' @description
#' The Many Centroids PSO algorithm is used to maximize a given cluster validity index
#' 
#' @param data_matrix a nxk matrix representing points to cluster
#' @param k a number of clusters
#' @param centroids_per_cluster numer of centroids for every cluster
#' @param labels an optional label vector from which the algorithm will start
#' @param control a list of control parameters given to PSO algorithm, most notably maxit
#' @param CVI_ptr optional cluster validity index calculator to optimize. When NULL, CalinskiHarabasz is created and used. 
#' 
#' @details
#' TBD
#' 
#' @return 
#' Returns the same value as psoptim.
#' 
#' @export
#' @family CVI
pso_many_centroids <- function(data_matrix, k, centroids_per_cluster=5, labels=NULL, control=list(maxit=10000, reltol=0.2), CVI_ptr=NULL)
{
  stopifnot(is.matrix(data_matrix))
  #stopifnot(is.integer(k) && length(k) == 1)
  stopifnot(is.numeric(k) && length(k) == 1 && k %% 1 == 0)
  #stopifnot(is.integer(centroids_per_cluster) && length(centroids_per_cluster) == 1)
  stopifnot(is.numeric(centroids_per_cluster) && length(centroids_per_cluster) == 1 && centroids_per_cluster %% 1 == 0)
  stopifnot(k>=2)
  stopifnot(centroids_per_cluster>=1)
  stopifnot(is.null(labels) || (length(labels)==nrow(data_matrix) && max(labels) == k))
  
  n <- nrow(data_matrix)
  d <- ncol(data_matrix)
  
  if(is.null(CVI_ptr))
  {
    cvi_ptr <- .CVI_create("CalinskiHarabasz", data_matrix, k)
  } else {
    cvi_ptr <- CVI_ptr
  }
  
  to_optimize <- function(centroids_coordinates_)
  {
    centroids_matrix <- matrix(centroids_coordinates_, byrow=TRUE, ncol=d)
    nearest_indices <- get.knnx(centroids_matrix, data_matrix, k=1)$nn.index[,1]
    cluster_indices <- as.integer((nearest_indices-1L)/centroids_per_cluster)+1L
    
    if(min(cluster_indices)!=1) { return(Inf) }
    if(max(cluster_indices)!=k) { return(Inf) }
    if(length(unique(cluster_indices))!=k) { return(Inf) }
    
    # NaN is possible when there is smaller number of clusters than k
    .CVI_set_labels(cvi_ptr, cluster_indices)
    w <- .CVI_compute(cvi_ptr)
    
    if(is.nan(w))
      return(Inf)
    else
      return(-w)
  }
  
  mins <- apply(data_matrix,2,min)
  maxs <- apply(data_matrix,2,max)
  
  if(is.null(labels))
  {
    centroids_coordinates = runif(k*centroids_per_cluster*d, min=mins, max=maxs)
  } else {
    kmeans_centroids <- vector("list", k)
    for(cluster_index in 1:k)
    {
      subset <- data_matrix[labels==cluster_index,]
      if(sum(labels==cluster_index) == 1)
      {
        subset <- matrix(subset, nrow=1)
      } else if(d == 1) {
        subset <- matrix(subset, ncol=1)
      }
      
      if(nrow(subset) <= centroids_per_cluster)
      {
        subset_additionalrows <- subset[sample(1:nrow(subset), centroids_per_cluster-nrow(subset), replace=TRUE), ]
        if(d == 1) {
          subset_additionalrows <- matrix(subset_additionalrows, ncol=1)
        }
        kmeans_result_centers <- rbind(subset, subset_additionalrows)
        kmeans_centroids[[cluster_index]] <- as.vector(t(kmeans_result_centers)) # by row
      } else {
        kmeans_result <- kmeans(subset, centroids_per_cluster)  
        kmeans_centroids[[cluster_index]] <- as.vector(t(kmeans_result$centers)) # by row
      }
    }
    centroids_coordinates <- unlist(kmeans_centroids)
  }
  
  res <- psoptim(par=centroids_coordinates, fn=to_optimize, lower=mins, upper=maxs, control=control)
  
  centroids_matrix <- matrix(res$par, byrow=TRUE, ncol=d)
  nearest_indices <- get.knnx(centroids_matrix, data_matrix, k=1)$nn.index[,1]
  cluster_indices <- as.integer((nearest_indices-1L)/centroids_per_cluster)+1L
  
  list(res=res, cluster_indices=cluster_indices)
}