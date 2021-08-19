## PSO algorithm
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
#' The PSO algorithm
#' 
#' @description
#' The PSO algorithm is used to maximize a given cluster validity index
#' 
#' @param data_matrix a nxk matrix representing points to cluster
#' @param k a number of clusters
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
pso <- function(data_matrix, k, labels=rep(NA, nrow(data_matrix)), control=list(maxit=10000, reltol=0.2), CVI_ptr=NULL)
{
  stopifnot(is.matrix(data_matrix))
  stopifnot(length(labels) == nrow(data_matrix))
  #stopifnot(is.integer(k) && length(k) == 1)
  stopifnot(is.numeric(k) && length(k) == 1 && k%%1==0)
  stopifnot(max(labels) == k || all(is.na(labels)))
  stopifnot(k>=2)
  
  n <- nrow(data_matrix)
  if(is.null(CVI_ptr))
  {
    cvi_ptr <- .CVI_create("CalinskiHarabasz", data_matrix, k)
  } else {
    cvi_ptr <- CVI_ptr
  }
  
  to_optimize <- function(labels_)
  {
    .CVI_set_labels(cvi_ptr, labels_)
    w <- .CVI_compute(cvi_ptr)
    
    if(is.nan(w))
      return(Inf)
    else
      return(-w)
  }
  
  res <- psoptim(par=labels, fn=to_optimize, lower=1, upper=k+0.999, control=control)
  res
}