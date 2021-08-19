## DEoptim algorithm
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
#' The DEoptim algorithm
#' 
#' @description
#' The DEoptim algorithm is used to maximize a given cluster validity index
#' 
#' @param data_matrix a nxk matrix representing points to cluster
#' @param k a number of clusters
#' @param labels an optional label vector from which the algorithm will start
#' @param population_size population size, default is as in DEoptim algorithm: 10 times length of a gene vector. It is only used if labels are not NULL. To set population size with labels equal to NULL use NP parameter in control parameter.
#' @param control a list of control parameters given to PSO algorithm, most notably itermax
#' @param CVI_ptr optional cluster validity index calculator to optimize. When NULL, CalinskiHarabasz is created and used. 
#' 
#' @details
#' TBD
#' 
#' @return 
#' Returns the same value as DEoptim
#' 
#' @export
#' @family CVI
deoptim <- function(data_matrix,
                    k,
                    labels=NULL, 
                    population_size=10*length(labels),
                    control=list(itermax=1000),
                    CVI_ptr=NULL)
{
  stopifnot(is.matrix(data_matrix))
  #stopifnot(is.integer(k) && length(k) == 1)
  stopifnot(is.numeric(k) && length(k) == 1 && k%%1==0)
  stopifnot(k>=2)
  stopifnot(is.null(labels) || (length(labels)==nrow(data_matrix) && max(labels) == k))
  stopifnot(is.numeric(population_size) && length(population_size) == 1)
  
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
  if(is.null(labels))
  {
    res <- DEoptim(fn=to_optimize, lower=rep(1, n), upper=rep(k+0.999, n), control=control)
  } else {
    initialpop=rbind(randomizeMatrix(t(replicate(population_size-1, labels)),null.model = "richness",iterations = 1000), labels)
    res <- DEoptim(fn=to_optimize, lower=rep(1, n), upper=rep(k+0.999, n), control=c(control, list(initialpop=initialpop)))
  }
  
  res
}