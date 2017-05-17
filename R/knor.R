# Copyright 2017 Neurodata (http://neurodata.io)
# Written by Disa Mhembere (disa@jhu.edu)
#
# This file is part of knorR.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#'  Perform k-means clustering on a data matrix.
#'
#' K-means provides `k` disjoint sets for a dataset.
#'
#' @param data Data file name on disk
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param k The desired number of clusters
#' @param max.iters Then maximum number of iteration of k-means to perform
#' @param nthread The number of parallel thread to run
#' @param centers Pre-computed centroids filename or in-memory data
#' @param init The type of initialization to use
#' @param tolerance The convergence tolerance
#' @param dist.type What dissimilarity metric to use
#' @param  omp Use (slower) OpenMP threads rather than pthreads (default: FALSE)
#'
#' @return A list containing the output of kmeans.
#'
#' @export
#' @name kmeans
#' @author Disa Mhembere <disa@@jhu.edu>
#' @rdname kmeans

kmeans <- function(data, centers, nrow=-1, ncol=-1,
                   max.iters=.Machine$integer.max, nthread=-1,
                   init=c("kmeanspp", "random", "forgy", "none"),
                   tolerance=1E-6, dist.type=c("eucl", "cos"), omp=FALSE) {

    if (class(data) == "character") {
        if (class(centers) == "numeric") {
            ret <- .Call("R_knor_kmeans", as.character(data),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol),
                         as.double(max.iters), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         PACKAGE="knorR")
        } else if (class(centers) == "matrix") {
            cat("ERROR: Precomputed centers not yet supported!\n")
            stopifnot(FALSE)
        }
    } else if (class(data) == "matrix") {
        cat("In memory Converting matrix in memory ...\n")
        stopifnot(FALSE)
    }
}
