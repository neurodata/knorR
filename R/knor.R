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
#' @param datafn Data file name on disk
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param k The desired number of clusters
#' @param max.iters Then maximum number of iteration of k-means to perform
#' @param nthread The number of parallel thread to run
#' @param centers Pre-computed centroids
#' @param init The type of initialization to use
#' @param tolerance The convergence tolerance
#' @param dist.type What dissimilarity metric to use
#' @param  centersfn Filename containing centers
#' @param  omp Use (slower) OpenMP threads rather than pthreads (default: FALSE)
#'
#' @return A list containing the output of kmeans.
#'
#' @name kmeans
#' @author Disa Mhembere <disa@@jhu.edu>
#' @rdname kmeans

kmeans <- function(datafn, nrow, ncol, k,
        max.iters=.Machine$integer.max, nthread=-1,
        centers=NULL, init=c("kmeanspp", "random", "forgy", "none"),
        tolerance=1E-6, dist.type=c("eucl", "cos"),
        centersfn="", omp=FALSE) {

	stopifnot(class(datafn) == "character")
	stopifnot(class(datafn) == "character")
	stopifnot(class(datafn) == "character")
	stopifnot(class(datafn) == "character")

	ret <- .Call("R_knor_kmeans", datafn, nrow, ncol, k,
                 max.iters, nthread, centers, init, tolerance,
                 dist.type, centersfn, omp, PACKAGE="knorR")
}
