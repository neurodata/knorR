# Copyright 2017 Neurodata (http://neurodata.io)
# Written by Disa Mhembere (disa@jhu.edu)
#
# This file is part of knor.
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

#' Perform k-means clustering on a data matrix.
#'
#' K-means provides \strong{k} disjoint sets for a dataset using a parallel and fast
#' NUMA optimized version of Lloyd's algorithm. The details of which are found
#' in this paper https://arxiv.org/pdf/1606.08905.pdf.
#'
#' @param data Data file name on disk or In memory data matrix
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param nthread The number of parallel thread to run
#' @param centers Either (i) The number of centers (i.e., k), or
#'  (ii) an In-memory data matrix, or (iii) A 2-Element \emph{list} with element 1
#'  being a filename for precomputed centers, and element 2
#'  the number of centroids.
#' @param init The type of initialization to use c("kmeanspp", "random",
#'  "forgy", "none")
#' @param tolerance The convergence tolerance
#' @param dist.type What dissimilarity metric to use
#' @param omp Use (slower) OpenMP threads rather than pthreads
#'
#' @return A list containing the attributes of the output of kmeans.
#'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
#'          which each point is allocated.
#'  centers: A matrix of cluster centres.
#'  size: The number of points in each cluster.
#'  iter: The number of (outer) iterations.
#'
#' @examples
#' iris.mat <- as.matrix(iris[,1:4])
#' k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
#' kms <- Kmeans(iris.mat, k)
#'
#' @export
#' @name Kmeans
#' @author Disa Mhembere <disa@@jhu.edu>
#' @rdname Kmeans

Kmeans <- function(data, centers, nrow=-1, ncol=-1,
                   iter.max=.Machine$integer.max, nthread=-1,
                   init=c("kmeanspp", "random", "forgy", "none"),
                   tolerance=1E-6, dist.type=c("eucl", "cos"),
                   omp=FALSE) {

    if (inherits(data, "character")) {
        if (inherits(centers, "numeric") ||
            inherits(centers, "integer")) {
            ret <- .Call("R_knor_kmeans", normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         PACKAGE="knor")
        } else if (inherits(centers, "matrix")) {
            ret <- .Call("R_knor_kmeans_centroids_im",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         PACKAGE="knor")
        } else if (inherits(centers, "list")) {
            ret <- .Call("R_knor_kmeans_data_centroids_em",
                         normalizePath(as.character(data)),
                         normalizePath(as.character(centers[1])),
                         as.integer(centers[2]),
                         as.double(nrow), as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         PACKAGE="knor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (inherits(data, "matrix")) {
        if (inherits(centers, "numeric") || inherits(centers, "integer")) {
            ret <- .Call("R_knor_kmeans_data_im", as.matrix(data),
                         as.integer(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         PACKAGE="knor")
        } else if (inherits(centers, "matrix")) {
            ret <- .Call("R_knor_kmeans_data_centroids_im", as.matrix(data),
                         as.matrix(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         PACKAGE="knor")
        } else if (inherits(centers, "character")) {
            ret <- .Call("R_knor_kmeans_data_im_centroids_em", as.matrix(data),
                         normalizePath(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         PACKAGE="knor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else {
        stop(paste("Cannot handle data of type", class(data), "\n"))
    }
}
