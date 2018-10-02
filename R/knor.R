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
#' @param numa.opt When passing \emph{data} as an in-memory data matrix you can
#'  optimize memory placement for Linux NUMA machines. \strong{NOTE:}
#'  performance may degrade with very large data & it requires
#'  2*memory of that without this.
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
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname Kmeans

Kmeans <- function(data, centers, nrow=-1, ncol=-1,
                   iter.max=.Machine$integer.max, nthread=-1,
                   init=c("kmeanspp", "random", "forgy", "none"),
                   tolerance=1E-6, dist.type=c("eucl", "cos", "taxi"),
                   omp=FALSE, numa.opt=FALSE) {

    if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_kmeans", normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         PACKAGE="knor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_kmeans_centroids_im",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         PACKAGE="knor")
        }
        else if (class(centers) == "list") {
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
    } else if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_kmeans_data_im", as.matrix(data),
                         as.integer(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         as.logical(numa.opt),
                         PACKAGE="knor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_kmeans_data_centroids_im", as.matrix(data),
                         as.matrix(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         as.logical(numa.opt),
                         PACKAGE="knor")
        } else if (class(centers) == "character") {
            ret <- .Call("R_knor_kmeans_data_im_centroids_em", as.matrix(data),
                         normalizePath(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         as.character(dist.type), as.logical(omp),
                         as.logical(numa.opt),
                         PACKAGE="knor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else {
        stop(paste("Cannot handle data of type", class(data), "\n"))
    }
}

#' Perform k-medoids clustering on a data matrix.
#' After initialization the k-medoids algorithm partitions data by testing which
#'  data member of a cluster Ci may make a better candidate as medoid (centroid)
#'  by reducing the sum of distance (usually taxi), then running a reclustering
#'  step with updated medoids.
#'
#' @param data Data file name on disk or In memory data matrix
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param nthread The number of parallel thread to run
#' @param centers Either (i) The number of centers (i.e., k), or
#'  (ii) an In-memory data matrix
#' @param init The type of initialization to use c("forgy", "none")
#' @param tolerance The convergence tolerance
#' @param dist.type What dissimilarity metric to use c("taxi", "eucl", "cos")
#'
#' @return A list containing the attributes of the output of kmedoids.
#'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
#'          which each point is allocated.
#'  centers: A matrix of cluster centres.
#'  size: The number of points in each cluster.
#'  iter: The number of (outer) iterations.
#'
#' @examples
#' iris.mat <- as.matrix(iris[,1:4])
#' k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
#' km <- Kmedoids(iris.mat, k)
#'
#' @export
#' @name Kmedoids
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname Kmedoids

Kmedoids <- function(data, centers, nrow=-1, ncol=-1,
                   iter.max=.Machine$integer.max, nthread=-1,
                   init=c("forgy", "none"), tolerance=1E-6,
                   dist.type=c("taxi", "eucl", "cos")) {

    if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_kmedoids_data_im", as.matrix(data),
                         as.integer(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         as.character(dist.type),
                         PACKAGE="knor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_kmedoids_data_centroids_im", as.matrix(data),
                         as.matrix(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         as.character(dist.type),
                         PACKAGE="knor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_kmedoids_data_em",
                         normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         as.character(dist.type),
                         PACKAGE="knor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_kmedoids_centroids_im(",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         as.character(dist.type),
                         PACKAGE="knor")
        }
    } else {
        stop(paste("Cannot handle data of type", class(data), "\n"))
    }
}

#' Perform spherical k-means clustering on a data matrix.
#' Similar to the k-means algorithm differing only in that data features are
#'  min-max normalized the dissimilarity metric is Cosine distance.
#'
#' @param data Data file name on disk (NUMA optmized) or In-memory data matrix
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param nthread The number of parallel thread to run
#' @param centers Either (i) The number of centers (i.e., k), or
#'  (ii) an In-memory data matrix
#' @param init The type of initialization to use c("kmeanspp",
#'  "random", "forgy", "none")
#' @param tolerance The convergence tolerance
#'
#' @return A list containing the attributes of the output of kmedoids.
#'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
#'          which each point is allocated.
#'  centers: A matrix of cluster centres.
#'  size: The number of points in each cluster.
#'  iter: The number of (outer) iterations.
#'
#' @examples
#' iris.mat <- as.matrix(iris[,1:4])
#' k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
#' km <- Skmeans(iris.mat, k)
#'
#' @export
#' @name Skmeans
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname Skmeans

Skmeans <- function(data, centers, nrow=-1, ncol=-1,
                   iter.max=.Machine$integer.max, nthread=-1,
                   init=c("kmeanspp", "random", "forgy", "none"),
                   tolerance=1E-6) {

    if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_skmeans_data_im", as.matrix(data),
                         as.integer(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         PACKAGE="knor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_skmeans_data_centroids_im", as.matrix(data),
                         as.matrix(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         PACKAGE="knor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_skmeans_data_em",
                         normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         PACKAGE="knor")
            stop(paste("Cannot handle data of type", class(data), "\n"))
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_skmeans_centroids_im(",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         PACKAGE="knor")
            stop(paste("Cannot handle data of type", class(data), "\n"))
        }
    } else {
        stop(paste("Cannot handle data of type", class(data), "\n"))
    }
}

#' Perform the k-means++ clustering algorithm on a data matrix.
#'
#' A parallel and scalable implementation of the algorithm described in
#' Ostrovsky, Rafail, et al. "The effectiveness of Lloyd-type methods for
#'  the k-means problem." Journal of the ACM (JACM) 59.6 (2012): 28.
#'
#' @param data Data file name on disk or In memory data matrix
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param centers The number of centers (i.e., k)
#' @param nstart The convergence tolerance
#' @param nthread The number of parallel thread to run
#'
#' @return A list containing the attributes of the output of kmedoids.
#'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
#'          which each point is allocated.
#'  centers: A matrix of cluster centres.
#'  size: The number of points in each cluster.
#'
#' @examples
#' iris.mat <- as.matrix(iris[,1:4])
#' k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
#' nstart <- 2
#' km <- KmeansPP(iris.mat, k, nstart)
#'
#' @export
#' @name KmeansPP
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname KmeansPP

KmeansPP <- function(data, centers, nrow=-1, ncol=-1, nstart=1, nthread=-1) {
    if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- NULL
            if (nstart == 1) {
                Kmeans(data, centers, nrow, ncol, iter.max=1, init=c("kmeanspp"))
            } else {
                ret <- Kmeans(data, centers, nrow, ncol,
                              iter.max=1, init=c("kmeanspp"))

                for (starts in 2:nstart) {
                    # Seed kmeans++
                    ret <- Kmeans(data, ret$centers, init=c("kmeanspp"))
                }
            }
            ret # return this
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            # TODO
            stop(paste("Cannot handle data of type", class(data), "\n"))
        }
    } else {
        stop(paste("Cannot handle data of type", class(data), "\n"))
    }
}

#' Perform parallel hierarchical clustering on a data matrix.
#'
#' A recursive (not acutally implemented as recursion) partitioning of data into
#'  two disjoint sets at every level as described in
#'  https://en.wikipedia.org/wiki/Hierarchical_clustering
#'
#' @param data Data file name on disk or In memory data matrix
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param level.max Maximmum
#' @param init The type of initialization to use c("kmeanspp", "random",
#'  "forgy", "none")
#' @param tolerance The convergence tolerance
#' @param dist.type What dissimilarity metric to use
#' @param algo What algorithm to use c("kmeans")
#' @param nthread The number of parallel thread to run
#'
#' @return A list of lists containing the attributes of the output of kmeans.
#'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
#'          which each point is allocated.
#'  centers: A matrix of cluster centres.
#'  size: The number of points in each cluster.
#'  iter: The number of (outer) iterations.
#'
#' @examples
#' iris.mat <- as.matrix(iris[,1:4])
#' k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
#' kms <- Hclust(iris.mat, k, level.max=2)
#'
#' @export
#' @name Hclust
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname Hclust

Hclust <- function(data, centers, nrow=-1, ncol=-1,
                   iter.max=.Machine$integer.max, nthread=-1,
                   init=c("kmeanspp", "random", "forgy", "none"),
                   tolerance=1E-6, dist.type=c("eucl", "cos", "taxi")) {

    if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            # TODO
        } else if (class(centers) == "matrix") {
            # TODO
        }
        else if (class(centers) == "list") {
            # TODO
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            # TODO
        } else if (class(centers) == "matrix") {
            # TODO
        } else if (class(centers) == "character") {
            # TODO
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else {
        stop(paste("Cannot handle data of type", class(data), "\n"))
    }
}
