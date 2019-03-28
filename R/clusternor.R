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
#' @param data Data file name on disk (NUMA optimized) or In memory data matrix
#' @param centers Either (i) The number of centers (i.e., k), or
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param nthread The number of parallel threads to run
#'  (ii) an In-memory data matrix, or (iii) A 2-Element \emph{list} with element 1
#'  being a filename for precomputed centers, and element 2
#'  the number of centroids.
#' @param init The type of initialization to use c("kmeanspp", "random",
#'  "forgy", "none")
#' @param tolerance The convergence tolerance
#' @param dist.type What dissimilarity metric to use
#'
#' @return A list containing the attributes of the output.
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
                   tolerance=1E-6, dist.type=c("eucl", "sqeucl", "cos", "taxi")) {

    if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_kmeans", normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol), as.double(iter.max),
                         as.integer(nthread), as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_kmeans_centroids_im",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        }
        else if (class(centers) == "list") {
            ret <- .Call("R_kmeans_data_centroids_em",
                         normalizePath(as.character(data)),
                         normalizePath(as.character(centers[1])),
                         as.integer(centers[2]),
                         as.double(nrow), as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_kmeans_data_im", as.matrix(data),
                         as.integer(centers), as.double(iter.max),
                         as.integer(nthread), as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_kmeans_data_centroids_im", as.matrix(data),
                         as.matrix(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else if (class(centers) == "character") {
            ret <- .Call("R_kmeans_data_im_centroids_em", as.matrix(data),
                         normalizePath(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else {
        stop(paste("Cannot handle data of type", class(data), "\n"))
    }
}

##' Perform k-medoids clustering on a data matrix.
##' After initialization the k-medoids algorithm partitions data by testing which
##'  data member of a cluster Ci may make a better candidate as medoid (centroid)
##'  by reducing the sum of distance (usually taxi), then running a reclustering
##'  step with updated medoids.
##'
##' @param data Data file name on disk or In memory data matrix
##' @param centers The number of centers (i.e., k)
##' @param nrow The number of samples in the dataset
##' @param ncol The number of features in the dataset
##' @param iter.max The maximum number of iteration of k-means to perform
##' @param nthread The number of parallel threads to run
##' @param init The type of initialization to use c("forgy")
##' @param tolerance The convergence tolerance
##' @param dist.type What dissimilarity metric to use
##'
##' @return A list containing the attributes of the output of kmedoids.
##'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
##'          which each point is allocated.
##'  centers: A matrix of cluster centres.
##'  size: The number of points in each cluster.
##'  iter: The number of (outer) iterations.
##'
##' @examples
##' iris.mat <- as.matrix(iris[,1:4])
##' k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
##' km <- Kmedoids(iris.mat, k)
##'
##' @export
##' @name Kmedoids
##' @author Disa Mhembere <disa@@cs.jhu.edu>
##' @rdname Kmedoids

#Kmedoids <- function(data, centers, nrow=-1, ncol=-1,
                   #iter.max=.Machine$integer.max, nthread=-1,
                   #init=c("forgy"), tolerance=1E-6,
                   #dist.type=c("taxi", "eucl", "cos")) {

    #if (class(data) == "matrix") {
        #if (class(centers) == "numeric" || class(centers) == "integer") {
            #ret <- .Call("R_kmedoids_data_im", as.matrix(data),
                         #as.integer(centers),
                         #as.double(iter.max), as.integer(nthread),
                         #as.character(init), as.double(tolerance),
                         #as.character(dist.type),
                         #PACKAGE="clusternor")
        #} else {
            #stop(paste("Cannot handle centers of type", class(centers), "\n"))
        #}
    #} else if (class(data) == "character") {
        #if (class(centers) == "numeric" || class(centers) == "integer") {
            #ret <- .Call("R_kmedoids_data_em",
                         #normalizePath(as.character(data)),
                         #as.integer(centers), as.double(nrow),
                         #as.double(ncol),
                         #as.double(iter.max), as.integer(nthread),
                         #as.character(init), as.double(tolerance),
                         #as.character(dist.type),
                         #PACKAGE="clusternor")
        #} else {
            #stop(paste("Cannot handle data of type", class(data), "\n"))
        #}
    #}
#}

#' Perform spherical k-means clustering on a data matrix.
#' Similar to the k-means algorithm differing only in that data features are
#'  min-max normalized the dissimilarity metric is Cosine distance.
#'
#' @param data Data file name on disk (NUMA optmized) or In-memory data matrix
#' @param centers Either (i) The number of centers (i.e., k), or
#'  (ii) an In-memory data matrix
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param nthread The number of parallel threads to run
#' @param init The type of initialization to use c("kmeanspp",
#'  "random", "forgy", "none")
#' @param tolerance The convergence tolerance
#'
#' @return A list containing the attributes of the output.
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
            ret <- .Call("R_skmeans_data_im", as.matrix(data),
                         as.integer(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         PACKAGE="clusternor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_skmeans_data_centroids_im", as.matrix(data),
                         as.matrix(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         PACKAGE="clusternor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_skmeans_data_em",
                         normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         PACKAGE="clusternor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_skmeans_centroids_im",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         PACKAGE="clusternor")
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
#' @param data Data file name on disk (NUMA optimized) or In memory data matrix
#' @param centers The number of centers (i.e., k)
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param nstart The number of iterations of kmeans++ to run
#' @param nthread The number of parallel threads to run
#' @param dist.type What dissimilarity metric to use c("taxi", "eucl", "cos")
#'
#' @return A list containing the attributes of the output.
#'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
#'          which each point is allocated.
#'  centers: A matrix of cluster centres.
#'  size: The number of points in each cluster.
#'  energy: The sum of distances for each sample from it's closest cluster.
#'  best.start: The sum of distances for each sample from it's closest cluster.
#'
#' @examples
#' iris.mat <- as.matrix(iris[,1:4])
#' k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
#' nstart <- 3
#' km <- KmeansPP(iris.mat, k, nstart=nstart)
#'
#' @export
#' @name KmeansPP
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname KmeansPP

KmeansPP <- function(data, centers, nrow=-1, ncol=-1,
                     nstart=1, nthread=-1,
                     dist.type=c("sqeucl", "eucl","cos", "taxi")) {
    if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_kmeanspp_data_im", as.matrix(data),
                          as.integer(centers), as.integer(nstart),
                          as.integer(nthread), as.character(dist.type),
                          PACKAGE="clusternor")
            ret$iters <- NULL
            ret
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_kmeanspp_data_em",
                         normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol), as.integer(nstart),
                         as.integer(nthread), as.character(dist.type),
                         PACKAGE="clusternor")
            ret$iters <- NULL
            ret
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else {
        stop(paste("Cannot handle data of type", class(data), "\n"))
    }
}

#' A randomized dataset sub-sample algorithm that approximates the k-means
#'  algorithm. See: https://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf
#   for details.
#'
#' @param data Data file name on disk (NUMA optimized) or In memory data matrix
#' @param centers Either (i) The number of centers (i.e., k), or
#'  (ii) an In-memory data matrix, or (iii) A 2-Element \emph{list} with element 1
#'  being a filename for precomputed centers, and element 2
#'  the number of centroids.
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param batch.size Size of the mini batches
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param nthread The number of parallel threads to run
#' @param init The type of initialization to use c("kmeanspp", "random",
#'          "forgy", "none")
#' @param tolerance The convergence tolerance
#' @param dist.type What dissimilarity metric to use
#' @param max.no.improvement Control early stopping based on the consecutive
#'      number of mini batches that does not yield an improvement on the
#'      smoothed inertia
#'
#' @return A list containing the attributes of the output.
#'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
#'          which each point is allocated.
#'  centers: A matrix of cluster centres.
#'  size: The number of points in each cluster.
#'  iter: The number of (outer) iterations.
#'
#' @examples
#' iris.mat <- as.matrix(iris[,1:4])
#' k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
#' kms <- MiniBatchKmeans(iris.mat, k, batch.size=5)
#'
#' @export
#' @name MiniBatchKmeans
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname MiniBatchKmeans

MiniBatchKmeans <- function(data, centers, nrow=-1, ncol=-1,
                            batch.size=100,
                   iter.max=.Machine$integer.max, nthread=-1,
                   init=c("kmeanspp", "random", "forgy", "none"),
                   tolerance=1E-2, dist.type=c("sqeucl", "eucl","cos", "taxi"),
                   max.no.improvement=3) {

    # TODO: Use a batch size of .2 if not provided
    if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_mbkmeans", normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol), as.integer(batch.size),
                         as.double(iter.max),
                         as.integer(nthread), as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_mbkmeans_centroids_im",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow),
                         as.integer(batch.size),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_mbkmeans_data_im", as.matrix(data),
                         as.integer(centers), as.integer(batch.size),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         as.character(dist.type),
                         PACKAGE="clusternor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_mbkmeans_data_centroids_im", as.matrix(data),
                         as.matrix(centers), as.integer(batch.size),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else {
        stop(paste("Cannot handle data of type", class(data), "\n"))
    }
}

#' Perform Fuzzy C-means clustering on a data matrix.
#' A soft variant of the kmeans algorithm where each data point are assigned a
#'  contribution weight to each cluster
#'
#' See: https://en.wikipedia.org/wiki/Fuzzy_clustering#Fuzzy_C-means_clustering
#'
#' @param data Data file name on disk (NUMA optimized) or In memory data matrix
#' @param centers Either (i) The number of centers (i.e., k), or
#'  (ii) an In-memory data matrix
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param nthread The number of parallel threads to run
#' @param fuzz.index The fuzziness coefficient/index (> 1 and < inf)
#' @param init The type of initialization to use c("forgy", "none")
#' @param tolerance The convergence tolerance
#' @param dist.type What dissimilarity metric to use
#'
#' @return A list containing the attributes of the output.
#'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
#'          which each point is allocated.
#'  centers: A matrix of cluster centres.
#'  size: The number of points in each cluster.
#'  iter: The number of (outer) iterations.
#'  contrib.mat: The data point to cluster contribution matrix
#'
#' @examples
#' iris.mat <- as.matrix(iris[,1:4])
#' k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
#' fcm <- FuzzyCMeans(iris.mat, k, iter.max=5)
#'
#' @export
#' @name FuzzyCMeans
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname FuzzyCMeans

FuzzyCMeans <- function(data, centers, nrow=-1, ncol=-1,
                   iter.max=.Machine$integer.max, nthread=-1,
                   fuzz.index=2, init=c("forgy", "none"), tolerance=1E-6,
                   dist.type=c("sqeucl", "eucl","cos", "taxi")) {

    if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_fcm_data_em", normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol), as.double(iter.max),
                         as.integer(nthread),
                         as.integer(fuzz.index), as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_fcm_data_em_centroids_im",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow), as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.integer(fuzz.index),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_fcm_data_im", as.matrix(data),
                         as.integer(centers), as.double(iter.max),
                         as.integer(nthread), as.integer(fuzz.index),
                         as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_fcm_data_centroids_im", as.matrix(data),
                         as.matrix(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.integer(fuzz.index),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="clusternor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
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
#' @param data Data file name on disk (NUMA optmized) or In memory data matrix
#' @param kmax The maximum number of centers
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param nthread The number of parallel threads to run
#' @param init The type of initialization to use c("forgy") or initial centers
#' @param tolerance The convergence tolerance for k-means at each
#'      hierarchical split
#' @param dist.type What dissimilarity metric to use
#' @param min.clust.size The minimum size of a cluster when it cannot be split
#'
#' @return A list of lists containing the attributes of the output.
#'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
#'          which each point is allocated.
#'  centers: A matrix of cluster centres.
#'  size: The number of points in each cluster.
#'  iter: The number of (outer) iterations.
#'
#' @examples
#' iris.mat <- as.matrix(iris[,1:4])
#' kmax <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
#' kms <- Hmeans(iris.mat, kmax)
#'
#' @export
#' @name Hmeans
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname Hmeans

Hmeans <- function(data, kmax, nrow=-1, ncol=-1, iter.max=20,
                   nthread=-1, init=c("forgy"), tolerance=1E-6,
                   dist.type=c("eucl", "cos", "sqeucl", "taxi"),
                   min.clust.size=1) {

    if (class(data) == "character") {
        if (class(init) == "character") {
            ret <- .Call("R_hmeans_data_em_init", as.character(data),
                         as.integer(kmax),
                         as.double(nrow), as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         as.character(dist.type), as.integer(min.clust.size),
                         PACKAGE="clusternor")
        } else if (class(init) == "matrix") {
            if (!(all(dim(init) == c(2, ncol), TRUE)))
                stop("init centers must have dim: `c(2, ncol)'")

            ret <- .Call("R_hmeans_data_em_centers", as.character(data),
                         as.integer(kmax),
                         as.double(nrow), as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.matrix(init), as.double(tolerance),
                         as.character(dist.type), as.integer(min.clust.size),
                         PACKAGE="clusternor")
        } else {
            stop(paste("Cannot handle init of type", class(init), "\n"))
        }
    } else if (class(data) == "matrix") {
        if (class(init) == "character") {

            ret <- .Call("R_hmeans_data_im_init", as.matrix(data),
                         as.integer(kmax), as.double(iter.max),
                         as.integer(nthread), as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         as.integer(min.clust.size),
                         PACKAGE="clusternor")
        } else if (class(init) == "matrix") {
            if (!(all(dim(init) == c(2, dim(data)[2]), TRUE)))
                stop("init centers must have dim: `c(2, dim(data)[1])'")

            ret <- .Call("R_hmeans_data_im_centers", as.matrix(data),
                         as.integer(kmax), as.double(iter.max),
                         as.integer(nthread), as.matrix(init),
                         as.double(tolerance), as.character(dist.type),
                         as.integer(min.clust.size),
                         PACKAGE="clusternor")
        } else {
            stop(paste("Cannot handle init of type", class(init), "\n"))
        }
    }
}

#' Perform a parallel hierarchical clustering using the x-means algorithm
#'
#' A recursive (not acutally implemented as recursion) partitioning of data into
#'  two disjoint sets at every level as described in:
#'  http://cs.uef.fi/~zhao/Courses/Clustering2012/Xmeans.pdf
#'
#' @param data Data file name on disk (NUMA optmized) or In memory data matrix
#' @param kmax The maximum number of centers
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param nthread The number of parallel threads to run
#' @param init The type of initialization to use c("forgy") or initial centers
#' @param tolerance The convergence tolerance for k-means at each hierarchical split
#' @param dist.type What dissimilarity metric to use
#' @param min.clust.size The minimum size of a cluster when it cannot be split
#'
#' @return A list of lists containing the attributes of the output.
#'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
#'          which each point is allocated.
#'  centers: A matrix of cluster centres.
#'  size: The number of points in each cluster.
#'  iter: The number of (outer) iterations.
#'
#' @examples
#' iris.mat <- as.matrix(iris[,1:4])
#' kmax <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
#' xms <- Xmeans(iris.mat, kmax)
#'
#' @export
#' @name Xmeans
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname Xmeans

Xmeans <- function(data, kmax, nrow=-1, ncol=-1, iter.max=20,
                   nthread=-1, init=c("forgy"), tolerance=1E-6,
                   dist.type=c("eucl", "cos", "taxi"), min.clust.size=1) {

    if (class(data) == "character") {
        if (class(init) == "character") {

            ret <- .Call("R_xmeans_data_em_init", as.character(data),
                         as.integer(kmax),
                         as.double(nrow), as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         as.character(dist.type), as.integer(min.clust.size),
                         PACKAGE="clusternor")
        } else if (class(init) == "matrix") {
            if (!(all(dim(init) == c(2, ncol), TRUE)))
                stop("init centers must have dim: `c(2, ncol)'")

            ret <- .Call("R_xmeans_data_em_centers", as.character(data),
                         as.integer(kmax),
                         as.double(nrow), as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.matrix(init), as.double(tolerance),
                         as.character(dist.type), as.integer(min.clust.size),
                         PACKAGE="clusternor")
        } else {
            stop(paste("Cannot handle init of type", class(init), "\n"))
        }
    } else if (class(data) == "matrix") {
        if (class(init) == "character") {

            ret <- .Call("R_xmeans_data_im_init", as.matrix(data),
                         as.integer(kmax), as.double(iter.max),
                         as.integer(nthread), as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         as.integer(min.clust.size),
                         PACKAGE="clusternor")
        } else if (class(init) == "matrix") {
            if (!(all(dim(init) == c(2, dim(data)[2]), TRUE)))
                stop("init centers must have dim: `c(2, dim(data)[1])'")

            ret <- .Call("R_xmeans_data_im_centers", as.matrix(data),
                         as.integer(kmax), as.double(iter.max),
                         as.integer(nthread), as.matrix(init),
                         as.double(tolerance), as.character(dist.type),
                         as.integer(min.clust.size),
                         PACKAGE="clusternor")
        } else {
            stop(paste("Cannot handle init of type", class(init), "\n"))
        }
    }
}

##' Perform a parallel hierarchical clustering using the g-means algorithm
##'
##' A hierarchical cluster algorithm that chooses the number of clusters based on
##'  the Anderson Darling statistic described in:
##'  http://papers.nips.cc/paper/2526-learning-the-k-in-k-means.pdf
##'
##' @param data Data file name on disk (NUMA optmized) or In memory data matrix
##' @param kmax The maximum number of centers
##' @param nrow The number of samples in the dataset
##' @param ncol The number of features in the dataset
##' @param iter.max The maximum number of iteration of k-means to perform
##' @param nthread The number of parallel threads to run
##' @param init The type of initialization to use c("forgy") or initial centers
##' @param tolerance The convergence tolerance for k-means at each
##'          hierarchical split
##' @param dist.type What dissimilarity metric to use
##' @param min.clust.size The minimum size of a cluster when it cannot be split
##' @param strictness The Anderson-Darling strictness level.
##'          Should be between 1 and 4 inclusive
##'
##' @return A list of lists containing the attributes of the output.
##'  cluster: A vector of integers (from 1:\strong{k}) indicating the cluster to
##'          which each point is allocated.
##'  centers: A matrix of cluster centres.
##'  size: The number of points in each cluster.
##'  iter: The number of (outer) iterations.
##'
##' @examples
##' iris.mat <- as.matrix(iris[,1:4])
##' kmax <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
##' xms <- Gmeans(iris.mat, kmax)
##'
##' @export
##' @name Gmeans
##' @author Disa Mhembere <disa@@cs.jhu.edu>
##' @rdname Gmeans

#Gmeans <- function(data, kmax, nrow=-1, ncol=-1, iter.max=20,
                   #nthread=-1, init=c("forgy"), tolerance=1E-6,
                   #dist.type=c("eucl", "cos", "taxi"), min.clust.size=1,
                   #strictness=4) {
    #if (strictness < 1 || strictness > 4)
        #stop("strictness must be between 1 and 4 (inclusive)")

    #if (class(data) == "character") {
        #if (class(init) == "character") {

            #ret <- .Call("R_gmeans_data_em_init", as.character(data),
                         #as.integer(kmax),
                         #as.double(nrow), as.double(ncol),
                         #as.double(iter.max), as.integer(nthread),
                         #as.character(init), as.double(tolerance),
                         #as.character(dist.type), as.integer(min.clust.size),
                         #as.integer(strictness),
                         #PACKAGE="clusternor")
        #} else if (class(init) == "matrix") {
            #if (!(all(dim(init) == c(2, ncol), TRUE)))
                #stop("init centers must have dim: `c(2, ncol)'")

            #ret <- .Call("R_gmeans_data_em_centers", as.character(data),
                         #as.integer(kmax),
                         #as.double(nrow), as.double(ncol),
                         #as.double(iter.max), as.integer(nthread),
                         #as.matrix(init), as.double(tolerance),
                         #as.character(dist.type), as.integer(min.clust.size),
                         #as.integer(strictness),
                         #PACKAGE="clusternor")
        #} else {
            #stop(paste("Cannot handle init of type", class(init), "\n"))
        #}
    #} else if (class(data) == "matrix") {
        #if (class(init) == "character") {

            #ret <- .Call("R_gmeans_data_im_init", as.matrix(data),
                         #as.integer(kmax), as.double(iter.max),
                         #as.integer(nthread), as.character(init),
                         #as.double(tolerance), as.character(dist.type),
                         #as.integer(min.clust.size), as.integer(strictness),
                         #PACKAGE="clusternor")
        #} else if (class(init) == "matrix") {
            #if (!(all(dim(init) == c(2, dim(data)[2]), TRUE)))
                #stop("init centers must have dim: `c(2, dim(data)[1])'")

            #ret <- .Call("R_gmeans_data_im_centers", as.matrix(data),
                         #as.integer(kmax), as.double(iter.max),
                         #as.integer(nthread), as.matrix(init),
                         #as.double(tolerance), as.character(dist.type),
                         #as.integer(min.clust.size), as.integer(strictness),
                         #PACKAGE="clusternor")
        #} else {
            #stop(paste("Cannot handle init of type", class(init), "\n"))
        #}
    #}
#}
