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
                   tolerance=1E-6, dist.type=c("eucl", "cos", "taxi")) {

    if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_kmeans", normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol), as.double(iter.max),
                         as.integer(nthread), as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_kmeans_centroids_im",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        }
        else if (class(centers) == "list") {
            ret <- .Call("R_knor_kmeans_data_centroids_em",
                         normalizePath(as.character(data)),
                         normalizePath(as.character(centers[1])),
                         as.integer(centers[2]),
                         as.double(nrow), as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_kmeans_data_im", as.matrix(data),
                         as.integer(centers), as.double(iter.max),
                         as.integer(nthread), as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_kmeans_data_centroids_im", as.matrix(data),
                         as.matrix(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else if (class(centers) == "character") {
            ret <- .Call("R_knor_kmeans_data_im_centroids_em", as.matrix(data),
                         normalizePath(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
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
                         as.double(tolerance),
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
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_skmeans_centroids_im(",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance),
                         PACKAGE="knor")
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
#' @param nthread The number of parallel thread to run
#' @param dist.type What dissimilarity metric to use c("taxi", "eucl", "cos")
#'
#' @return A list containing the attributes of the output of kmedoids.
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
                     nstart=1, nthread=-1, dist.type=c("taxi", "eucl", "cos")) {
    if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_kmeanspp_data_im", as.matrix(data),
                          as.integer(centers), as.integer(nstart),
                          as.integer(nthread), as.character(dist.type),
                          PACKAGE="knor")
            ret$iters <- NULL
            ret
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_kmeanspp_data_em",
                         normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol), as.integer(nstart),
                         as.integer(nthread), as.character(dist.type),
                         PACKAGE="knor")
            ret$iters <- NULL
            ret
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
#' @param centers The number of centers (i.e., k)
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param dist.type What dissimilarity metric to use
#' @param algo What splitting algorithm to use c("kmeans")
#' @param min.clust.size The minimum size of a cluster when it cannot be split
#' @param init The type of initialization to use c("kmeanspp", "random",
#'  "forgy", "none")
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
#' kms <- Hclust(iris.mat, k)
#'
#' @export
#' @name Hclust
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname Hclust

Hclust <- function(data, centers, nrow=-1, ncol=-1,
                   dist.type=c("eucl", "cos", "taxi"), algo=c("kmeans"),
                   min.clust.size=1,
                   init=c("kmeanspp", "random", "forgy", "none"), nthread=-1) {

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
#' @param nthread The number of parallel thread to run
#' @param init The type of initialization to use c("kmeanspp", "random",
#'  "forgy", "none")
#' @param tolerance The convergence tolerance
#' @param dist.type What dissimilarity metric to use
#' @param max.no.improvement Control early stopping based on the consecutive
#'      number of mini batches that does not yield an improvement on the
#'      smoothed inertia
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
                   tolerance=1E-2, dist.type=c("eucl", "cos", "taxi"),
                   max.no.improvement=3) {

    # TODO: Use a batch size of .2 if not provided
    if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_mbkmeans", normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol), as.integer(batch.size),
                         as.double(iter.max),
                         as.integer(nthread), as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_mbkmeans_centroids_im",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow),
                         as.integer(batch.size),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_mbkmeans_data_im", as.matrix(data),
                         as.integer(centers), as.integer(batch.size),
                         as.double(iter.max), as.integer(nthread),
                         as.character(init), as.double(tolerance),
                         as.character(dist.type),
                         PACKAGE="knor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_mbkmeans_data_centroids_im", as.matrix(data),
                         as.matrix(centers), as.integer(batch.size),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else if (class(centers) == "character") {
            ret <- .Call("R_knor_mbkmeans_data_im_centroids_em", as.matrix(data),
                         normalizePath(centers), as.integer(batch.size),
                         as.double(iter.max), as.integer(nthread),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else {
        stop(paste("Cannot handle data of type", class(data), "\n"))
    }
}

#' Fit a mixture of Gaussians distribution to the data.
#'
#'
#' @param data Data file name on disk (NUMA optimized) or In memory data matrix
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param k The number of gaussians to estimate
#' @param reg.covar The covariance matrix regularization term
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param nthread The number of parallel thread to run
#' @param tolerance The convergence tolerance
#' @param cov.type is the type of covariance matrix. It can be
#'       one of {"full", "tied", "diag", "spherical"}.
#'       \itemize{
#'       \item{"full"}{each component has its own general covariance matrix.}
#'       \item{"tied"}{all components share the same general covariance matrix.}
#'       \item{"diag"}{each component has its own diagonal covariance matrix.}
#'       \item{"spherical"}{each component has its own single variance.}
#'       }
#'
#' @return A list containing the attributes of the output of GMM.
#'       \itemize{
#'        \item{loglik}{a n x k matrix, whose \code{[i, k]}th entry is
#'                the conditional probability of the ith observation
#'                belonging to the kth component of the mixture.}
#'        \item{iter}{the number of iterations}
#'        \item{parameters}{parameters of the mixture of Gaussian distribution.
#'             \itemize{
#'             \item{weights}{a vector of k elements. Each element is
#'              the weight of the Gaussian distribution in the mixture.}
#'             \item{means}{a k x d matrix. Each row is the mean of a
#'               Gaussian distribution.}
#'             \item{covs}{a list of matrices, a matrix or a vector,
#'               depending on \code{cov.type}}}
#'        }
#'      }
#'
#' @examples
#' iris.mat <- as.matrix(iris[,1:4])
#' k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
#' kms <- GMM(iris.mat, k, cov=c("full"))
#'
#' @export
#' @name GMM
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname  GMM

GMM <- function(data, nrow=-1, ncol=-1, k, reg.covar=1E-6,
                iter.max=100, nthread=-1, tolerance=1E-6,
                cov.type=c("full", "tied", "diag", "spherical")) {

}

#' Perform Fuzzy C-means clustering on a data matrix.
#' A soft variant of the kmeans algorithm where each data point are assigned a
#'  contribution weight to each cluster
#'
#' See: https://en.wikipedia.org/wiki/Fuzzy_clustering#Fuzzy_C-means_clustering
#'
#' @param data Data file name on disk (NUMA optimized) or In memory data matrix
#' @param nrow The number of samples in the dataset
#' @param ncol The number of features in the dataset
#' @param iter.max The maximum number of iteration of k-means to perform
#' @param nthread The number of parallel thread to run
#' @param centers Either (i) The number of centers (i.e., k), or
#'  (ii) an In-memory data matrix
#' @param init The type of initialization to use c("kmeanspp", "random",
#'  "forgy", "none")
#' @param fuzziness.index The fuzziness coefficient/index (> 1 and < inf)
#' @param tolerance The convergence tolerance
#' @param dist.type What dissimilarity metric to use
#'
#' @return A list containing the attributes of the output of kmeans.
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
#' fcm <- FuzzyCMeans(iris.mat, k)
#'
#' @export
#' @name FuzzyCMeans
#' @author Disa Mhembere <disa@@cs.jhu.edu>
#' @rdname FuzzyCMeans

FuzzyCMeans <- function(data, centers, nrow=-1, ncol=-1,
                   iter.max=.Machine$integer.max, nthread=-1,
                   fuzz.index=2, init=c("forgy", "kmeanspp", "none"),
                   tolerance=1E-6, dist.type=c("eucl", "cos", "taxi")) {

    if (class(data) == "character") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_fcm_data_em", normalizePath(as.character(data)),
                         as.integer(centers), as.double(nrow),
                         as.double(ncol), as.double(iter.max),
                         as.integer(nthread),
                         as.integer(fuzz.index), as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_fcm_data_em_centroids_im",
                         normalizePath(as.character(data)),
                         as.matrix(centers), as.double(nrow), as.double(ncol),
                         as.double(iter.max), as.integer(nthread),
                         as.integer(fuzz.index),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else if (class(data) == "matrix") {
        if (class(centers) == "numeric" || class(centers) == "integer") {
            ret <- .Call("R_knor_fcm_data_im", as.matrix(data),
                         as.integer(centers), as.double(iter.max),
                         as.integer(nthread), as.integer(fuzz.index),
                         as.character(init),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else if (class(centers) == "matrix") {
            ret <- .Call("R_knor_fcm_data_centroids_im", as.matrix(data),
                         as.matrix(centers),
                         as.double(iter.max), as.integer(nthread),
                         as.integer(fuzz.index),
                         as.double(tolerance), as.character(dist.type),
                         PACKAGE="knor")
        } else {
            stop(paste("Cannot handle centers of type", class(centers), "\n"))
        }
    } else {
        stop(paste("Cannot handle data of type", class(data), "\n"))
    }
}
