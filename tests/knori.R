# Copyright 2016 neurodata (http://neurodata.io/)
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

require(knorR)

fn <- "../data/matrix_r50_c5_rrw.bin"
centroidfn <- "../data/init_clusters_k8_c5.bin"
k <- 8
nrow <- 50
ncol <- 5
nthread <- 2


# Data in memory, compute centroids
test.data.in.mem <- function() {
    cat("Data ==> memory, centroids ==> compute\n\n")
    #data <- replicate(ncol, rnorm(nrow))
    data <- matrix(readBin(fn, "double", (nrow*ncol), 8),
                        nrow=nrow, ncol=ncol)
    print(kmeans(data, k, nrow, ncol, nthread=nthread))
}

# Data on disk, compute centroids
test.data.ex.mem <- function() {
    cat("Data ==> disk, centroids ==> compute\n\n")
    print(kmeans(fn, k, nrow, ncol, nthread=nthread))
}

# Data on disk, centroids in memory
test.centroids.in.mem <- function() {
    cat("Data ==> disk, centroids ==> memory\n\n")
    centroids <- matrix(readBin(centroidfn, "double", (k*ncol), 8),
                        nrow=k, ncol=ncol)
    print(kmeans(fn, centroids, nrow, nthread=nthread))
}

# Data in memory, centroids in memory
test.data.centroids.in.mem <- function() {
    cat("Data ==> memory, centroids ==> memory\n\n")
    data <- matrix(readBin(fn, "double", (nrow*ncol), 8),
                        nrow=nrow, ncol=ncol)
    centroids <- matrix(readBin(centroidfn, "double", (k*ncol), 8),
                        nrow=k, ncol=ncol)

    print(kmeans(data, centroids, nthread=nthread))
}

# Data in memory, centroids on disk
test.data.in.mem.centroids.em <- function() {
    cat("Data ==> memory, centroids ==> disk\n\n")
    data <- matrix(readBin(fn, "double", (nrow*ncol), 8),
                        nrow=nrow, ncol=ncol)
    print(kmeans(data, centroidfn, nthread=nthread))
}

# Data on disk, centroids on disk
test.data.centroids.em <- function() {
    cat("Data ==> disk, centroids ==> disk\n\n")
    print(kmeans(fn, list(centroidfn, k), nrow=nrow,
                 ncol=ncol,nthread=nthread))
}

# Main
test.data.in.mem()
test.data.in.mem()
test.data.ex.mem()
test.centroids.in.mem()
test.data.centroids.in.mem()
test.data.in.mem.centroids.em()
test.data.centroids.em()
