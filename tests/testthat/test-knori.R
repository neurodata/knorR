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

require(knor)
require(testthat)

fn <- "test.data.bin"
centroidfn <- "test.centroids.bin"

k <- 8
nrow <- 50
ncol <- 5
nthread <- 2

# Data in memory, compute centroids
test.data.in.mem <- function() {
    cat("Data ==> memory, centroids ==> compute\n\n")
    print(Kmeans(test_data, k, nrow, ncol, nthread=nthread))
}

# Data on disk, compute centroids
test.data.ex.mem <- function() {
    cat("Data ==> disk, centroids ==> compute\n\n")
    print(Kmeans(fn, k, nrow, ncol, nthread=nthread))
}

# Data on disk, centroids in memory
test.centroids.in.mem <- function() {
    cat("Data ==> disk, centroids ==> memory\n\n")
    print(Kmeans(fn, test_centroids, nrow, nthread=nthread))
}

# Data in memory, centroids in memory
test.data.centroids.in.mem <- function() {
    cat("Data ==> memory, centroids ==> memory\n\n")
    Kmeans(test_data, test_centroids, nthread=nthread)
}

# Data in memory, centroids on disk
test.data.in.mem.centroids.em <- function() {
    cat("Data ==> memory, centroids ==> disk\n\n")
    Kmeans(test_data, centroidfn, nthread=nthread)
}

# Data on disk, centroids on disk
test.data.centroids.em <- function() {
    cat("Data ==> disk, centroids ==> disk\n\n")
    print(Kmeans(fn, list(centroidfn, k), nrow=nrow,
                 ncol=ncol,nthread=nthread))
}

# Data in memory, centroids in memory, numa reorg
test.data.centroids.in.mem.numa.reorg <- function() {
    cat("Data ==> memory, centroids ==> memory, NUMA reorg\n\n")
    ret <- Kmeans(test_data, test_centroids, nthread=nthread, numa.opt=TRUE)
}

# Data in memory, centroids on disk, numa reorg
test.data.in.mem.centroids.em.numa.reorg <- function() {
    cat("Data ==> memory, centroids ==> disk, NUMA reorg\n\n")
    Kmeans(test_data, centroidfn, nthread=nthread, numa.opt=TRUE)
}

# Main
test.data.in.mem()
test.data.ex.mem()

test.centroids.in.mem()
test.data.centroids.em()

ret1 <- test.data.centroids.in.mem()
ret2 <- test.data.centroids.in.mem.numa.reorg()
test_that("Data in-mem compared to numa reorg", {
              expect_identical(ret1, ret2)
})

ret1 <- test.data.in.mem.centroids.em()
ret2 <- test.data.in.mem.centroids.em.numa.reorg()
test_that("centroids EM compared to numa reorg", {
              expect_identical(ret1, ret2)
})

source("verify-correctness.R")
test.iris()
