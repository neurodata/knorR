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

require(clusternor)
require(testthat)

fn <- "test.data.bin"
centroidfn <- "test.centroids.bin"

k <- 8
nrow <- 50
ncol <- 5
nthread <- 2
mb.size <- 20
iter.max <- 100

# Data in memory, compute centroids
test.data.in.mem <- function() {
    cat("Data ==> memory, centroids ==> compute\n\n")
    print(MiniBatchKmeans(test_data, k, nrow, ncol, mb.size,
                          iter.max=iter.max, nthread=nthread, init="kmeanspp"))
}

# Data on disk, compute centroids
test.data.ex.mem <- function() {
    cat("Data ==> disk, centroids ==> compute\n\n")
    print(MiniBatchKmeans(fn, k, nrow, ncol, mb.size, nthread=nthread,
                          iter.max=iter.max, init="random"))
}

# Data on disk, centroids in memory
test.centroids.in.mem <- function() {
    cat("Data ==> disk, centroids ==> memory\n\n")
    print(MiniBatchKmeans(fn, test_centroids, nrow,
                          iter.max=iter.max, batch.size, nthread=nthread))
}

# Data in memory, centroids in memory
test.data.centroids.in.mem <- function() {
    cat("Data ==> memory, centroids ==> memory\n\n")
    MiniBatchKmeans(test_data, test_centroids, mb.size,
                    iter.max=iter.max, nthread=nthread)
}

# Main

cat("\n\n***Running test for kmeans***\n\n")
test.centroids.in.mem()
ret1 <- test.data.centroids.in.mem()
ret2 <- test.data.centroids.in.mem()
test_that("Data in-mem compared to same", {
              expect_identical(ret1, ret2)
})
ret1 <- test.data.in.mem()
ret2 <- test.data.in.mem()
test_that("data IM compared to same", {
              expect_identical(ret1, ret2)
})

ret1 <- test.data.ex.mem()
ret2 <- test.data.ex.mem()
test_that("data EM compared to same", {
              expect_identical(ret1, ret2)
})
