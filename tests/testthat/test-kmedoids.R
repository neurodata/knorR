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

# Data in memory, compute centroids
test.data.in.mem <- function() {
    cat("Data ==> memory, centroids ==> compute\n\n")
    print(Kmedoids(test_data, k, nrow, ncol, nthread=nthread))
    cat("Computed successfully\n")
}

# Data on disk, compute centroids
test.data.ex.mem <- function() {
    cat("Data ==> disk, centroids ==> compute\n\n")
    print(Kmedoids(fn, k, nrow, ncol, nthread=nthread))
    cat("Computed successfully\n")
}

cat("\n\n***Running test for kmedoids***\n\n")
ret1 <- test.data.in.mem()
ret2 <- test.data.in.mem()
test_that("Data IM other compared to same", {
              expect_identical(ret1, ret2)
})

ret1 <- test.data.ex.mem()
ret2 <- test.data.ex.mem()
test_that("Data EM compared to each other", {
              expect_identical(ret1, ret2)
})
