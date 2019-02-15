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
    print(KmeansPP(test_data, k, nrow, ncol, nstart=10, nthread=nthread))
}

# Data on disk, compute centroids
test.data.ex.mem <- function() {
    cat("Data ==> disk, centroids ==> compute\n\n")
    print(KmeansPP(fn, k, nrow, ncol, nstart=10, nthread=nthread))
}

# Main
test.data.in.mem()
test.data.ex.mem()

ret1 <- test.data.in.mem()
ret2 <- test.data.in.mem()
test_that("Data in-mem compared to same", {
              expect_identical(ret1, ret2)
})

ret1 <- test.data.ex.mem()
ret2 <- test.data.ex.mem()
test_that("centroids EM compared to same", {
              expect_identical(ret1, ret2)
})
