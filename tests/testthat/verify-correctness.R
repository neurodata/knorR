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

verify.equivalent <- function(data, k) {
    suppressWarnings(centroids <-
        stats::kmeans(data, k, iter.max=1, algorithm="Lloyd")$centers)

    kms <- stats::kmeans(data, centroids, iter.max=10, algorithm="Lloyd")
    clusternor.kms <- clusternor::Kmeans(data, centroids, iter.max=10, nthread=4)
    test_that("Verify equivalent test",{
                  expect_equivalent(clusternor.kms$centers, kms$centers)
        })
}

test.iris <- function() {
    iris.mat <- as.matrix(iris[,1:4])
    k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
    verify.equivalent(iris.mat, k)
}
