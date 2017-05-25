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

verify.equivalent <- function(data, nclust) {
    centroids <- kmeans(iris.mat, k, iter.max=1, algorithm="Lloyd")$centers

    kms <- stats::kmeans(iris.mat, centroids, iter.max=10, algorithm="Lloyd")
    knor.kms <- knorR::kmeans(iris.mat, centroids, max.iters=10, nthread=4)

    # We shouldn't test attributes
    if (all.equal(knor.kms$centers, kms$centers, check.attributes=FALSE)) {
        cat("\nstat::kmeans == knorR::kmeans SUCCESS!\n")
        return TRUE
    } else {
        cat("[ERROR]: stat::kmeans == knorR::kmeans SUCCESS!\n")
        return FALSE
    }
}

test.iris() <- function() {
    iris.mat <- as.matrix(iris[,1:4])
    k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
    stopifnot(verify.equivalent(iris.mat, k))
}
