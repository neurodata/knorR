# Copyright 2017 Neurodata (http://neurodata.io)
# Written by Disa Mhembere (disa@jhu.edu)
#
# This file is part of knorR.
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

#' A small dataset of dim: (50,5) used as for micro-benchmarks of the knorR
#' package. The data are randomly generated hence a clear number of clusters
#' will be hard to find.
#'
#' @name test_data
#' @docType data
#' @usage data(test_data)
#' @format An object of class \code{"matrix"}
#' @keywords datasets
#' @examples
#' ncenters <- 8
#' kms <- kmeans(test_data, ncenters)
"test_data"

#' A small example of centroids of dim: (8,5) used as for micro-benchmarks
#' of the knorR package. The data are randomly generated.
#'
#' @name test_centroids
#' @docType data
#' @usage data(test_centroids)
#' @format An object of class \code{"matrix"}
#' @keywords datasets
#' @examples
#' data(test_centroids)
#' kms <- kmeans(test_data, test_centroids)
"test_centroids"
