/*
 * Copyright 2017 neurodata (http://neurodata.io/)
 * Written by Disa Mhembere (disa@jhu.edu)
 *
 * This file is part of knorR
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY CURRENT_KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <numa.h>
#include <Rcpp.h>

#include <unordered_map>
#include "binding/knori.hpp"

RcppExport SEXP R_knor_kmeans(SEXP rdatafn, SEXP rk,
        SEXP rnrow, SEXP rncol,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance,
        SEXP rdist_type, SEXP romp) {

    std::string datafn = CHAR(STRING_ELT(rdatafn,0));
	unsigned k = INTEGER(rk)[0];
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
    bool omp = INTEGER(romp)[0];

	if (nthread == -1) {
        nthread = kpmeans::base::get_num_omp_threads();
        std::cout << "Running on " << nthread << " threads!\n";
    }

    srand(1234);

    kpmeans::base::kmeans_t kret = kpmeans::base::kmeans(datafn,
            nrow, ncol, k, max_iters, numa_num_task_nodes(), nthread, NULL,
            init, tolerance, dist_type, "", omp);

	Rcpp::List ret;
    ret["nrow"] = kret.nrow;
    ret["ncol"] = kret.ncol;
    ret["iters"] = kret.iters;
    ret["k"] = kret.k;

	Rcpp::NumericMatrix centers = Rcpp::NumericMatrix(k, ncol);
#pragma omp parallel for shared(centers)
    for (unsigned row = 0; row < k; row++)
        for (size_t col = 0; col < ncol; col++)
            centers(row, col) = kret.centroids[row*ncol + col];
	ret["centers"] = centers;

	Rcpp::IntegerVector assignment(nrow);
	for (size_t rid = 0; rid < nrow; rid++)
		assignment[rid] = kret.assignments[rid] + 1; // 1-based indexing
	ret["cluster"] = assignment;

	Rcpp::IntegerVector size(k);
	for (unsigned i = 0; i < k; i++)
		size[i] = kret.assignment_count[i];
    ret["size"] = size;

	return ret;
}
