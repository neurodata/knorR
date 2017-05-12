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

#include <unordered_map>
#include <Rcpp.h>
#include "binding/knori.hpp"

RcppExport SEXP R_kmeans(SEXP rdatafn, SEXP rnrow, SEXP rncol, SEXP rk,
        SEXP rmax_iters, SEXP rnnodes, SEXP rnthread,
        SEXP rp_centers, SEXP rinit,
        SEXP rtolerance, SEXP rdist_type,
        SEXP rcentersfn, SEXP romp) {

    std::string datafn = CHAR(STRING_ELT(rdatafn,0));
	size_t nrow = INTEGER(rnrow)[0];
	size_t ncol = INTEGER(rncol)[0];
	unsigned k = INTEGER(rk)[0];
	size_t max_iters = INTEGER(rmax_iters)[0];
	unsigned nnodes = INTEGER(rnnodes)[0];
	int nthread = INTEGER(rnthread)[0];
	double* p_centers = REAL(rp_centers);
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	std::string centersfn = CHAR(STRING_ELT(rcentersfn,0));
    bool omp = INTEGER(romp)[0];

	if (rnthread < 1) {
		fprintf(stderr, "# threads must be >= 1");
		return R_NilValue;
	}

    srand(1234);

    kpmeans::base::kmeans_t kret = kpmeans::base::kmeans(datafn,
            nrow, ncol, k, max_iters, nnodes, nthread, p_centers,
            init, tolerance, dist_type, centersfn, omp);

	Rcpp::List ret;
    ret["nrow"] = kret.nrow;
    ret["ncol"] = kret.ncol;
    ret["iters"] = kret.iters;
    ret["k"] = kret.k;

	Rcpp::NumericMatrix centers = Rcpp::NumericMatrix(k, ncol);
#pragma omp parallel for firstprivate(p_centers) shared(centers)
    for (unsigned row = 0; row < k; row++)
        for (size_t col = 0; col < ncol; col++)
            centers(row, col) = kret.centers[row*ncol + col];
	ret["centers"] = centers;

	Rcpp::IntegerVector assignment(nrow);
	for (size_t rid = 0; rid < nrow; rid++)
		assignment[rid] = kret.assignment[rid] + 1; // 1-based indexing
	ret["cluster"] = assignment;

	Rcpp::IntegerVector size(k);
	for (unsigned i = 0; i < k; i++)
		size[i] = kret.assignment_count[i];
    ret["size"] = size;

	return ret;
}
