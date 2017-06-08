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

#ifdef LINUX
#include <numa.h>
#endif
#include <Rcpp.h>

#include <unordered_map>
#include "binding/knori.hpp"
#include "libkcommon/io.hpp"

/**
  * Transform the C output to R
  **/
static void marshall_c_to_r(const kpmbase::kmeans_t& kret,
        Rcpp::List& ret) {

    ret["nrow"] = kret.nrow;
    ret["ncol"] = kret.ncol;
    ret["iters"] = kret.iters;
    ret["k"] = kret.k;

	Rcpp::NumericMatrix centers = Rcpp::NumericMatrix(kret.k, kret.ncol);
#pragma omp parallel for shared(centers)
    for (unsigned row = 0; row < kret.k; row++)
        for (size_t col = 0; col < kret.ncol; col++)
            centers(row, col) = kret.centroids[row*kret.ncol + col];
	ret["centers"] = centers;

	Rcpp::IntegerVector assignment(kret.nrow);
	for (size_t rid = 0; rid < kret.nrow; rid++)
		assignment[rid] = kret.assignments[rid] + 1; // 1-based indexing
	ret["cluster"] = assignment;

	Rcpp::IntegerVector size(kret.k);
	for (unsigned i = 0; i < kret.k; i++)
        size[i] = kret.assignment_count[i];
    ret["size"] = size;
}

/**
  * Data and centroids in-memory
  **/
RcppExport SEXP R_knor_kmeans_data_centroids_im(SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rtolerance,
        SEXP rdist_type, SEXP romp, SEXP rnuma_opt) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);

	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
    bool omp = INTEGER(romp)[0];
	unsigned k = centroids.nrow();
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    bool numa_opt = INTEGER(rnuma_opt)[0];

    std::vector<double> cdata(nrow*ncol);
    std::vector<double> ccentroids(k*ncol);

    if (nthread == -1)
        nthread = kpmeans::base::get_num_omp_threads();

#ifdef LINUX
    unsigned nnodes = numa_num_task_nodes();
#else
    unsigned nnodes = 1;
#endif

// TODO: Slow transpose
#pragma omp parallel for firstprivate(data) shared (cdata)
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    kpmeans::base::kmeans_t kret = kpmeans::base::kmeans(&cdata[0],
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type, omp, numa_opt);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data only in-memory
  **/
RcppExport SEXP R_knor_kmeans_data_im(SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance,
        SEXP rdist_type, SEXP romp) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned k = INTEGER(rk)[0];
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
    bool omp = INTEGER(romp)[0];
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    std::vector<double> cdata(nrow*ncol);

    if (nthread == -1)
        nthread = kpmeans::base::get_num_omp_threads();

#ifdef LINUX
    unsigned nnodes = numa_num_task_nodes();
#else
    unsigned nnodes = 1;
#endif

// TODO: Slow transpose
#pragma omp parallel for firstprivate(data) shared (cdata)
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    kpmeans::base::kmeans_t kret = kpmeans::base::kmeans(&cdata[0],
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance, dist_type, omp);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Centroids only in-memory
  */
RcppExport SEXP R_knor_kmeans_centroids_im(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rmax_iters, SEXP rnthread, SEXP rtolerance,
        SEXP rdist_type, SEXP romp) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
    bool omp = INTEGER(romp)[0];

    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);
	unsigned k = centroids.nrow();
	const size_t ncol = centroids.ncol();
    std::vector<double> ccentroids(k*ncol);

    if (nthread == -1)
        nthread = kpmeans::base::get_num_omp_threads();

#ifdef LINUX
    unsigned nnodes = numa_num_task_nodes();
#else
    unsigned nnodes = 1;
#endif

#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    kpmeans::base::kmeans_t kret = kpmeans::base::kmeans(data,
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type, omp);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data on disk, centroids computed by init method
  */
RcppExport SEXP R_knor_kmeans(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rncol,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance,
        SEXP rdist_type, SEXP romp) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	unsigned k = INTEGER(rk)[0];
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
    bool omp = INTEGER(romp)[0];

    if (nthread == -1)
        nthread = kpmeans::base::get_num_omp_threads();

#ifdef LINUX
    unsigned nnodes = numa_num_task_nodes();
#else
    unsigned nnodes = 1;
#endif

    kpmeans::base::kmeans_t kret = kpmeans::base::kmeans(data,
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance, dist_type, omp);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data in-memory and centroids on disk
  **/
RcppExport SEXP R_knor_kmeans_data_im_centroids_em(
        SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rtolerance,
        SEXP rdist_type, SEXP romp) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
    bool omp = INTEGER(romp)[0];
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    std::vector<double> cdata(nrow*ncol);

    // Figure out k and load centroids
    std::string centroidfn = CHAR(STRING_ELT(rk,0));
    std::ifstream in(centroidfn, std::ifstream::ate | std::ifstream::binary);
    size_t k = in.tellg()/(sizeof(double)*ncol);
    std::vector<double> ccentroids(k*ncol);
    kpmbase::bin_rm_reader<double> br(centroidfn);
    br.read(ccentroids);

    if (nthread == -1)
        nthread = kpmeans::base::get_num_omp_threads();

#ifdef LINUX
    unsigned nnodes = numa_num_task_nodes();
#else
    unsigned nnodes = 1;
#endif

// TODO: Slow transpose
#pragma omp parallel for firstprivate(data) shared (cdata)
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    kpmeans::base::kmeans_t kret = kpmeans::base::kmeans(&cdata[0],
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type, omp);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data on disk and centroids on disk
  **/
RcppExport SEXP R_knor_kmeans_data_centroids_em(
        SEXP rdata, SEXP rcentroids, SEXP rk,
        SEXP rnrow, SEXP rncol,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rtolerance,
        SEXP rdist_type, SEXP romp) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);

	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
    bool omp = INTEGER(romp)[0];

    std::vector<double> cdata(nrow*ncol);

    // Figure out k and load centroids
    std::string centroidfn = CHAR(STRING_ELT(rcentroids,0));
	unsigned k = INTEGER(rk)[0];

    std::vector<double> ccentroids(k*ncol);
    kpmbase::bin_rm_reader<double> br(centroidfn);
    br.read(ccentroids);

    if (nthread == -1)
        nthread = kpmeans::base::get_num_omp_threads();

#ifdef LINUX
    unsigned nnodes = numa_num_task_nodes();
#else
    unsigned nnodes = 1;
#endif

    kpmeans::base::kmeans_t kret = kpmeans::base::kmeans(data,
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type, omp);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}
