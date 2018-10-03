/*
 * Copyright 2017 neurodata (http://neurodata.io/)
 * Written by Disa Mhembere (disa@jhu.edu)
 *
 * This file is part of knor
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

#include <Rcpp.h>

#include <unordered_map>
#include "cknor/binding/knori.hpp"
#include "cknor/libkcommon/io.hpp"
#include "cknor/binding/kmedoids.hpp"
#include "cknor/binding/skmeans.hpp"

/**
  * Transform the C output to R
  **/

static void marshall_c_to_r(const knor::base::cluster_t& kret,
        Rcpp::List& ret) {
    ret["nrow"] = kret.nrow;
    ret["ncol"] = kret.ncol;
    ret["iters"] = kret.iters;
    ret["k"] = kret.k;

	Rcpp::NumericMatrix centers = Rcpp::NumericMatrix(kret.k, kret.ncol);

#ifdef _OPENMP
#pragma omp parallel for shared(centers)
#endif
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
        SEXP rmax_iters, SEXP rnthread, SEXP rtolerance, SEXP rdist_type) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);

	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned k = centroids.nrow();
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();

    std::vector<double> cdata(nrow*ncol);
    std::vector<double> ccentroids(k*ncol);

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

// TODO: Slow transpose
#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    knor::base::cluster_t kret = knor::base::kmeans(&cdata[0],
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type);

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
        SEXP rdist_type) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned k = INTEGER(rk)[0];
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    std::vector<double> cdata(nrow*ncol);

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

// TODO: Slow transpose
#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    knor::base::cluster_t kret = knor::base::kmeans(&cdata[0],
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance, dist_type);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Centroids only in-memory
  */
RcppExport SEXP R_knor_kmeans_centroids_im(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rmax_iters, SEXP rnthread, SEXP rtolerance,
        SEXP rdist_type) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));

    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);
	unsigned k = centroids.nrow();
	const size_t ncol = centroids.ncol();
    std::vector<double> ccentroids(k*ncol);

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    knor::base::cluster_t kret = knor::base::kmeans(data,
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type);

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
        SEXP rdist_type) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	unsigned k = INTEGER(rk)[0];
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

    knor::base::cluster_t kret = knor::base::kmeans(data,
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance, dist_type);

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
        SEXP rdist_type) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    std::vector<double> cdata(nrow*ncol);

    // Figure out k and load centroids
    std::string centroidfn = CHAR(STRING_ELT(rk,0));
    std::ifstream in(centroidfn, std::ifstream::ate | std::ifstream::binary);
    size_t k = in.tellg()/(sizeof(double)*ncol);
    std::vector<double> ccentroids(k*ncol);
    knor::base::bin_rm_reader<double> br(centroidfn);
    br.read(ccentroids);

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

// TODO: Slow transpose
#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    knor::base::cluster_t kret = knor::base::kmeans(&cdata[0],
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type);

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
        SEXP rtolerance, SEXP rdist_type) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);

	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));

    std::vector<double> cdata(nrow*ncol);

    // Figure out k and load centroids
    std::string centroidfn = CHAR(STRING_ELT(rcentroids,0));
	unsigned k = INTEGER(rk)[0];

    std::vector<double> ccentroids(k*ncol);
    knor::base::bin_rm_reader<double> br(centroidfn);
    br.read(ccentroids);

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

    knor::base::cluster_t kret = knor::base::kmeans(data,
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

////////////////////////////////// KMEDOIDS ////////////////////////////////////
/**
  * Data in memory
**/
RcppExport SEXP R_knor_kmedoids_data_im(SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance,
        SEXP rdist_type) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned k = INTEGER(rk)[0];
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    std::vector<double> cdata(nrow*ncol);

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    knor::base::cluster_t kret = knor::base::kmedoids(&cdata[0],
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance, dist_type);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data and centroids in-memory
**/
RcppExport SEXP R_knor_kmedoids_data_centroids_im(SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rtolerance, SEXP rdist_type) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);

	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned k = centroids.nrow();
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();

    std::vector<double> cdata(nrow*ncol);
    std::vector<double> ccentroids(k*ncol);

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    knor::base::cluster_t kret = knor::base::kmedoids(&cdata[0],
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Centroids only in-memory
  */
RcppExport SEXP R_knor_kmedoids_centroids_im(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rmax_iters, SEXP rnthread, SEXP rtolerance,
        SEXP rdist_type) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));

    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);
	unsigned k = centroids.nrow();
	const size_t ncol = centroids.ncol();
    std::vector<double> ccentroids(k*ncol);

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);
knor::base::cluster_t kret = knor::base::kmedoids(data, nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data on disk, centroids computed by init method
  */
RcppExport SEXP R_knor_kmedoids_data_em(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rncol,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance,
        SEXP rdist_type) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	unsigned k = INTEGER(rk)[0];
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

    knor::base::cluster_t kret = knor::base::kmedoids(data,
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance, dist_type);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}
//////////////////////////////// END KMEDOIDS //////////////////////////////////

////////////////////////////////// SKMEANS /////////////////////////////////////
/**
  * Data in memory
**/
RcppExport SEXP R_knor_skmeans_data_im(SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned k = INTEGER(rk)[0];
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    std::vector<double> cdata(nrow*ncol);

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    knor::base::cluster_t kret = knor::base::skmeans(&cdata[0],
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data and centroids in-memory
  **/
RcppExport SEXP R_knor_skmeans_data_centroids_im(SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rtolerance) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);

	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	unsigned k = centroids.nrow();
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();

    std::vector<double> cdata(nrow*ncol);
    std::vector<double> ccentroids(k*ncol);

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    knor::base::cluster_t kret = knor::base::skmeans(&cdata[0],
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data on disk, centroids computed by init method
  */
RcppExport SEXP R_knor_skmeans_data_em(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rncol,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	unsigned k = INTEGER(rk)[0];
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

    knor::base::cluster_t kret = knor::base::skmeans(data,
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Centroids only in-memory
  */
RcppExport SEXP R_knor_skmeans_centroids_im(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rmax_iters, SEXP rnthread, SEXP rtolerance) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];

    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);
	unsigned k = centroids.nrow();
	const size_t ncol = centroids.ncol();
    std::vector<double> ccentroids(k*ncol);

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    unsigned nnodes = knor::base::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    knor::base::cluster_t kret =
        knor::base::skmeans(data, nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}
/////////////////////////////// END SKMEANS ////////////////////////////////////

///////////////////////////////// KMEANSPP /////////////////////////////////////
/**
  * Data on disk
  */
RcppExport SEXP R_knor_kmeanspp_data_em(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rncol, SEXP rnstart, SEXP rnthread, SEXP rdist_type) {

    std::string data = CHAR(STRING_ELT(rdata,0));
    unsigned k = INTEGER(rk)[0];
    size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
    size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
    unsigned nstart = INTEGER(rnstart)[0];
    int nthread = INTEGER(rnthread)[0];
    std::string dist_type = CHAR(STRING_ELT(rdist_type,0));

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    std::pair<std::pair<unsigned, double>, knor::base::cluster_t> kret =
        knor::base::kmeansPP(data,
                nrow, ncol, k, nstart, nthread, dist_type);

    Rcpp::List ret;
    marshall_c_to_r(kret.second, ret);
    ret["best.start"] = kret.first.first;
    ret["energy"] = kret.first.second;
    return ret;
}

/**
  * Data in memory
  */
RcppExport SEXP R_knor_kmeanspp_data_im(SEXP rdata, SEXP rk,
        SEXP rnstart, SEXP rnthread, SEXP rdist_type) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned k = INTEGER(rk)[0];
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
	const unsigned nstart = INTEGER(rnstart)[0];
	int nthread = INTEGER(rnthread)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));

    if (nthread == -1)
        nthread = knor::base::get_num_omp_threads();

    // Convert R matrix to vector
    std::vector<double> cdata(nrow*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    std::pair<std::pair<unsigned, double>, knor::base::cluster_t> kret =
        knor::base::kmeansPP(
                &cdata[0], nrow, ncol, k, nstart, nthread, dist_type);

	Rcpp::List ret;
    marshall_c_to_r(kret.second, ret);
    ret["best.start"] = kret.first.first;
    ret["energy"] = kret.first.second;
	return ret;
}
/////////////////////////////// END KMEANSPP ///////////////////////////////////
