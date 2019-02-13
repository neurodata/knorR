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
#include "cknor/libkcommon/io.hpp"
#include "cknor/libkcommon/types.hpp"
#include "cknor/binding/kmeanspp.hpp"

#include "cknor/libman/kmeans_task_coordinator.hpp"
#include "cknor/libman/kmeans_coordinator.hpp"
#include "cknor/libman/skmeans_coordinator.hpp"
#include "cknor/libman/fcm_coordinator.hpp"
#include "cknor/libman/medoid_coordinator.hpp"
#include "cknor/libman/hclust_coordinator.hpp"
#include "cknor/libman/xmeans_coordinator.hpp"
#include "cknor/libman/gmeans_coordinator.hpp"

/**
  * Transform the C output to R
  **/
namespace kprune = knor::prune;
namespace kbase = knor::base;

static void marshall_c_to_r(const kbase::cluster_t& kret,
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
  * BEGIN ALGORITHMS
  */

//////////////////////////////////// KMEANS ////////////////////////////////////
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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

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

    kbase::cluster_t kret = kprune::kmeans_task_coordinator::create(
            "", nrow, ncol, k, max_iters, nnodes,
            nthread, &ccentroids[0], "none", tolerance,
            dist_type)->run(&cdata[0]);

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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

// TODO: Slow transpose
#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    kbase::cluster_t kret = kprune::kmeans_task_coordinator::create(
            "", nrow, ncol, k, max_iters, nnodes,
            nthread, NULL, init, tolerance, dist_type)->run(&cdata[0]);

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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    kbase::cluster_t kret = kprune::kmeans_task_coordinator::create(
            data, nrow, ncol, k, max_iters, nnodes,
            nthread, &ccentroids[0], "none", tolerance, dist_type)->run();

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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

    kbase::cluster_t kret = kprune::kmeans_task_coordinator::create(
            data, nrow, ncol, k, max_iters, nnodes,
            nthread, NULL, init, tolerance, dist_type)->run();

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
    kbase::bin_rm_reader<double> br(centroidfn);
    br.read(ccentroids);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    kbase::cluster_t kret = kprune::kmeans_task_coordinator::create("",
            nrow, ncol, k, max_iters, nnodes, nthread, &ccentroids[0],
            "none", tolerance, dist_type)->run(&cdata[0]);

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
    kbase::bin_rm_reader<double> br(centroidfn);
    br.read(ccentroids);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

    kbase::cluster_t kret = kprune::kmeans_task_coordinator::create(data,
            nrow, ncol, k, max_iters, nnodes,
            nthread, &ccentroids[0], "none", tolerance, dist_type)->run();

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}
//////////////////////////////// END KMEANS ////////////////////////////////////

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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    auto kret = knor::medoid_coordinator::create("",
                nrow, ncol, k, max_iters, nnodes, nthread, NULL,
                init, tolerance, dist_type)->run(&cdata[0]);

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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

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

    kbase::cluster_t kret = knor::medoid_coordinator::create("",
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type)->run(&cdata[0]);

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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    kbase::cluster_t kret = knor::medoid_coordinator::create(
            data, nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance, dist_type)->run();

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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

    kbase::cluster_t kret = knor::medoid_coordinator::create(data,
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance, dist_type)->run();

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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    kbase::cluster_t kret = knor::skmeans_coordinator::create("",
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance)->run(&cdata[0]);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data and centroids in-memory
  **/
RcppExport SEXP R_knor_skmeans_data_centroids_im(SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread, SEXP rtolerance) {

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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

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

    kbase::cluster_t kret = knor::skmeans_coordinator::create("",
            nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance)->run(&cdata[0]);

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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

    kbase::cluster_t kret = knor::skmeans_coordinator::create(data,
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance)->run();

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
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    kbase::cluster_t kret = knor::skmeans_coordinator::create(
            data, nrow, ncol, k, max_iters, nnodes, nthread,
            &ccentroids[0], "none", tolerance)->run();

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
        nthread = kbase::get_num_omp_threads();

    std::pair<std::pair<unsigned, double>, kbase::cluster_t> kret =
        kbase::kmeansPP(data,
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
        nthread = kbase::get_num_omp_threads();

    // Convert R matrix to vector
    std::vector<double> cdata(nrow*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    std::pair<std::pair<unsigned, double>, kbase::cluster_t> kret =
        kbase::kmeansPP(
                &cdata[0], nrow, ncol, k, nstart, nthread, dist_type);

	Rcpp::List ret;
    marshall_c_to_r(kret.second, ret);
    ret["best.start"] = kret.first.first;
    ret["energy"] = kret.first.second;
	return ret;
}
/////////////////////////////// END KMEANSPP ///////////////////////////////////

/////////////////////////////// MINIBATCH KMEANS ///////////////////////////////
/**
  * Data on disk, centroids computed by init method
  */
RcppExport SEXP R_knor_mbkmeans(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rncol, SEXP rmb_size,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance,
        SEXP rdist_type) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	unsigned k = INTEGER(rk)[0];
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	unsigned mb_size = INTEGER(rmb_size)[0];
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

    knor::coordinator::ptr coord = kprune::kmeans_task_coordinator::create(
            data, nrow, ncol, k, max_iters, nnodes,
            nthread, NULL, init, tolerance, dist_type);
    auto kc = std::static_pointer_cast<kprune::kmeans_task_coordinator>(coord);
    kc->set_mini_batch_size(mb_size);
    kbase::cluster_t kret = kc->mb_run();

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Centroids only in-memory
  */
RcppExport SEXP R_knor_mbkmeans_centroids_im(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rmb_size, SEXP rmax_iters, SEXP rnthread,
        SEXP rtolerance, SEXP rdist_type) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned mb_size = INTEGER(rmb_size)[0];

    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);
	unsigned k = centroids.nrow();
	const size_t ncol = centroids.ncol();
    std::vector<double> ccentroids(k*ncol);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    knor::coordinator::ptr coord = kprune::kmeans_task_coordinator::create(
            data, nrow, ncol, k, max_iters, nnodes,
            nthread, &ccentroids[0], "none", tolerance, dist_type);
    auto kc = std::static_pointer_cast<kprune::kmeans_task_coordinator>(coord);
    kc->set_mini_batch_size(mb_size);
    kbase::cluster_t kret = kc->mb_run();

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

RcppExport SEXP R_knor_mbkmeans_data_centroids_im(SEXP rdata, SEXP rk,
        SEXP rmb_size, SEXP rmax_iters, SEXP rnthread, SEXP rtolerance,
        SEXP rdist_type) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);

	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned k = centroids.nrow();
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
	unsigned mb_size = INTEGER(rmb_size)[0];

    std::vector<double> cdata(nrow*ncol);
    std::vector<double> ccentroids(k*ncol);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

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

    knor::coordinator::ptr coord = kprune::kmeans_task_coordinator::create(
            "", nrow, ncol, k, max_iters, nnodes,
            nthread, &ccentroids[0], "none", tolerance, dist_type);
    auto kc = std::static_pointer_cast<kprune::kmeans_task_coordinator>(coord);
    kc->set_mini_batch_size(mb_size);

    kbase::cluster_t kret = kc->mb_run(&cdata[0]);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data only in-memory
  **/
RcppExport SEXP R_knor_mbkmeans_data_im(SEXP rdata, SEXP rk,
        SEXP rmb_size, SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance, SEXP rdist_type) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned k = INTEGER(rk)[0];
	unsigned mb_size = INTEGER(rmb_size)[0];
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    std::vector<double> cdata(nrow*ncol);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    knor::coordinator::ptr coord = kprune::kmeans_task_coordinator::create(
            "", nrow, ncol, k, max_iters, nnodes,
            nthread, NULL, init, tolerance, dist_type);
    auto kc = std::static_pointer_cast<kprune::kmeans_task_coordinator>(coord);
    kc->set_mini_batch_size(mb_size);

    kbase::cluster_t kret = kc->mb_run(&cdata[0]);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}
///////////////////////////// END MINIBATCH KMEANS /////////////////////////////

////////////////////////////// FUZZY C-MEANS ///////////////////////////////////
/**
  * Data in memory
**/
RcppExport SEXP R_knor_fcm_data_im(SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread, SEXP rfuzzindex,
        SEXP rinit, SEXP rtolerance, SEXP rdist_type) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned k = INTEGER(rk)[0];
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	const unsigned fuzzindex = INTEGER(rfuzzindex)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    std::vector<double> cdata(nrow*ncol);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

     auto kret  = knor::fcm_coordinator::create("",
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance, dist_type, fuzzindex)->run(&cdata[0]);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data and centroids in-memory
  **/
RcppExport SEXP R_knor_fcm_data_centroids_im(SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread, SEXP rfuzzindex,
        SEXP rtolerance, SEXP rdist_type) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);

	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	unsigned k = centroids.nrow();
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
	const unsigned fuzzindex = INTEGER(rfuzzindex)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));

    std::vector<double> cdata(nrow*ncol);
    std::vector<double> ccentroids(k*ncol);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

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

     auto kret = knor::fcm_coordinator::create("",
            nrow, ncol, k, max_iters, nnodes, nthread, &ccentroids[0],
            "none", tolerance, dist_type, fuzzindex)->run(&cdata[0]);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data on disk, centroids computed by init method
  */
RcppExport SEXP R_knor_fcm_data_em(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rncol,
        SEXP rmax_iters, SEXP rnthread, SEXP rfuzzindex,
        SEXP rinit, SEXP rtolerance, SEXP rdist_type) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	unsigned k = INTEGER(rk)[0];
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	const unsigned fuzzindex = INTEGER(rfuzzindex)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

    auto kret = knor::fcm_coordinator::create(data,
            nrow, ncol, k, max_iters, nnodes, nthread, NULL,
            init, tolerance, dist_type, fuzzindex)->run();

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data on disk, centroids in memory
  */
RcppExport SEXP R_knor_fcm_data_em_centroids_im(SEXP rdata, SEXP rk,
        SEXP rnrow, SEXP rncol,
        SEXP rmax_iters, SEXP rnthread, SEXP rfuzzindex,
        SEXP rtolerance, SEXP rdist_type) {

    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk);
    std::string data = CHAR(STRING_ELT(rdata,0));
	unsigned k = centroids.nrow();
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	const unsigned fuzzindex = INTEGER(rfuzzindex)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

    std::vector<double> ccentroids(k*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    auto kret = knor::fcm_coordinator::create(data,
            nrow, ncol, k, max_iters, nnodes, nthread, &ccentroids[0],
            "none", tolerance, dist_type, fuzzindex)->run();

	Rcpp::List ret; marshall_c_to_r(kret, ret);
	return ret;
}
////////////////////////////// END FUZZY C-MEANS ///////////////////////////////

////////////////////////////////// HMEANS //////////////////////////////////////

/**
  * Data on disk provided k
  */
RcppExport SEXP R_knor_hmeans_data_em_k(SEXP rdata, SEXP rnrow, SEXP rncol,
        SEXP rk, SEXP rmax_iters, SEXP rnthread, SEXP rinit,
        SEXP rtolerance, SEXP rdist_type, SEXP rmin_clust_size) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	unsigned k = INTEGER(rk)[0]; // FIXME: kmax
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();
    unsigned nnodes = kbase::get_num_nodes();

    kbase::cluster_t kret =
                knor::hclust_coordinator::create(data, nrow, ncol, k,
                max_iters, nnodes, nthread, NULL,
                init, tolerance, dist_type, min_clust_size)->run();

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data on disk centers provided
  */
RcppExport SEXP R_knor_hmeans_data_em_centers(SEXP rdata, SEXP rnrow,
        SEXP rncol, SEXP rk, SEXP rmax_iters, SEXP rnthread,
        SEXP rtolerance, SEXP rdist_type, SEXP rmin_clust_size) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);

    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk); // FIXME
	unsigned k = INTEGER(rk)[0]; // FIXME: kmax

	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];

    std::vector<double> ccentroids(centroids.nrow()*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < centroids.nrow(); row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();
    unsigned nnodes = kbase::get_num_nodes();

    kbase::cluster_t kret =
        knor::hclust_coordinator::create(data, nrow, ncol, k,
                max_iters, nnodes, nthread, &ccentroids[0], "none",
                tolerance, dist_type, min_clust_size)->run();

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data in memory and provided k
**/
RcppExport SEXP R_knor_hmeans_data_im_k(SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance,
        SEXP rdist_type, SEXP rmin_clust_size) {
    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned k = INTEGER(rk)[0]; // FIXME: kmax
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];

	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    std::vector<double> cdata(nrow*ncol);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    kbase::cluster_t kret =
        knor::hclust_coordinator::create("", nrow, ncol, k,
                max_iters, nnodes, nthread, NULL, init, tolerance,
                dist_type, min_clust_size)->run(&cdata[0]);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data in memory and provided centers
**/
RcppExport SEXP R_knor_hmeans_data_im_centers(SEXP rdata, SEXP rk,
        SEXP rmax_iters, SEXP rnthread, SEXP rtolerance,
        SEXP rdist_type, SEXP rmin_clust_size) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rk); // FIXME
	unsigned k = INTEGER(rk)[0]; // FIXME: kmax
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];
	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();
    unsigned nnodes = kbase::get_num_nodes();

    std::vector<double> cdata(nrow*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    std::vector<double> ccentroids(centroids.nrow()*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < k; row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    kbase::cluster_t kret =
        knor::hclust_coordinator::create("", nrow, ncol, k,
                max_iters, nnodes, nthread, &ccentroids[0], "none", tolerance,
                dist_type, min_clust_size)->run(&cdata[0]);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}
////////////////////////////// END HMEANS //////////////////////////////////////

////////////////////////////////// XMEANS //////////////////////////////////////
/**
  * Data on disk str init
  */
RcppExport SEXP R_knor_xmeans_data_em_init(SEXP rdata, SEXP rkmax,
        SEXP rnrow, SEXP rncol, SEXP rmax_iters, SEXP rnthread, SEXP rinit,
        SEXP rtolerance, SEXP rdist_type, SEXP rmin_clust_size) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	unsigned kmax = INTEGER(rkmax)[0];
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();
    unsigned nnodes = kbase::get_num_nodes();

    kbase::cluster_t kret =
                knor::xmeans_coordinator::create(data, nrow, ncol, kmax,
                max_iters, nnodes, nthread, NULL,
                init, tolerance, dist_type, min_clust_size)->run();

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data on disk centers provided
  */
RcppExport SEXP R_knor_xmeans_data_em_centers(SEXP rdata, SEXP rkmax,
        SEXP rnrow, SEXP rncol, SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance, SEXP rdist_type, SEXP rmin_clust_size) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	int kmax = INTEGER(rkmax)[0];
    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rinit);
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];

    std::vector<double> ccentroids(centroids.nrow()*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < centroids.nrow(); row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();
    unsigned nnodes = kbase::get_num_nodes();

    kbase::cluster_t kret =
        knor::xmeans_coordinator::create(data, nrow, ncol, kmax,
                max_iters, nnodes, nthread, &ccentroids[0], "none",
                tolerance, dist_type, min_clust_size)->run();

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data in memory and str init
**/
RcppExport SEXP R_knor_xmeans_data_im_init(SEXP rdata, SEXP rkmax,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance,
        SEXP rdist_type, SEXP rmin_clust_size) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned kmax = INTEGER(rkmax)[0];
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];

	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    std::vector<double> cdata(nrow*ncol);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    kbase::cluster_t kret =
        knor::xmeans_coordinator::create("", nrow, ncol, kmax,
                max_iters, nnodes, nthread, NULL, init, tolerance,
                dist_type, min_clust_size)->run(&cdata[0]);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data in memory and provided init centers
**/
RcppExport SEXP R_knor_xmeans_data_im_centers(SEXP rdata, SEXP rkmax,
        SEXP rmax_iters, SEXP rnthread, SEXP rinit, SEXP rtolerance,
        SEXP rdist_type, SEXP rmin_clust_size) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned kmax = INTEGER(rkmax)[0];
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rinit);
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];

	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();
    unsigned nnodes = kbase::get_num_nodes();

    std::vector<double> cdata(nrow*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    std::vector<double> ccentroids(centroids.nrow()*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < centroids.nrow(); row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    kbase::cluster_t kret =
        knor::xmeans_coordinator::create("", nrow, ncol, kmax,
                max_iters, nnodes, nthread, &ccentroids[0], "none", tolerance,
                dist_type, min_clust_size)->run(&cdata[0]);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}
////////////////////////////// END XMEANS //////////////////////////////////////

////////////////////////////////// GMEANS //////////////////////////////////////
/**
  * Data on disk str init
  */
RcppExport SEXP R_knor_gmeans_data_em_init(SEXP rdata, SEXP rkmax,
        SEXP rnrow, SEXP rncol, SEXP rmax_iters, SEXP rnthread, SEXP rinit,
        SEXP rtolerance, SEXP rdist_type, SEXP rmin_clust_size,
        SEXP rstrictness) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	unsigned kmax = INTEGER(rkmax)[0];
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];
	unsigned strictness = INTEGER(rstrictness)[0];

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();
    unsigned nnodes = kbase::get_num_nodes();

    kbase::cluster_t kret =
                knor::gmeans_coordinator::create(data, nrow, ncol, kmax,
                max_iters, nnodes, nthread, NULL,
                init, tolerance, dist_type, min_clust_size, strictness)->run();

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data on disk centers provided
  */
RcppExport SEXP R_knor_gmeans_data_em_centers(SEXP rdata, SEXP rkmax,
        SEXP rnrow, SEXP rncol, SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance, SEXP rdist_type,
        SEXP rmin_clust_size, SEXP rstrictness) {

    std::string data = CHAR(STRING_ELT(rdata,0));
	size_t nrow = static_cast<size_t>(REAL(rnrow)[0]);
	size_t ncol = static_cast<size_t>(REAL(rncol)[0]);
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	int kmax = INTEGER(rkmax)[0];
    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rinit);
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];
	unsigned strictness = INTEGER(rstrictness)[0];

    std::vector<double> ccentroids(centroids.nrow()*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < centroids.nrow(); row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();
    unsigned nnodes = kbase::get_num_nodes();

    kbase::cluster_t kret =
        knor::gmeans_coordinator::create(data, nrow, ncol, kmax,
                max_iters, nnodes, nthread, &ccentroids[0], "none",
                tolerance, dist_type, min_clust_size, strictness)->run();

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data in memory and str init
**/
RcppExport SEXP R_knor_gmeans_data_im_init(SEXP rdata, SEXP rkmax,
        SEXP rmax_iters, SEXP rnthread,
        SEXP rinit, SEXP rtolerance,
        SEXP rdist_type, SEXP rmin_clust_size, SEXP rstrictness) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned kmax = INTEGER(rkmax)[0];
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
	std::string init = CHAR(STRING_ELT(rinit,0));
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];
	unsigned strictness = INTEGER(rstrictness)[0];

	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();
    std::vector<double> cdata(nrow*ncol);

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();

    unsigned nnodes = kbase::get_num_nodes();

#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    kbase::cluster_t kret =
        knor::gmeans_coordinator::create("", nrow, ncol, kmax,
                max_iters, nnodes, nthread, NULL, init, tolerance,
                dist_type, min_clust_size, strictness)->run(&cdata[0]);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}

/**
  * Data in memory and provided init centers
**/
RcppExport SEXP R_knor_gmeans_data_im_centers(SEXP rdata, SEXP rkmax,
        SEXP rmax_iters, SEXP rnthread, SEXP rinit, SEXP rtolerance,
        SEXP rdist_type, SEXP rmin_clust_size, SEXP rstrictness) {

    Rcpp::NumericMatrix data = Rcpp::NumericMatrix(rdata);
	unsigned kmax = INTEGER(rkmax)[0];
	size_t max_iters = static_cast<size_t>(REAL(rmax_iters)[0]);
	int nthread = INTEGER(rnthread)[0];
    Rcpp::NumericMatrix centroids = Rcpp::NumericMatrix(rinit);
	double tolerance = REAL(rtolerance)[0];
	std::string dist_type = CHAR(STRING_ELT(rdist_type,0));
	unsigned min_clust_size = INTEGER(rmin_clust_size)[0];
	unsigned strictness = INTEGER(rstrictness)[0];

	const size_t nrow = data.nrow();
	const size_t ncol = data.ncol();

    if (nthread == -1)
        nthread = kbase::get_num_omp_threads();
    unsigned nnodes = kbase::get_num_nodes();

    std::vector<double> cdata(nrow*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(data) shared (cdata)
#endif
	for (size_t row = 0; row < nrow; row++)
		for (size_t col = 0; col < ncol; col++)
			cdata[row*ncol + col] = data(row, col);

    std::vector<double> ccentroids(centroids.nrow()*ncol);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(centroids) shared (ccentroids)
#endif
	for (size_t row = 0; row < centroids.nrow(); row++)
		for (size_t col = 0; col < ncol; col++)
			ccentroids[row*ncol + col] = centroids(row, col);

    kbase::cluster_t kret =
        knor::gmeans_coordinator::create("", nrow, ncol, kmax,
                max_iters, nnodes, nthread, &ccentroids[0], "none", tolerance,
                dist_type, min_clust_size, strictness)->run(&cdata[0]);

	Rcpp::List ret;
    marshall_c_to_r(kret, ret);
	return ret;
}
////////////////////////////// END GMEANS //////////////////////////////////////
