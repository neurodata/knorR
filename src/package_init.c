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

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP R_knor_kmeans(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_knor_kmeans_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_knor_kmeans_data_centroids_em(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_knor_kmeans_data_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_knor_kmeans_data_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_knor_kmeans_data_im_centroids_em(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// kmedoids
extern SEXP R_knor_kmedoids_data_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_knor_kmedoids_data_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_knor_kmedoids_data_em(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_knor_kmedoids_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"R_knor_kmeans",                      (DL_FUNC) &R_knor_kmeans,                      10},
    {"R_knor_kmeans_centroids_im",         (DL_FUNC) &R_knor_kmeans_centroids_im,          8},
    {"R_knor_kmeans_data_centroids_em",    (DL_FUNC) &R_knor_kmeans_data_centroids_em,    10},
    {"R_knor_kmeans_data_centroids_im",    (DL_FUNC) &R_knor_kmeans_data_centroids_im,     8},
    {"R_knor_kmeans_data_im",              (DL_FUNC) &R_knor_kmeans_data_im,               9},
    {"R_knor_kmeans_data_im_centroids_em", (DL_FUNC) &R_knor_kmeans_data_im_centroids_em,  8},

    {"R_knor_kmedoids_data_im",             (DL_FUNC) &R_knor_kmedoids_data_im,            7},
    {"R_knor_kmedoids_data_centroids_im",   (DL_FUNC) &R_knor_kmedoids_data_centroids_im,  6},
    {"R_knor_kmedoids_data_em",             (DL_FUNC) &R_knor_kmedoids_data_em,            9},
    {"R_knor_kmedoids_centroids_im",   (DL_FUNC) &R_knor_kmedoids_centroids_im,       7},
    {NULL, NULL, 0}
};

void R_init_knor(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
