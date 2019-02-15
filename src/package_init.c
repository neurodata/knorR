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
extern SEXP R_kmeans(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_kmeans_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_kmeans_data_centroids_em(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_kmeans_data_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_kmeans_data_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_kmeans_data_im_centroids_em(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// kmedoids
extern SEXP R_kmedoids_data_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_kmedoids_data_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_kmedoids_data_em(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_kmedoids_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// skmeans
extern SEXP R_skmeans_data_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_skmeans_data_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_skmeans_data_em(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_skmeans_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// kmeansPP
extern SEXP R_kmeanspp_data_em(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_kmeanspp_data_im(SEXP, SEXP, SEXP, SEXP, SEXP);

// Mini-batch
extern SEXP R_mbkmeans_data_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_mbkmeans(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_mbkmeans_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_mbkmeans_data_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// Fuzzy C-means
extern SEXP R_fcm_data_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_fcm_data_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_fcm_data_em_centroids_im(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_fcm_data_em(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// Hmeans
extern SEXP R_hmeans_data_em_init(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_hmeans_data_em_centers(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_hmeans_data_im_init(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_hmeans_data_im_centers(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// Xmeans
extern SEXP R_xmeans_data_em_init(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_xmeans_data_em_centers(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_xmeans_data_im_init(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_xmeans_data_im_centers(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// Gmeans
extern SEXP R_gmeans_data_em_init(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_gmeans_data_em_centers(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_gmeans_data_im_init(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_gmeans_data_im_centers(
        SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"R_kmeans",                      (DL_FUNC) &R_kmeans,                       9},
    {"R_kmeans_centroids_im",         (DL_FUNC) &R_kmeans_centroids_im,          7},
    {"R_kmeans_data_centroids_em",    (DL_FUNC) &R_kmeans_data_centroids_em,     9},
    {"R_kmeans_data_centroids_im",    (DL_FUNC) &R_kmeans_data_centroids_im,     6},
    {"R_kmeans_data_im",              (DL_FUNC) &R_kmeans_data_im,               7},
    {"R_kmeans_data_im_centroids_em", (DL_FUNC) &R_kmeans_data_im_centroids_em,  6},

    {"R_kmedoids_data_im",            (DL_FUNC) &R_kmedoids_data_im,             7},
    {"R_kmedoids_data_centroids_im",  (DL_FUNC) &R_kmedoids_data_centroids_im,   6},
    {"R_kmedoids_data_em",            (DL_FUNC) &R_kmedoids_data_em,             9},
    {"R_kmedoids_centroids_im",       (DL_FUNC) &R_kmedoids_centroids_im,        7},

    {"R_skmeans_data_im",             (DL_FUNC) &R_skmeans_data_im,              6},
    {"R_skmeans_data_centroids_im",   (DL_FUNC) &R_skmeans_data_centroids_im,    5},
    {"R_skmeans_data_em",             (DL_FUNC) &R_skmeans_data_em,              8},
    {"R_skmeans_centroids_im",        (DL_FUNC) &R_skmeans_centroids_im,         6},

    {"R_kmeanspp_data_em",            (DL_FUNC) &R_kmeanspp_data_em,             7},
    {"R_kmeanspp_data_im",            (DL_FUNC) &R_kmeanspp_data_im,             5},

    {"R_mbkmeans_data_im",            (DL_FUNC) &R_mbkmeans_data_im,             8},
    {"R_mbkmeans",                    (DL_FUNC) &R_mbkmeans,                    10},
    {"R_mbkmeans_centroids_im",       (DL_FUNC) &R_mbkmeans_centroids_im,        8},
    {"R_mbkmeans_data_centroids_im",  (DL_FUNC) &R_mbkmeans_data_centroids_im,   7},

    {"R_fcm_data_im",                 (DL_FUNC) &R_fcm_data_im,                  8},
    {"R_fcm_data_centroids_im",       (DL_FUNC) &R_fcm_data_centroids_im,        7},
    {"R_fcm_data_em_centroids_im",    (DL_FUNC) &R_fcm_data_em_centroids_im,     9},
    {"R_fcm_data_em",                 (DL_FUNC) &R_fcm_data_em,                 10},

    {"R_hmeans_data_em_init",         (DL_FUNC) &R_hmeans_data_em_init,          10},
    {"R_hmeans_data_em_centers",      (DL_FUNC) &R_hmeans_data_em_centers,       10},
    {"R_hmeans_data_im_init",         (DL_FUNC) &R_hmeans_data_im_init,           8},
    {"R_hmeans_data_im_centers",      (DL_FUNC) &R_hmeans_data_im_centers,        8},

    {"R_xmeans_data_em_init",         (DL_FUNC) &R_xmeans_data_em_init,          10},
    {"R_xmeans_data_em_centers",      (DL_FUNC) &R_xmeans_data_em_centers,       10},
    {"R_xmeans_data_im_init",         (DL_FUNC) &R_xmeans_data_im_init,           8},
    {"R_xmeans_data_im_centers",      (DL_FUNC) &R_xmeans_data_im_centers,        8},

    {"R_gmeans_data_em_init",         (DL_FUNC) &R_gmeans_data_em_init,          11},
    {"R_gmeans_data_em_centers",      (DL_FUNC) &R_gmeans_data_em_centers,       11},
    {"R_gmeans_data_im_init",         (DL_FUNC) &R_gmeans_data_im_init,           9},
    {"R_gmeans_data_im_centers",      (DL_FUNC) &R_gmeans_data_im_centers,        9},

    {NULL, NULL, 0}
};

void R_init(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
