/*  This file is part of the 'quitefastmst' package.
 *
 *  Copyleft (C) 2025-2025, Marek Gagolewski <https://www.gagolewski.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License
 *  Version 3, 19 November 2007, published by the Free Software Foundation.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Affero General Public License Version 3 for more details.
 *  You should have received a copy of the License along with this program.
 *  If this is not the case, refer to <https://www.gnu.org/licenses/>.
 */


#include "c_fastmst.h"
#include "c_common.h"
#include "c_mst_triple.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>


#define MST_OMP_CHUNK_SIZE 1024



/*! Order the n-1 edges of a spanning tree of n points in place,
 * w.r.t. the weights increasingly, resolving ties if needed based on
 * the points' IDs.
 *
 * @param n
 * @param mst_dist [in/out] size m
 * @param mst_ind [in/out] size m*2
 */
template <class FLOAT>
void Ctree_order(Py_ssize_t m, FLOAT* tree_dist, Py_ssize_t* tree_ind)
{
    QUITEFASTMST_PROFILER_USE

    std::vector< CMstTriple<FLOAT> > mst(m);

    for (Py_ssize_t i=0; i<m; ++i) {
        mst[i] = CMstTriple<FLOAT>(tree_ind[2*i+0], tree_ind[2*i+1], tree_dist[i]);
    }

    QUITEFASTMST_PROFILER_START
    std::sort(mst.begin(), mst.end());
    QUITEFASTMST_PROFILER_STOP("mst sort");

    for (Py_ssize_t i=0; i<m; ++i) {
        tree_dist[i]    = mst[i].d;
        tree_ind[2*i+0] = mst[i].i1;  // i1 < i2
        tree_ind[2*i+1] = mst[i].i2;
    }
}



/*! Reconnect leaves of the minimum spanning tree w.r.t. the
 *  M-mutual reachability distances without affecting the minimality condition:
 *  such that, if possible, they are adjacent to inner vertices that have
 *  them amongst their M nearest neighbours;
 *  prefer vertices of the smallest core distances
 *
 *  NOTE This procedure should be run a few times until no more leaves
 *  change their parents.
 *
 * @param n number of rows
 * @param m number of edges in the tree
 * @param M the degree of the "core" distance, M > 0
 * @param tree_dist [in/out] size m - edge weights
 * @param tree_ind [in/out] size m*2 - edges of the tree
 * @param nn_dist [out] n*M Euclidean distances
 *        to the n points' M nearest neighbours
 * @param nn_ind [out] n*M indexes of the n points' M nearest neighbours
 *
 * @return the number of leaves reconnected
 */
template <class FLOAT>
Py_ssize_t Cleaves_reconnect_dcore_min(
    Py_ssize_t m, Py_ssize_t n, Py_ssize_t M,
    FLOAT* tree_dist, Py_ssize_t* tree_ind,
    FLOAT* nn_dist, Py_ssize_t* nn_ind
) {
    std::vector<Py_ssize_t> degrees(n, 0);
    for (Py_ssize_t i=0; i<2*m; ++i) {
        QUITEFASTMST_ASSERT(tree_ind[i] >= 0 && tree_ind[i] < n);
        degrees[tree_ind[i]]++;
    }

    std::vector<Py_ssize_t> closest_inlier(n, -1);
    for (Py_ssize_t v=0; v<n; ++v) {
        QUITEFASTMST_ASSERT(degrees[v] > 0);
        if (degrees[v] == 1) continue;  // a leaf

        FLOAT dcore_v = nn_dist[v*M+(M-1)];
        for (Py_ssize_t j=0; j<M; ++j) {
            Py_ssize_t u  = nn_ind[v*M+j];
            FLOAT dcore_u = nn_dist[u*M+(M-1)];
            if (dcore_v > dcore_u) continue;   // v cannot become adjacent to u (minimality condition!)
            if (closest_inlier[u] < 0 || dcore_v < nn_dist[closest_inlier[u]*M+(M-1)])
                closest_inlier[u] = v;

            // choose v if u is amongst M NNs of v and v itself has "small" core distance
        }
    }

    Py_ssize_t num_changes = 0;
    for (Py_ssize_t i=0; i<m; ++i) {
        for (Py_ssize_t j=0; j<=1; ++j) {
            Py_ssize_t u = tree_ind[i*2+j];
            if (degrees[u] > 1) continue;  // we want u to be a leaf

            Py_ssize_t v = tree_ind[i*2+(1-j)];
            QUITEFASTMST_ASSERT(degrees[v] > 1);  // v is a non-leaf

            Py_ssize_t w = closest_inlier[u];
            if (w >= 0 && w != v) {
                // w will now be the vertex adjacent to u
                num_changes++;
                degrees[v]--;
                degrees[w]++;
                tree_ind[i*2+(1-j)] = w;
            }
        }
    }

    return num_changes;
}


/*! A Jarník (Prim/Dijkstra)-like algorithm for determining
 *  a(*) Euclidean minimum spanning tree (MST) or
 *  one w.r.t. an M-mutual reachability distance.
 *
 *  If `M>1`, the spanning tree is the smallest w.r.t. the degree-`M`
 *  mutual reachability distance [9]_ given by
 *  :math:`d_M(i, j)=\\max\\{ c_M(i), c_M(j), d(i, j)\\}`, where :math:`d(i,j)`
 *  is the Euclidean distance between the `i`-th and the `j`-th point,
 *  and :math:`c_M(i)` is the `i`-th `M`-core distance defined as the distance
 *  between the `i`-th point and its `M`-th nearest neighbour
 *  (not including the query points themselves).
 *
 *  Note that [9]_ defines the core distance as the distance to the (M-1)-th NN.
 *
 *  (*) We note that if there are many pairs of equidistant points,
 *  there can be many minimum spanning trees. In particular, it is likely
 *  that there are point pairs with the same mutual reachability distances.
 *  To make the definition less ambiguous (albeit with no guarantees),
 *  internally, we resolve ties as follows.
 *  The `mutreach_ties` argument indicates the preference towards
 *  connecting to farther(-1)/closer(1) points with respect to the original
 *  metric or having smaller(-2)/larger(2) core distances.
 *
 *  Time complexity: O(n^2). It is assumed that M is rather small
 *  (say, M <= 20). If M>1, all pairwise the distances are computed twice
 *  (first for the neighbours/core distance, then to determine the tree).
 *
 *
 *  References:
 *  ----------
 *
 *  [1] V. Jarník, O jistém problému minimálním,
 *  Práce Moravské Přírodovědecké Společnosti 6, 1930, 57–63.
 *
 *  [2] C.F. Olson, Parallel algorithms for hierarchical clustering,
 *  Parallel Computing 21(8), 1995, 1313–1325.
 *
 *  [3] R. Prim, Shortest connection networks and some generalizations,
 *  The Bell System Technical Journal 36(6), 1957, 1389–1401.
 *
 *  [9] R.J.G.B. Campello, D. Moulavi, J. Sander, Density-based clustering based
 *  on hierarchical density estimates, Lecture Notes in Computer Science 7819,
 *  2013, 160–172. DOI: 10.1007/978-3-642-37456-2_14.
 *
 *
 * @param X [destroyable] a C-contiguous data matrix, shape n*d
 * @param n number of rows
 * @param d number of columns
 * @param M the degree of the "core" distance if M > 0
 * @param mst_dist [out] vector of length n-1, gives weights of the
 *        resulting MST edges in nondecreasing order
 * @param mst_ind [out] vector of length 2*(n-1), representing
 *        a c_contiguous array of shape (n-1,2), defining the edges
 *        corresponding to mst_d, with mst_i[j,0] < mst_i[j,1] for all j
 * @param nn_dist [out] NULL for M==0 or the n*M Euclidean distances
 *        to the n points' M nearest neighbours
 * @param nn_ind [out] NULL for M==0 or the n*M indexes of the n points'
 *        M nearest neighbours
 * @param mutreach_ties adjustment for mutual reachability distance ambiguity
 *        (for M>1): -2 and 2 prefer connecting to points with,
 *        respectively, smaller and larger core distance; -1 and 1 prefer,
 *        respectively, farther and closer nearest neighbours
 * @param verbose should we output diagnostic/progress messages?
 */
template <class FLOAT>
void Cmst_euclid_brute(
    FLOAT* X, Py_ssize_t n, Py_ssize_t d, Py_ssize_t M,
    FLOAT* mst_dist, Py_ssize_t* mst_ind,
    FLOAT* nn_dist, Py_ssize_t* nn_ind,
    Py_ssize_t mutreach_ties,
    bool verbose
) {
    if (n <= 0)   throw std::domain_error("n <= 0");
    if (d <= 0)   throw std::domain_error("d <= 0");
    if (M <  0)   throw std::domain_error("M <  0");
    if (M >= n)   throw std::domain_error("M >= n");
    QUITEFASTMST_ASSERT(mst_dist);
    QUITEFASTMST_ASSERT(mst_ind);

    bool mutreach_adj_via_dcore = (std::abs(mutreach_ties) >= 2);
    FLOAT mutreach_adj = ((mutreach_ties<0)?-1:1);
    FLOAT mutreach_adj_factr = 0.00000011920928955078125;  // 2**-23

    std::vector<FLOAT> d_core;
    if (M > 1) {
        d_core.resize(n);
        QUITEFASTMST_ASSERT(nn_dist);
        QUITEFASTMST_ASSERT(nn_ind);
        Cknn1_euclid_brute(X, n, d, M, nn_dist, nn_ind,
                           /*squared=*/true, verbose);
        for (Py_ssize_t i=0; i<n; ++i) d_core[i] = nn_dist[i*M+M-1];

        // for M==1, we can fetch d_core, nn_dist, nn_ind from the resulting MST,
        // as nearest neighbours are connected by an edge (see below)
    }

    if (verbose) QUITEFASTMST_PRINT("[quitefastmst] Computing the MST... %3d%%", 0);


    // ncl_ind[j] is the vertex from the current tree closest to vertex j
    std::vector<Py_ssize_t> ncl_ind(n);
    std::vector<FLOAT> ncl_dist(n, INFINITY);        // ncl_dist[j] = d_M(j, ncl_ind[j])

    std::vector<FLOAT> ncl_dist_adj;   // ncl_dist[j] = adjustment for d_M(j, ncl_ind[j])'s ambiguity
    if (M > 1) ncl_dist_adj.resize(n, INFINITY);

    std::vector<Py_ssize_t> remaining_ind(n);  // a.k.a. perm
    for (Py_ssize_t i=0; i<n; ++i) remaining_ind[i] = i;

    std::vector< CMstTriple<FLOAT> > mst(n-1);

    for (Py_ssize_t i=1; i<n; ++i) {
        // i-1 is the vertex most recently added to the tree
        // i, i+1, ..., n-1 - vertices not yet in the tree

        FLOAT* x_cur = X+(i-1)*d;

        // compute the distances
        // between the (i-1)-th vertex and all j=i, i+1, ..., n-1:
#if 0
        // NOTE two-stage Euclidean distance computation: slower -> removed
#else
        if (M <= 1) {
            #if OPENMP_IS_ENABLED
            #pragma omp parallel for schedule(static,MST_OMP_CHUNK_SIZE)  /* chunks get smaller and smaller... */
            #endif
            for (Py_ssize_t j=i; j<n; ++j) {
                FLOAT dd = 0.0;
                for (Py_ssize_t u=0; u<d; ++u)
                    dd += square(x_cur[u]-X[j*d+u]);

                if (dd < ncl_dist[j]) {
                    ncl_dist[j] = dd;
                    ncl_ind[j] = i-1;
                }
            }
        }
        else
        {
            #if OPENMP_IS_ENABLED
            #pragma omp parallel for schedule(static,MST_OMP_CHUNK_SIZE)
            #endif
            for (Py_ssize_t j=i; j<n; ++j) {
                if (ncl_dist[j] < d_core[i-1]) continue;
                if (ncl_dist[j] < d_core[j]) continue;

                FLOAT dd = 0.0;
                for (Py_ssize_t u=0; u<d; ++u)
                    dd += square(x_cur[u]-X[j*d+u]);

                if (ncl_dist[j] < dd) continue;  // nothing to do

                FLOAT dd_orig = dd;
                FLOAT dd_adj;
                FLOAT d_core_min;
                FLOAT d_core_max;

                if (d_core[i-1] <= d_core[j]) {
                    d_core_min = d_core[i-1];
                    d_core_max = d_core[j];
                }
                else {
                    d_core_min = d_core[j];
                    d_core_max = d_core[i-1];
                }

                if (dd > d_core_max) {
                    dd_adj = 0.0;
                }
                else {
                    dd = d_core_max;
                    if (mutreach_adj_via_dcore)
                        dd_adj = mutreach_adj*(-d_core_min+mutreach_adj_factr*dd_orig);
                    else
                        dd_adj = mutreach_adj*dd_orig;
                }

                if (dd < ncl_dist[j] || (dd == ncl_dist[j] && dd_adj < ncl_dist_adj[j]))  {
                    ncl_dist[j] = dd;
                    ncl_dist_adj[j] = dd_adj;
                    ncl_ind[j] = i-1;
                }
            }
        }
#endif

        // we want to include the vertex that is closest to
        // the vertices of the tree constructed so far
        Py_ssize_t best_j = i;
        for (Py_ssize_t j=i+1; j<n; ++j) {
            if (ncl_dist[j] < ncl_dist[best_j] || (ncl_dist[j] == ncl_dist[best_j] && ncl_dist_adj[j] < ncl_dist_adj[best_j]))
                best_j = j;
        }

        // with swapping we get better locality of reference
        std::swap(remaining_ind[best_j], remaining_ind[i]);
        std::swap(ncl_dist[best_j], ncl_dist[i]);
        std::swap(ncl_ind[best_j], ncl_ind[i]);
        for (Py_ssize_t u=0; u<d; ++u) std::swap(X[best_j*d+u], X[i*d+u]);


        if (M > 1) {
            std::swap(d_core[best_j], d_core[i]);
            std::swap(ncl_dist_adj[best_j], ncl_dist_adj[i]);
        }

        // don't visit i again - it's being added to the tree

        // connect best_remaining_ind with the tree: add a new edge {best_remaining_ind, ncl_ind[best_remaining_ind]}
        QUITEFASTMST_ASSERT(ncl_ind[i] < i);
        mst[i-1] = CMstTriple<FLOAT>(remaining_ind[ncl_ind[i]], remaining_ind[i], ncl_dist[i], /*order=*/true);


        if (verbose) QUITEFASTMST_PRINT("\b\b\b\b%3d%%", (int)((n-1+n-i-1)*(i+1)*100/n/(n-1)));

        if (i % MST_OMP_CHUNK_SIZE == MST_OMP_CHUNK_SIZE-1) {
            #if QUITEFASTMST_R
            Rcpp::checkUserInterrupt();  // throws an exception, not a longjmp
            #elif QUITEFASTMST_PYTHON
            if (PyErr_CheckSignals() != 0) throw std::runtime_error("signal caught");
            #endif
        }
    }

    // sort the resulting MST edges in increasing order w.r.t. d
    std::sort(mst.begin(), mst.end());

    for (Py_ssize_t i=0; i<n-1; ++i) {
        mst_dist[i]    = sqrt(mst[i].d);
        mst_ind[2*i+0] = mst[i].i1; // i1 < i2
        mst_ind[2*i+1] = mst[i].i2;
    }

    // set up nn_dist (and nn_ind) - only for M>=1
    if (M > 1) {
        for (Py_ssize_t i=0; i<n*M; ++i)
            nn_dist[i] = sqrt(nn_dist[i]);
    }
    else if (M == 1) {
        // for M==1 we just need the nearest neighbours,
        // as the MST connects them with each other
        for (Py_ssize_t i=0; i<n; ++i)
            nn_dist[i] = INFINITY;

        for (Py_ssize_t i=0; i<n-1; ++i) {
            if (nn_dist[mst_ind[2*i+0]] > mst_dist[i]) {
                nn_dist[mst_ind[2*i+0]] = mst_dist[i];
                nn_ind[mst_ind[2*i+0]] = mst_ind[2*i+1];
            }
            if (nn_dist[mst_ind[2*i+1]] > mst_dist[i]) {
                nn_dist[mst_ind[2*i+1]] = mst_dist[i];
                nn_ind[mst_ind[2*i+1]] = mst_ind[2*i+0];
            }
        }
    }

    if (verbose) QUITEFASTMST_PRINT("\b\b\b\bdone.\n");
}


// instantiate:
template void Ctree_order<float>(Py_ssize_t m, float* tree_dist, Py_ssize_t* tree_ind);

template void Ctree_order<double>(Py_ssize_t m, double* tree_dist, Py_ssize_t* tree_ind);

template Py_ssize_t Cleaves_reconnect_dcore_min<float>(
    Py_ssize_t m, Py_ssize_t n, Py_ssize_t M,
    float* tree_dist, Py_ssize_t* tree_ind,
    float* nn_dist, Py_ssize_t* nn_ind
);

template Py_ssize_t Cleaves_reconnect_dcore_min<double>(
    Py_ssize_t m, Py_ssize_t n, Py_ssize_t M,
    double* tree_dist, Py_ssize_t* tree_ind,
    double* nn_dist, Py_ssize_t* nn_ind
);

template void Cmst_euclid_brute<float>(
    float* X, Py_ssize_t n, Py_ssize_t d, Py_ssize_t M,
    float* mst_dist, Py_ssize_t* mst_ind,
    float* nn_dist, Py_ssize_t* nn_ind,
    Py_ssize_t mutreach_ties,
    bool verbose
);

template void Cmst_euclid_brute<double>(
    double* X, Py_ssize_t n, Py_ssize_t d, Py_ssize_t M,
    double* mst_dist, Py_ssize_t* mst_ind,
    double* nn_dist, Py_ssize_t* nn_ind,
    Py_ssize_t mutreach_ties,
    bool verbose
);
