/*  Borůvka-type algorithms for finding minimum spanning trees
 *  wrt the Euclidean metric or the thereon-based mutual reachability distance.
 *
 *  The dual-tree Borůvka version is, in principle, based on
 *  "Fast Euclidean Minimum Spanning Tree: Algorithm, Analysis,
 *  and Applications" by W.B. March, P. Ram, A.G. Gray published
 *  in ACM SIGKDD 2010.  As far as our implementation
 *  is concerned, the dual-tree approach is only faster in 2- and
 *  3-dimensional spaces, for M <= 2, and in a single-threaded setting
 *  (in the current implementation, only the first iteration is parallelised
 *  anyway).
 *
 *  The single-tree version (iteratively find each point's nearest neighbour
 *  outside its own cluster, i.e., nearest alien) is naively parallelisable.
 *
 *  The "sesqui-tree" variant (by the current author) is a mixture of the two
 *  approaches:  it compares leaves against the full tree.  It is usually
 *  faster than the single- and dual-tree methods in very low dimensional
 *  spaces and usually not much slower than the single-tree variant otherwise.
 *
 *  For more details on our implementation of K-d trees, see
 *  the source file defining the base class.
 *
 *
 *
 *  Copyleft (C) 2025, Marek Gagolewski <https://www.gagolewski.com>
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


#ifndef __c_kdtree_boruvka_h
#define __c_kdtree_boruvka_h

#include "c_common.h"
#include "c_kdtree.h"
#include "c_disjoint_sets.h"


namespace quitefastkdtree {

template <typename FLOAT, Py_ssize_t D>
struct kdtree_node_clusterable : public kdtree_node_base<FLOAT, D>
{
    kdtree_node_clusterable* left;
    kdtree_node_clusterable* right;

    Py_ssize_t cluster_repr;  //< representative point index if all descendants are in the same cluster, -1 otherwise

    struct t_dtb_data { FLOAT cluster_max_dist; FLOAT min_dcore /* M>2 */; };
    struct t_qtb_data { FLOAT lastbest_dist; Py_ssize_t lastbest_ind; Py_ssize_t lastbest_from; };

    union {
        t_dtb_data dtb_data;
        t_qtb_data qtb_data;
    };


    kdtree_node_clusterable()
    {
        left = nullptr;
        // right = nullptr;
    }

    inline bool is_leaf() const
    {
        return left == nullptr /*&& right == nullptr*/; // either both null or none
    }
};


template <
    typename FLOAT,
    Py_ssize_t D,
    typename DISTANCE=kdtree_distance_sqeuclid<FLOAT,D>,
    typename NODE=kdtree_node_clusterable<FLOAT, D>
>
struct kdtree_node_orderer {
    NODE* nearer_node;
    NODE* farther_node;
    FLOAT nearer_dist;
    FLOAT farther_dist;

    kdtree_node_orderer(const FLOAT* x, NODE* to1, NODE* to2)  // QTB, STB
    {
        nearer_dist  = DISTANCE::point_node(
            x, to1->bbox_min.data(),  to1->bbox_max.data()
        );

        farther_dist = DISTANCE::point_node(
            x, to2->bbox_min.data(),  to2->bbox_max.data()
        );

        if (nearer_dist <= farther_dist) {
            nearer_node  = to1;
            farther_node = to2;
        }
        else {
            std::swap(nearer_dist, farther_dist);
            nearer_node  = to2;
            farther_node = to1;
        }
    }

    kdtree_node_orderer(NODE* from, NODE* to1, NODE* to2, bool use_min_dcore=false)  // DTB
    {
        nearer_dist  = DISTANCE::node_node(
            from->bbox_min.data(), from->bbox_max.data(),
            to1->bbox_min.data(),  to1->bbox_max.data()
        );

        farther_dist = DISTANCE::node_node(
            from->bbox_min.data(), from->bbox_max.data(),
            to2->bbox_min.data(),  to2->bbox_max.data()
        );

        if (use_min_dcore)
        {
            nearer_dist  = max3(nearer_dist,  from->dtb_data.min_dcore, to1->dtb_data.min_dcore);
            farther_dist = max3(farther_dist, from->dtb_data.min_dcore, to2->dtb_data.min_dcore);
        }

        if (nearer_dist <= farther_dist) {
            nearer_node  = to1;
            farther_node = to2;
        }
        else {
            std::swap(nearer_dist, farther_dist);
            nearer_node  = to2;
            farther_node = to1;
        }
    }
};



/** A class enabling searching for the nearest neighbour
 *  outside of the current point's cluster;
 *  (for the "sesqui-tree" and "single-tree" Borůvka algo); it is thread-safe
 */
template <
    typename FLOAT,
    Py_ssize_t D,
    typename DISTANCE=kdtree_distance_sqeuclid<FLOAT,D>,
    typename NODE=kdtree_node_clusterable<FLOAT, D>
>
class kdtree_nearest_outsider
{
private:
    const FLOAT* data;      ///< the dataset
    const FLOAT* dcore;     ///< the "core" distances
    Py_ssize_t M;

    const Py_ssize_t* ds_par;  ///< points' cluster IDs (par[i]==ds.find(i)!)

    FLOAT nn_dist;          ///< shortest distance
    Py_ssize_t nn_ind;      ///< index of the nn
    Py_ssize_t nn_from;

    const FLOAT* x;     ///< the point itself (shortcut) / first point
    NODE* curleaf;      ///< nullptr or a whole leaf
    Py_ssize_t which;   ///< for which point are we getting the nns / first point index
    Py_ssize_t cluster; ///< the point's / points' cluster



    template <bool USE_DCORE>
    inline void point_vs_points(Py_ssize_t idx_from, Py_ssize_t idx_to)
    {
        const FLOAT* y = data+D*idx_from;
        for (Py_ssize_t j=idx_from; j<idx_to; ++j, y+=D) {
            if (cluster == ds_par[j]) continue;

            if (USE_DCORE && dcore[j] >= nn_dist) continue;
            FLOAT dd = DISTANCE::point_point(x, y);
            if (USE_DCORE) dd = max3(dd, dcore[which], dcore[j]);
            if (dd < nn_dist) {
                nn_dist = dd;
                nn_ind = j;
            }
        }
    }


    template <bool USE_DCORE>
    void find_nn_single(const NODE* root)
    {
        if (root->cluster_repr == cluster) {
            // nothing to do - all are members of the x's cluster
            return;
        }

        if (root->is_leaf()/* || root->idx_to-root->idx_from <= max_brute_size*/) {
            if (which < root->idx_from || which >= root->idx_to)
                point_vs_points<USE_DCORE>(root->idx_from, root->idx_to);
            else {
                point_vs_points<USE_DCORE>(root->idx_from, which);
                point_vs_points<USE_DCORE>(which+1, root->idx_to);
            }
            return;
        }


        kdtree_node_orderer<FLOAT, D, DISTANCE, NODE> sel(
            x, root->left, root->right
        );

        if (sel.nearer_dist < nn_dist) {
            find_nn_single<USE_DCORE>(sel.nearer_node);

            if (sel.farther_dist < nn_dist)
                find_nn_single<USE_DCORE>(sel.farther_node);
        }
    }


    template <bool USE_DCORE>
    void find_nn_multi(const NODE* root)
    {
        if (root->cluster_repr == curleaf->cluster_repr) {
            // nothing to do - all are members of the x's cluster
            return;
        }

        if (root->is_leaf()) {
            const FLOAT* _y = data+D*root->idx_from;
            for (Py_ssize_t j=root->idx_from; j<root->idx_to; ++j, _y+=D) {
                if (curleaf->cluster_repr == ds_par[j]) continue;
                if (USE_DCORE && dcore[j] >= nn_dist) continue;

                const FLOAT* _x = x;
                for (Py_ssize_t i=curleaf->idx_from; i<curleaf->idx_to; ++i, _x+=D) {
                    if (USE_DCORE && dcore[i] >= nn_dist) continue;
                    FLOAT dd = DISTANCE::point_point(_x, _y);
                    if (USE_DCORE) dd = max3(dd, dcore[i], dcore[j]);
                    if (dd < nn_dist) {
                        nn_dist = dd;
                        nn_ind = j;
                        nn_from = i;
                    }
                }
            }
            return;
        }


        kdtree_node_orderer<FLOAT, D, DISTANCE, NODE> sel(
            curleaf, root->left, root->right
        );

        if (sel.nearer_dist < nn_dist) {
            find_nn_multi<USE_DCORE>(sel.nearer_node);

            if (sel.farther_dist < nn_dist)
                find_nn_multi<USE_DCORE>(sel.farther_node);
        }
    }


public:
    kdtree_nearest_outsider(
        const FLOAT* data,
        FLOAT* dcore,
        Py_ssize_t M,
        const Py_ssize_t* ds_par
    ) :
        data(data), dcore(dcore), M(M), ds_par(ds_par)
    {
        ;
    }


    /**
     *  @param curleaf
     *  @param root
     *  @param nn_dist best nn_dist found so far for the current cluster
     */
    void find_multi(NODE* curleaf, const NODE* root, FLOAT nn_dist=INFINITY)
    {
        this->nn_dist = nn_dist;
        this->nn_ind  = -1;
        this->nn_from = -1;

        this->curleaf = curleaf;
        this->which = curleaf->idx_from;
        this->x = data+D*this->which;
        this->cluster = curleaf->cluster_repr;

        if (M>2) find_nn_multi<true>(root);
        else find_nn_multi<false>(root);
    }


    /**
     *  @param which
     *  @param root
     *  @param nn_dist best nn_dist found so far for the current cluster
     */
    void find_single(Py_ssize_t which, const NODE* root, FLOAT nn_dist=INFINITY)
    {
        this->nn_dist = nn_dist;
        this->nn_ind  = -1;
        this->nn_from = which;

        this->curleaf = nullptr;
        this->which = which;
        this->x = data+D*this->which;
        this->cluster = ds_par[this->which];

        if (M>2) find_nn_single<true>(root);
        else find_nn_single<false>(root);
    }


    inline FLOAT get_nn_dist()      { return nn_dist; }
    inline Py_ssize_t get_nn_ind()  { return nn_ind; }
    inline Py_ssize_t get_nn_from() { return nn_from; }
};




template <
    typename FLOAT,
    Py_ssize_t D,
    typename DISTANCE=kdtree_distance_sqeuclid<FLOAT,D>,
    typename NODE=kdtree_node_clusterable<FLOAT, D>
>
class kdtree_boruvka : public kdtree<FLOAT, D, DISTANCE, NODE>
{
protected:
    FLOAT*  tree_dist;     ///< size n-1
    Py_ssize_t* tree_ind;  ///< size 2*(n-1)
    Py_ssize_t  tree_edges;  /// number of MST edges already found
    Py_ssize_t tree_iter;
    CDisjointSets ds;

    std::vector<FLOAT>      ncl_dist;  // ncl_dist[find(i)] - distance to i's nn
    std::vector<Py_ssize_t> ncl_ind;   // ncl_ind[find(i)]  - index of i's nn
    std::vector<Py_ssize_t> ncl_from;  // ncl_from[find(i)] - the relevant member of i

    const Py_ssize_t first_pass_max_brute_size;  // used in the first iter (finding 1-nns)

    enum BORUVKA_TYPE { BORUVKA_STB, BORUVKA_QTB, BORUVKA_DTB };
    BORUVKA_TYPE boruvka_variant;
    bool reset_nns;

    const FLOAT mutreach_adj;  // M>2 only

    std::vector<FLOAT> lastbest_dist;   // !use_dtb only
    std::vector<Py_ssize_t> lastbest_ind;   // !use_dtb only

    const Py_ssize_t M;              // mutual reachability distance - "smoothing factor"
    std::vector<FLOAT> dcore;        // distances to the (M-1)-th nns of each point if M>1 or 1-NN for M==1
    std::vector<FLOAT> Mnn_dist;     // M-1 nearest neighbours of each point if M>1
    std::vector<Py_ssize_t> Mnn_ind;

    #if OPENMP_IS_ENABLED
    omp_lock_t omp_lock;
    #endif
    int omp_nthreads;


    std::vector<NODE*> leaves; // sesquitree only


    inline void tree_add(Py_ssize_t i, Py_ssize_t j, FLOAT d)
    {
        tree_ind[tree_edges*2+0] = i;
        tree_ind[tree_edges*2+1] = j;
        tree_dist[tree_edges] = d;
        ds.merge(i, j);
        tree_edges++;
    }


    void setup_leaves()
    {
        QUITEFASTMST_ASSERT(boruvka_variant == BORUVKA_QTB);

        leaves.resize(this->nleaves);

        Py_ssize_t _leafnum = 0;
        for (auto curnode = this->nodes.begin(); curnode != this->nodes.end(); ++curnode) {
            if (curnode->is_leaf()) {
                leaves[_leafnum++] = &(*curnode);
                curnode->qtb_data.lastbest_dist = 0.0;
                curnode->qtb_data.lastbest_ind  = -1;
                curnode->qtb_data.lastbest_from = -1;
            }
        }
        QUITEFASTMST_ASSERT(_leafnum == this->nleaves);
    }


    void setup_min_dcore()
    {
        QUITEFASTMST_ASSERT(M>=2);
        QUITEFASTMST_ASSERT(boruvka_variant == BORUVKA_DTB);

        for (auto curnode = this->nodes.rbegin(); curnode != this->nodes.rend(); ++curnode)
        {
            if (curnode->is_leaf()) {
                curnode->dtb_data.min_dcore = dcore[curnode->idx_from];
                for (Py_ssize_t i=curnode->idx_from+1; i<curnode->idx_to; ++i) {
                    if (dcore[i] < curnode->dtb_data.min_dcore)
                        curnode->dtb_data.min_dcore = dcore[i];
                }
            }
            else {
                // all descendants have already been processed as children in `nodes` occur after their parents
                curnode->dtb_data.min_dcore = std::min(
                    curnode->left->dtb_data.min_dcore,
                    curnode->right->dtb_data.min_dcore
                );
            }
        }
    }


    void update_node_data()
    {
        // Performed in each iteration

        // ds.find(i) == ds.get_parent(i) for all i

        // nodes is a deque...
        for (auto curnode = this->nodes.rbegin(); curnode != this->nodes.rend(); ++curnode)
        {
            if (boruvka_variant == BORUVKA_DTB) curnode->dtb_data.cluster_max_dist = INFINITY;  // for DTB

            if (curnode->cluster_repr >= 0) {
                curnode->cluster_repr = ds.get_parent(curnode->cluster_repr);
                continue;
            }

            if (curnode->is_leaf()) {
                curnode->cluster_repr = ds.get_parent(curnode->idx_from);
                for (Py_ssize_t j=curnode->idx_from+1; j<curnode->idx_to; ++j) {
                    if (curnode->cluster_repr != ds.get_parent(j)) {
                        curnode->cluster_repr = -1;  // not all are members of the same cluster
                        break;
                    }
                }

                if (curnode->cluster_repr >= 0 && boruvka_variant == BORUVKA_QTB) {
                    Py_ssize_t i=curnode->idx_from;
                    curnode->qtb_data.lastbest_dist = lastbest_dist[i];
                    curnode->qtb_data.lastbest_ind  = lastbest_ind[i];
                    curnode->qtb_data.lastbest_from = i;
                    for (++i; i<curnode->idx_to; ++i) {
                        if (curnode->qtb_data.lastbest_dist > lastbest_dist[i]) {
                            curnode->qtb_data.lastbest_dist = lastbest_dist[i];
                            curnode->qtb_data.lastbest_ind  = lastbest_ind[i];
                            curnode->qtb_data.lastbest_from = i;
                        }
                    }
                }
            }
            else {
                // all descendants have already been processed as children in `nodes` occur after their parents
                if (curnode->left->cluster_repr >= 0) {
                    // if both children only feature members of the same cluster, update the cluster repr for the current node;
                    if (curnode->left->cluster_repr == curnode->right->cluster_repr)
                        curnode->cluster_repr = curnode->left->cluster_repr;
                }
                // else curnode->cluster_repr = -1;  // it already is
            }
        }
    }


    void update_nn_data()
    {
        if (boruvka_variant != BORUVKA_DTB && tree_iter > 1) {
            // if tree_iter == 1, then all lastbest_ind[i] == -1;
            // we don't get access to individual NNs in DTB, except in the 1st iter

            for (Py_ssize_t i=0; i<this->n; ++i) {
                if (lastbest_ind[i] < 0) continue;

                Py_ssize_t ds_find_i = ds.get_parent(i);
                Py_ssize_t ds_find_j = ds.get_parent(lastbest_ind[i]);
                if (ds_find_i == ds_find_j) {
                    lastbest_ind[i] = -1;
                    continue;
                }

                if (ncl_dist[ds_find_i] > lastbest_dist[i]) {
                    ncl_dist[ds_find_i] = lastbest_dist[i];
                    ncl_ind[ds_find_i]  = lastbest_ind[i];
                    ncl_from[ds_find_i] = i;
                }

                // ok even if nthreads>1
                if (ncl_dist[ds_find_j] > lastbest_dist[i]) {
                    ncl_dist[ds_find_j] = lastbest_dist[i];
                    ncl_ind[ds_find_j]  = i;
                    ncl_from[ds_find_j] = lastbest_ind[i];
                }
            }
        }

        if (M > 2) {
            // reuse M-1 NNs if d==dcore[i] as an initialiser to ncl_ind/dist/from;
            // good speed-up sometimes (we'll be happy with any match; leaves
            // are formed in the 1st iteration of the algorithm)
            const Py_ssize_t k = M-1;
            for (Py_ssize_t i=0; i<this->n; ++i) {
                Py_ssize_t ds_find_i = ds.get_parent(i);
                if (ncl_dist[ds_find_i] <= lastbest_dist[i] || lastbest_dist[i] > dcore[i]) continue;
                for (Py_ssize_t v=0; v<k; ++v)
                {
                    Py_ssize_t j = Mnn_ind[i*(M-1)+((mutreach_adj<0.0)?(k-1-v):(v))];
                    if (ds_find_i == ds.get_parent(j) || dcore[i] < dcore[j]) continue;

                    ncl_dist[ds_find_i] = dcore[i];
                    ncl_ind[ds_find_i]  = j;
                    ncl_from[ds_find_i] = i;

                    lastbest_dist[i] = dcore[i];  // actually unchanged
                    lastbest_ind[i] = j;

                    if (boruvka_variant == BORUVKA_DTB || omp_nthreads == 1) {
                        Py_ssize_t ds_find_j = ds.get_parent(j);
                        if (ncl_dist[ds_find_j] > dcore[i]) {
                            ncl_dist[ds_find_j] = dcore[i];
                            ncl_ind[ds_find_j]  = i;
                            ncl_from[ds_find_j] = j;
                        }
                    }

                    break;  // other candidates have d_M >= dcore[i] anyway
                }
            }
        }
    }



    void find_mst_first_1()
    {
        QUITEFASTMST_ASSERT(M <= 2);
        const Py_ssize_t k = 1;

        for (Py_ssize_t i=0; i<this->n; ++i) ncl_dist[i] = INFINITY;
        for (Py_ssize_t i=0; i<this->n; ++i) ncl_ind[i] = -1;

        // find 1-nns of each point using max_brute_size,
        // preferably with max_brute_size>max_leaf_size
        #if OPENMP_IS_ENABLED
        #pragma omp parallel for schedule(static)
        #endif
        for (Py_ssize_t i=0; i<this->n; ++i) {
            kdtree_kneighbours<FLOAT, D, DISTANCE, NODE> nn(
                this->data, nullptr, i, &ncl_dist[i], &ncl_ind[i], k,
                first_pass_max_brute_size
            );
            nn.find(&this->nodes[0], /*reset=*/false);

            if (omp_nthreads == 1 && ncl_dist[i] < ncl_dist[ncl_ind[i]]) {
                // the speed up is rather small...
                ncl_dist[ncl_ind[i]] = ncl_dist[i];
                ncl_ind[ncl_ind[i]] = i;
            }

            lastbest_ind[i] = -1;  // inactive
            lastbest_dist[i] = ncl_dist[i];

            if (M > 1) {
                dcore[i]    = ncl_dist[i];
                Mnn_dist[i] = ncl_dist[i];
                Mnn_ind[i]  = ncl_ind[i];
            }

        }

        // connect nearest neighbours with each other
        for (Py_ssize_t i=0; i<this->n; ++i) {
            if (ds.find(i) != ds.find(ncl_ind[i])) {
                tree_add(i, ncl_ind[i], ncl_dist[i]);
            }
        }
    }


    void find_mst_first_M()
    {
        QUITEFASTMST_ASSERT(M>1);
        const Py_ssize_t k = M-1;
        // find (M-1)-nns of each point

        for (size_t i=0; i<Mnn_dist.size(); ++i) Mnn_dist[i] = INFINITY;
        for (size_t i=0; i<Mnn_ind.size(); ++i)  Mnn_ind[i] = -1;

        #if OPENMP_IS_ENABLED
        #pragma omp parallel for schedule(static)
        #endif
        for (Py_ssize_t i=0; i<this->n; ++i) {
            kdtree_kneighbours<FLOAT, D, DISTANCE, NODE> nn(
                this->data, nullptr, i, Mnn_dist.data()+k*i, Mnn_ind.data()+k*i, k,
                first_pass_max_brute_size
            );
            nn.find(&this->nodes[0], /*reset=*/false);
            dcore[i] = Mnn_dist[i*k+(k-1)];

            lastbest_dist[i] = dcore[i];  // merely a lower bound
            lastbest_ind[i] = -M;
        }


        // k-nns wrt Euclidean distances are not necessarily k-nns wrt M-mutreach
        // k-nns have d_M >= d_core

        // dcore[i] is definitely the smallest possible d_M(i, *); i!=*
        // we can only be sure that j is a NN if d_M(i, j) == dcore[i]

        // but NNs wrt d_m might be ambiguous - we might want to pick,
        // e.g., the farthest or the closest one wrt the original dist

        // the correction for ambiguity is only applied at this stage!

        if (mutreach_adj <= -1 || mutreach_adj >= 1) {
            for (Py_ssize_t i=0; i<this->n; ++i) {
                // mutreach_adj <= -1 - connect with j whose dcore[j] is the smallest
                // mutreach_adj >=  1 - connect with j whose dcore[j] is the largest

                Py_ssize_t bestj = -1;
                FLOAT bestdcorej = (mutreach_adj <= -1)?INFINITY:(-INFINITY);
                for (Py_ssize_t v=0; v<k; ++v)
                {
                    Py_ssize_t j = Mnn_ind[i*k+v];
                    if (dcore[i] >= dcore[j] && ds.find(i) != ds.find(j)) {
                        if (
                            (mutreach_adj <= -1 && bestdcorej >= dcore[j]) ||
                            (mutreach_adj >=  1 && bestdcorej <  dcore[j])
                        ) {
                            bestj = j;
                            bestdcorej = dcore[j];
                        }
                    }
                }
                if (bestj >= 0) tree_add(i, bestj, dcore[i]);
            }
        }
        else {
            for (Py_ssize_t i=0; i<this->n; ++i) {
                // connect with j whose d(i,j) is the smallest (1>mutreach_adj>0) or largest (-1<mutreach_adj<0)
                // stops searching early, because the original distances are sorted

                for (Py_ssize_t v=0; v<k; ++v)
                {
                    Py_ssize_t j = Mnn_ind[i*k+((mutreach_adj<0.0)?(k-1-v):(v))];
                    if (dcore[i] >= dcore[j] && ds.find(i) != ds.find(j)) {
                        // j is the nearest neighbour of i wrt mutreach dist.
                        tree_add(i, j, dcore[i]);
                        break;  // other candidates have d_M >= dcore[i] anyway
                    }
                }
            }
        }
    }


    void find_mst_first()
    {
        // the 1st iteration: connect nearest neighbours with each other
        if (M <= 2) find_mst_first_1();
        else        find_mst_first_M();
    }


    template <bool USE_DCORE>
    inline void leaf_vs_leaf_dtb(NODE* roota, NODE* rootb)
    {
        // assumes ds.find(i) == ds.get_parent(i) for all i!
        const FLOAT* _x = this->data + roota->idx_from*D;
        for (Py_ssize_t i=roota->idx_from; i<roota->idx_to; ++i, _x += D)
        {
            Py_ssize_t ds_find_i = ds.get_parent(i);
            if (USE_DCORE && dcore[i] >= ncl_dist[ds_find_i]) continue;

            for (Py_ssize_t j=rootb->idx_from; j<rootb->idx_to; ++j)
            {
                Py_ssize_t ds_find_j = ds.get_parent(j);
                if (ds_find_i == ds_find_j) continue;
                if (USE_DCORE && dcore[j] >= ncl_dist[ds_find_i]) continue;

                FLOAT dij = DISTANCE::point_point(_x, this->data+j*D);

                if (USE_DCORE) dij = max3(dij, dcore[i], dcore[j]);

                if (dij < ncl_dist[ds_find_i]) {
                    ncl_dist[ds_find_i] = dij;
                    ncl_ind[ds_find_i]  = j;
                    ncl_from[ds_find_i] = i;
                }
            }
        }
    }


    void find_mst_next_dtb(NODE* roota, NODE* rootb)
    {
        // we have ds.find(i) == ds.get_parent(i) for all i!

        if (roota->cluster_repr >= 0 && roota->cluster_repr == rootb->cluster_repr) {
            // both consist of members of the same cluster - nothing to do
            return;
        }

        if (roota->is_leaf()) {
            if (rootb->is_leaf()) {

                if (M>2) leaf_vs_leaf_dtb<true>(roota, rootb);
                else     leaf_vs_leaf_dtb<false>(roota, rootb);

                if (roota->cluster_repr >= 0) {  // all points are in the same cluster
                    roota->dtb_data.cluster_max_dist = ncl_dist[roota->cluster_repr];
                }
                else {
                    roota->dtb_data.cluster_max_dist = ncl_dist[ds.get_parent(roota->idx_from)];
                    for (Py_ssize_t i=roota->idx_from+1; i<roota->idx_to; ++i) {
                        FLOAT dist_cur = ncl_dist[ds.get_parent(i)];
                        if (dist_cur > roota->dtb_data.cluster_max_dist)
                            roota->dtb_data.cluster_max_dist = dist_cur;
                    }
                }
            }
            else {
                // nearer node first -> faster!
                kdtree_node_orderer<FLOAT, D, DISTANCE, NODE> sel(roota, rootb->left, rootb->right, (M>2));

                // prune nodes too far away if we have better candidates
                if (roota->dtb_data.cluster_max_dist > sel.nearer_dist) {
                    find_mst_next_dtb(roota, sel.nearer_node);
                    if (roota->dtb_data.cluster_max_dist > sel.farther_dist)
                        find_mst_next_dtb(roota, sel.farther_node);
                }


                // roota->dtb_data.cluster_max_dist updated above
            }
        }
        else {  // roota is not a leaf
            if (rootb->is_leaf()) {
                kdtree_node_orderer<FLOAT, D, DISTANCE, NODE> sel(rootb, roota->left, roota->right, (M>2));
                if (sel.nearer_node->dtb_data.cluster_max_dist > sel.nearer_dist)
                    find_mst_next_dtb(sel.nearer_node, rootb);
                if (sel.farther_node->dtb_data.cluster_max_dist > sel.farther_dist)  // separate if!
                    find_mst_next_dtb(sel.farther_node, rootb);
            }
            else {
                kdtree_node_orderer<FLOAT, D, DISTANCE, NODE> sel(roota->left, rootb->left, rootb->right, (M>2));
                if (roota->left->dtb_data.cluster_max_dist > sel.nearer_dist) {
                    find_mst_next_dtb(roota->left, sel.nearer_node);
                    if (roota->left->dtb_data.cluster_max_dist > sel.farther_dist)
                        find_mst_next_dtb(roota->left, sel.farther_node);
                }

                sel = kdtree_node_orderer<FLOAT, D, DISTANCE, NODE>(roota->right, rootb->left, rootb->right, (M>2));
                if (roota->right->dtb_data.cluster_max_dist > sel.nearer_dist) {
                    find_mst_next_dtb(roota->right, sel.nearer_node);
                    if (roota->right->dtb_data.cluster_max_dist > sel.farther_dist)
                        find_mst_next_dtb(roota->right, sel.farther_node);
                }
            }

            roota->dtb_data.cluster_max_dist = std::max(
                roota->left->dtb_data.cluster_max_dist,
                roota->right->dtb_data.cluster_max_dist
            );
        }
    }


    void find_mst_next_dtb()
    {
        find_mst_next_dtb(&this->nodes[0], &this->nodes[0]);
    }


    void find_nn_next_multi(NODE* curleaf)  // QTB
    {
    QUITEFASTMST_ASSERT(curleaf->cluster_repr == ds.get_parent(curleaf->idx_from));
        Py_ssize_t ds_find_i = curleaf->cluster_repr;

        // NOTE: assumption: no race condition/atomic read...
        FLOAT ncl_dist_cur = ncl_dist[ds_find_i];

        if (ncl_dist_cur <= curleaf->qtb_data.lastbest_dist) return;

        if (curleaf->qtb_data.lastbest_ind >= 0) {
            Py_ssize_t ds_find_j = ds.get_parent(curleaf->qtb_data.lastbest_ind);
            if (ds_find_i == ds_find_j)
                curleaf->qtb_data.lastbest_ind = -1;
        }

        if (curleaf->qtb_data.lastbest_ind < 0) {
            kdtree_nearest_outsider<FLOAT, D, DISTANCE, NODE> nn(
                this->data, (M>2)?(this->dcore.data()):NULL,
                M, ds.get_parents()
            );
            nn.find_multi(curleaf, &this->nodes[0], reset_nns?INFINITY:ncl_dist_cur);
            if (nn.get_nn_ind() >= 0) {
                curleaf->qtb_data.lastbest_ind  = nn.get_nn_ind();
                curleaf->qtb_data.lastbest_dist = nn.get_nn_dist();
                curleaf->qtb_data.lastbest_from = nn.get_nn_from();
            }
        }

        if (curleaf->qtb_data.lastbest_ind < 0) return;

        #if OPENMP_IS_ENABLED
        if (omp_nthreads > 1) omp_set_lock(&omp_lock);
        #endif
        if (curleaf->qtb_data.lastbest_dist < ncl_dist[ds_find_i]) {
            ncl_dist[ds_find_i] = curleaf->qtb_data.lastbest_dist;
            ncl_ind[ds_find_i]  = curleaf->qtb_data.lastbest_ind;
            ncl_from[ds_find_i] = curleaf->qtb_data.lastbest_from;
        }

        if (omp_nthreads == 1) {  // otherwise slightly worse performance...
            Py_ssize_t ds_find_j = ds.get_parent(curleaf->qtb_data.lastbest_ind);
            QUITEFASTMST_ASSERT(ds_find_i != ds_find_j);
            if (curleaf->qtb_data.lastbest_dist < ncl_dist[ds_find_j]) {
                ncl_dist[ds_find_j] = curleaf->qtb_data.lastbest_dist;
                ncl_ind[ds_find_j]  = curleaf->qtb_data.lastbest_from;
                ncl_from[ds_find_j] = curleaf->qtb_data.lastbest_ind;
            }
        }
        #if OPENMP_IS_ENABLED
        if (omp_nthreads > 1) omp_unset_lock(&omp_lock);
        #endif
    }


    void find_nn_next_single(Py_ssize_t i)  // STB, QTB
    {
        // Py_ssize_t i = (M<2)?u:ptperm[u];
        Py_ssize_t ds_find_i = ds.get_parent(i);

        // NOTE: assumption: no race condition/atomic read...
        FLOAT ncl_dist_cur = ncl_dist[ds_find_i];

        if (ncl_dist_cur <= lastbest_dist[i]) return;  // speeds up even for M==1

        if (lastbest_ind[i] < 0) {
            kdtree_nearest_outsider<FLOAT, D, DISTANCE, NODE> nn(
                this->data, (M>2)?(this->dcore.data()):NULL,
                M, ds.get_parents()
            );
            nn.find_single(i, &this->nodes[0], reset_nns?INFINITY:ncl_dist_cur);
            lastbest_ind[i] = nn.get_nn_ind();  // can be negative if best found >= ncl_dist_cur
            if (lastbest_ind[i] >= 0)
                lastbest_dist[i] = nn.get_nn_dist();
        }

        if (lastbest_ind[i] < 0) return;

        #if OPENMP_IS_ENABLED
        if (omp_nthreads > 1) omp_set_lock(&omp_lock);
        #endif

        if (lastbest_dist[i] < ncl_dist[ds_find_i]) {
            ncl_dist[ds_find_i] = lastbest_dist[i];
            ncl_ind[ds_find_i]  = lastbest_ind[i];
            ncl_from[ds_find_i] = i;
        }

        if (omp_nthreads == 1) {  // otherwise slightly worse performance...
            Py_ssize_t ds_find_j = ds.get_parent(lastbest_ind[i]);
            QUITEFASTMST_ASSERT(ds_find_i != ds_find_j);
            if (lastbest_dist[i] < ncl_dist[ds_find_j]) {
                ncl_dist[ds_find_j] = lastbest_dist[i];
                ncl_ind[ds_find_j]  = i;
                ncl_from[ds_find_j] = lastbest_ind[i];
            }
        }

        #if OPENMP_IS_ENABLED
        if (omp_nthreads > 1) omp_unset_lock(&omp_lock);
        #endif
    }


    void find_mst_next_qtb()
    {
        // find the point from another cluster that is closest to the i-th point
        // i.e., the nearest "alien"
        // go leaf-by-leaf
        #if OPENMP_IS_ENABLED
        #pragma omp parallel for schedule(static)
        #endif
        for (Py_ssize_t l=0; l<this->nleaves; ++l) {
            NODE* curleaf = leaves[l];

            if (curleaf->cluster_repr >= 0 && curleaf->idx_to - curleaf->idx_from > 1)  // all elems in the same cluster
            {
                find_nn_next_multi(curleaf);
            }
            else {
                for (Py_ssize_t i=curleaf->idx_from; i<curleaf->idx_to; ++i) {
                    find_nn_next_single(i);  // updates lastbest_dist[i] and ncl_dist[ds_find_i] if necessary
                }
            }
        }
    }


    void find_mst_next_stb()
    {
        // find the point from another cluster that is closest to the i-th point
        // i.e., the nearest "alien"
        #if OPENMP_IS_ENABLED
        #pragma omp parallel for schedule(static)
        #endif
        for (Py_ssize_t i=0; i<this->n; ++i) {
            find_nn_next_single(i);  // updates lastbest_dist[i] and ncl_dist[ds_find_i] if necessary
        }
    }


    void find_mst()
    {
        QUITEFASTMST_PROFILER_USE

        QUITEFASTMST_PROFILER_START
        // the 1st iteration: connect nearest neighbours with each other
        find_mst_first();
        QUITEFASTMST_PROFILER_STOP("find_mst_first")

        if (boruvka_variant == BORUVKA_DTB && M>2) {
            QUITEFASTMST_PROFILER_START
            setup_min_dcore();
            QUITEFASTMST_PROFILER_STOP("setup_min_dcore")
        }

        if (boruvka_variant == BORUVKA_QTB) {
            QUITEFASTMST_PROFILER_START
            setup_leaves();
            QUITEFASTMST_PROFILER_STOP("setup_leaves")
        }

        std::vector<Py_ssize_t> ds_parents(this->n);
        Py_ssize_t ds_k;

        while (tree_edges < this->n-1) {
            #if QUITEFASTMST_R
            Rcpp::checkUserInterrupt();  // throws an exception, not a longjmp
            #elif QUITEFASTMST_PYTHON
            if (PyErr_CheckSignals() != 0) throw std::runtime_error("signal caught");
            #endif

            tree_iter++;
            QUITEFASTMST_PROFILER_START

            ds_k = 0;
            for (Py_ssize_t i=0; i<this->n; ++i) {
                if (i == this->ds.find(i)) {
                    ncl_dist[i] = INFINITY;
                    ncl_ind[i]  = -1;
                    ncl_from[i] = -1;
                    ds_parents[ds_k++] = i;
                }
            }
            // now ds.find(i) == ds.get_parent(i) for all i

            update_nn_data();  // update lastbest_dist etc.

            update_node_data();  // reset cluster_max_dist and set up cluster_repr


            if (boruvka_variant == BORUVKA_DTB)
                find_mst_next_dtb();
            else if (boruvka_variant == BORUVKA_QTB) // TODO
                find_mst_next_qtb();
            else
                find_mst_next_stb();

            for (Py_ssize_t j=0; j<ds_k; ++j) {
                Py_ssize_t i = ds_parents[j];
                QUITEFASTMST_ASSERT(ncl_ind[i] >= 0 && ncl_ind[i] < this->n);
                if (ds.find(i) != ds.find(ncl_ind[i])) {
                    QUITEFASTMST_ASSERT(ncl_from[i] >= 0 && ncl_from[i] < this->n);
                    QUITEFASTMST_ASSERT(ds.find(i) == ds.find(ncl_from[i]));
                    tree_add(ncl_from[i], ncl_ind[i], ncl_dist[i]);
                }
            }

            QUITEFASTMST_PROFILER_STOP("find_mst iter #%d (tree_edges=%d)", (int)tree_iter, tree_edges)
        }
    }


public:
    kdtree_boruvka()
        : kdtree<FLOAT, D, DISTANCE, NODE>()
    {
        omp_nthreads = -1;
    }


    /**!
     * see fastmst.h for the description of the parameters,
     * no need to repeat that here
     */
    kdtree_boruvka(
        FLOAT* data, const Py_ssize_t n, const Py_ssize_t M=1,
        const Py_ssize_t max_leaf_size=16,
        const Py_ssize_t first_pass_max_brute_size=16,
        const FLOAT boruvka_variant=1.5,
        const FLOAT mutreach_adj=-INFINITY
    ) :
        kdtree<FLOAT, D, DISTANCE, NODE>(data, n, max_leaf_size), tree_edges(0), tree_iter(0),
        ds(n), ncl_dist(n), ncl_ind(n), ncl_from(n),
        first_pass_max_brute_size(first_pass_max_brute_size),
        mutreach_adj(mutreach_adj), M(M)
    {
        QUITEFASTMST_ASSERT(M>0);

        if (M >= 2) {
            dcore.resize(n);
            Mnn_dist.resize(n*(M-1));
            Mnn_ind.resize(n*(M-1));
        }

        lastbest_dist.resize(n);
        lastbest_ind.resize(n);

        if (boruvka_variant == 2.0)
            this->boruvka_variant = BORUVKA_DTB;
        else if (boruvka_variant == 1.0)
            this->boruvka_variant = BORUVKA_STB;
        else
            this->boruvka_variant = BORUVKA_QTB;  // 1.5 ;)

        reset_nns = (M<=2);  // plain Euclidean MST benefits from this


        #if OPENMP_IS_ENABLED
        omp_nthreads = Comp_get_max_threads();
        if (omp_nthreads > 1) omp_init_lock(&omp_lock);
        #else
        omp_nthreads = 1;
        #endif
    }


    ~kdtree_boruvka()
    {
        #if OPENMP_IS_ENABLED
        if (omp_nthreads > 1) omp_destroy_lock(&omp_lock);
        #endif
    }


    void mst(FLOAT* tree_dist, Py_ssize_t* tree_ind)
    {
        this->tree_dist = tree_dist;
        this->tree_ind  = tree_ind;

        if (ds.get_k() != (Py_ssize_t)this->n) ds.reset();
        tree_edges = 0;
        tree_iter = 0;

        for (Py_ssize_t i=0; i<this->n-1; ++i)     tree_dist[i] = INFINITY;
        for (Py_ssize_t i=0; i<2*(this->n-1); ++i) tree_ind[i]  = -1;

        // nodes is a deque...
        for (auto curnode = this->nodes.rbegin(); curnode != this->nodes.rend(); ++curnode)
            curnode->cluster_repr = -1;

        find_mst();
    }


    inline const FLOAT* get_Mnn_dist() const
    {
        QUITEFASTMST_ASSERT(M>1);
        return this->Mnn_dist.data();
    }

    inline const Py_ssize_t* get_Mnn_ind() const {
        QUITEFASTMST_ASSERT(M>1);
        return this->Mnn_ind.data();
    }

    inline const FLOAT* get_dcore() const {
        QUITEFASTMST_ASSERT(M>1);
        return this->dcore.data();
    }

    inline Py_ssize_t get_M() const { return this->M; }
};



/*!
 * Find a minimum spanning tree of X (in the tree)
 *
 * see _mst_euclid_kdtree
 *
 * @param tree a pre-built K-d tree containing n points
 * @param tree_dist [out] size n*k
 * @param tree_ind [out] size n*k
 * @param nn_dist [out] distances to M-1 nns of each point
 * @param nn_ind  [out] indexes of M-1 nns of each point
 */
template <typename FLOAT, Py_ssize_t D, typename DISTANCE, typename TREE>
void mst(
    TREE& tree,
    FLOAT* tree_dist,           // size n-1
    Py_ssize_t* tree_ind,       // size 2*(n-1),
    FLOAT* nn_dist=nullptr,     // size n*(M-1)
    Py_ssize_t* nn_ind=nullptr  // size n*(M-1)
) {
    tree.mst(tree_dist, tree_ind);

    Py_ssize_t n = tree.get_n();
    Py_ssize_t M = tree.get_M();
    const Py_ssize_t* perm = tree.get_perm();

    if (M > 1) {
        QUITEFASTMST_ASSERT(nn_dist);
        QUITEFASTMST_ASSERT(nn_ind);
        const FLOAT*      _nn_dist = tree.get_Mnn_dist();
        const Py_ssize_t* _nn_ind  = tree.get_Mnn_ind();

        for (Py_ssize_t i=0; i<n; ++i) {
            for (Py_ssize_t j=0; j<M-1; ++j) {
                nn_dist[perm[i]*(M-1)+j] = _nn_dist[i*(M-1)+j];
                nn_ind[perm[i]*(M-1)+j]  = perm[_nn_ind[i*(M-1)+j]];
            }
        }
    }

    for (Py_ssize_t i=0; i<n-1; ++i) {
        Py_ssize_t i1 = tree_ind[2*i+0];
        Py_ssize_t i2 = tree_ind[2*i+1];
        QUITEFASTMST_ASSERT(i1 != i2);
        QUITEFASTMST_ASSERT(i1 >= 0 && i1 < n);
        QUITEFASTMST_ASSERT(i2 >= 0 && i2 < n);
        tree_ind[2*i+0] = perm[i1];
        tree_ind[2*i+1] = perm[i2];
    }

    // the edges are not ordered, use Cmst_order
}



};  // namespace

#endif
