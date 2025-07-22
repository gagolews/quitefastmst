# distutils: language=c++
# cython: boundscheck=False
# cython: cdivision=True
# cython: nonecheck=False
# cython: wraparound=False
# cython: language_level=3


"""
k-nearest neighbours and minimum spanning trees with respect to the Euclidean
metric or the thereon-based mutual reachability distances. The module provides
access to a quite fast implementation of K-d trees.

For best speed, consider building the package from sources
using, e.g., ``-O3 -march=native`` compiler flags and with OpenMP support on.
"""


# ############################################################################ #
#                                                                              #
#   Copyleft (C) 2025-2025, Marek Gagolewski <https://www.gagolewski.com>      #
#                                                                              #
#                                                                              #
#   This program is free software: you can redistribute it and/or modify       #
#   it under the terms of the GNU Affero General Public License                #
#   Version 3, 19 November 2007, published by the Free Software Foundation.    #
#   This program is distributed in the hope that it will be useful,            #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               #
#   GNU Affero General Public License Version 3 for more details.              #
#   You should have received a copy of the License along with this program.    #
#   If this is not the case, refer to <https://www.gnu.org/licenses/>.         #
#                                                                              #
# ############################################################################ #


import numpy as np
cimport numpy as np
np.import_array()
import os
import warnings

cimport libc.math
from libcpp cimport bool
from libcpp.vector cimport vector
from libc.math cimport INFINITY


ctypedef fused T:
    int
    long
    long long
    Py_ssize_t
    float
    double

ctypedef fused floatT:
    float
    double




cdef extern from "../src/c_common.h":

    int Comp_set_num_threads(int n_threads)
    int Comp_get_max_threads()


cdef extern from "../src/c_fastmst.h":

    void Ctree_order[T](Py_ssize_t n, T* tree_dist, Py_ssize_t* tree_ind)

    void Cknn1_euclid_kdtree[T](
        T* X, Py_ssize_t n, Py_ssize_t d, Py_ssize_t k,
        T* nn_dist, Py_ssize_t* nn_ind,
        Py_ssize_t max_leaf_size, bint squared,
        bint verbose
    ) except +

    void Cknn1_euclid_brute[T](
        T* X, Py_ssize_t n, Py_ssize_t d, Py_ssize_t k,
        T* nn_dist, Py_ssize_t* nn_ind, bint squared,
        bint verbose
    ) except +

    void Cknn2_euclid_kdtree[T](
        T* X, Py_ssize_t n, T* Y, Py_ssize_t m, Py_ssize_t d, Py_ssize_t k,
        T* nn_dist, Py_ssize_t* nn_ind,
        Py_ssize_t max_leaf_size, bint squared,
        bint verbose
    ) except +

    void Cknn2_euclid_brute[T](
        T* X, Py_ssize_t n, T* Y, Py_ssize_t m, Py_ssize_t d, Py_ssize_t k,
        T* nn_dist, Py_ssize_t* nn_ind, bint squared,
        bint verbose
    ) except +

    void Cmst_euclid_kdtree[T](
        T* X, Py_ssize_t n, Py_ssize_t d, Py_ssize_t M,
        T* mst_dist, Py_ssize_t* mst_ind,
        T* nn_dist, Py_ssize_t* nn_ind,
        Py_ssize_t max_leaf_size, Py_ssize_t first_pass_max_brute_size,
        T boruvka_variant, T mutreach_adj, bint verbose
    ) except +

    void Cmst_euclid_brute[T](
        T* X, Py_ssize_t n, Py_ssize_t d, Py_ssize_t M,
        T* mst_dist, Py_ssize_t* mst_ind,
        T* nn_dist, Py_ssize_t* nn_ind, T mutreach_adj, bint verbose
    ) except +


################################################################################

# cdef void _openmp_set_num_threads():
#     c_omp.Comp_set_num_threads(int(os.getenv("OMP_NUM_THREADS", -1)))

cpdef int omp_set_num_threads(int n_threads):
    """
    quitefastmst.omp_set_num_threads(n_threads)

    Sets the number of threads to be used in the subsequent calls to parallel
    code regions managed by OpenMP.

    The function has no effect if there is no built-in support for OpenMP.


    Parameters
    ----------

    n_threads : int


    Returns
    -------

    A single integer: the previous value of ``max_threads``
    or 1 if OpenMP is disabled.

    """
    return Comp_set_num_threads(n_threads)


cpdef int omp_get_max_threads():
    """
    quitefastmst.omp_get_max_threads()

    The function's name is confusing: it returns the maximal number
    of threads that will be used during the next call to a parallelised
    function, not the maximal number of threads possibly available.

    It there is no built-in support for OpenMP, 1 is always returned.


    Returns
    -------

    A single integer.


    """
    return Comp_get_max_threads()


################################################################################



cpdef tuple tree_order(const floatT[:] tree_dist, const Py_ssize_t[:,:] tree_ind):
    """
    quitefastmst.tree_order(tree_dist, tree_ind)


    Orders the edges of a graph (e.g., a spanning tree) wrt the weights
    increasingly, resolving ties if needed based on the points' IDs,
    i.e., the triples (dist, ind1, ind2) are sorted lexicographically.


    Parameters
    ----------

    tree_dist : c_contiguous ndarray, shape (m,)
        The m edges' weights

    tree_ind : c_contiguous ndarray, shape (m,2)
        The corresponding pairs of IDs of the incident nodes


    Returns
    -------

    pair : tuple
        A pair (tree_dist, tree_ind) after the ordering.

    """
    cdef Py_ssize_t m = tree_dist.shape[0]

    if tree_ind.shape[0] != m or tree_ind.shape[1] != 2:
        raise ValueError("incorrect shape of tree_ind")

    cdef np.ndarray[floatT] tree_dist_ret = np.asarray(tree_dist, order="C", copy=True)
    cdef np.ndarray[Py_ssize_t,ndim=2] tree_ind_ret = np.asarray(tree_ind, order="C", copy=True)

    Ctree_order(m, &tree_dist_ret[0], &tree_ind_ret[0,0])

    return tree_dist_ret, tree_ind_ret




cpdef tuple knn_euclid(
    const floatT[:,:] X,
    Py_ssize_t k=1,
    const floatT[:,:] Y=None,
    str algorithm="auto",
    int max_leaf_size=0,
    bint squared=False,
    bint verbose=False
):
    """
    quitefastmst.knn_euclid(X, k=1, Y=None, algorithm="auto", max_leaf_size=0, squared=False, verbose=False)

    Euclidean Nearest Neighbours

    If `Y` is ``None``, then the function determines the first `k` nearest
    neighbours of each point in `X` with respect to the Euclidean distance.
    It is assumed that each query point is not its own neighbour.

    Otherwise, for each point in `Y`, this function determines the `k` nearest
    points thereto from `X`.


    Notes
    -----

    The implemented algorithms, see the ``algorithm`` parameter, assume that
    `k` is rather small; say, `k ≤ 20`.

    Our implementation of K-d trees [1]_ has been quite optimised; amongst
    others, it has good locality of reference (at the cost of making a
    copy of the input dataset), features the sliding midpoint (midrange) rule
    suggested in [2]_, node pruning strategies inspired by some ideas
    from [3]_, and a couple of further tuneups proposed by the current author.
    Still, it is well-known that K-d trees perform well only in spaces of low
    intrinsic dimensionality.  Thus, due to the so-called curse of
    dimensionality, for high `d`, the brute-force algorithm is preferred.

    The number of threads used is controlled via the ``OMP_NUM_THREADS``
    environment variable or via ``quitefastmst.omp_set_num_threads``
    at runtime.  For best speed, consider building the package from sources
    using, e.g., ``-O3 -march=native`` compiler flags.


    References
    ----------

    .. [1]
        J.L. Bentley, Multidimensional binary search trees used for associative
        searching, Communications of the ACM 18(9), 509–517, 1975,
        https://doi.org/10.1145/361002.361007

    .. [2]
        S. Maneewongvatana, D.M. Mount, It's okay to be skinny, if your friends
        are fat, 4th CGC Workshop on Computational Geometry, 1999

    .. [3]
        N. Sample, M. Haines, M. Arnold, T. Purcell, Optimizing search
        strategies in K-d Trees, 5th WSES/IEEE Conf. on Circuits, Systems,
        Communications & Computers (CSCC'01), 2001


    Parameters
    ----------

    X : matrix, shape `(n,d)`
        the `n` input points in :math:`\\mathbb{R}^d` (the "database")
    k : int `< n`
        requested number of nearest neighbours
        (should be rather small, say, `≤ 20`)
    Y : None or an ndarray, shape `(m,d)`
        the "query points"; note that setting ``Y=X``, contrary to ``Y=None``,
        will include the query points themselves amongst their own neighbours
    algorithm : {``"auto"``, ``"kd_tree"``, ``"brute"``}, default ``"auto"``
        K-d trees can only be used for `d` between 2 and 20.
        ``"auto"`` selects ``"kd_tree"`` in low-dimensional spaces
    max_leaf_size : int
        maximal number of points in the K-d tree leaves;
        smaller leaves use more memory, yet are not necessarily faster;
        use ``0`` to select the default value, currently set to 32.
    squared : False
        whether the output ``dist`` should use the squared Euclidean distance
    verbose: bool
        whether to print diagnostic messages


    Returns
    -------

    pair : tuple
        A pair ``(dist, ind)`` representing the `k`-nearest neighbour graph, where:

            dist : a c_contiguous ndarray, shape `(n,k)` or `(m,k)`
                ``dist[i,:]`` is sorted nondecreasingly for all `i`,
                ``dist[i,j]`` gives the weight of the edge ``{i, ind[i,j]}``,
                i.e., the distance between the `i`-th point and its `j`-th NN

            ind : a c_contiguous ndarray of the same shape
                ``ind[i,j]`` is the index (between `0` and `n-1`)
                of the `j`-th nearest neighbour of `i`
    """
    cdef Py_ssize_t n = X.shape[0]
    cdef Py_ssize_t d = X.shape[1]
    cdef Py_ssize_t m

    if n < 1 or d <= 1: raise ValueError("X is ill-shaped");
    if k < 1: raise ValueError("k must be >= 1")

    if algorithm == "auto":
        algorithm = "kd_tree" if 2 <= d <= 20 else "brute"

    if algorithm == "kd_tree":
        if not 2 <= d <= 20:
            raise ValueError("K-d trees can only be used for 2 <= d <= 20")

        if max_leaf_size == 0: max_leaf_size = 32  # the current default

        if max_leaf_size <= 0:
            raise ValueError("max_leaf_size must be positive")

        use_kdtree = True

    elif algorithm == "brute":
        use_kdtree = False
    else:
        raise ValueError("invalid 'algorithm'")

    cdef np.ndarray[Py_ssize_t,ndim=2] ind
    cdef np.ndarray[floatT,ndim=2]     dist

    cdef np.ndarray[floatT,ndim=2] X2
    X2 = np.asarray(X, order="C", copy=True)  # destroyable

    cdef np.ndarray[floatT,ndim=2] Y2

    if Y is None:
        if k >= n:
            raise ValueError("too many neighbours requested")

        ind  = np.empty((n, k), dtype=np.intp)
        dist = np.empty((n, k), dtype=np.float32 if floatT is float else np.float64)

        if use_kdtree:
            Cknn1_euclid_kdtree(&X2[0,0], n, d, k, &dist[0,0], &ind[0,0], max_leaf_size, squared, verbose)
        else:
            Cknn1_euclid_brute(&X2[0,0], n, d, k, &dist[0,0], &ind[0,0], squared, verbose)
    else:
        if k > n:
            raise ValueError("too many neighbours requested")

        Y2 = np.asarray(Y, order="C", copy=True, dtype=np.float32 if floatT is float else np.float64)

        if Y2.ndim == 1: Y2 = Y2.reshape(1, -1)

        if Y2.ndim != 2 or Y2.shape[1] != d:
            raise ValueError("Y's dimensionality does not match that of X")

        m = Y2.shape[0]

        ind  = np.empty((m, k), dtype=np.intp)
        dist = np.empty((m, k), dtype=np.float32 if floatT is float else np.float64)

        if use_kdtree:
            Cknn2_euclid_kdtree(&X2[0,0], n, &Y2[0,0], m, d, k, &dist[0,0], &ind[0,0], max_leaf_size, squared, verbose)
        else:
            Cknn2_euclid_brute(&X2[0,0], n, &Y2[0,0], m, d, k, &dist[0,0], &ind[0,0], squared, verbose)

    return dist, ind



cpdef tuple mst_euclid(
    const floatT[:,:] X,
    Py_ssize_t M=1,
    str algorithm="auto",
    int max_leaf_size=0,
    int first_pass_max_brute_size=0,
    floatT mutreach_adj=-1.00000011920928955078125,
    bint verbose=False
):
    """
    quitefastmst.mst_euclid(
        X, M=1, algorithm="auto", max_leaf_size=0,
        first_pass_max_brute_size=0, mutreach_adj=-1.00000011920928955078125,
        verbose=False
    )

    Euclidean and Mutual Reachability Minimum Spanning Trees


    The function determines the/a(\*) minimum spanning tree (MST) of a set
    of `n` points, i.e., an acyclic undirected connected graph whose:
    vertices represent the points,
    edges are weighted by the distances between point pairs,
    and have minimal total weight.

    MSTs have many uses in, amongst others, topological data analysis
    (clustering, density estimation, dimensionality reduction,
    outlier detection, etc.).

    For `M ≤ 2`, we get a spanning tree that minimises the sum of Euclidean
    distances between the points, i.e., the classic Euclidean minimum
    spanning tree (EMST).  If `M = 2`, the function additionally returns
    the distance to each point's nearest neighbour.

    If `M > 2`, the spanning tree is the smallest wrt the degree-`M`
    mutual reachability distance [9]_ given by
    :math:`d_M(i, j)=\\max\\{ c_M(i), c_M(j), d(i, j)\\}`, where :math:`d(i,j)`
    is the Euclidean distance between the `i`-th and the `j`-th point,
    and :math:`c_M(i)` is the `i`-th `M`-core distance defined as the distance
    between the `i`-th point and its `(M-1)`-th nearest neighbour
    (not including the query point itself).



    Notes
    -----

    (\*) We note that if there are many pairs of equidistant points,
    there can be many minimum spanning trees. In particular, it is likely
    that there are point pairs with the same mutual reachability distances.

    To make the definition less ambiguous (albeit with no guarantees),
    internally, the brute-force algorithm relies on the adjusted distance:
    :math:`d_M(i, j)=\\max\\{c_M(i), c_M(j), d(i, j)\\}+\\varepsilon d(i, j)` or
    :math:`d_M(i, j)=\\max\\{c_M(i), c_M(j), d(i, j)\\}-\\varepsilon \\min\\{c_M(i), c_M(j)\\}`,
    where ε is close to 0.
    `|mutreach_adj| < 1` selects the former formula (`ε=mutreach_adj`)
    whilst `1 < |mutreach_adj| < 2` chooses the latter (`ε=mutreach_adj±1`).

    For the K-d tree-based methods, on the other hand, `mutreach_adj`
    indicates the preference towards connecting to farther/closer
    points wrt the original metric or having smaller/larger core distances
    if a point `i` has multiple nearest-neighbour candidates `j'`, `j''` with
    :math:`c_M(i) \geq \\max\\{d(i, j'),  c_M(j')\\}` and
    :math:`c_M(i) \geq \\max\\{d(i, j''), c_M(j'')\\}`.
    Generally, the smaller the `mutreach_adj`, the more leaves there should
    be in the tree (note that there are only four types of adjustments, though).

    The implemented algorithms, see the `algorithm` parameter, assume that
    `M` is rather small; say, `M ≤ 20`.

    Our implementation of K-d trees [6]_ has been quite optimised; amongst
    others, it has good locality of reference (at the cost of making a
    copy of the input dataset), features the sliding midpoint (midrange) rule
    suggested in [7]_, node pruning strategies inspired by some ideas
    from [8]_, and a couple of further tuneups proposed by the current author.

    The "single-tree" version of the Borůvka algorithm is naively
    parallelisable: in every iteration, it seeks each point's nearest "alien",
    i.e., the nearest point thereto from another cluster.
    The "dual-tree" Borůvka version of the algorithm is, in principle, based
    on [5]_.  As far as our implementation is concerned, the dual-tree approach
    is often only faster in 2- and 3-dimensional spaces, for `M ≤ 2`, and in
    a single-threaded setting.  For another (approximate) adaptation
    of the dual-tree algorithm to the mutual reachability distance, see [11]_.

    The "sesqui-tree" variant (by the current author) is a mixture of the two
    approaches:  it compares leaves against the full tree.  It is usually
    faster than the single- and dual-tree methods in very low dimensional
    spaces and usually not much slower than the single-tree variant otherwise.

    Nevertheless, it is well-known that K-d trees perform well only in spaces
    of low intrinsic dimensionality (a.k.a. the "curse").  For high `d`,
    the "brute-force" algorithm is recommended.  Here, we provided a
    parallelised [2]_ version of the Jarník [1]_ (a.k.a.
    Prim [3]_ or Dijkstra) algorithm, where the distances are computed
    on the fly (only once for `M ≤ 2`).

    The number of threads is controlled via the ``OMP_NUM_THREADS``
    environment variable or via ``quitefastmst.omp_set_num_threads``
    at runtime. For best speed, consider building the package from sources
    using, e.g., ``-O3 -march=native`` compiler flags.



    References
    ----------

    .. [1]
        V. Jarník, O jistém problému minimálním,
        Práce Moravské Přírodovědecké Společnosti 6, 1930, 57–63

    .. [2]
        C.F. Olson, Parallel algorithms for hierarchical clustering,
        Parallel Computing 21(8), 1995, 1313–1325,
        https://doi.org/10.1016/0167-8191(95)00017-I

    .. [3]
        R. Prim, Shortest connection networks and some generalizations,
        The Bell System Technical Journal 36(6), 1957, 1389–1401,
        https://doi.org/10.1002/j.1538-7305.1957.tb01515.x

    .. [4]
        O. Borůvka, O jistém problému minimálním,
        Práce Moravské Přírodovědecké Společnosti 3, 1926, 37–58

    .. [5]
        W.B. March, R. Parikshit, A.G. Gray, Fast Euclidean minimum spanning
        tree: Algorithm, analysis, and applications, Proc. 16th ACM SIGKDD Intl.
        Conf. Knowledge Discovery and Data Mining (KDD '10), 2010, 603–612

    .. [6]
        J.L. Bentley, Multidimensional binary search trees used for associative
        searching, Communications of the ACM 18(9), 509–517, 1975,
        https://doi.org/10.1145/361002.361007

    .. [7]
        S. Maneewongvatana, D.M. Mount, It's okay to be skinny, if your friends
        are fat, The 4th CGC Workshop on Computational Geometry, 1999

    .. [8]
        N. Sample, M. Haines, M. Arnold, T. Purcell, Optimizing search
        strategies in K-d Trees, 5th WSES/IEEE Conf. on Circuits, Systems,
        Communications & Computers (CSCC'01), 2001

    .. [9]
        R.J.G.B. Campello, D. Moulavi, J. Sander, Density-based clustering based
        on hierarchical density estimates, Lecture Notes in Computer Science 7819,
        2013, 160–172, https://doi.org/10.1007/978-3-642-37456-2_14

    .. [10]
        R.J.G.B. Campello, D. Moulavi, A. Zimek, J. Sander, Hierarchical density
        estimates for data clustering, visualization, and outlier detection,
        ACM Transactions on Knowledge Discovery from Data (TKDD) 10(1),
        2015, 1–51, https://doi.org/10.1145/2733381

    .. [11]
        L. McInnes, J. Healy, Accelerated hierarchical density-based
        clustering, IEEE Intl. Conf. Data Mining Workshops (ICMDW), 2017, 33–42,
        https://doi.org/10.1109/ICDMW.2017.12


    Parameters
    ----------

    X : matrix, shape `(n,d)`
        the `n` input points in :math:`\\mathbb{R}^d`
    M : int `< n`
        the degree of the mutual reachability distance (should be rather small,
        say, `≤ 20`). `M ≤ 2` denotes the ordinary Euclidean distance
    algorithm : {``"auto"``, ``"single_kd_tree"``, ``"sesqui_kd_tree"``, ``"dual_kd_tree"``, ``"brute"``}, default ``"auto"``
        K-d trees can only be used for `d` between 2 and 20 only.
        ``"auto"`` selects ``"sesqui_kd_tree"`` for `d ≤ 20`.
        ``"brute"`` is used otherwise
    max_leaf_size : int
        maximal number of points in the K-d tree leaves;
        smaller leaves use more memory, yet are not necessarily faster;
        use ``0`` to select the default value, currently set to 32 for the
        single-tree and sesqui-tree and 8 for the dual-tree Borůvka algorithm
    first_pass_max_brute_size : int
        minimal number of points in a node to treat it as a leaf (unless
        it actually is a leaf) in the first iteration of the algorithm;
        use ``0`` to select the default value, currently set to 32
    mutreach_adj : float
        adjustment for mutual reachability distance ambiguity (for M>2)
        whose fractional part should be close to 0:
        values in `(-1,0)` prefer connecting to farther NNs,
        values in `(0, 1)` fall for closer NNs (which is what many other
        implementations provide), values in `(-2,-1)` prefer connecting to
        points with smaller core distances, values in `(1, 2)` favour larger
        core distances; see above for more details
    verbose: bool
        whether to print diagnostic messages


    Returns
    -------

    tuple
        If `M = 1`, a pair ``(mst_dist, mst_index)`` defining the `n-1` edges
        of the computed spanning tree is returned:

            mst_dist : an array of length `(n-1)`
                ``mst_dist[i]`` gives the weight of the `i`-th edge

            mst_index : a matrix with `n-1` rows and `2` columns
                ``{mst_index[i,0], mst_index[i,1]}`` defines the `i`-th edge
                of the tree

        The tree edges are ordered w.r.t. weights nondecreasingly, and then by
        the indexes (lexicographic ordering of the ``(weight, index1, index2)``
        triples).  For each `i`, it holds ``mst_index[i,0]<mst_index[i,1]``.

        For `M > 1`, we additionally get:

            nn_dist : an `n` by `M-1` matrix
                it gives the distances between
                each point and its `M-1` nearest neighbours

            nn_index : a matrix of the same shape
                it provides the corresponding indexes of the neighbours
    """

    cdef Py_ssize_t n = X.shape[0]
    cdef Py_ssize_t d = X.shape[1]

    if n < 1 or d <= 1: raise ValueError("X is ill-shaped");
    if M < 1 or M > n-1: raise ValueError("incorrect M")

    cdef floatT boruvka_variant = 1.5
    cdef bool use_kdtree = True

    if algorithm == "auto":
        if 2 <= d <= 20:
            #if d <= 3:
            algorithm = "sesqui_kd_tree"
            #else:
            #    algorithm = "single_kd_tree"
        else:
            algorithm = "brute"

    if algorithm in ("single_kd_tree", "sesqui_kd_tree", "dual_kd_tree"):
        if not 2 <= d <= 20:
            raise ValueError("K-d trees can only be used for 2 <= d <= 20")

        use_kdtree = True

        if algorithm == "single_kd_tree":
            if max_leaf_size == 0:
                max_leaf_size = 32  # the current default
            if first_pass_max_brute_size == 0:
                first_pass_max_brute_size = 32  # the current default
            boruvka_variant = 1.0
        elif algorithm == "sesqui_kd_tree":
            if max_leaf_size == 0:
                max_leaf_size = 32  # the current default
            if first_pass_max_brute_size == 0:
                first_pass_max_brute_size = 32  # the current default
            boruvka_variant = 1.5
        elif algorithm == "dual_kd_tree":
            if max_leaf_size == 0:
                max_leaf_size = 8  # the current default
            if first_pass_max_brute_size == 0:
                first_pass_max_brute_size = 32  # the current default
            boruvka_variant = 2.0

        if max_leaf_size <= 0:
            raise ValueError("max_leaf_size must be positive")
        if first_pass_max_brute_size <= 0:
            raise ValueError("first_pass_max_brute_size must be positive")

    elif algorithm == "brute":
        use_kdtree = False
    else:
        raise ValueError("invalid 'algorithm'")

    cdef np.ndarray[Py_ssize_t,ndim=2] mst_ind  = np.empty((n-1, 2),
        dtype=np.intp)
    cdef np.ndarray[floatT] mst_dist = np.empty(n-1,
        dtype=np.float32 if floatT is float else np.float64)

    cdef np.ndarray[floatT,ndim=2] X2
    X2 = np.asarray(X, order="C", copy=True)  # destroyable

    cdef np.ndarray[Py_ssize_t,ndim=2] nn_ind
    cdef np.ndarray[floatT,ndim=2]     nn_dist
    if M > 1:
        nn_ind  = np.empty((n, M-1), dtype=np.intp)
        nn_dist = np.empty((n, M-1), dtype=np.float32 if floatT is float else np.float64)

    if use_kdtree:
        Cmst_euclid_kdtree(
            &X2[0,0], n, d, M,
            &mst_dist[0], &mst_ind[0,0],
            <floatT*>(0) if M==1 else &nn_dist[0,0],
            <Py_ssize_t*>(0) if M==1 else &nn_ind[0,0],
            max_leaf_size, first_pass_max_brute_size,
            boruvka_variant, mutreach_adj, verbose
        )
    else:
        Cmst_euclid_brute(
            &X2[0,0], n, d, M,
            &mst_dist[0], &mst_ind[0,0],
            <floatT*>(0) if M==1 else &nn_dist[0,0],
            <Py_ssize_t*>(0) if M==1 else &nn_ind[0,0],
            mutreach_adj,
            verbose
        )

    if M == 1:
        return mst_dist, mst_ind
    else:
        return mst_dist, mst_ind, nn_dist, nn_ind
