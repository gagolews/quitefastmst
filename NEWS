# Changelog

## To Do

*   [HELP NEEDED] [Python] Set up OpenMP on macOS.

*   Parallelise the K-d tree building procedure.

*   In the Borůvka algorithm based on K-d trees, apply the correction
    for ambiguity of mutual reachability distances (`mutreach_adj`)
    also when considering non-M first neighbours.

*   Extend the online documentation: Tutorials, benchmarks, definitions.


## 0.9.1 (2026-02-11)

*   [NEW FEATURE]  The `mutreach_leaves` argument to `mst_euclid`
    controls the postprocessing of tree leaves.

*   [BACKWARD INCOMPATIBILITY]  The definition of the mutual reachability
    distance has changed (for notational prudence).  Unlike in
    Campello et al.'s 2013 paper, now the core distance is the distance
    to the M-th nearest neighbour, not the (M-1)-th one (not including self).

*   [BACKWARD INCOMPATIBILITY]  The `mutreach_adj` argument to `mst_euclid`
    was removed.  Instead, the `mutreach_ties` argument is now available.
    It defaults to `"dist_min"` for (rough) compatibility with other packages.

*   [BUGFIX]  #3: SIGSEGV on duplicated inputs in `mst_euclid`
    with `algorithm="brute"` was fixed.


## 0.9.0 (2025-07-22)

*   [R]  Initial CRAN release.

*   [Python]  Initial PyPI release.
