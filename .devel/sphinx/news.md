# Changelog

## To Do

*   Parallelise the K-d tree building procedure.

*   In the Boruvka algorithm based on K-d trees, apply the correction
    for ambiguity of mutual reachability distances (`mutreach_adj`)
    also when considering non-M first neighbours

*   [Python] Set up OpenMP on macOS. (help needed)

*   Extend the online documentation: Tutorials, benchmarks, definitions.

*   PyPI build for Python 3.14




## 0.9.1 (2025-10-21)

*   [BACKWARD INCOMPATIBILITY]  The definition of the mutual reachability
    distance has changed.  Unlike in Campello et al.'s 2013 paper,
    now the core distance is the distance to the M-th nearest neighbour,
    not the (M-1)-th one (not including self).

*   [BACKWARD INCOMPATIBILITY]  The `mutreach_adj` argument to `mst_euclid`
    was removed.  Instead, the `mutreach_ties` argument is now available.


## 0.9.0 (2025-07-22)

*   [R] Initial CRAN release.

*   [Python] Initial PyPI release.
