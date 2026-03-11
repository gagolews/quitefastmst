# *quitefastmst*: Euclidean and Mutual Reachability Minimum Spanning Trees


**Keywords**: Euclidean minimum spanning tree, MST, EMST,
mutual reachability distance, nearest neighbours, k-nn, k-d tree,
Borůvka, Prim, Jarník, Kruskal, Genie, HDBSCAN\*, DBSCAN,
clustering, outlier detection.


Package **features**:

* [Euclidean Minimum Spanning Trees](https://en.wikipedia.org/wiki/Euclidean_minimum_spanning_tree)
    using single-, sesqui-, and dual-tree Borůvka algorithms, which are quite fast
    in spaces of low intrinsic dimensionality,

* minimum spanning trees with respect to mutual reachability distances based
    on the Euclidean metric (used in the definition of the HDBSCAN\* algorithm;
    see {cite}`hdbscan`),

* Euclidean nearest neighbours with nicely-optimised K-d trees,

* relatively fast fallback algorithms for spaces of higher dimensionality,

* supports multiprocessing via OpenMP (on selected platforms).


Possible applications in topological data analysis:
clustering ([HDBSCAN\*](https://hdbscan.readthedocs.io/en/latest/index.html),
[Lumbermark](https://lumbermark.gagolewski.com/),
[Genie](https://genieclust.gagolewski.com/), Single linkage, etc.),
outlier detection ([Deadwood](https://deadwood.gagolewski.com/)),
density estimation, dimensionality reduction,
and many more.


## Contributing

**quitefastmst** is distributed under the open source GNU AGPL v3 license.
Its source code can be downloaded from [GitHub](https://github.com/gagolews/quitefastmst).

The Python version is available from [PyPI](https://pypi.org/project/quitefastmst).
The R version can be fetched from [CRAN](https://CRAN.R-project.org/package=quitefastmst).

The core functionality is implemented in the form of a C++ library.
It can thus be easily adapted for use in other environments.
New contributions are welcome, e.g., Julia, Matlab/GNU Octave wrappers.


**Author and Maintainer**: [Marek Gagolewski](https://www.gagolewski.com/)



::::{toctree}
:maxdepth: 2
:caption: quitefastmst
:hidden:

About <self>
Author <https://www.gagolewski.com/>
Source Code (GitHub) <https://github.com/gagolews/quitefastmst>
Bug Tracker and Feature Suggestions <https://github.com/gagolews/quitefastmst/issues>
PyPI Entry <https://pypi.org/project/quitefastmst>
CRAN Entry <https://CRAN.R-project.org/package=quitefastmst>
::::



::::{toctree}
:maxdepth: 1
:caption: Documentation
:hidden:

weave/definitions
weave/benchmarks
weave/python
pythonapi
weave/r
rapi
::::


::::{toctree}
:maxdepth: 1
:caption: Other
:hidden:

deadwood <https://deadwood.gagolewski.com/>
lumbermark <https://lumbermark.gagolewski.com/>
genieclust <https://genieclust.gagolewski.com/>
Clustering Benchmarks <https://clustering-benchmarks.gagolewski.com/>
Minimalist Data Wrangling in Python <https://datawranglingpy.gagolewski.com/>
Deep R Programming <https://deepr.gagolewski.com/>
news
z_bibliography
::::


<!--
Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
-->
