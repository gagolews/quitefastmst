# *quitefastmst*: Euclidean and Mutual Reachability Minimum Spanning Trees


**Keywords**: Euclidean minimum spanning tree, MST, EMST,
mutual reachability distance, nearest neighbours, k-nn, k-d tree,
Boruvka, Prim, Jarnik, Kruskal, Genie, HDBSCAN\*, DBSCAN,
clustering, outlier detection.


Package **features**:

* [Euclidean Minimum Spanning Trees](https://en.wikipedia.org/wiki/Euclidean_minimum_spanning_tree)
    using single-, sesqui-, and dual-tree Borůvka algorithms – quite fast
    in spaces of low intrinsic dimensionality,

* Minimum spanning trees with respect to mutual reachability distances based
    on the Euclidean metric (used in the definition of the HDBSCAN\* algorithm;
    see {cite}`hdbscan`),

* Euclidean nearest neighbours with nicely-optimised K-d trees,

* relatively fast fallback algorithms for spaces of higher dimensionality,

* supports multiprocessing via OpenMP (on selected platforms).


Possible applications in data analysis:
clustering ([HDBSCAN\*](https://hdbscan.readthedocs.io/en/latest/index.html),
[Genie](https://genieclust.gagolewski.com/), Lumbermark, Single linkage, etc.),
classification and regression (k-nearest neighbours),
outlier and noise point detection, and many more.




## How to Install

### Python Version

To install from [PyPI](https://pypi.org/project/quitefastmst), call:

```bash
pip3 install quitefastmst  # python3 -m pip install quitefastmst
```

*To learn more about Python, check out my open-access textbook*
[Minimalist Data Wrangling in Python](https://datawranglingpy.gagolewski.com/).



For best performance, advanced users will benefit from compiling the package
from sources:

```bash
CPPFLAGS="-O3 -march=native" pip3 install quitefastmst --force --no-binary="quitefastmst"
```

🚧 TO DO (help needed): How to enable OpenMP support on macOS/Darwin in `setup.py`?



### R Version

To install from [CRAN](https://CRAN.R-project.org/package=quitefastmst), call:

```r
install.packages("quitefastmst")
```

*To learn more about R, check out my open-access textbook*
[Deep R Programming](https://deepr.gagolewski.com/).




For best performance, advanced users will benefit from compiling the package
from sources:

```r
Sys.setenv(CXX_DEFS="-O3 -march=native")  # for gcc and clang
install.packages("quitefastmst", type="source")
```




## Contributing

*quitefastmst* is distributed under the open source GNU AGPL v3 license
and can be downloaded from [GitHub](https://github.com/gagolews/quitefastmst).
The core functionality is implemented in the form of a C++ library,
so it may be adapted to new projects relatively easily:
any valuable contributions are welcome (Julia or Matlab bindings, etc.).


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
:caption: Reference Manual
:hidden:

pythonapi
rapi
::::


::::{toctree}
:maxdepth: 1
:caption: Other
:hidden:

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
