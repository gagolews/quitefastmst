<a href="https://quitefastmst.gagolewski.com"><img src="https://www.gagolewski.com/_static/img/quitefastmst.png" align="right" height="128" width="128" /></a>
# [*quitefastmst*](https://quitefastmst.gagolewski.com/) Package for R and Python

## Euclidean and Mutual Reachability Minimum Spanning Trees


![quitefastmst for Python](https://github.com/gagolews/quitefastmst/workflows/quitefastmst%20for%20Python/badge.svg)
![quitefastmst for R](https://github.com/gagolews/quitefastmst/workflows/quitefastmst%20for%20R/badge.svg)

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
    see Campello, Moulavi, Sander, 2013),

* Euclidean nearest neighbours with nicely-optimised K-d trees,

* relatively fast fallback algorithms for spaces of higher dimensionality,

* supports multiprocessing via OpenMP (on selected platforms).


Refer to the package **homepage** at <https://quitefastmst.gagolewski.com/>
for the reference manual, tutorials, examples, and benchmarks.

**Author and maintainer**: [Marek Gagolewski](https://www.gagolewski.com/)


Possible applications in topological data analysis:
clustering ([HDBSCAN\*](https://hdbscan.readthedocs.io/en/latest/index.html),
[Genie](https://genieclust.gagolewski.com/), Lumbermark, Single linkage, etc.),
density estimation, dimensionality reduction,
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




### Other

The core functionality is implemented in the form of a C++ library.
It can thus be easily adapted for use in other projects.

New contributions are welcome, e.g., Julia, Matlab/GNU Octave wrappers.



## License

Copyright (C) 2025–2025 Marek Gagolewski <https://www.gagolewski.com/>

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU Affero General Public License Version 3,
19 November 2007, published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero
General Public License Version 3 for more details. You should have
received a copy of the License along with this program. If not, see
<https://www.gnu.org/licenses/>.



## References

Borůvka, O., O jistém problému minimálním,
*Práce Moravské Přírodovědecké Společnosti* **3**, 1926, 37–58.

Bentley, J.L., Multidimensional binary search trees used for associative
searching, *Communications of the ACM* **18**(9), 509–517, 1975,
[DOI: 10.1145/361002.361007](https://doi.org/10.1145/361002.361007).

Campello, R.J.G.B., Moulavi, D., Zimek, A., Sander, J., Hierarchical
density estimates for data clustering, visualization, and outlier detection,
*ACM Transactions on Knowledge Discovery from Data (TKDD)* **10**(1),
2015, 1–51, [DOI: 10.1145/2733381](https://doi.org/10.1145/2733381).

Campello, R.J.G.B., Moulavi, D., Sander, J.,
Density-based clustering based on hierarchical density estimates,
*Lecture Notes in Computer Science* **7819**, 2013, 160–172.
[DOI: 10.1007/978-3-642-37456-2_14](https://doi.org/10.1007/978-3-642-37456-2_14).

Gagolewski, M., Cena, A., Bartoszuk, M., Brzozowski, L.,
Clustering with minimum spanning trees: How good can it be?,
*Journal of Classification* **42**, 2025, 90–112.
[DOI: 10.1007/s00357-024-09483-1](https://doi.org/10.1007/s00357-024-09483-1).

Gagolewski, M., A framework for benchmarking clustering algorithms,
*SoftwareX* **20**, 2022, 101270.
[DOI: 10.1016/j.softx.2022.101270](https://doi.org/10.1016/j.softx.2022.101270).
<https://clustering-benchmarks.gagolewski.com/>.

Jarník, V., O jistém problému minimálním,
*Práce Moravské Přírodovědecké Společnosti* **6**, 1930, 57–63.

Maneewongvatana, S., Mount, D.M., It's okay to be skinny, if your friends
are fat, *The 4th CGC Workshop on Computational Geometry*, 1999.

March, W.B., Parikshit, R., Gray, A.G., Fast Euclidean minimum spanning
tree: Algorithm, analysis, and applications,
*Proc. 16th ACM SIGKDD Intl. Conf. Knowledge Discovery and Data Mining (KDD '10)*,
2010, 603–612.

Olson C.F., Parallel algorithms for hierarchical clustering,
*Parallel Computing* **21**(8), 1995, 1313–1325.

McInnes, L., Healy, J., Accelerated hierarchical density-based
clustering, *IEEE Intl. Conf. Data Mining Workshops (ICMDW)*, 2017, 33–42,
[DOI: 10.1109/ICDMW.2017.12](https://doi.org/10.1109/ICDMW.2017.12).

Prim, R., Shortest connection networks and some generalizations,
*The Bell System Technical Journal* **36**(6), 1957, 1389–1401.

Sample, N., Haines, M., Arnold, M., Purcell, T.,
Optimizing search strategies in K-d Trees,
*5th WSES/IEEE Conf. on Circuits, Systems, Communications & Computers (CSCC'01)*,
2001.


See **quitefastmst**'s [homepage](https://quitefastmst.gagolewski.com/)
for more references.
