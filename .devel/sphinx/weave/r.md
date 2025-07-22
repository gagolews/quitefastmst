



# R Examples

*To learn more about R, check out my open-access textbook*
[Deep R Programming](https://deepr.gagolewski.com/).





## How to Install


To install the package from [CRAN](https://CRAN.R-project.org/package=quitefastmst), call:


```r
install.packages("quitefastmst")
```



For best performance, advanced users will benefit from compiling the package
from sources:


```r
Sys.setenv(CXX_DEFS="-O3 -march=native")  # for gcc and clang
install.packages("quitefastmst", type="source")
```



## Basic Use


::::{note}
This section is a work in progress.
In the meantime, take a look at the documentation of the
[mst_euclid](../rapi/mst_euclid) and
[knn_euclid](../rapi/knn_euclid) functions.
::::

