



# Benchmarks

::::{note}
This section is a work in progress.
Below we present some preliminary results.
::::



## Platform

All timings were performed on a mid-tier laptop with an
Intel Core i7-1355U (AVX2-VNNI) (12 threads, 2 performance cores, 8 efficient cores) CPU with 12MB cache
and 32 GB RAM, running GNU/Linux 6.15.5-1-default #1 SMP PREEMPT_DYNAMIC (478c062) x86_64 openSUSE Tumbleweed.

Python 3.13.5 [gcc] OpenBLAS 0.3.29;
packages:
**quitefastmst** 0.9.0,
**numba** 0.61.2,
**cython** 3.1.2,
**numpy** 2.2.6,
**scipy** 1.15.3,
**sklearn** 1.7.0,
**[hdbscan](https://github.com/scikit-learn-contrib/hdbscan)** {cite}`hdbscanpkg` 0.8.40,
**[fast_hdbscan](https://github.com/TutteInstitute/fast_hdbscan/)**+ {cite}`fasthdbscan` 0.2.2 (compiled via numba),
**rpy2** 3.6.1.

R version 4.5.1 (2025-06-13) Rblas;
R packages (via **rpy2**):
**[mlpack](https://github.com/mlpack/mlpack)**+ {cite}`mlpack` 4.6.2.

Other:
**[wangyiqiu_hdbscan](https://github.com/wangyiqiu/hdbscan)**+ e120dcf.

*+* means that a package was compiled with `-O3 -march=native`
using gcc 15.1.1 20250626



## Experiment Setup

For different $n$ and $d$, $n$ points in $\mathbb{R}^d$ were generated,
with each coordinate sampled independently from the standard normal distribution.
This is a pretty "difficult" setting for the spatial search data structures
as data are very "unstructured".  Most real-world datasets are more well-behaving.

We tested smoothing parameters $M=1$ (standard Euclidean MSTs) and $M=10$.

Number of threads was set to 1 or 10.  Note that **mlpack** does not support
parallel processing and only implements standard EMSTs.
**[wangyiqiu_hdbscan](https://github.com/wangyiqiu/hdbscan)**
required some trickery to be run from the terminal – it does not work out
of the box.

For each $n$, $d$, and $M$, three runs were performed and the smallest
timing was recorded.  Moreover, there was a pre-run on a subset of the
input dataset to ensure **fast_hdbscan** got compiled on-the-fly by **numba**.
**hdbscan_kdtree** was not run for larger $n$s.

See the relevant Python scripts:
[`perf_mst_202506.py`](https://github.com/gagolews/quitefastmst/blob/0b29075589475223f4eb43a16204b82df3a82cde/.devel/benchmarks/perf_mst_202506.py),
[`perf_mst_202506_defs`.py](https://github.com/gagolews/quitefastmst/blob/0b29075589475223f4eb43a16204b82df3a82cde/.devel/benchmarks/perf_mst_202506_defs.py)



## Results

All timings are given in seconds.
Less is better.





### 1 thread

|       n |   d |   M |   quitefastmst |   wangyiqiu |   fasthdbscan |   mlpack |   hdbscan |
|--------:|----:|----:|---------------:|------------:|--------------:|---------:|----------:|
|  250000 |   2 |   1 |           0.26 |        0.64 |          3.90 |     1.15 |      6.06 |
|  250000 |   2 |  10 |           0.39 |        1.41 |          2.26 |       |      6.22 |
|  250000 |   5 |   1 |           1.26 |        6.61 |         19.80 |    10.77 |     37.28 |
|  250000 |   5 |  10 |           1.50 |       19.29 |          5.67 |       |     26.70 |
| 2500000 |   2 |   1 |           2.98 |        8.33 |         83.73 |       |        |
| 2500000 |   2 |  10 |           4.35 |       17.93 |         45.79 |       |        |
| 2500000 |   5 |   1 |          14.21 |       67.87 |        460.17 |       |        |
| 2500000 |   5 |  10 |          16.59 |      217.54 |        141.18 |       |        | 



### 10 threads

|       n |   d |   M |   quitefastmst |   wangyiqiu |   fasthdbscan |   hdbscan |
|--------:|----:|----:|---------------:|------------:|--------------:|----------:|
|  250000 |   2 |   1 |           0.19 |        0.20 |          1.02 |      1.73 |
|  250000 |   2 |  10 |           0.18 |        0.39 |          0.60 |      1.52 |
|  250000 |   5 |   1 |           0.31 |        1.32 |          4.78 |     28.86 |
|  250000 |   5 |  10 |           0.38 |        4.41 |          1.52 |     14.74 |
| 2500000 |   2 |   1 |           2.08 |        2.65 |         25.18 |        |
| 2500000 |   2 |  10 |           2.03 |        4.92 |         14.36 |        |
| 2500000 |   5 |   1 |           4.41 |       15.81 |        149.27 |        |
| 2500000 |   5 |  10 |           5.19 |       61.80 |         50.12 |        | 


Thanks.
