import numpy as np
import sklearn.neighbors
import scipy.spatial.distance
import time
import gc
import quitefastmst
import sklearn.datasets
import pandas as pd
import sys

import os
if os.path.exists(".devel/benchmark_data"):
    path = ".devel/benchmark_data"
elif os.path.exists("benchmark_data"):
    path = "benchmark_data"
else:
    path = "../benchmark_data"

# NOTE: more tests in the genieclust package



def mst_check(X, **kwargs):
    n = X.shape[0]
    d = X.shape[1]

    for M in [0, 1, 2, 3, 5, n-1]:
        res = []
        for algo in ["brute", "single_kd_tree", "sesqui_kd_tree", "dual_kd_tree", "auto"]:
            if d > 20 and algo in ["single_kd_tree", "sesqui_kd_tree", "dual_kd_tree"]:
                continue

            print("    fastmst_%s (M=%d)   " % (algo, M), end="", file=sys.stderr, flush=True)
            t0 = time.time()
            res.append(quitefastmst.mst_euclid(X, algorithm=algo, M=M))
            print("%10.3fs" % (time.time()-t0,), file=sys.stderr)

            assert np.allclose(res[0][0].sum(), res[-1][0].sum())
            assert np.allclose(res[0][0], res[-1][0])

            assert np.all(np.bincount(res[-1][1].ravel(), minlength=X.shape[0]) >= 1)

            if M == 0:
                assert len(res[-1]) == 2
                assert np.all(res[0][1] == res[-1][1])  # is unique
            else:
                assert len(res[-1]) == 4
                assert np.allclose(res[0][2].sum(), res[-1][2].sum())
                assert np.allclose(res[0][2], res[-1][2])
                assert np.all(res[0][3] == res[-1][3])

    return True



def test_MST():
    for dataset in ["pathbased", "big_one", "h2mg_64_50"]:
        if dataset == "big_one":
            X = np.random.rand(1000, 20)
        else:
            X = np.loadtxt("%s/%s.data.gz" % (path,dataset), ndmin=2)

        # center X + scale (NOT: standardize!)
        X = (X-X.mean(axis=0))/X.std(axis=None, ddof=1)
        X += np.random.normal(0, 0.0001, X.shape)

        print(dataset, file=sys.stderr)
        mst_check(X)
        gc.collect()


def test_borderline():
    # test if does not crash

    np.random.seed(123)
    n = 3

    X = np.random.rand(n, 2)
    X = np.vstack((X,X))  # repeated

    Y = np.zeros((n, 2))

    for algo in ["brute", "single_kd_tree", "sesqui_kd_tree", "dual_kd_tree"]:
        for M in [0, 1, 2]:
            print("repeated", algo, M, file=sys.stderr)
            res = quitefastmst.mst_euclid(X, algorithm=algo, M=M)
            assert np.all(np.bincount(res[1].ravel(), minlength=X.shape[0]) >= 1)

            print("zeros", algo, M, file=sys.stderr)
            res = quitefastmst.mst_euclid(Y, algorithm=algo, M=M)
            assert np.all(np.bincount(res[1].ravel(), minlength=Y.shape[0]) >= 1)


if __name__ == "__main__":
    test_borderline()
    test_MST()
