import numpy as np
import sklearn.neighbors
import scipy.spatial.distance
import time
import gc
import quitefastmst

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

    for M in [1, 2, 3, 5, n-1]:
        res = []
        for algo in ["auto", "brute", "single_kd_tree", "sesqui_kd_tree", "dual_kd_tree"]:
            if d > 20 and algo in ["single_kd_tree", "sesqui_kd_tree", "dual_kd_tree"]:
                continue

            t0 = time.time()
            res.append(quitefastmst.mst_euclid(X, algorithm=algo, M=M))
            print("    fastmst_%s (M=%d)   %10.3fs" % (algo, M, time.time()-t0,))

            assert np.allclose(res[0][0].sum(), res[-1][0].sum())
            assert np.allclose(res[0][0], res[-1][0])

            if M == 1:
                assert len(res[-1]) == 2
                assert np.all(res[0][1] == res[-1][1])  # is unique
            else:
                assert len(res[-1]) == 4
                assert np.allclose(res[0][2].sum(), res[-1][2].sum())
                assert np.allclose(res[0][2], res[-1][2])
                assert np.all(res[0][3] == res[-1][3])

    return True



def test_MST():
    for dataset in ["big_one", "pathbased", "h2mg_64_50"]:
        if dataset == "big_one":
            X = np.random.rand(1000, 20)
        else:
            X = np.loadtxt("%s/%s.data.gz" % (path,dataset), ndmin=2)

        # center X + scale (NOT: standardize!)
        X = (X-X.mean(axis=0))/X.std(axis=None, ddof=1)
        X += np.random.normal(0, 0.0001, X.shape)

        print(dataset)
        mst_check(X)
        gc.collect()


if __name__ == "__main__":
    test_MST()
