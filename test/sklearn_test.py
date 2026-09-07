#!/bin/env python

import numpy as np
import time
import sys

from sklearn.neighbors import NearestNeighbors

def print_memory():
    import os
    pid = os.getpid()
    pidfile = f"/proc/{pid}/status"
    with open(pidfile) as f:
        for line in f:
            if "VmPeak" in line:
                print(line.strip())

def test(n, k, ndim=3):
    rng = np.random.default_rng()
    min_bounds = np.zeros(ndim)
    max_bounds = 1024*np.ones(ndim)
    X = rng.uniform(min_bounds, max_bounds, size=(n, ndim))
    t0 = time.perf_counter()
    nbrs = NearestNeighbors(n_neighbors=k, algorithm='kd_tree').fit(X)
    t1 = time.perf_counter()
    # all vs all
    distances, indices = nbrs.kneighbors(X)
    t2 = time.perf_counter()
    t_construct = 1000*(t1-t0)
    t_scan = 1000*(t2-t1)
    t_total = t_construct + t_scan
    print(f"| sklearn | {n} | {t_construct:.3f} | {t_scan:.3f} | {t_total:.3f} |")

if __name__ == '__main__':

    ndim = 3
    if len(sys.argv) > 1:
        ndim = int(sys.argv[1])

    print(f'ndim={ndim}')
    print("| method |       N | t_construct [ms] | t_query [ms] | t_total [ms] |")
    print("|--------|--------:|-----------------:|-------------:|-------------:|")
    k = 5
    n = 128
    while n < 2**21:
        n *= 2
        test(n, k, ndim=ndim)

    print(f"n = {n}, k = {k}")

    print_memory()
