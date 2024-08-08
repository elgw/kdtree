from sklearn.neighbors import NearestNeighbors
import numpy as np
import time

n = 5000
k = 5

rng = np.random.default_rng()
X = rng.uniform([0, 0, 0], [1024, 1024, 1024], size=(n, 3))


t0 = time.perf_counter()
nbrs = NearestNeighbors(n_neighbors=k, algorithm='kd_tree').fit(X)
t1 = time.perf_counter()
# all vs all
distances, indices = nbrs.kneighbors(X)
t2 = time.perf_counter()
print(f"Took {t2-t0:02f} s using kd_tree. Construction {t1-t0:02f} Query: {t2-t1:02f}")


import os
pid = os.getpid()
pidfile = f"/proc/{pid}/status"
with open(pidfile) as f:
    for line in f:
        if "VmPeak" in line:
            print(line)
