# k-d tree but only 3D

[K-d trees](https://en.wikipedia.org/wiki/K-d_tree) are fun data
structures. Useful for finding k-nearest neighbours and neighbours
within some distance in point clouds.

## Usage

``` C++
#include <kdtree.h>
...
/* X: N 3D points [3 x N] */
kdtree_t * T = kdtree_new(X, N, 20);

/* Find the k nearest neighbours to Q [3 x 1] */
size_t * knn = kdtree_query_knn(T, Q, k);

/* Find any point within a distance of radius to */
size_t * idx = kdtree_query_radius(T, Q, radius, &n);

/* When done */
kdtree_free(T);
```

See `kdtree.h` for the complete function signatures and some
documentation. Look in `kdtree_ut.c` for usage examples.

## Details:
- Partitioning the data using Hoare's scheme.

- The memory layout of the nodes has a big impact on performance and
  memory usage. The code in this repo use the
  [Eytzinger](https://arxiv.org/abs/1509.05053) layout, which is the
  same as used in [binary
  heaps](https://en.wikipedia.org/wiki/Binary_heap). On the positive
  side this give a good memory locality and fast queries. The major
  downside is that we need to decide upfront how deep the tree should
  be, which means that the memory usage (for the tree, excluding the
  data points) will grow in steps of approximately 2.


Supported operations:
- k nearest neighbours
- all points within some radius
- KDE estimator with a Gaussian kernel.


- [ ] My own median routine using quickselect



## Maybe
- [Implicit](https://en.wikipedia.org/wiki/Implicit_k-d_tree)
- Option to pass a list of pointers instead of just returning indexes (when I need it).

- [ ] write some test image [https://github.com/skeeto/bmp/blob/master/test.c]


## Performance

Finding the 5 nearest neighbours for each point among N=10,000,000:


| Software | Tree construction | Query | Total time |  VmPeak |
| -------- | ----------------- | ----- | ---------- | ------- |
| this     |               3 s |  34 s |       37 s |  795 MB |
| sklearn  |              13 s |  80 s |       93 s | 3419 MB |


For the python code, see `test_python.py`.
