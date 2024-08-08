# k-d tree but only 3D

[K-d trees](https://en.wikipedia.org/wiki/K-d_tree) are fun data
structures that are useful for finding k-nearest neighbours and
neighbours within some distance in point clouds. The can of course
also be used to construct kernel density estimators. Although they can
be implemented for high dimensions, the benefit compared to brute
force drops quickly with the number of dimensions.

This repo is used for some of my other projects and currently supports
exactly what I need and nothing more. I'd be happy if anyone else finds
it useful and can send me a bug report now and then, or even a pull
request :)

Not bullet tested, so please use something else for what is important!

## Usage
Below are examples of the supported methods:

``` C++
#include <kdtree.h>
...
/* X: N 3D points [3 x N] */
kdtree_t * T = kdtree_new(X, N, 20);

/* Find the k nearest neighbours to Q [3 x 1] */
size_t * knn = kdtree_query_knn(T, Q, k);

/* Find any point within a distance of radius to Q */
size_t * idx = kdtree_query_radius(T, Q, radius, &n);

/* Evaluate the point density under the point, using
 * an isotropic Gaussian controlled by sigma.        */
double v = kdtree_kde(T, Q, sigma);

/* When done */
kdtree_free(T);
```

see `kdtree.h` for the complete function signatures and some
documentation. Look in `kdtree_ut.c` for complete usage examples.

## Details
- Data partitioning using Hoare's scheme, typically used in
  quicksort and quickselect.
- For finding the k nearest neighbours the candidates are put in a
  priority queue, implemented by a binary heap.
- The memory layout of the nodes has a big impact on performance and
  memory usage. The code in this repo use the
  [Eytzinger](https://arxiv.org/abs/1509.05053) layout, which is the
  same as used in [binary
  heaps](https://en.wikipedia.org/wiki/Binary_heap). On the positive
  side this give a good memory locality and fast queries. The major
  downside is that we need to decide upfront how deep the tree should
  be, which means that the memory usage (for the tree, excluding the
  data points) will grow in steps of approximately 2 when the number
  of points passes some boundaries.
- There is no parallel code for the tree construction at the moment,
  although that would be possible to do. `kdtree_query_radius` and
  `kdtree_kde` are thread safe. `kdtree_query_knn` is not tread safe
  but `kdtree_query_knn_multi` can be used to query multiple points in
  parallel.
- Uses quickselect for median finding. Optionally GSL can be used for
  this. Performance seems equal but I would trust the GSL library more
  than my code in this case.

## Performance

Finding the 5 nearest neighbours for each point among N=1,000. Only
using one thread. In this case "this" means `kdtree_query_knn` from
this repo.

| Software | Tree construction |  Query | Total time |  VmPeak |
| -------- | ----------------- | ------ | ---------- | ------- |
| this     |            0.2 ms | 0.8 s  |     1.3 ms |    4 MB |
| sklearn  |            0.5 ms | 1.8 ms |     2.4 ms | 1502 MB |

N=5000

| Software | Tree construction |  Query | Total time |  VmPeak |
| -------- | ----------------- | ------ | ---------- | ------- |
| this     |            1.3 ms | 4.9 s  |     6.2 ms |    4 MB |
| sklearn  |            1.7 ms | 9.9 ms |    11.5 ms | 1504 MB |

N=1,000,000:

| Software | Tree construction | Query | Total time |  VmPeak |
| -------- | ----------------- | ----- | ---------- | ------- |
| this     |             0.3 s | 2.0 s |      2.3 s |   87 MB |
| sklearn  |             0.8 s | 5.8 s |      6.5 s | 1695 MB |

N=10,000,000:

| Software | Tree construction | Query | Total time |  VmPeak |
| -------- | ----------------- | ----- | ---------- | ------- |
| this     |               3 s |  37 s |       40 s |  792 MB |
| sklearn  |              13 s |  80 s |       93 s | 3419 MB |

N = 100,000,000

| Software | Tree construction | Query | Total time |   VmPeak |
| -------- | ----------------- | ----- | ---------- | -------- |
| this     |              40 s | 483 s |      524 s |  8878 MB |
| sklearn  |             187 s | 968 s |     1155 s | 20580 MB |



For the python code, see `test_python.py`. Sklearn is proably built
upon
[ckdtree](https://github.com/scipy/scipy/tree/main/scipy/spatial/ckdtree/src)
but that is just a hypothesis, not a fact.

## TODO
- [ ] For N dimensions. Remove `#define KDTREE_DIM 3` and make it a
      parameter. Write tests and do the small adjustments needed.
- [ ] Remove `<pthread.h>` from the main code for portability.
- [ ] Also fix so that `kdtree_query_knn` is thread safe.

## Maybe
- [Implicit](https://en.wikipedia.org/wiki/Implicit_k-d_tree)
- Option to pass a list of pointers instead of just returning indexes (when I need it).

- [ ] write some test image [https://github.com/skeeto/bmp/blob/master/test.c]
