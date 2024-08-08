# k-d tree (but only the 3D case)

The [K-d tree](https://en.wikipedia.org/wiki/K-d_tree) is a fun data
structure, useful for finding k-nearest neighbours and neighbours
within some distance in point clouds. The benefits it provides
compared to brute force drops quickly with the number of dimensions as
you can read on the Wiki page.

This repo supports exactly what I need and nothing more, so the
functionality is quite minimal. I'd be happy if anyone else finds it
useful and can send me a bug report now and then, or even a pull
request :)

Usual warnings applies, use with caution!

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

- There is no parallel code for the tree construction at the moment
  although that would be possible to do. `kdtree_query_radius` and
  `kdtree_kde` should be thread safe. `kdtree_query_knn` is not tread safe
  but `kdtree_query_knn_multi` can be used to query multiple points in
  parallel.

- Can use GSL (`gsl_stats_median`) to find the pivot or the provided
  quick select implementation.

## Performance hints

Only using one thread. In this case "this" means `kdtree_query_knn`
from this repo.

For reference, there are also results from
`sklearn.neighbors.NearestNeighbors`, see `test_python.py` for the
test code. Sklearn is probably an interface to
[ckdtree](https://github.com/scipy/scipy/tree/main/scipy/spatial/ckdtree/src)
but that is just a hypothesis, not a fact. The comparison is not fair,
as comparisons seldom are, since the Python code stores the full
results (an Nxk matrix) while that is discarded in my code. For
sklearn, the memory measurement includes the whole Python environment,
not just the algorithm and the associated data.

Finding the k=5 nearest neighbours for each point among
N=1,000.


| Software | Tree construction |  Query | Total time |  VmPeak |
| -------- | ----------------- | ------ | ---------- | ------- |
| this     |            0.2 ms | 0.8 s  |     1.3 ms |    4 MB |
| sklearn  |            0.5 ms | 1.8 ms |     2.4 ms | 1502 MB |

N=5000, k = 5

| Software | Tree construction |  Query | Total time |  VmPeak |
| -------- | ----------------- | ------ | ---------- | ------- |
| this     |            1.3 ms | 4.9 s  |     6.2 ms |    4 MB |
| sklearn  |            1.7 ms | 9.9 ms |    11.5 ms | 1504 MB |

N=1,000,000, k = 5

| Software | Tree construction | Query | Total time |  VmPeak |
| -------- | ----------------- | ----- | ---------- | ------- |
| this     |             0.3 s | 2.0 s |      2.3 s |   87 MB |
| sklearn  |             0.8 s | 5.8 s |      6.5 s | 1695 MB |


N=1,000,000, **k = 100**

| Software | Tree construction | Query | Total time |  VmPeak |
| -------- | ----------------- | ----- | ---------- | ------- |
| this     |            0.3 s  |17.2 s |     17.5 s |   87 MB |
| sklearn  |            0.8 s  |35.8 s |     36.6 s | 4664 MB |

N = 100,000,000

| Software | Tree construction | Query | Total time |   VmPeak |
| -------- | ----------------- | ----- | ---------- | -------- |
| this     |              40 s | 483 s |      524 s |  8878 MB |
| sklearn  |             187 s | 968 s |     1155 s | 20580 MB |


## Current validation steps

More tests should be written. Especially to cover corner cases. The
current test pack includes:

- Should compile with zero warnings using `gcc -Wall -Wextra -pedantic
  -std=gnu11 -g3 -fanalyzer`.
- Passes the few tests in `kdtree_ut.c`, some of them are comparisons
  to brute force calculations.
- **valgrind** finds no issues.

## To Do
- [ ] For N dimensions. Remove `#define KDTREE_DIM 3` and make it a
      parameter. Write tests and do the small adjustments needed.
- [ ] Remove `<pthread.h>` from the main code for portability.
- [ ] Also fix so that `kdtree_query_knn` is thread safe (via a query
      object containing the per-thread data).
- [ ] More validation.
