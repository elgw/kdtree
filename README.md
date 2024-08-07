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
kdtree_free(&T);
assert(T = NULL);
```

See `kdtree.h` for the complete function signatures and some
documentation. Look in `kdtree_ut.c` for usage examples.

## Details:
- Partitioning the data using Hoare's scheme.
- Using [Eytzinger](https://arxiv.org/abs/1509.05053) for the node
  layout in memory which is the same as used in [binary
  heaps](https://en.wikipedia.org/wiki/Binary_heap)


Supported operations:
- k nearest neighbours
- all points within some radius
- KDE estimator with a Gaussian kernel.


- [ ] My own median routine using quickselect



## Maybe
- [Implicit](https://en.wikipedia.org/wiki/Implicit_k-d_tree)
- Option to pass a list of pointers instead of just returning indexes (when I need it).

- [ ] write some test image [https://github.com/skeeto/bmp/blob/master/test.c]
