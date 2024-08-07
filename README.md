# k-d tree

[K-d trees](https://en.wikipedia.org/wiki/K-d_tree) are fun data structures. Useful for finding k-nearest neighbours and neighbours within some distance in point clouds.

## Details:
- Using Hoare's partition scheme.
- Node layout a la [Eytzinger](https://arxiv.org/abs/1509.05053) which is the same as used in [binary heaps](https://en.wikipedia.org/wiki/Binary_heap)


For up to 4 dimensions (more if you change `KDTREE_MAXDIM` in `kdtree.h`).

Supported operations:
- k nearest neighbours
- all points within some radius
- KDE estimator with a Gaussian kernel.

## TODO
- [ ] Start with interface and write from scratch again :)

- [ ] Dense copy of the input array, which is later just sorted as in qsort. Three allocs/frees in total.

- [ ] My own median routine using quickselect



## Maybe
- [Implicit](https://en.wikipedia.org/wiki/Implicit_k-d_tree)
- Option to pass a list of pointers instead of just returning indexes (when I need it).
