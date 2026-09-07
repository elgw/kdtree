# K-d tree algorithm/data structure

Version 0.2.0 2026-09-17.

The [K-d tree](https://en.wikipedia.org/wiki/K-d_tree) is a fun data
structure, useful for finding k-nearest neighbours and neighbours
within some distance in point clouds. The benefits it provides
compared to brute force drops quickly with the number of dimensions as
you can read on the Wiki page.

This repo supports exactly what I need and nothing more, so the
functionality is quite minimal. I'd be happy if anyone else finds it
useful and can send me a bug report now and then, or even a pull
request :)

Notes:

- C99, no dependencies.
- Builds with gcc, clang and musl-gcc and zig cc under linux.
- The size of `libkdtree.so` is 26 KB.
- Usual warnings applies, use with caution!

## Usage
Below are examples of the supported methods:

``` C
#include <kdtree.h>
...
// X: N k-D points [k x N]
kdtree_t * T = kdtree_new(X, N, k, 20);

// Find the k nearest neighbours to Q [k x 1]
size_t * knn = kdtree_query_knn(T, Q, k);

// Find any point within a distance of radius to Q
size_t * idx = kdtree_query_radius(T, Q, radius, &n);

// Evaluate the point density under the point, using
// an isotropic Gaussian controlled by sigma.
double v = kdtree_kde(T, Q, sigma);

// When done
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
  although that would be possible to do. The query methods are not
  thread safe (but that is of course doable as well).

- Can use GSL (`gsl_stats_median`) to find the pivot or the provided
  quick select implementation.

## Performance hints

Finding the k=5 nearest neighbours for each point among
N 3D points (`./kdtree_ut --table1`) gives:

<details><summary>kdtree 2D</summary>

| method |       N | t_construct [ms] | t_query [ms] | t_total [ms] |
|--------|--------:|-----------------:|-------------:|-------------:|
| kdtree |     128 |            0.019 |        0.076 |        0.095 |
| kdtree |     256 |            0.023 |        0.193 |        0.215 |
| kdtree |     512 |            0.051 |        0.340 |        0.391 |
| kdtree |    1024 |            0.122 |        0.718 |        0.840 |
| kdtree |    2048 |            0.289 |        1.496 |        1.786 |
| kdtree |    4096 |            0.650 |        3.094 |        3.744 |
| kdtree |    8192 |            1.429 |        6.275 |        7.703 |
| kdtree |   16384 |            3.086 |       12.917 |       16.003 |
| kdtree |   32768 |            6.897 |       27.305 |       34.203 |
| kdtree |   65536 |           15.170 |       56.631 |       71.800 |
| kdtree |  131072 |           33.080 |      116.275 |      149.355 |
| kdtree |  262144 |           70.644 |      242.413 |      313.057 |
| kdtree |  524288 |          151.798 |      515.187 |      666.986 |
| kdtree | 1048576 |          314.865 |     1228.529 |     1543.394 |
| kdtree | 2097152 |          677.605 |     2556.722 |     3234.326 |

</details>

<details><summary>kdtree 3D</summary>

| method |       N | t_construct [ms] | t_query [ms] | t_total [ms] |
|--------|--------:|-----------------:|-------------:|-------------:|
| kdtree |     128 |            0.024 |        0.099 |        0.122 |
| kdtree |     256 |            0.024 |        0.218 |        0.242 |
| kdtree |     512 |            0.053 |        0.500 |        0.553 |
| kdtree |    1024 |            0.140 |        1.163 |        1.303 |
| kdtree |    2048 |            0.297 |        2.401 |        2.698 |
| kdtree |    4096 |            0.653 |        5.394 |        6.047 |
| kdtree |    8192 |            1.493 |       11.178 |       12.670 |
| kdtree |   16384 |            3.232 |       23.700 |       26.932 |
| kdtree |   32768 |            7.202 |       52.745 |       59.947 |
| kdtree |   65536 |           16.522 |      107.931 |      124.454 |
| kdtree |  131072 |           33.604 |      218.450 |      252.054 |
| kdtree |  262144 |           72.223 |      461.946 |      534.169 |
| kdtree |  524288 |          154.701 |     1067.012 |     1221.714 |
| kdtree | 1048576 |          337.677 |     2401.210 |     2738.887 |
| kdtree | 2097152 |          688.714 |     5695.521 |     6384.235 |

</details>

<details><summary>kdtree 7D</summary>

| method |       N | t_construct [ms] | t_query [ms] | t_total [ms] |
|--------|--------:|-----------------:|-------------:|-------------:|
| kdtree |     128 |            0.022 |        0.150 |        0.171 |
| kdtree |     256 |            0.029 |        0.473 |        0.502 |
| kdtree |     512 |            0.074 |        1.481 |        1.555 |
| kdtree |    1024 |            0.161 |        4.137 |        4.298 |
| kdtree |    2048 |            0.355 |       11.261 |       11.616 |
| kdtree |    4096 |            0.834 |       30.004 |       30.838 |
| kdtree |    8192 |            1.760 |       75.012 |       76.772 |
| kdtree |   16384 |            3.820 |      183.927 |      187.748 |
| kdtree |   32768 |            8.241 |      460.892 |      469.133 |
| kdtree |   65536 |           17.459 |     1027.477 |     1044.936 |
| kdtree |  131072 |           31.084 |     2189.055 |     2220.139 |
| kdtree |  262144 |           68.500 |     5616.061 |     5684.561 |
| kdtree |  524288 |          177.351 |    15701.070 |    15878.421 |
| kdtree | 1048576 |          317.223 |    39659.393 |    39976.617 |
| kdtree | 2097152 |          684.292 |    94619.149 |    95303.441 |

</details>

For reference, results from
`sklearn.neighbors.NearestNeighbors` are also given (see `test_python.py` for the
test code). Sklearn is probably an interface to
[ckdtree](https://github.com/scipy/scipy/tree/main/scipy/spatial/ckdtree/src)
but that is just a hypothesis, not a fact. The comparison is not fair,
comparisons seldom are, since the Python code stores the full
result (an Nxk matrix) at the end.

In short, these lines were used:

``` Python
NearestNeighbors(n_neighbors=k, algorithm='kd_tree').fit(X)
distances, indices = nbrs.kneighbors(X)
```

which gave:


<details><summary>sklearn 2D</summary>

``` shell
$ python test/sklearn_test.py 2
```

| method  |       N | t_construct [ms] | t_query [ms] | t_total [ms] |
|---------|--------:|-----------------:|-------------:|-------------:|
| sklearn |     256 |            0.420 |        0.694 |        1.114 |
| sklearn |     512 |            0.382 |        0.891 |        1.273 |
| sklearn |    1024 |            0.421 |        1.636 |        2.057 |
| sklearn |    2048 |            0.836 |        4.167 |        5.003 |
| sklearn |    4096 |            1.488 |        6.346 |        7.834 |
| sklearn |    8192 |            2.602 |       12.380 |       14.982 |
| sklearn |   16384 |            5.588 |       28.196 |       33.785 |
| sklearn |   32768 |           11.697 |       59.825 |       71.522 |
| sklearn |   65536 |           25.911 |      126.584 |      152.495 |
| sklearn |  131072 |           57.091 |      257.201 |      314.292 |
| sklearn |  262144 |          123.887 |      531.859 |      655.746 |
| sklearn |  524288 |          256.934 |     1268.977 |     1525.911 |
| sklearn | 1048576 |          647.956 |     3601.681 |     4249.637 |
| sklearn | 2097152 |         1633.978 |     8429.846 |    10063.824 |

</details>

<details><summary>sklearn 3D</summary>

``` shell
$ python test/sklearn_test.py 3
```

| method  |       N | t_construct [ms] | t_query [ms] | t_total [ms] |
|---------|--------:|-----------------:|-------------:|-------------:|
| sklearn |     256 |            0.398 |        0.745 |        1.143 |
| sklearn |     512 |            0.376 |        1.053 |        1.429 |
| sklearn |    1024 |            0.472 |        2.196 |        2.668 |
| sklearn |    2048 |            0.923 |        4.448 |        5.372 |
| sklearn |    4096 |            1.687 |        9.338 |       11.026 |
| sklearn |    8192 |            3.280 |       20.254 |       23.533 |
| sklearn |   16384 |            7.097 |       44.846 |       51.943 |
| sklearn |   32768 |           15.880 |       98.863 |      114.742 |
| sklearn |   65536 |           34.257 |      212.100 |      246.357 |
| sklearn |  131072 |           74.622 |      429.024 |      503.645 |
| sklearn |  262144 |          140.516 |      778.479 |      918.995 |
| sklearn |  524288 |          290.906 |     2555.740 |     2846.646 |
| sklearn | 1048576 |          842.556 |     6875.754 |     7718.310 |
| sklearn | 2097152 |         2110.700 |    15457.541 |    17568.241 |

</details>

<details><summary>sklearn 7D</summary>

``` shell
$ python test/sklearn_test.py 7
```

| method  |       N | t_construct [ms] | t_query [ms] | t_total [ms] |
|---------|--------:|-----------------:|-------------:|-------------:|
| sklearn |     256 |            0.456 |        1.278 |        1.734 |
| sklearn |     512 |            0.471 |        3.292 |        3.763 |
| sklearn |    1024 |            0.791 |        8.524 |        9.316 |
| sklearn |    2048 |            1.414 |       22.159 |       23.573 |
| sklearn |    4096 |            2.694 |       56.367 |       59.061 |
| sklearn |    8192 |            6.035 |      144.327 |      150.362 |
| sklearn |   16384 |           13.113 |      368.192 |      381.305 |
| sklearn |   32768 |           28.875 |      844.882 |      873.757 |
| sklearn |   65536 |           52.246 |     1837.930 |     1890.176 |
| sklearn |  131072 |          113.787 |     4639.780 |     4753.567 |
| sklearn |  262144 |          280.968 |    15538.331 |    15819.299 |
| sklearn |  524288 |          777.966 |    43915.266 |    44693.232 |
| sklearn | 1048576 |         1888.079 |   105688.990 |   107577.068 |
| sklearn | 2097152 |         4468.667 |   249584.014 |   254052.681 |

</details>


## Current validation steps

More tests should be written. Especially to cover corner cases. The
current test pack includes:

- Compile with zero warnings using `gcc -Wall -Wextra -pedantic
  -std=gnu11 -g3 -Og -fanalyzer`.
- `scan-build make -B` reports 0 issues.
- Passes the few tests in `kdtree_ut.c`, some of them are comparisons
  to brute force calculations.
- The tests suite can be run with either valgrind or
  `-fsantize=address` without any issues.
- If you plan to use it, please add some extra tests!

## Build/Install

Use the makefile to generate the test program.

To build and install the library, please use

``` shell
mkdir build
cd build
cmake ..
make
sudo make install
```

To include in your project with CMake, you could copy the files to a
subfolder, for example named `modules/`, and then add something like
this to your `CMakeLists.txt`:

``` CMake
  add_subdirectory("modules/kdtree/")
  target_include_directories(myAPP PUBLIC "modules/kdtree/include/")
  target_link_directories(myAPP PUBLIC "modules/kdtree/")
  target_link_libraries(myAPP kdtree)
```

## To Do

- [ ] Make it thread safe.
- [ ] More validation.

## History

0.2.0, 2026-09-07

- Improvement: should work with any number of dimensions
- New: Callback interface for detecting collisions.
- Bug fix: The box sphere intersection routine looked weird and is replaced.

0.1.1, 2025-04-05,

-  Added `kdtree_kde_mean` for mean shift algorithms.

0.1.0 2024-08-09.

- Hello world

## References

To read:

- [K-d tree page on Wikipedia](https://en.wikipedia.org/wiki/K-d_tree)
- [Blog post: Color quantization, minimizing variance, and k-d trees](https://www.crisluengo.net/archives/932/)

Implementations:

- Python: [sklearn.neighbors.KDTree](https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KDTree.html)
- Matlab: [knnsearch](https://se.mathworks.com/help/stats/knnsearch.html)

See also:

- [collide3](https://github.com/elgw/collide3) which provides
  collision detection for 3D points in a much more time-efficient way
  when the points are quite uniformly distributed in $`[-1,1]^3`$.
