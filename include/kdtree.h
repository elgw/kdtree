#pragma once

// A space partitioning K-d tree, somewhat like described on
// https://en.wikipedia.org/wiki/K-d_tree
// Erik Wernersson 2025-2026

#include <stdint.h>
#include <stddef.h>

#define KDTREE_VERSION_MAJOR 0
#define KDTREE_VERSION_MINOR 2
#define KDTREE_VERSION_PATCH 1

#ifdef WIN32
#define PUB __declspec(dllexport)
#else
#define PUB __attribute__((visibility("default")))
#endif

struct kdtree_node_struct;
typedef struct kdtree_node_struct kdtree_node_t;

struct kdtree_struct;
typedef struct kdtree_struct kdtree_t;

// Construct a new tree based on the N points stored in X
//
// binsize (or max_leaf_size) is an algorithmic parameter. According
// to [1] bin size a of 4-32 elements is optimal regardless of the
// number of dimensions
PUB kdtree_t *
kdtree_new(const double * X,
           uint32_t N, uint32_t ndim,
           int binsize);


// Frees all resources associated with a tree
PUB void kdtree_free(kdtree_t * T);

// Query one point for its k nearest neighbours.  The returned array
// contains the index of k points, sorted according to the distance of
// the points, with the closest point first.
//
// Important: The returned array is owned by the tree and should not
// be freed. It will be re-used with the next call to kdtree_query_
PUB uint32_t *
kdtree_query_knn(kdtree_t * T,
                 const double * Q,
                 uint32_t k);

// Find all points within some radius of Q
//
// Returns the number of found points, their indexes can be found in
// the results array.
//
// The logic of the result array is the same as with getline, i.e.,
// the result array grows as needed and might change address.  Free
// after final use. Can be NULL, however result_capacity must be a
// valid pointer.
PUB uint32_t
kdtree_query_radius(const kdtree_t * T,
                    const double * Q,
                    const double radius,
                    uint32_t ** result,
                    uint32_t * result_capacity);

//
// Function below are not equally well tested; and might be removed in
// the future
//

// Estimate the local density using non-normalized Gaussian symmetric
// kernel with a fixed sigma.
//
// G(x, sigma) = exp(-x^2 / (2*sigma^2))
//
// For a custom kernel please use kdtree_query_radius and calculate
// the KDE based on the found points.
//
// Cutoff: will use points up to sigma*cutoff away from the query
// point Q. A default value will be used if cutoff == -1.
PUB double
kdtree_kde(const kdtree_t * T,
           const double * Q,
           double sigma,
           double cutoff);

// Calculate the weighted mean position of a Gaussian KDE at the point
// Q as the following sum over all the neighbors
//
// m(x) = \frac{ \sum K(x_i-x) x_i }{\sum K(x_i -x)}
//
// returns: m -- the weighted mean position. Will be Q if no neighbours found.
//
// This is an ingredient of the mean shift algorithm. Please not that
// using a Gaussian kernel a quite large radius contributes to the
// kde, set cutoff to you liking. A cutoff of 0 means that it will be
// set automagically.
PUB void
kdtree_kde_mean(const kdtree_t *,
                const double * Q,
                double sigma,
                double cutoff,
                double * mean);

// Wanted: Expectation Maximization (EM) with Gaussian Mixture Model (GMM)
// void kdtree_emgmm(const kdtree_t * T,
// const gaussian ** G0,
// gaussian ** Gfinal);

// Make a shallow copy of a kd-tree for usage by another thread
PUB kdtree_t * kdtree_copy_shallow(kdtree_t * );

// Free a tree returned from kdtree_copy_shallow
PUB void kdtree_free_shallow(kdtree_t * T);

// Performs an all-vs-all collision test to detect points that
// are withing radius distance from each other.
//
// In the callback function, u and v refer to the points of
// the array used to construct the tree, i.e.,
// double * point_u = X + 3*u;
// double * point_v = X + 3*v;
// The callback is only called once per pair, i.e., if it
// is called with pair of indexes, (u,v), then it will not be called with
// (v, u). It will not be called for self collisions, i.e. u!=v.
typedef void (*kdtree_collide_cb)(uint32_t u, uint32_t v, void * data);
PUB void kdtree_collide(const kdtree_t * T, double radius,
                        kdtree_collide_cb cb, void * cb_data);
