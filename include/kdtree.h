#pragma once

// A space partitioning K-d tree, somewhat like described on
// https://en.wikipedia.org/wiki/K-d_tree
// Erik Wernersson 2025-2026

#include <stdint.h>
#include <stddef.h>

#define KDTREE_VERSION_MAJOR 0
#define KDTREE_VERSION_MINOR 2
#define KDTREE_VERSION_PATCH 0

#ifdef WIN32
#define PUB __declspec(dllexport)
#else
#define PUB __attribute__((visibility("default")))
#endif

struct pqheap;

// A node, called a leaf it is does not have any children.
// the nodes use a binary heap layout so there is no need to store
// pointer to the children nodes.
typedef struct {
    // numerical identifier of the node, i.e. where in T->nodes it can be found
    // todo: Remove this since it can be calculated based on the address if needed
    size_t id;
    // Location in T->X and T->id where the elements of the node are stored
    size_t offset;
    // Number of points of this node
    uint32_t n_points;
    // Dimension to split on, or set to ndim if to indicate a leaf node
    uint8_t split_dim;
    double pivot; // the location of the split along the split dimension
} kdtree_node_t;

typedef struct{
    uint32_t ndim; // Number of dimensions

    //
    // Per node / region data.
    //

    // Node k is stored in T->nodes[k],
    // and the corresponding bbx at T->boxes[k*T->ndim]
    size_t n_nodes_alloc; // Total number of nodes
    kdtree_node_t * nodes; /* Array of nodes, nodes[0] is the root */
    // Bounding boxes [minx, maxx,  miny, maxy,  ... ]
    double * boxes;

    //
    // Per point data
    //

    // Storage for coordinates, note: the order will be scrambled
    // along the tree construction, however orig_idx keeps track
    // of the original id of the points
    // The points for a given node N are stored in
    // T->X node->offset * T->ndim
    // and forwards
    // And the original location of those points at
    // T->OID + node->offset
    double * X;
    uint32_t * OID; // original id

    // Maximum number of points per leaf (i.e. end node)
    size_t max_leaf_size;
    size_t n_points; // Number of supplied points


    // Temporary buffer used during tree construction
    double * median_buffer;

    // at least ndim*sizeof(double) large
    double * point_buffer;

    //
    // State variables for querying the closest points, should not be
    // here really. Must have had a lazy day. TODO
    //
    struct pqheap * pq; // used for k-nearest queries
    int direct_path;
    // The latest query is stored internally to avoid an abundant
    // number of malloc/free. Can of course be copied by the caller. */
    size_t * result; // KN for storing idx of K neighbours
    size_t result_alloc; /* number of elements allocated for result */
} kdtree_t;

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
PUB size_t *
kdtree_query_knn(kdtree_t * T,
                 const double * Q,
                 size_t k);

// Find all points within some radius of Q
// Returns a newly allocated array of indexes of length nfound
// On failure: Returns NULL and sets nfound to 0
PUB size_t *
kdtree_query_radius(const kdtree_t * T,
                    const double * Q,
                    const double radius,
                    size_t * nfound);

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

// Find the index of the closest point
PUB size_t kdtree_query_closest(kdtree_t * T, double * X);

PUB void node_print_bbx(const kdtree_t * T, const kdtree_node_t * N);

// Make a shallow copy of a kd-tree for usage by another thread
PUB kdtree_t * kdtree_copy_shallow(kdtree_t * );

// Free a tree returned from kdtree_copy_shallow
PUB void kdtree_free_shallow(kdtree_t * T);

// Run some self-tests
PUB void kdtree_validate(kdtree_t * T);

PUB void kdtree_print_info(kdtree_t * T);


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
