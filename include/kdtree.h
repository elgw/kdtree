#pragma once
#include <stdint.h>

/* Will only work when KDTREE_DIM=3 unless modified */
#define KDTREE_DIM 3

struct pqheap;

/* A tree or a leaf */
typedef struct {
    size_t id; // number; // node number TODO not needed

    /* Bounding box, [minx, maxx, miny, maxy, minz, maxz] */
    double bbx[2*KDTREE_DIM];

    /* Tells where in XID that the data for this node can be found */
    size_t data_idx;

    /* The number of points associated with the node
     * which are found in T->XID */
    uint32_t n_points;

    /* If not a leaf, this tells how this node was split */
    uint8_t split_dim; /* is set to KDTREE_DIM for leafs */
    double pivot;
} kdtree_node_t;

typedef struct{
    /* Node allocation  */
    kdtree_node_t * nodes; /* Array of nodes, nodes[0] is the root */
    size_t n_nodes_alloc; // Total number of nodes

    /* Storage for coordinates and their indexes, [(x, y, z), id] */
    double * XID;

    /* Maximum number of points per leaf (i.e. end node) */
    size_t max_leaf_size;
    size_t n_points; // Number of supplied points

    /* Temporary buffer used during tree construction */
    double * median_buffer;

    /* State variables for queries */
    struct pqheap * pq; // needs to be here?
    int direct_path;
    /* The latest query is stored internally to avoid an abundant
       number of malloc/free. Can of course be copied by the caller. */
    size_t * result; // KN for storing idx of K neighbours
    size_t result_alloc; /* number of elements allocated for result */
} kdtree_t;

/* Construct a new tree based on the N points in X
 *
 * According to [1] a binsize (max_leaf_size) of 4-32 elements is optimal
 * regardless of k and the number of dimensions */
kdtree_t * kdtree_new(const double * X,
                      size_t N, size_t ndim, int binsize);


/* Frees all resources associated with a tree and sets the pointer to NULL */
void kdtree_free(kdtree_t ** _T);

/* Query one point for its k nearest neighbours.  The returned array
 * contains the index of k points, sorted according to the distance of
 * the points, with the closest point first.
 *
 * Important: The returned array is owned by the tree and should not
 * be freed. It will be re-used with the next call to kdtree_query_*
 */
size_t * kdtree_query_knn(kdtree_t * T, const double * Q, size_t k);

/* TODO
 * Find all points within some radius of Q
 * Returns a newly allocated array of indexes of length nfound
 * On failure: Returns NULL and sets nfound to 0
*/
size_t *
kdtree_query_radius(const kdtree_t * T,
                    const double * Q,
                    const double radius,
                    size_t * nfound);

/* TODO: Estimate the local density using a Gaussian symmetric kernel
   with sigma. Esentially this calls kdtree_query_radius and applies
   the KDE to the found points. */
double
kdtree_kde(const kdtree_t * T,
           const double * Q,
           double sigma);

/* Query the nQ points in Q for the k nearest neighbors
 *
 * The returned matrix is kxN elements large and should be freed by
 * the caller.
 */
size_t * kdtree_query_knn_multi(kdtree_t * T,
                                const double * Q, size_t nQ,
                                int k, int ntheads);

size_t kdtree_query_closest(kdtree_t * T, double * X);

void node_print_bbx(const kdtree_node_t * N);

/* Run some self-tests */
void kdtree_validate(kdtree_t * T);
