#ifndef __kdtree_h__
#define __kdtree_h__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <gsl/gsl_statistics_double.h>

#include <pthread.h>

#define KDTREE_MAXDIM 4

struct pqheap;

typedef struct {
    size_t number; // node number
    size_t n_points; // number of points
    double geometry[2*KDTREE_MAXDIM]; // TODO: replaces xmin, xmax, ...
    double xmin; // domain of points
    double xmax;
    double ymin;
    double ymax;
    int var; // split variable
    double pivot; //
    int node_left; // children
    int node_right;
    size_t * idx; // Pointer to the local idx
    double * X; // Pointer to the local coordinates

} kdtree_node_t;

typedef struct{
    kdtree_node_t * root;
    size_t n_nodes; // Total number of nodes
    size_t next; // What node to write to during construction
    size_t ndim; // Number of dimensions
    size_t * idx_local; // For node idx storage
    double * X_local; // For node coordinate storage
    size_t next_idx; // Where to start writing idx
    size_t binsize;
    double * X; // pointer to supplied data
    size_t N; // Number of supplied points
    int direct_path;
    // For queries
    int k;
    struct pqheap * pq;
    size_t * KN; // for storing idx of K neighbours
    // FILE * log;
    // Temporary buffers
    double * median_buffer;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
} kdtree_t;

/* According to [1] a binsize of 4-32 elements is optimal
 * regardless of k and the number of dimensions */
kdtree_t * kdtree_new(const double * X,
                      size_t N, int binsize);


void kdtree_free(kdtree_t ** _T);
// Query one point. The returned array is k long and should not be freed
size_t * kdtree_query_knn(kdtree_t * T, const double * Q, int k);

/* TODO
 * Find all points within some radius of Q
 * Returns a newly allocated array of indexes of length nfound
 * On failure: Returns NULL and sets nfound to 0
*/
size_t *
kdtree_query_radius(const kdtree_t * T,
                    const double * Q,
                    const double radius,
                    int * nfound);

/* TODO:
   Estimate the local density using a Gaussian symmetric kernel
   with sigma. */
double
kdtree_kde(const kdtree_t * T,
           const double * Q,
           double sigma);

// Query nQ points for the k nearest neighbors
// The returned matrix is kxN elements large and should be freed by the caller.
size_t * kdtree_query_knn_multi(kdtree_t * T, const double * Q, size_t nQ, int k, int ntheads);

size_t kdtree_query_closest(kdtree_t * T, double * X);

#endif
