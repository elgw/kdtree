#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "kdtree.h"
#include "pqheap.h"
#include "quickselect.h"

// To enable costly checks
// #define KDTREE_DEBUG

typedef int64_t i64;
typedef uint32_t u32;

// The type of the index arrays as well as the pointers
// to them
typedef uint32_t kdtree_index;

// A node, called a leaf it is does not have any children.
// the layout of the nodes is of binary heap type so there is no need to store
// pointers to the children nodes.
struct kdtree_node_struct {
    // numerical identifier of the node, i.e. where in T->nodes it can be found
    // Not strictly needed since it can be calculated based on the address
    kdtree_index id;
    // Location in T->X and T->id where the elements of the node are stored
    kdtree_index point_offset;
    // Number of points of this node
    kdtree_index n_point;
    // Dimension to split on, or set to ndim if to indicate a leaf node
    u32 split_dim;
    double pivot; // the location of the split along the split dimension
};

struct kdtree_struct {
    u32 ndim; // Number of dimensions

    //
    // Per node / region data.
    //

    // Node k is stored in T->nodes[k],
    // and the corresponding bbx at T->boxes[k*T->ndim]
    kdtree_index n_nodes_alloc; // Total number of nodes
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
    // T->X node->point_offset * T->ndim
    // and forwards
    // And the original location of those points at
    // T->OID + node->point_offset
    double * X;
    kdtree_index * OID; // original id

    // Maximum number of points per leaf (i.e. end node)
    kdtree_index max_leaf_size;
    kdtree_index n_point; // Number of supplied points

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
    kdtree_index * result; // KN for storing idx of K neighbours
    kdtree_index result_alloc; /* number of elements allocated for result */
};

// Resolve the index of the right child based on the index of the
// parent, following a Eytzinger scheme
static kdtree_index
node_right_child_id(kdtree_index node_id)
{
    return 2*node_id+2;
}

static kdtree_index
node_left_child_id(kdtree_index node_id)
{
    return 2*node_id+1;
}

static kdtree_index
sizeof_bbx(kdtree_index ndim) {
    return 2*ndim*sizeof(double);
}

// A node is marked as final/leaf when the split_dim is impossibly high
static void
node_set_final(const kdtree_t * T, kdtree_node_t * node)
{
    node->split_dim = T->ndim;
    return;
}

static int
node_is_final(const kdtree_t * T, const kdtree_node_t * node)
{
    return (node->split_dim == T->ndim);
}

static void
swap_doubles(double * restrict X, double * restrict Y, kdtree_index ndim)
{
    double T[ndim];
    memcpy(T, X, // T = X
           ndim*sizeof(double));
    memcpy(X, Y, // X = Y
           ndim*sizeof(double));
    memcpy(Y, T, // Y = T (= copy of input X)
           ndim*sizeof(double));
    return;
}

/*  Hoare's partition scheme for vectors where the partitioning is
 *  performed based on a single dimension or index.  Adopted from
 *  arch/24/03/11_quickselect
 *
 * X : the data of size [ndim x n] which will be partitioned.
 * ID : the data of size [n] will be partitioned in the same way
 *
 * vim: the index or dimension of the pivot
 *
 * Returns:
 * nLow: the number of vectors where v[vdim] <= pivot
 * nHigh: the number of vectors where v[vdim] > pivot
 */
static void
partition_vectors(double * restrict X,
                  kdtree_index * ID,
                  const kdtree_index n, /* Number of points */
                  const kdtree_index ndim, // number of dimensions
                  const kdtree_index vdim, /* Dimension to take value from */
                  const double pivot,
                  kdtree_index * nLow, kdtree_index * nHigh)
{
    int64_t low = -1;
    int64_t high = n;
    int64_t n2 = n;

    while(1)
    {
        do { low++; } while ( (low < n2)  && X[low*ndim + vdim] <= pivot );

        do { high--; } while ( (high > 0) && X[high*ndim + vdim] > pivot );

        if(low >= high)
        { *nLow = low;  *nHigh = n-*nLow;
#ifdef KDTREE_DEBUG
            assert(*nLow + *nHigh == n );
            for(int64_t kk = 0; kk < low; kk++)
            {
                //printf("Pivot = %f\n", pivot);
                //print_XID(X, low);
                assert(X[kk*ndim + vdim] <= pivot);
            }
            for(int64_t kk = low; kk < n2; kk++)
            {
                assert(X[kk*ndim + vdim] > pivot);
            }
#endif
            return;
        }
        swap_doubles(X + low*ndim, X + high*ndim, ndim);
        kdtree_index t = ID[low]; ID[low] = ID[high]; ID[high] = t;
    }
    return;
}

#ifdef KDTREE_DEBUG
void kdtree_validate(kdtree_t * T)
{
    // Check that all points are within bounds
    for(kdtree_index n = 0 ; n < T->n_nodes_alloc; n++){
        kdtree_node_t * node = T->nodes + n;
        if(node->n_point > 0) {
            if(node_is_final(T, node)){
                for(kdtree_index pp = 0 ; pp < node->n_point; pp++){
                    const double * X =  T->X + T->ndim*(node->point_offset+pp);
                    for(kdtree_index dd = 0; dd < T->ndim; dd ++){
                        const double * bbx = T->boxes + T->ndim*2*node->id;
                        if(X[dd] < bbx[2*dd]){
                            fprintf(stderr, "Point outside of box\n");
                            exit(EXIT_FAILURE);
                        }
                        if(X[dd] > bbx[2*dd+1]){
                            fprintf(stderr, "Point outside of box\n");
                            exit(EXIT_FAILURE);
                        }
                    }
                }
            }
        }
    }
}
#endif


void kdtree_free(kdtree_t * T)
{
    if(T == NULL) {
        return;
    }
    free(T->boxes);
    free(T->OID);
    free(T->X);
    free(T->nodes);
    if(T->pq != NULL) {
        pqheap_free(&T->pq);
    }
    free(T->result);
    free(T->point_buffer);
    free(T);
    return;
}

kdtree_t * kdtree_copy_shallow(kdtree_t * _T)
{
    assert(_T != NULL);
    if(_T == NULL)
        return NULL;
    kdtree_t * T = calloc(1, sizeof(kdtree_t));
    assert(T != NULL);
    memcpy(T, _T, sizeof(kdtree_t));
    T->result = NULL;
    T->result_alloc = 0;
    T->pq = NULL;
    return T;
}

void kdtree_free_shallow(kdtree_t * T)
{
    pqheap_free(&T->pq);
    free(T->result);
    free(T);
    return;
}

// Euclidean distance squared
static double
eudist_sq(const double * A, const double * B, const u32 ndim)
{
    double sum = 0;
    for(kdtree_index ii = 0; ii < ndim; ii++) {
        sum+=pow(A[ii]-B[ii], 2);
    }
    return sum;
}

double get_median_from_strided(const double * X, // data
                               kdtree_index N, // number of points
                               double * T, // temp buffer
                               kdtree_index stride) // stride
{
    // T is a temporary buffer, should be N elements large
    // https://www.gnu.org/software/gsl/doc/html/statistics.html
    // quickselect
    for(kdtree_index kk = 0; kk < N; kk++)
    {
        T[kk] = X[stride*kk];
        //printf("(%f) ", T[kk]);
    }
    //printf("\n");
    //printf("N=%zu, N/2=%zu\n", N, N/2);
#ifdef GSL
    double median = gsl_stats_median(T, 1, N);
#else
    double median = quickselect(T, N, N/2);
#endif
    return median;
}


void bounding_box(const double * restrict X,
                  const kdtree_index N, const kdtree_index ndim,
                  double * restrict bbx)
{
    for(kdtree_index dd = 0 ; dd < ndim; dd++)
    {
        bbx[2*dd] = X[dd]; // Min along dimension dd
        bbx[2*dd+1] = X[dd]; // Max along dimensions dd
    }
    for(kdtree_index nn = 0; nn < N; nn++)
    {
        for(kdtree_index dd = 0 ; dd < ndim; dd++)
        {
            X[ndim*nn + dd] < bbx[2*dd + 0] ? bbx[2*dd + 0] = X[ndim*nn + dd] : 0;
            X[ndim*nn + dd] > bbx[2*dd + 1] ? bbx[2*dd + 1] = X[ndim*nn + dd] : 0;
        }
    }
    return;
}


static void
print_bbx(const double * bbx, int ndim)
{
    printf("bbx=[");
    for(int kk = 0; kk < ndim; kk++)
    {
        printf("[%f, %f]", bbx[2*kk], bbx[2*kk+1]);
        if(kk + 1 == ndim){
            printf("]\n");
        } else {
            printf(", ");
        }
    }
}


/* Recursive splitting  */
void
kdtree_split(kdtree_t * T,
             kdtree_index node_id)
{
    kdtree_node_t * node = T->nodes + node_id;

    /* Possible to append children without running out of nodes? */
    if(2*node_id+2 >= T->n_nodes_alloc) {
        goto final;
    }

    if(node->n_point < (kdtree_index) T->max_leaf_size)
    {
    final: ; // Construct a "final" node without children
        node_set_final(T, node);
        return;
    }

    // Decide along which dimension to split using the bbx from the parent
    double * bbx = T->boxes + node_id*2*T->ndim;
    kdtree_index split_dim = T->ndim; // dimension or variable to split on
    {
        double max_size = bbx[1] - bbx[0];
        for(kdtree_index dd = 0; dd < T->ndim; dd++)
        {
            double t = bbx[2*dd+1] - bbx[2*dd];
            if(t >= max_size)
            {
                split_dim = dd;
                max_size = t;
            }
        }
        assert(max_size >= 0.0);
    }
    if(split_dim == T->ndim){
        print_bbx(bbx, T->ndim);
        fprintf(stderr, "Could not find a dimension to split on\n");
        exit(EXIT_FAILURE);
    }
    node->split_dim = split_dim;

    double pivot =
        get_median_from_strided( // coordinate split_dim of the first point that
                                 // belongs to the node
                                 T->X + node->point_offset*T->ndim + split_dim,
                                 node->n_point,
                                 T->median_buffer,
                                 T->ndim);
    //printf("split_dim = %u, pivot = %f, from %u points\n", split_dim, pivot, node->n_point);
    node->pivot = pivot;
    //printf("[%f   (pivot=%f)   %f]\n", node->bbx[2*split_dim], pivot, node->bbx[2*split_dim+1]);
    //assert(node->bbx[2*split_dim] <= pivot);
    //assert(node->bbx[2*split_dim+1] >= pivot);


    /* Avoid infinite recursion.
       This happens when points are flat in the splitting dimension */
    if(pivot == bbx[2*split_dim]
       || pivot == bbx[2*split_dim+1])
    {
        goto final;
    }

    /* Partition the data  */
    kdtree_index nLow = 0;
    kdtree_index nHigh = 0;

    double * node_X = T->X + T->ndim*node->point_offset;
    kdtree_index * node_OID = T->OID + node->point_offset;
    partition_vectors(node_X, node_OID,
                      node->n_point, T->ndim,
                      split_dim, pivot, &nLow, &nHigh);
    double * bbx_node = T->boxes + node->id*2*T->ndim;
    {
        kdtree_index left_id = node_left_child_id(node_id);

        assert(left_id < T->n_nodes_alloc);
        kdtree_node_t * node_left = T->nodes+left_id;
        assert(node_left->id == 0);
        node_left->id = left_id;
        double * bbx_left = T->boxes + left_id*2*T->ndim;
        memcpy(bbx_left,
               bbx_node,
               sizeof_bbx(T->ndim));
        bbx_left[2*split_dim + 1] = pivot;
        assert(bbx_left[2*split_dim] < pivot);
        node_left->n_point = nLow;
        node_left->point_offset = node->point_offset;
        kdtree_split(T, left_id);
    }

    {
        kdtree_index right_id = node_right_child_id(node_id);
        assert(right_id < T->n_nodes_alloc);
        kdtree_node_t * node_right = T->nodes+right_id;
        assert(node_right->id == 0); /* Unused? */
        node_right->id = right_id;
        double * bbx_right = T->boxes + right_id*2*T->ndim;
        memcpy(bbx_right, bbx_node, sizeof_bbx(T->ndim));
        bbx_right[2*split_dim] = pivot;
        assert(bbx_right[2*split_dim+1] > pivot);
        node_right->n_point = nHigh;
        node_right->point_offset = node->point_offset + nLow;
        kdtree_split(T, right_id);
    }

    return;
}

kdtree_t *
kdtree_new(const double * X,
           u32 N, u32 ndim,
           int max_leaf_size)
{

    if(max_leaf_size < 1){
        printf("kdtree_new: invalid bin size, use for example 10\n");
        return NULL;
    }

    if(N < 1){
        printf("kdtree_new: At least one data point needed\n");
        return NULL;
    }

    kdtree_t * T = calloc(1, sizeof(kdtree_t));
    if(T == NULL) { return NULL; }
    T->ndim = ndim;
    T->max_leaf_size = max_leaf_size;
    T->n_point = N;

    // Allocate storage for the nodes. We allocate enough
    // nodes for a complete binary tree up to some depth.
    // I.e. we will have 1, 3, 7, 15, ... (2^(L+1)-1) nodes where L
    // is the number of leafs.

    {
        /* If each leaf is 50% full we will have approximately */
        double n_leafs0 = 2.0 * (double) N / (double) max_leaf_size;
        /* Since it has to be a power of two we pick */
        double n_leafs = pow(2.0, ceil(log2(n_leafs0)));
        /* Then the number of nodes needed is */
        T->n_nodes_alloc = n_leafs*2 - 1;
        T->n_nodes_alloc < 3 ? T->n_nodes_alloc = 3 : 0;
    }

    // Per node data
    T->nodes = calloc(T->n_nodes_alloc, sizeof(kdtree_node_t));
    if(T->nodes == NULL) { goto failTree; }
    T->boxes = malloc(T->n_nodes_alloc*T->ndim*2*sizeof(double));

    // Per point data
    T->X = malloc(N*T->ndim*sizeof(double));
    if(T->X == NULL) { goto failTree; }
    memcpy(T->X, X, N*T->ndim*sizeof(double));
    T->OID = malloc(T->n_point*sizeof(kdtree_index));
    if(T->OID == NULL) { goto failTree; }
    for(u32 kk = 0; kk < T->n_point; kk++){
        T->OID[kk] = kk;
    }

    T->point_buffer = malloc(ndim*sizeof(double));
    if(T->point_buffer == NULL) { goto failTree; }
    T->median_buffer = calloc(N, sizeof(double));
    if(T->median_buffer == NULL) { goto failTree; }

    // Create the root node
    kdtree_node_t * node = T->nodes;
    double * bbx = T->boxes + node->id*2*T->ndim;
    bounding_box(X, N, T->ndim, bbx);
    node->n_point = N;
    node->point_offset = 0;

    // Recursive construction
    kdtree_split(T, // Tree
                 0); // node_id (location in array)

    free(T->median_buffer);
    T->median_buffer = NULL;
#ifdef KDTREE_DEBUG
    kdtree_validate(T);
#endif
    return T;

 failTree:
    kdtree_free(T);
    return NULL;
}


// Return 1 if the disk centered at Q
// with radius r is FULLY inside the node bounding box
// else 0
static int
within_bounds(const double * bbx,
              const kdtree_index ndim,
              const double * Q, const double r)
{
    for(u32 kk = 0; kk < ndim; kk++){
        if(Q[kk] + r > bbx[2*kk+1] || Q[kk] - r < bbx[2*kk]) {
            return 0;
        }
    }
    return 1;
}

// See if the axis aligned bounding box bbx overlaps the sphere centered at S
// and with a squared radius of r2.
static int
aa_box_hit_sphere_test(const double * restrict bbx,
                       kdtree_index ndim,
                       const double * restrict S,
                       const double r2,
                       double * restrict B)
{
    // Will eventually be the point in the bbx which is closest to
    // the sphere
    memcpy(B, S, ndim*sizeof(double));
    for(u32 dd = 0; dd < ndim; dd++){
        B[dd] < bbx[2*dd+0] ? B[dd] = bbx[2*dd+0] : 0;
        B[dd] > bbx[2*dd+1] ? B[dd] = bbx[2*dd+1] : 0;
    }
    double d = eudist_sq(S, B, ndim) <= r2;
    return d;
}

static int
bounds_overlap_ball(const kdtree_t * T,
                    const kdtree_node_t * node,
                    const double * Q)
{

    const double rmax = pqheap_get_max_value(T->pq);
    return aa_box_hit_sphere_test(T->boxes + T->ndim*2*node->id,
                                  T->ndim, Q, rmax,
                                  T->point_buffer);
}

// Recursive search until no more points can be found
//  Return 1 if we are done
//  Return 0 else
static int
kdtree_search_knn(kdtree_t * T, const kdtree_node_t * node, const double * Q)
{
    pqheap_t * pq = T->pq;

    if(node_is_final(T, node))
    {
        T->direct_path = 0;

        // Add all points
        for(kdtree_index kk = 0; kk<node->n_point; kk++)
        {
            const double * point_X = T->X + T->ndim*(node->point_offset + kk);
            kdtree_index point_ID = T->OID[node->point_offset + kk];
            double d2 = eudist_sq(point_X, Q, T->ndim);
            pqheap_insert(pq, d2, point_ID);
        }

        // Check if the most distal point in the priority queue
        // is confined within the bounding box of this leaf
        // what if the leaf contain less than the wanted number of points? TODO
        double rmax = sqrt(pqheap_get_max_value(pq));
        int done = within_bounds(T->boxes + 2*T->ndim*node->id,
                                  T->ndim, Q, rmax);
        return done;
    }

    // Descend depending on pivot
    // First take the path that gets us closer to the query point
    int split_dim = node->split_dim;
    int done = 0;
    if(Q[split_dim] > node->pivot)
    {
        // correct direction
        if(T->direct_path || bounds_overlap_ball(T, T->nodes + node_right_child_id(node->id), Q))
        {
            done = kdtree_search_knn(T, T->nodes + node_right_child_id(node->id), Q);
            if(done == 1)
            {
                return done;
            }
        }
        // "wrong direction"
        if(bounds_overlap_ball(T, T->nodes + node_left_child_id(node->id), Q))
        {
            done = kdtree_search_knn(T, T->nodes + node_left_child_id(node->id), Q);
            if(done == 1)
            {
                return done;
            }
        }
    } else {
        // "correct" direction
        if(T->direct_path || bounds_overlap_ball(T, T->nodes + node_left_child_id(node->id), Q))
        {
            done = kdtree_search_knn(T, T->nodes + node_left_child_id(node->id), Q);
            if(done)
            {
                return 1;
            }
        }

        // "wrong" direction
        if(bounds_overlap_ball(T, T->nodes + node_right_child_id(node->id), Q))
        {
            done = kdtree_search_knn(T, T->nodes + node_right_child_id(node->id), Q);
        }
        if(done)
        {
            return 1;
        }
    }
    // Now we have added all sub regions so we can check if the ball falls
    // within this non-end-node as well

    double rmax = sqrt(pqheap_get_max_value(pq));

    if(within_bounds(T->boxes + T->ndim*2*node->id,
                     T->ndim, Q, rmax))
    {
        return 1;
    }
    return 0;
}

kdtree_index * kdtree_query_knn(kdtree_t * T, const double * Q, kdtree_index k)
{
    if(k > T->n_point) {
        fprintf(stderr,
                "kdtree_query_knn error: Impossible to call for %u points when\n"
                "there are only %u in the tree\n", k, T->n_point);
        return NULL;
    }

    // If k changed from the last query, update:
    if(T->result_alloc != k) {
        if(T->pq != NULL)
        {
            pqheap_free(&T->pq);
        }
        T->pq = NULL;
        free(T->result);
        T->result = NULL;
        T->result_alloc = 0;
    }

    // Set up priority queue
    if(T->pq == NULL){
        T->pq = pqheap_new(k);
    }
    pqheap_t * pq = T->pq;
    pq->n = 0;
    pqheap_insert(pq, 1e99, 0);


    if(T->result == NULL){
        T->result = calloc(k, sizeof(kdtree_index));
        assert(T->result != NULL);
    }

    // Traverse the tree
    T->direct_path = 1;
    kdtree_search_knn(T, T->nodes, Q);

    // Move resulting indices from pq to array
    for(kdtree_index kk = 0; kk<k; kk++){
        double val = 0;
        uint64_t idx = 0;
        pqheap_pop(pq, &val, &idx);
        //printf("Popped: %lu, d = %f\n", idx, val);
        T->result[k-kk-1] = idx;
    }

    return T->result;
}


// Dynamic array to store the results during kdtree_query_radius
struct darray {
    kdtree_index * data;
    kdtree_index n_used;
    kdtree_index n_alloc;
};

// grow dynamic array
static void darray_n_more(struct darray * A, kdtree_index nmore)
{
    if(A->n_used + nmore >= A->n_alloc)
    {
        kdtree_index new_size = A->n_alloc + nmore;
        if(new_size < 1.2 *A->n_alloc)
        {
            new_size = 1.2*A->n_alloc;
        }
        A->data = realloc(A->data, new_size*sizeof(kdtree_index));
        assert(A->data != NULL);
        A->n_alloc = new_size;
    }
}

static void
_kdtree_query_radius(const kdtree_t * T,
                     const double * Q,
                     kdtree_index node_id,
                     const double r,
                     const double r2,
                     struct darray * res)
{
    kdtree_node_t * node = T->nodes + node_id;
    if( ! aa_box_hit_sphere_test(T->boxes + T->ndim*2*node->id,
                                 T->ndim, Q, r2,
                                 T->point_buffer) ) {
        return;
    }

    // If we reached a leaf see what points match the criteria
    if(node_is_final(T, node)){
        double * node_X = T->X + T->ndim*node->point_offset;
        kdtree_index * node_ID = T->OID + node->point_offset;
        darray_n_more(res, node->n_point);
        for(kdtree_index kk = 0; kk < node->n_point; kk++){
            if(eudist_sq(node_X + kk*T->ndim, Q, T->ndim) < r2){
                res->data[res->n_used] = node_ID[kk];
                res->n_used++;
            }
        }
        return;
    }
    // If not in a leaf, we see what children it makes sense to traverse
    // Some linear algebra: If Q + (mid-Q)/||mid-Q||*r crosses any of the 6 faces
    // we need to check. Simpler way to determine ?
    // if (x_intersect) if (y_intersect) if (z_intersect) then traverse ...
    _kdtree_query_radius(T, Q,
                         node_left_child_id(node_id),
                         r, r2, res);
    _kdtree_query_radius(T, Q,
                         node_right_child_id(node_id),
                         r, r2, res);
    return;
}

kdtree_index
kdtree_query_radius(const kdtree_t * T,
                    const double * Q,
                    const double radius,
                    kdtree_index ** result,
                    uint32_t * result_capacity)
{
    struct darray * res = calloc(1, sizeof(struct darray));
    kdtree_index * data = result[0];
    if(data == NULL){
        res->n_alloc = 100;
        res->data = calloc(res->n_alloc, sizeof(kdtree_index));
    } else {
        res->n_alloc = *result_capacity;
        res->data = result[0];
    }

    _kdtree_query_radius(T, Q, 0, radius, pow(radius, 2), res);

    *result = res->data;
    u32 nfound = res->n_used;
    *result_capacity = res->n_alloc;
    result[0] = res->data;
    free(res);
    return nfound;
}

static double gaussian(double d2, double sigma22)
{
    // sigma22 = 2*sigma^2
    // d2 = d^2
    return exp(-d2/sigma22);
}

static void
_kdtree_kde_mean(const kdtree_t * T,
                 const double * Q,
                 kdtree_index node_id,
                 const double r2,
                 const double sigma22,
                 double * xmeank,
                 double * meank,
                 kdtree_index * npoint)
{

    kdtree_node_t * node = T->nodes + node_id;
    /* Termination condition */
    if( ! aa_box_hit_sphere_test(T->boxes + T->ndim*2*node->id,
                                 T->ndim, Q, r2, T->point_buffer) )
    {
        return;
    }

    if(node_is_final(T, node)){
        double * node_X = T->X + node->point_offset*T->ndim;
        *npoint += node->n_point;

        for(kdtree_index kk = 0; kk < node->n_point; kk++){

            double * X = node_X + kk*T->ndim;
            double d2 = eudist_sq(X, Q, T->ndim);
            // Possibly check r2 criteria here
            double kde = gaussian(d2, sigma22);

            *meank += kde;
            for(int ll = 0; ll < 3; ll++){
                xmeank[ll] += kde*X[ll];
            }
        }
        return;
    }


    /* If not in a leaf, we see what children it makes sense to traverse */
    // Some linear algebra: If Q + (mid-Q)/||mid-Q||*r crosses any of the 6 faces
    // we need to check. Simpler way to determine ?
    // if (x_intersect) if (y_intersect) if (z_intersect) then traverse ...

    _kdtree_kde_mean(T, Q,
                     node_left_child_id(node_id),
                     r2, sigma22,
                     xmeank, meank, npoint);

    _kdtree_kde_mean(T, Q,
                     node_right_child_id(node_id),
                     r2, sigma22,
                     xmeank, meank, npoint);

    return;
}

static double
_kdtree_kde(const kdtree_t * T,
            const double * Q,
            kdtree_index node_id,
            const double r2,
            const double sigma22)
{
    kdtree_node_t * node = T->nodes + node_id;
    if( ! aa_box_hit_sphere_test(T->boxes + T->ndim*2*node->id,
                                 T->ndim, Q, r2, T->point_buffer) )
    {
        return 0;
    }

    double kde = 0;
    /* If we reached a leaf see what points match the criteria */
    if(node_is_final(T, node))
    {
        double * node_X = T->X + T->ndim*node->point_offset;

        for(kdtree_index kk = 0; kk < node->n_point; kk++)
        {
            double d2 = eudist_sq(node_X + kk*T->ndim, Q, T->ndim);
            // Possibly check r2 criteria here
            kde += gaussian(d2, sigma22);
        }
        return kde;
    }


    /* If not in a leaf, we see what children it makes sense to traverse */
    // Some linear algebra: If Q + (mid-Q)/||mid-Q||*r crosses any of the 6 faces
    // we need to check. Simpler way to determine ?
    // if (x_intersect) if (y_intersect) if (z_intersect) then traverse ...

    kde += _kdtree_kde(T, Q,
                       node_left_child_id(node_id),
                       r2, sigma22);

    kde += _kdtree_kde(T, Q,
                       node_right_child_id(node_id),
                       r2, sigma22);

    return kde;
}


double kdtree_kde(const kdtree_t * T,
                  const double * Q,
                  const double sigma,
                  const double cutoff)
{
    /* How distant points are of interest for the given sigma value ?
     */
    double r = 2.5*sigma;
    if(cutoff > 0.0)
    {
        r = cutoff*sigma;
    }
    return _kdtree_kde(T, Q, 0, pow(r,2.0), 2.0*pow(sigma, 2.0));
}

void
kdtree_kde_mean(const kdtree_t * T,
                const double * Q,
                double sigma,
                double cutoff,
                double * mean)
{
    double r = 4.0*2.5*sigma;
    if(cutoff > 0.0)
    {
        r = cutoff*sigma;
    }
    double xmeank[3] = {0};
    double meank = 0;
    kdtree_index npoint = 0;

    _kdtree_kde_mean(T, Q,
                     0, // start node == root
                     pow(r,2.0), // radius squared
                     2.0*pow(sigma, 2.0), // divisor for exponent
                     xmeank, // accumulate positions
                     &meank, // accumulate kernel values
                     &npoint); // number of points

    if(meank < 1e-9)
    {
        memcpy(mean, Q, 3*sizeof(double));
    } else {
        for(int kk = 0; kk < 3; kk++)
        {
            mean[kk] = xmeank[kk] / meank;
        }
    }
    //printf("npoint = %zu, meank = %F\n", npoint, meank);
    return;

}

static void
internal_kdtree_collide(const kdtree_t * T,
                        i64 u,
                        kdtree_index node_id,
                        const double r2,
                        kdtree_collide_cb cb_fun,
                        void * cb_data)

{
    const double * Q = T->X + T->ndim*u;
    kdtree_node_t * node = T->nodes + node_id;
    if( ! aa_box_hit_sphere_test(T->boxes + T->ndim*2*node->id,
                                 T->ndim, Q, r2, T->point_buffer) )
    {
        return;
    }

    if(node_is_final(T, node))
    {
        double * X = T->X + T->ndim*node->point_offset;
        for(kdtree_index kk = 0; kk < node->n_point; kk++)
        {
            if(eudist_sq(X + kk*T->ndim, Q, T->ndim) < r2)
            {
                i64 v = kk + node->point_offset;
                i64 u_orig = T->OID[u];
                i64 v_orig = T->OID[v];
                if(u < v){
                    cb_fun(u_orig, v_orig, cb_data);
                }
            }
        }
        return;
    }

    internal_kdtree_collide(T, u,
                            node_left_child_id(node_id),
                            r2, cb_fun, cb_data);
    internal_kdtree_collide(T, u,
                            node_right_child_id(node_id),
                            r2, cb_fun, cb_data);
    return;
}


void
kdtree_collide(const kdtree_t * T,
               double radius,
               kdtree_collide_cb cb_fun,
               void * cb_data)
{
    for(u32 u = 0; u < T->n_point; u++) {
        internal_kdtree_collide(T,
                                u,
                                0, // root
                                radius*radius,
                                cb_fun,
                                cb_data);
    }
    return;
}
