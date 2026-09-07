#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "kdtree.h"
#include "pqheap.h"
#include "quickselect.h"

typedef int64_t i64;
typedef uint32_t u32;

// Resolve the index of the right child based on the index of the
// parent, following a Eytzinger scheme
static size_t
node_right_child_id(size_t node_id)
{
    return 2*node_id+2;
}

static size_t
node_left_child_id(size_t node_id)
{
    return 2*node_id+1;
}

static size_t
sizeof_bbx(size_t ndim) {
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
swap_doubles(double * restrict X, double * restrict Y, size_t ndim)
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
                  u32 * ID,
                  const size_t n, /* Number of points */
                  const size_t ndim, // number of dimensions
                  const size_t vdim, /* Dimension to take value from */
                  const double pivot,
                  size_t * nLow, size_t * nHigh)
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
#ifndef NDEBUG
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
        u32 t = ID[low]; ID[low] = ID[high]; ID[high] = t;
    }
    return;
}

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
    for(size_t ii = 0; ii < ndim; ii++) {
        sum+=pow(A[ii]-B[ii], 2);
    }
    return sum;
}

double get_median_from_strided(const double * X, // data
                               size_t N, // number of points
                               double * T, // temp buffer
                               size_t stride) // stride
{
    // T is a temporary buffer, should be N elements large
    // https://www.gnu.org/software/gsl/doc/html/statistics.html
    // quickselect
    for(size_t kk = 0; kk < N; kk++)
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
                  const size_t N, const size_t ndim,
                  double * restrict bbx)
{
    for(size_t dd = 0 ; dd < ndim; dd++)
    {
        bbx[2*dd] = X[dd]; // Min along dimension dd
        bbx[2*dd+1] = X[dd]; // Max along dimensions dd
    }
    for(size_t nn = 0; nn < N; nn++)
    {
        for(size_t dd = 0 ; dd < ndim; dd++)
        {
            X[ndim*nn + dd] < bbx[2*dd + 0] ? bbx[2*dd + 0] = X[ndim*nn + dd] : 0;
            X[ndim*nn + dd] > bbx[2*dd + 1] ? bbx[2*dd + 1] = X[ndim*nn + dd] : 0;
        }
    }
    return;
}


/* Recursive splitting  */
void
kdtree_split(kdtree_t * T,
             size_t node_id)
{
    kdtree_node_t * node = T->nodes + node_id;

    /* Possible to append children without running out of nodes? */
    if(2*node_id+2 >= T->n_nodes_alloc) {
        goto final;
    }

    if(node->n_points < (size_t) T->max_leaf_size)
    {
    final: ; // Construct a "final" node without children
        node_set_final(T, node);
        return;
    }

    /* Decide along which dimension to split */
    double * bbx = T->boxes + node_id*2*T->ndim;
    size_t split_dim = 0; // dimension or variable to split on
    {
        double max_size = bbx[1] - bbx[0];
        for(size_t dd = 0; dd < T->ndim; dd++)
        {
            double t = bbx[2*dd+1] - bbx[2*dd];
            assert(t >= 0);
            if(t > max_size)
            {
                split_dim = dd;
                max_size = t;
            }
        }
    }
    node->split_dim = split_dim;

    double pivot =
        get_median_from_strided( // coordinate split_dim of the first point that
                                 // belongs to the node
                                T->X + node->offset*T->ndim + split_dim,
                                node->n_points,
                                T->median_buffer,
                                T->ndim);

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
    size_t nLow = 0;
    size_t nHigh = 0;

    double * node_X = T->X + T->ndim*node->offset;
    u32 * node_OID = T->OID + node->offset;
    partition_vectors(node_X, node_OID,
                      node->n_points, T->ndim,
                      split_dim, pivot, &nLow, &nHigh);
    double * bbx_node = T->boxes + node->id*2*T->ndim;
    {
        size_t left_id = node_left_child_id(node_id);

        assert(left_id < T->n_nodes_alloc);
        kdtree_node_t * node_left = T->nodes+left_id;
        assert(node_left->id == 0);
        node_left->id = left_id;
        double * bbx_left = T->boxes + left_id*2*T->ndim;
        memcpy(bbx_left,
               bbx_node,
               sizeof_bbx(T->ndim));
        bbx_left[2*split_dim + 1] = pivot;
        node_left->n_points = nLow;
        node_left->offset = node->offset;
        kdtree_split(T, left_id);
    }

    {
        size_t right_id = node_right_child_id(node_id);
        assert(right_id < T->n_nodes_alloc);
        kdtree_node_t * node_right = T->nodes+right_id;
        assert(node_right->id == 0); /* Unused? */
        node_right->id = right_id;
        double * bbx_right = T->boxes + right_id*2*T->ndim;
        memcpy(bbx_right, bbx_node, 2*T->ndim*sizeof(size_t));
        bbx_right[2*split_dim] = pivot;
        node_right->n_points = nHigh;

        node_right->offset = node->offset + nLow;

        kdtree_split(T, right_id);
    }

    return;
}

kdtree_t *
kdtree_new(const double * X,
           u32 N, u32 ndim,
           int max_leaf_size)
{

    if(max_leaf_size < 1)
    {
        printf("kdtree_new: invalid bin size, use for example 10\n");
        return NULL;
    }

    if(N < 1)
    {
        printf("kdtree_new: At least one data point needed\n");
        return NULL;
    }

#ifndef NDEBUG
    printf("kdtree warning: Not compiled with -DNDEBUG. Performance will be restrained. \n");
#endif

    /* Set up the tree and the basic settings*/
    kdtree_t * T = calloc(1, sizeof(kdtree_t));
    if(T == NULL) {
        return NULL;
    }
    T->ndim = ndim;
    T->max_leaf_size = max_leaf_size;
    T->n_points = N;

    /* Allocate storage for the nodes. We allocate enough
     * nodes for a complete binary tree up to some depth.
     * I.e. we will have 1, 3, 7, 15, ... (2^(L+1)-1) nodes where L
     * is the number of leafs.
     */

    {
        /* If each leaf is 50% full we will have approximately */
        double n_leafs0 = 2.0* (double) N / (double) max_leaf_size;
        /* Since it has to be a power of two we pick */
        double n_leafs = pow(2.0, ceil(log2(n_leafs0)));
        /* Then the number of nodes needed is */
        T->n_nodes_alloc = n_leafs*2-1;
        T->n_nodes_alloc < 3 ? T->n_nodes_alloc = 3 : 0;
    }
    //T->n_nodes_alloc = 2*N;
    T->nodes = calloc(T->n_nodes_alloc, sizeof(kdtree_node_t));
    T->X = malloc(N*T->ndim*sizeof(double));
    memcpy(T->X, X, N*T->ndim*sizeof(double));
    T->OID = malloc(T->n_points*sizeof(u32));
    T->boxes = malloc(T->n_nodes_alloc*T->ndim*2*sizeof(double));
    for(u32 kk = 0; kk < T->n_points; kk++){
        T->OID[kk] = kk;
    }
    assert(T->nodes != NULL);
    if(T->nodes == NULL)
    {
        printf("kdtree_new: Memory allocation failed. Tried to allocate for %zu nodes\n"
               "            but couldn't get it from the system\n",
               T->n_nodes_alloc);
        kdtree_free(T);
        return NULL;
    }

    T->median_buffer = calloc(N, sizeof(double));
    assert(T->median_buffer != NULL);

    // Create the root node
    kdtree_node_t * node = T->nodes;
    double * bbx = T->boxes + node->id*2*T->ndim;
    bounding_box(X, N, T->ndim, bbx);
    node->n_points = N;
    node->offset = 0;

    // Recursive construction
    kdtree_split(T, // Tree
                 0); // node_id (location in array)

    free(T->median_buffer);
    T->median_buffer = NULL;
    T->point_buffer = malloc(ndim*sizeof(double));
    return T;
}

size_t kdtree_query_closest(kdtree_t * T, double * X)
{
    kdtree_node_t * N = T->nodes;
    while( ! node_is_final(T, N) ){
        int split_dim = N->split_dim;
        if(X[split_dim] > N->pivot){
            N = T->nodes + node_right_child_id(N->id);
        } else {
            N = T->nodes+ + node_left_child_id(N->id);
        }
    }

    double * node_X = T->X + T->ndim*N->offset;
    u32 * node_OID = T->OID + N->offset;
    double dmin2 = eudist_sq(X, node_X, T->ndim);
    u32 imin = node_OID[0]; // N->idx[0];
    for(size_t kk = 0; kk<N->n_points; kk++){
        double d2 = eudist_sq(X, node_X + kk*T->ndim, T->ndim);
        if(d2 < dmin2){
            imin = node_OID[kk];
            dmin2 = d2;
        }
    }
    assert(dmin2 < 1e-9);
    return imin;
}

// Return 1 if the disk centered at Q
// with radius r is FULLY inside the node bounding box
// else 0
static int
within_bounds(const double * bbx,
              const size_t ndim,
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
                       size_t ndim,
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

/* Recursive search until no more points can be found
   Return 1 if we are done
   Return 0 else
*/

static int kdtree_search(kdtree_t * T, const kdtree_node_t * node, const double * Q)
{
    pqheap_t * pq = T->pq;

    if(node_is_final(T, node))
    {
        T->direct_path = 0;

        // Add all points
        for(size_t kk = 0; kk<node->n_points; kk++)
        {
            const double * point_X = T->X + T->ndim*(node->offset + kk);
            u32 point_ID = T->OID[node->offset + kk];
            double d2 = eudist_sq(point_X, Q, T->ndim);
            pqheap_insert(pq, d2, point_ID);
        }

        double rmax = sqrt(pqheap_get_max_value(pq));

        int done =  within_bounds(T->boxes + 2*T->ndim*node->id,
                                  T->ndim, Q, rmax);
        //printf("rmax = %f, done = %d\n", rmax, done);
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
            done = kdtree_search(T, T->nodes + node_right_child_id(node->id), Q);
            if(done == 1)
            {
                return done;
            }
        }
        // "wrong direction"
        if(bounds_overlap_ball(T, T->nodes + node_left_child_id(node->id), Q))
        {
            done = kdtree_search(T, T->nodes + node_left_child_id(node->id), Q);
            if(done == 1)
            {
                return done;
            }
        }
    } else {
        // "correct" direction
        if(T->direct_path || bounds_overlap_ball(T, T->nodes + node_left_child_id(node->id), Q))
        {
            done = kdtree_search(T, T->nodes + node_left_child_id(node->id), Q);
            if(done)
            {
                return 1;
            }
        }

        // "wrong" direction
        if(bounds_overlap_ball(T, T->nodes + node_right_child_id(node->id), Q))
        {
            done = kdtree_search(T, T->nodes + node_right_child_id(node->id), Q);
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

size_t * kdtree_query_knn(kdtree_t * T, const double * Q, size_t k)
{
    if(k > T->n_points)
    {
        fprintf(stderr,
                "kdtree_query_knn error: Impossible to call for %zu points when\n"
                "there are only %zu in the tree\n", k, T->n_points);
        return NULL;
    }
    //    printf("-> Q = (%f, %f)\n", Q[0], Q[1]);

    // If k changed from the last query, update:
    if(T->result_alloc != k)
    {
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
    if(T->pq == NULL)
    {
        T->pq = pqheap_new(k);
    }
    pqheap_t * pq = T->pq;
    pq->n = 0;
    pqheap_insert(pq, 1e99, 0);


    if(T->result == NULL)
    {
        T->result = calloc(k, sizeof(size_t));
        assert(T->result != NULL);
    }


    // Traverse the tree
    T->direct_path = 1;
    kdtree_search(T, T->nodes, Q);

    // If we don't need an ordered answer we could just traverse
    // the pq and extract the elements as we go.

    for(size_t kk = 0; kk<k; kk++)
    {
        double val = 0;
        uint64_t idx = 0;
        pqheap_pop(pq, &val, &idx);
        //printf("Popped: %lu, d = %f\n", idx, val);
        T->result[k-kk-1] = idx;
    }

    return T->result;
}


void kdtree_validate(kdtree_t * T)
{
#ifdef NDEBUG
    printf("kdtree_validate does not work when NDEBUG is defined\n");
    if(T == NULL)
    {
        printf("T is null\n");
    }
#else
    printf("kdtree_validate()\n");
    assert(T != NULL);
    assert(sizeof(double) == sizeof(size_t));
    assert(T->n_nodes_alloc > 0);

    // Check that all points are within bounds
    for(size_t n = 0 ; n < T->n_nodes_alloc; n++)
    {
        kdtree_node_t * node = T->nodes + n;
        if(node->n_points > 0)
        {


            //printf("Node id: %zu (Left %d, Right %d), %zu points\n", n,
            //       node->node_left, node->node_right, node->n_points*(node->node_left == -1));
            // TODO: check XID that the points are within bounds ...
            // Looks like things are wrong. Too many points in some leafs.

            //node_print_bbx(node);


            if(node_is_final(T, node))
            {
                //print_XID(XID, node->n_points);
                for(size_t pp = 0 ; pp < node->n_points; pp++)
                {
                    const double * X =  T->X + T->ndim*(node->offset+pp);
                    for(size_t dd = 0; dd < T->ndim; dd ++)
                    {
                        const double * bbx = T->boxes + T->ndim*2*node->id;
                        assert(X[dd] >= bbx[2*dd]);
                        assert(X[dd] <= bbx[2*dd+1]);
                    }
                }
            }
        }
    }
    printf("done\n");
#endif
}

struct darray {
    size_t * data;
    size_t n_used;
    size_t n_alloc;
};

static void darray_n_more(struct darray * A, size_t nmore)
{
    if(A->n_used + nmore >= A->n_alloc)
    {
        size_t new_size = A->n_alloc + nmore;
        if(new_size < 1.2 *A->n_alloc)
        {
            new_size = 1.2*A->n_alloc;
        }
        A->data = realloc(A->data, new_size*sizeof(size_t));
        assert(A->data != NULL);
        A->n_alloc = new_size;
    }
}

static void _kdtree_query_radius(const kdtree_t * T,
                                 const double * Q,
                                 size_t node_id,
                                 const double r,
                                 const double r2,
                                 struct darray * res)
{
    kdtree_node_t * node = T->nodes + node_id;
    if( ! aa_box_hit_sphere_test(T->boxes + T->ndim*2*node->id,
                                 T->ndim, Q, r2,
            T->point_buffer) )
    {
        return;
    }

    /* If we reached a leaf see what points match the criteria */
    if(node_is_final(T, node))
    {
        double * node_X = T->X + T->ndim*node->offset;
        u32 * node_ID = T->OID + node->offset;
        darray_n_more(res, node->n_points);
        for(size_t kk = 0; kk < node->n_points; kk++)
        {
            if(eudist_sq(node_X + kk*T->ndim, Q, T->ndim) < r2)
            {
                res->data[res->n_used] = node_ID[kk];
                res->n_used++;
            }
        }
        return;
    }
    /* If not in a leaf, we see what children it makes sense to traverse */
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

size_t *
kdtree_query_radius(const kdtree_t * T,
                    const double * Q,
                    const double radius,
                    size_t * nfound)
{
    struct darray * res = calloc(1, sizeof(struct darray));
    res->n_alloc = 100;
    res->data = calloc(res->n_alloc, sizeof(size_t));

    _kdtree_query_radius(T, Q, 0, radius, pow(radius, 2), res);

    size_t * result = res->data;
    *nfound = res->n_used;
    free(res);
    return result;
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
                 size_t node_id,
                 const double r2,
                 const double sigma22,
                 double * xmeank,
                 double * meank,
                 size_t * npoint)
{

    kdtree_node_t * node = T->nodes + node_id;
    /* Termination condition */
    if( ! aa_box_hit_sphere_test(T->boxes + T->ndim*2*node->id,
                                 T->ndim, Q, r2, T->point_buffer) )
    {
        return;
    }

    if(node_is_final(T, node))
    {
        double * node_X = T->X + node->offset*T->ndim;
        *npoint += node->n_points;

        for(size_t kk = 0; kk < node->n_points; kk++)
        {
            double * X = node_X + kk*T->ndim;
            double d2 = eudist_sq(X, Q, T->ndim);
            // Possibly check r2 criteria here
            double kde = gaussian(d2, sigma22);

            *meank += kde;
            for(int ll = 0; ll < 3; ll++)
            {
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
            size_t node_id,
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
        double * node_X = T->X + T->ndim*node->offset;

        for(size_t kk = 0; kk < node->n_points; kk++)
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
    size_t npoint = 0;

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
                        size_t node_id,
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
        double * X = T->X + T->ndim*node->offset;
        for(size_t kk = 0; kk < node->n_points; kk++)
        {
            if(eudist_sq(X + kk*T->ndim, Q, T->ndim) < r2)
            {
                i64 v = kk + node->offset;
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
    for(u32 u = 0; u < T->n_points; u++) {
        internal_kdtree_collide(T,
                                u,
                                0, // root
                                radius*radius,
                                cb_fun,
                                cb_data);
    }
    return;
}
