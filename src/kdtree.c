#include "../include/kdtree.h"
#include "pqheap.h"
#include <string.h>

#define XID_STRIDE (KDTREE_DIM + 1)

// Make a shallow copy of a kd-tree for usage by another thread.
kdtree_t * kdtree_copy(kdtree_t * );

static void print_XID(double * XID, size_t N)
{
    size_t * ID = (size_t * ) XID;
    for(size_t kk = 0 ; kk< N; kk++)
    {
        printf("  [%6zu] ", ID[kk*(1+KDTREE_DIM) + KDTREE_DIM]);
        for(size_t dd = 0; dd < KDTREE_DIM; dd ++)
        {
            printf("%4.2f ", XID[kk * (1+KDTREE_DIM) + dd]);
        }
        printf("\n");


    }
}

/*  Hoare's partition scheme
    Adopted from arch/24/03/11_quickselect */
static void
partition(double * restrict X,
          const size_t n, /* Number of points */
          const size_t vdim, /* Dimension to take value from */
          const double pivot,
          size_t * nLow, size_t * nHigh)
{
    int64_t low = -1;
    int64_t high = n;
    int64_t n2 = n;

    while(1)
    {
        do { low++; } while ( X[low*XID_STRIDE+vdim] <= pivot && low < n2 );

        do { high--; } while ( X[high*XID_STRIDE+vdim] > pivot && high > 0);

        if(low >= high)
        { *nLow = low;  *nHigh = n-*nLow;
#ifndef NDEBUG
            assert(*nLow + *nHigh == n );
            for(int64_t kk = 0; kk < low; kk++)
            {
                //printf("Pivot = %f\n", pivot);
                //print_XID(X, low);
                assert(X[kk*XID_STRIDE+vdim] <= pivot);
            }
            for(int64_t kk = low; kk < n2; kk++)
            {
                assert(X[kk*XID_STRIDE+vdim] > pivot);
            }
#endif
            return;
        }

        /* Swap data */
        double t[XID_STRIDE];
        memcpy(t,
               X + low*XID_STRIDE,
               XID_STRIDE*sizeof(double));
        memcpy(X + low*XID_STRIDE,
               X + high*XID_STRIDE,
               XID_STRIDE*sizeof(double));
        memcpy(X + high*(KDTREE_DIM+1),
               t,
               XID_STRIDE*sizeof(double));
    }

    return;
}


void
node_print_bbx(const kdtree_node_t * N)
{

    printf("#%zu: ", N->id);
    for(size_t kk = 0; kk < KDTREE_DIM; kk++)
    {
        printf("[%f -- %f]", N->bbx[2*kk], N->bbx[2*kk+1]);
        if(kk + 1 < KDTREE_DIM)
        {
            printf("x");
        }
    }
    printf("\n");
    return;
}

void kdtree_free(kdtree_t ** _T)
{
    kdtree_t * T = _T[0];

    free(T->XID);
    free(T->nodes);
    if(T->pq != NULL)
    {
        pqheap_free(&T->pq);
    }
    if(T->KN != NULL)
    {
        free(T->KN);
    }
    free(T);
    _T[0] = NULL;
}

kdtree_t * kdtree_copy(kdtree_t * _T)
{
    assert(_T != NULL);
    if(_T == NULL)
        return NULL;
    kdtree_t * T = malloc(sizeof(kdtree_t));
    memcpy(T, _T, sizeof(kdtree_t));
    T->k = -1;
    T->pq = NULL;
    T->KN = NULL;
    return T;
}

/* Euclidean distance squared */
static double eudist_sq(const double * A, const double * B)
{
    double sum = 0;
    for(size_t ii = 0; ii<KDTREE_DIM; ii++)
    {
        sum+=pow(A[ii]-B[ii], 2);
    }
    return sum;
}


static double eudist(const double * A, const double * B)
{
    return sqrt(eudist_sq(A, B));
}


double get_median_from_strided(const double * X, // data
                               size_t N, // number of points
                               double * T, // temp buffer
                               size_t stride) // stride
{
    // T is a temporary buffer, should be N elements large
    // https://www.gnu.org/software/gsl/doc/html/statistics.html
    // quickselect
    assert(stride == KDTREE_DIM + 1);
    for(size_t kk = 0; kk < N; kk++)
    {
        T[kk] = X[stride*kk];
        //printf("(%f) ", T[kk]);
    }
    //printf("\n");
    double median = gsl_stats_median(T, 1, N);
    return median;
}


void bounding_box(const double * restrict X,
                  const size_t N, const size_t D,
                  double * restrict bbx)
{
    for(size_t dd = 0 ; dd < D; dd++)
    {
        bbx[2*dd] = X[dd]; // Min along dimension dd
        bbx[2*dd+1] = X[dd]; // Max along dimensions dd
    }
    for(size_t nn = 0; nn < N; nn++)
    {
        for(size_t dd = 0 ; dd < D; dd++)
        {
            X[D*nn + dd] < bbx[2*dd] ? bbx[2*dd] = X[D*nn + dd] : 0;
            X[D*nn + dd] > bbx[2*dd + 1] ? bbx[2*dd + 1] = X[D*nn + dd] : 0;
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
    assert(node->data_idx % XID_STRIDE == 0);

    node_print_bbx(node);

    /* Possible to append children without running out of nodes? */
    if(2*node_id+2 >= T->n_nodes_alloc) {
        goto final;
    }

    if(node->n_points < (size_t) T->binsize)
    {
    final: ; // Construct a "final" node without children
        node->node_left = -1;
        node->node_right = -1;
        return;
    }

    /* Decide along which dimension to split */
    size_t split_dim = 0; // dimension or variable to split on
    {
        double max_size = node->bbx[1] - node->bbx[0];
        for(size_t dd = 0; dd < KDTREE_DIM; dd++)
        {
            double t = node->bbx[2*dd+1] - node->bbx[2*dd];
            assert(t >= 0);
            if(t > max_size)
            {
                split_dim = dd;
                max_size = t;
            }
        }
    }
    node->split_dim = split_dim;

    double * XID = T->XID + node->data_idx;
    assert(node->data_idx + node->n_points <= T->N*XID_STRIDE);
    double pivot =
        get_median_from_strided(XID+split_dim,
                                node->n_points,
                                T->median_buffer,
                                XID_STRIDE);

    node->pivot = pivot;
    printf("[%f   (pivot=%f)   %f]\n", node->bbx[2*split_dim], pivot, node->bbx[2*split_dim+1]);
    assert(node->bbx[2*split_dim] <= pivot);
    assert(node->bbx[2*split_dim+1] >= pivot);


    /* Avoid infinite recursion.
     This happens when points are flat in the splitting dimension */
    if(pivot == node->bbx[2*split_dim]
       || pivot == node->bbx[2*split_dim+1])
    {
        goto final;
    }

    /* Partition the data  */
    size_t nLow = 0;
    size_t nHigh = 0;

    size_t * ID = (size_t *) XID;
    assert(XID[KDTREE_DIM] < T->N);
    partition(XID, node->n_points, split_dim, pivot, &nLow, &nHigh);



    {
        size_t left_id = 2*node_id+1; /* according to Etyzinger / Binary tree */
        node->node_left = left_id;
        assert(left_id < T->n_nodes_alloc);
        kdtree_node_t * node_left = T->nodes+left_id;
        assert(node_left->id == 0);
        node_left->id = left_id;
        memcpy(node_left->bbx, node->bbx, 2*KDTREE_DIM*sizeof(size_t));
        node_left->bbx[2*split_dim + 1] = pivot;
        node_left->n_points = nLow;

        node_left->data_idx = node->data_idx;

        kdtree_split(T, left_id);
    }
    {
        size_t right_id = 2*node_id+2;
        node->node_right = right_id;
        assert(right_id < T->n_nodes_alloc);
        kdtree_node_t * node_right = T->nodes+right_id;
        assert(node_right->id == 0); /* Unused? */
        node_right->id = right_id;
        memcpy(node_right->bbx, node->bbx, 2*KDTREE_DIM*sizeof(size_t));
        node_right->bbx[2*split_dim] = pivot;
        node_right->n_points = nHigh;

        node_right->data_idx = node->data_idx + nLow*XID_STRIDE;

        if(node_right->data_idx + node_right->n_points > T->N)
        {
            printf("node->data_idx = %zu, node->n_points = %zu, T->N = %zu\n", node->data_idx, node->n_points, T->N);
        }
        kdtree_split(T, right_id);
    }

    return;
}

kdtree_t *
kdtree_new(const double * X,
           size_t N, size_t ndim,
           int binsize)
{
    if(ndim !=  KDTREE_DIM)
    {
        printf("kdtree_new: Invalid number of dimensions\n");
        return NULL;
    }

    if(binsize < 1)
    {
        printf("kdtree_new: invalid bin size\n");
        return NULL;
    }

    if(N < 1)
    {
        printf("kdtree_new: At least one data point needed\n");
        return NULL;
    }

    /* Set up the tree and the basic settings*/
    kdtree_t * T = calloc(1, sizeof(kdtree_t));
    assert( T!= NULL);
    T->binsize = binsize;
    T->k = -1; // not set up for queries of any size
    T->N = N;

    /* Allocate storage for the nodes */
    /* This use too much memory, we only need approximately
     * N/binsize nodes. We could try 2*N/binsize and grow
     * later on if needed */
    T->n_nodes_alloc = 2*N;
    T->nodes = calloc(T->n_nodes_alloc, sizeof(kdtree_node_t));

    /* Copy the data points and attach an index */
    T->XID = malloc((KDTREE_DIM+1)*N*sizeof(double));
    for(size_t kk = 0; kk<N; kk++)
    {
        memcpy(T->XID+kk*(KDTREE_DIM+1),
               X+kk*(KDTREE_DIM),
               KDTREE_DIM*sizeof(double));
        size_t * ID = (size_t *) T->XID + (KDTREE_DIM+1)*kk + KDTREE_DIM;
        *ID = kk;
    }

    //print_XID(T->XID, T->N);

    // Expand the region for correct ball-within-region checking
    //double dx = xmax-xmin;
    //double dy = ymax-ymin;
    //xmin -= 2*dx; xmax += 2*dx;
    //ymin -= 2*dy; ymax += 2*dy;

    T->median_buffer = calloc(N, sizeof(double));

    kdtree_node_t * node = T->nodes;
    bounding_box(X, N, KDTREE_DIM, node->bbx);
    node->n_points = N;
    node->data_idx = 0;

    /* Recursive construction */
    kdtree_split(T, // Tree
                 0); // node_id (location in array)

    free(T->median_buffer);
    T->median_buffer = NULL;

    //print_XID(T->XID, T->N);
    return T;
}

size_t kdtree_query_closest(kdtree_t * T, double * X)
{

    kdtree_node_t * N = T->nodes;
    while(N->node_left > 0)
    {
        int split_dim = N->split_dim;
        if(X[split_dim] > N->pivot)
        {
            N = T->nodes+N->node_right;
        } else {
            N = T->nodes+N->node_left;
        }
    }
    assert(N->node_left == -1);
    assert(N->node_left == -1);
    double * NX = T->XID + N->data_idx;
    double dmin = eudist(X, NX);
    size_t imin = NX[KDTREE_DIM]; // N->idx[0];
    for(size_t kk = 0; kk<N->n_points; kk++)
    {
        double d = eudist(X, NX + kk*(KDTREE_DIM+1));
        if(d < dmin)
        {
            imin = NX[kk*(KDTREE_DIM+1) + KDTREE_DIM]; // ->idx[kk];
            dmin = d;
        }
    }
    assert(dmin < 1e-9);
    return imin;
}

int within_bounds(const kdtree_node_t * node, const double * Q, const double r)
{
    // Return 1 if the disk centered at Q
    // with radius r is inside the node
    // else 0
    const double x = Q[0];
    const double y = Q[1];
    const double z = Q[2];
    if(x + r > node->bbx[1] || x - r < node->bbx[0] ||
       y + r > node->bbx[3] || y - r < node->bbx[2] ||
       z + r > node->bbx[5] || z - r < node->bbx[4])
    {
        return 0;
    } else {
        //   printf("(%f, %f), r=%f is within [%f, %f, %f, %f]\n", x, y, r,
        //     node->xmin, node->xmax, node->ymin, node->ymax);
        return 1;
    }
}

static inline double min(double a, double b)
{
    if(a<b)
    {
        return a;
    }
    return b;
}

int bounds_overlap_ball(const kdtree_t * T, const kdtree_node_t * node, const double * Q)
{

    const double rmax = pqheap_get_max_value(T->pq);


    // find nearest x and nearest y
    double nx = 0;
    double ny = 0;
    double nz = 0;
    const double x = Q[0];
    const double y = Q[1];
    const double z = Q[2];

    if( x < node->bbx[0])
    {
        nx = node->bbx[0] - x;
    } else if (x > node->bbx[1])
    {
        nx = x - node->bbx[1];
    }

    if( y < node->bbx[2])
    {
        nx = node->bbx[2] - y;
    } else if (y > node->bbx[3])
    {
        ny = y - node->bbx[3];
    }

    if( z < node->bbx[4])
    {
        nz = node->bbx[4] - z;
    } else if (z > node->bbx[5])
    {
        nz = z - node->bbx[5];
    }

    if(pow(nx, 2) + pow(ny, 2) + pow(nz, 2) < rmax)
    {
        return 1;
    }

    return 0;

}

int kdtree_search(kdtree_t * T, const kdtree_node_t * node, const double * Q)
{
    /* Recursive search until no more points can be found */
    // Return 1 if we are done
    // Return 0 else
    if(0)
    {
        printf("Visiting node %6zu #=%7zu",
               node->id, node->n_points);
        node_print_bbx(node);
    }

    pqheap_t * pq = T->pq;

    // If final node
    if(node->node_left == -1)
    {
        T->direct_path = 0;
        const double * NX = T->XID + node->data_idx;
        const size_t * NID = (size_t *) T->XID + node->data_idx;

        // Add all points
        for(size_t kk = 0; kk<node->n_points; kk++)
        {
            double d2 = eudist_sq(NX + kk*XID_STRIDE, Q);
            //X += 2;
            //printf("Adding %f\n", sqrt(d2));
            //fprintf(T->log, "%f ", d2);

            pqheap_insert(pq, d2, NID[kk*XID_STRIDE + KDTREE_DIM]);
        }

        double rmax = sqrt(pqheap_get_max_value(pq));

        int done =  within_bounds(node, Q, rmax);
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
        if(T->direct_path || bounds_overlap_ball(T, T->nodes + node->node_right, Q))
        {
        done = kdtree_search(T, T->nodes + node->node_right, Q);
        if(done == 1)
        {
            return done;
        }
        }
        // "wrong direction"
        if(bounds_overlap_ball(T, T->nodes + node->node_left, Q))
        {
            done = kdtree_search(T, T->nodes + node->node_left, Q);
            if(done == 1)
            {
                return done;
            }
        }
    } else {
        // "correct" direction
        if(T->direct_path || bounds_overlap_ball(T, T->nodes + node->node_left, Q))
        {
        done = kdtree_search(T, T->nodes + node->node_left, Q);
        if(done)
        {
            return 1;
        }
        }

        // "wrong" direction
        if(bounds_overlap_ball(T, T->nodes + node->node_right, Q))
        {
          done = kdtree_search(T, T->nodes + node->node_right, Q);
        }
        if(done)
        {
            return 1;
        }
    }
    // Now we have added all sub regions so we can check if the ball falls
    // within this non-end-node as well

    double rmax = sqrt(pqheap_get_max_value(pq));

    if(within_bounds(node, Q, rmax))
    {
        return 1;
    }
    return 0;
}

size_t * kdtree_query_knn(kdtree_t * T, const double * Q, int k)
{
    //    printf("-> Q = (%f, %f)\n", Q[0], Q[1]);

    // If k changed from the last query, update:
    if(k != T->k)
    {
        free(T->pq);
        T->pq = NULL;
        free(T->KN);
        T->KN = NULL;
    }
    T->k = k;

    // Set up priority queue
    if(T->pq == NULL)
    {
        T->pq = pqheap_new(k);
    }
    pqheap_t * pq = T->pq;
    pq->n = 0;
    pqheap_insert(pq, 1e99, 0);


    if(T->KN == NULL)
    {
        T->KN = calloc(k, sizeof(size_t));
    }


    // Traverse the tree
    T->direct_path = 1;
    kdtree_search(T, T->nodes, Q);

    // If we don't need an ordered answer we could just traverse
    // the pq and extract the elements as we go.

    for(int kk = 0; kk<k; kk++)
    {
        double val = 0;
        uint64_t idx = 0;
        pqheap_pop(pq, &val, &idx);
        //printf("Popped: %lu, d = %f\n", idx, val);
        T->KN[k-kk-1] = idx;
    }

    return T->KN;
}

// Struct for parallel queries
typedef struct{
    kdtree_t * T;
    const double * Q;
     size_t nQ;
    int k;
    int thread;
    int nthreads;
    size_t * KNN;
} _p_query_t;

void * _p_query(void * _config)
{
    _p_query_t * config = (_p_query_t *) _config;
    kdtree_t * T = config->T;

    const double * Q = config->Q;
    const size_t nQ = config->nQ;
    const int k = config->k;
    const int thread = config->thread;
    const int nthreads = config->nthreads;
    size_t * KNN = config->KNN;
    for(size_t kk = thread; kk<nQ; kk+=nthreads)
    {
        size_t * knn = kdtree_query_knn(T, Q+2*kk, k);
        memcpy(KNN + k*kk, knn, k*sizeof(size_t));
    }
    return NULL;
}

size_t * kdtree_query_knn_multi(kdtree_t * T, const double * Q, size_t nQ, int k, int nthreads)
{
    size_t * KNN = malloc(nQ*k*sizeof(size_t));

    if(nthreads == 1)
    {
        for(size_t kk = 0; kk<nQ; kk++)
        {
            size_t * knn = kdtree_query_knn(T, Q+2*kk, k);
            memcpy(KNN + k*kk, knn, k*sizeof(size_t));
        }
        return KNN;
    } else {
        pthread_t * threads = malloc(nthreads*sizeof(pthread_t));
        _p_query_t * confs = malloc(nthreads*sizeof(_p_query_t));

        for(int kk = 0; kk<nthreads; kk++)
        {
            confs[kk].T = kdtree_copy(T);
            confs[kk].thread = kk;
            confs[kk].nthreads = nthreads;
            confs[kk].KNN = KNN;
            confs[kk].k = k;
            confs[kk].Q = Q;
            confs[kk].nQ = nQ;
            pthread_create(&threads[kk], NULL, _p_query, (void *) &confs[kk]);
        }

        for(int kk = 0; kk<nthreads; kk++)
        {
            pthread_join(threads[kk], NULL);
            // TODO: free some stuff in T
            pqheap_free(&confs[kk].T->pq);
            free(confs[kk].T->KN);
            free(confs[kk].T);
        }
        free(confs);
        free(threads);
    }

    return KNN;
}

void kdtree_validate(kdtree_t * T)
{
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

            double * XID = T->XID + node->data_idx;
            if(node->node_left == -1)
            {
                //print_XID(XID, node->n_points);
                for(size_t pp = 0 ; pp < node->n_points; pp++)
                {
                    for(size_t dd = 0; dd < KDTREE_DIM; dd ++)
                    {
                        assert(XID[pp*XID_STRIDE + dd] >= node->bbx[2*dd]);
                        assert(XID[pp*XID_STRIDE + dd] <= node->bbx[2*dd+1]);
                    }
            }
            }
        }
    }
    printf("done\n");
}
