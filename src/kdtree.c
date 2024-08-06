#include "../include/kdtree.h"
#include "pqheap.h"

// Make a shallow copy of a kd-tree for usage by another thread.
kdtree_t * kdtree_copy(kdtree_t * );

void node_print_bbx(const kdtree_t * T, const kdtree_node_t * N)
{

    for(size_t kk = 0; kk < T->ndim; kk++)
    {
        printf("[%f -- %f]", N->bbx[2*kk], N->bbx[2*kk+1]);
        if(kk + 1 < T->ndim)
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

    free(T->idx_local);
    free(T->X_local);
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
static double eudist_sq(const double * A, const double * B, size_t ndim)
{
    double sum = 0;
    for(size_t ii = 0; ii<ndim; ii++)
    {
        sum+=pow(A[ii]-B[ii], 2);
    }
    return sum;
}


static double eudist(const double * A, const double * B, size_t ndim)
{
    return sqrt(eudist_sq(A, B, ndim));
}


double get_median_from_strided(const double * X, size_t N, double * T, size_t stride)
{
    // T is a temporary buffer, should be N elements large
    // https://www.gnu.org/software/gsl/doc/html/statistics.html
    // quickselect
    for(size_t kk = 0; kk < N; kk++)
    {
        T[kk] = X[stride*kk];
    }
    double median =  gsl_stats_median(T, 1, N);
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
kdtree_split(kdtree_t * T, // Add kdtree_node_t * parent
             const double * X, size_t * idx, size_t n_points)
{
    kdtree_node_t * node = T->nodes + T->next;
    node->id = T->next;
    T->next++;

    node->n_points = n_points; // number of elements to split on

    // This should be given by the parent
    // No need to repeat all the time.
    // only if parent == NULL
    bounding_box(X, n_points, T->ndim, node->bbx);

    if(n_points < (size_t) T->binsize)
    {
    final: ; // Construct a "final" node without children
        node->node_left = -1;
        node->node_right = -1;

        node->X =  T->X_local + 2*T->next_idx;
        node->idx = T->idx_local + T->next_idx;
        //printf("T->next_idx = %zu, idx[0]=%zu\n", T->next_idx, idx[0]);
        //getchar();
        for(size_t kk = 0; kk < n_points; kk++)
        {
            node->idx[kk] = idx[kk];
            for(size_t ii = 0; ii < T->ndim; ii++)
            {
                node->X[T->ndim*kk + ii] =     X[T->ndim*kk+ii];
            }
        }
        //node_shuffle(node->X, node->idx, kk);

        T->next_idx += n_points;
        return;
    }

    /* Decide along which dimension to split */
    {
        size_t split_dim = 0; // dimension or variable to split on
        double max_size = node->bbx[1] - node->bbx[0];
        for(size_t dd = 0; dd < T->ndim; dd++)
        {
            double t = node->bbx[2*dd+1]-node->bbx[2*dd];
            assert(t >= 0);
            if(t > max_size)
            {
                split_dim = dd;
                max_size = t;
            }
        }
        node->split_dim = split_dim;
    }


    node->pivot =
        get_median_from_strided(X+node->split_dim,
                                n_points, T->median_buffer,
                                T->ndim);

    /* Avoid infinite recursion. Needed? */
    if(node->pivot == node->bbx[2*node->split_dim]
       || node->pivot == node->bbx[2*node->split_dim+1])
    {
        goto final;
    }

    size_t n_high = 0;
    for(size_t kk = 0; kk < n_points; kk++)
    {
        if(X[T->ndim*kk+node->split_dim] > node->pivot)
        {
            n_high++;
        }
    }


    // TODO: Looks like these use twice as much memory as needed together */
    // Not needed at all. We split the list of coordinates (and idx)
    // in the same way as we do with qsort

    size_t n_left = n_points - n_high;
    size_t n_right = n_high;
    printf("n_left: %zu, n_right: %zu\n", n_left, n_right);
    assert(n_left != 0);
    assert(n_right != 0);
    double * Xleft = malloc(T->ndim*n_left*sizeof(double));
    size_t * idx_left = malloc(n_left*sizeof(double));

    double * Xright = malloc(T->ndim*n_right*sizeof(double));
    size_t * idx_right = malloc(n_right*sizeof(double));

    size_t rpos = 0;
    size_t lpos = 0;
    for(size_t kk = 0; kk < n_points; kk++)
    {
        if(X[T->ndim*kk + node->split_dim] > node->pivot)
        { // Right=High=Above pivot
            for(size_t ii = 0; ii < T->ndim; ii++)
            {
                Xright[T->ndim*rpos+ii] = X[T->ndim*kk+ii];
            }
            idx_right[rpos] = idx[kk];
            rpos++;
        } else {
            for(size_t ii = 0; ii < T->ndim; ii++)
            {
                Xleft[T->ndim*lpos+ii] = X[T->ndim*kk];
            }
            idx_left[lpos] = idx[kk];
            lpos++;
        }
    }

    assert(lpos > 0); // if this is the case, this should be
    assert(rpos > 0); // made a final node (duplicate points...)
    assert(lpos == n_left);
    assert(rpos == n_right);

    assert(lpos+rpos == n_points);

    node->node_right = T->next;
    kdtree_split(T, Xright, idx_right, rpos);
    free(Xright);
    free(idx_right);


    node->node_left = T->next;;
    kdtree_split(T, Xleft, idx_left, lpos);

    free(Xleft);
    free(idx_left);
    return;
}

kdtree_t *
kdtree_new(const double * X,
           size_t N, size_t ndim,
           int binsize)
{
    if(ndim < 2 || ndim > KDTREE_MAXDIM)
    {
        printf("kdtree_new: Invalid number of dimensions\n");
        return NULL;
    }
    if(binsize < 1)
    {
        printf("kdtree_new: invalid bin size\n");
        return NULL;
    }

    /// Set up the tree
    kdtree_t * T = malloc(sizeof(kdtree_t));
    T->ndim = ndim;
    T->k = -1; // not set up for queries of any size
    /* This use too much memory, we only need approximately
     * N/binsize nodes. We could try 2*N/binsize and grow
     * later on if needed */
    T->nodes = malloc(N*sizeof(kdtree_node_t));

    T->n_nodes = 0;
    T->idx_local = malloc(N*sizeof(size_t));
    T->X_local = malloc(T->ndim*N*sizeof(size_t));
    T->pq = NULL;
    T->KN = NULL;

    T->next = 0; // next node to write to
    T->next_idx = 0; // next idx_local to write to
    T->binsize = binsize;

    size_t * idx_list = malloc(N*sizeof(size_t));
    for(size_t kk = 0; kk<N; kk++)
    {
        idx_list[kk] = kk;
    }

    /// Set up the root node
    //double bbx[KDTREE_MAXDIM*2] = {0};

    // Expand the region for correct ball-within-region checking
    //double dx = xmax-xmin;
    //double dy = ymax-ymin;
    //xmin -= 2*dx; xmax += 2*dx;
    //ymin -= 2*dy; ymax += 2*dy;

    T->median_buffer = malloc(N*sizeof(double));

    /* Recursive construction */
    kdtree_split(T, X, idx_list, N);

    free(T->median_buffer);
    T->median_buffer = NULL;
    T->n_nodes = T->next;
    /* TODO: did we allocate more than needed? */
    free(idx_list);

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
    double dmin = eudist(X, N->X, T->ndim);
    size_t imin = N->idx[0];
    for(size_t kk = 0; kk<N->n_points; kk++)
    {
        double d = eudist(X, N->X + 2*kk, T->ndim);
        if(0){
        printf("%zu (%f, %f), %f\n", N->idx[kk],
               N->X[2*kk], N->X[2*kk+1], d);
        }
        if(d < dmin)
        {
            imin = N->idx[kk];
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
    if(x + r > node->bbx[1] || x - r < node->bbx[0] ||
       y + r > node->bbx[3] || y - r < node->bbx[2])
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
    const double x = Q[0];
    const double y = Q[1];
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

    if(pow(nx, 2) + pow(ny, 2) < rmax)
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
        node_print_bbx(T, node);
    }

    pqheap_t * pq = T->pq;

    // If final node
    if(node->node_left == -1)
    {
        T->direct_path = 0;
        const double * X = node->X;
        const size_t * idx = node->idx;

        // Add all points
        for(size_t kk = 0; kk<node->n_points; kk++)
        {
            double d2 = eudist_sq(X+2*kk, Q, T->ndim);
            //X += 2;
            //printf("Adding %f\n", sqrt(d2));
            //fprintf(T->log, "%f ", d2);

            pqheap_insert(pq, d2, idx[kk]);
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
        T->KN = malloc(k*sizeof(size_t));
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
    // Check that all points are within bounds
    for(size_t n = 0 ; n < T->n_nodes; n++)
    {
        kdtree_node_t * node = T->nodes + n;
        printf("Node id: %zu (%d, %d)\n", n,
               node->node_left, node->node_right);

        node_print_bbx(T, node);
    }
}
