#include "../include/kdtree.h"
#include "pqheap.h"

// Make a shallow copy of a kd-tree for usage by another thread.
kdtree_t * kdtree_copy(kdtree_t * );


void kdtree_free(kdtree_t ** _T)
{
    kdtree_t * T = _T[0];

    free(T->idx_local);
    free(T->X_local);
    free(T->root);
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

static inline double eudist(const double * A, const double * B)
{
    return sqrt(pow(A[0]-B[0],2) + pow(A[1]-B[1], 2));
}

static inline double eudist2(const double * A, const double * B)
{
    return pow(A[0]-B[0],2) + pow(A[1]-B[1], 2);
}

double get_median(const double * X, size_t N, double * T)
{
    // T is a temporary buffer, should be N elements large
    // https://www.gnu.org/software/gsl/doc/html/statistics.html
    // quickselect
    for(size_t kk = 0; kk < N; kk++)
    {
        T[kk] = X[2*kk];
    }
    double median =  gsl_stats_median(T, 1, N);
    return median;
}



void get_min_max(const double * restrict X, const size_t N,
                 double * restrict xa, double * restrict xb,
                 double * restrict ya, double * restrict yb)
{
    double _xa = X[0];
    double _xb = X[0];
    double _ya = X[1];
    double _yb = X[1];
    for(size_t kk = 0; kk<N; kk++)
    {
        double x = X[2*kk];
        if(x < _xa)
        {
            _xa = x;
        }
        if(x > _xb)
        {
            _xb = x;
        }
        double y = X[2*kk+1];
        if(y < _ya)
        {
            _ya = y;
        }
        if(y > _yb)
        {
            _yb = y;
        }
    }
    xa[0] = _xa;
    xb[0] = _xb;
    ya[0] = _ya;
    yb[0] = _yb;
    return;
}

#if 0
void node_shuffle(double * X, size_t * idx, size_t kk)
{
}
#endif

// Sort the elements in the final nodes radially ?
// or unsort them?
void kdtree_split(kdtree_t * T, double * X, size_t * idx, size_t n_points, double xmin, double xmax, double ymin, double ymax)
{
    kdtree_node_t * node = T->root + T->next;

    node->number = T->next;
    node->n_points = n_points; // number of elements to split on
    node->xmin = xmin;
    node->xmax = xmax;
    node->ymin = ymin;
    node->ymax = ymax;

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
            node->X[2*kk] =     X[2*kk];
            node->X[2*kk + 1] = X[2*kk + 1];
        }
        //node_shuffle(node->X, node->idx, kk);

        T->next_idx += n_points;
        return;
    }
    //printf("%zu %zu points, [%f, %f, %f, %f]\n", next[0], N, xmin, xmax, ymin, ymax);

    double dx = xmax-xmin;
    assert(dx > 0);
    double dy = ymax-ymin;
    assert(dy > 0);
    int var = 0; // variable to split on
    if(dy > dx)
    {
        var = 1;
    }
    node->var = var;

    double pivot, rxmin, rxmax, rymin, rymax, lxmin, lxmax, lymin, lymax = 0;

    if(var == 0)
    {
        pivot = get_median(X, n_points, T->median_buffer);
        //pivot = xmin + 0.5*dx;
        rxmin = pivot; rxmax = xmax; rymin = ymin; rymax = ymax;
        lxmin = xmin; lxmax = pivot; lymin = ymin; lymax = ymax;
        if(pivot == xmin || pivot == xmax)
        { // Avoid infinite recursion
            goto final;
        }
    } else {
        //pivot = ymin + 0.5*dy;
        pivot = get_median(X+1, n_points, T->median_buffer);
        rxmin = xmin; rxmax = xmax; rymin = pivot; rymax = ymax;
        lxmin = xmin; lxmax = xmax; lymin = ymin; lymax = pivot;
        if(pivot == ymin || pivot == ymax)
        {
            goto final;
        }
    }
    node->pivot = pivot;

    double * Xleft = malloc(2*n_points*sizeof(double));
    size_t * idx_left = malloc(n_points*sizeof(double));
    double * Xright = malloc(2*n_points*sizeof(double));
    size_t * idx_right = malloc(n_points*sizeof(double));

    size_t rpos = 0;
    size_t lpos = 0;
    for(size_t kk = 0; kk < n_points; kk++)
    {
        if(X[2*kk+var] > pivot)
        {
            Xright[2*rpos] = X[2*kk];
            Xright[2*rpos+1] = X[2*kk+1];
            idx_right[rpos] = idx[kk];
            rpos++;
        } else {
            Xleft[2*lpos] = X[2*kk];
            Xleft[2*lpos+1] = X[2*kk+1];
            idx_left[lpos] = idx[kk];
            lpos++;
        }
    }
    assert(lpos > 0); // if this is the case, this should be
    assert(rpos > 0); // made a final node (duplicate points...)

    assert(lpos+rpos == n_points);
    T->next++;
    node->node_right = T->next;;
    kdtree_split(T, Xright, idx_right, rpos,
                 rxmin, rxmax, rymin, rymax);
    free(Xright);
    free(idx_right);

    T->next++;
    node->node_left = T->next;;
    kdtree_split(T, Xleft, idx_left, lpos,
                 lxmin, lxmax, lymin, lymax);
    free(Xleft);
    free(idx_left);
    return;
}

kdtree_t *
kdtree_new(const double * X, size_t N, int binsize)
{
    /// Set up the tree
    kdtree_t * T = malloc(sizeof(kdtree_t));
    T->k = -1; // not set up for queries of any size
    kdtree_node_t * nodes = malloc(N*sizeof(kdtree_node_t));
    T->root = nodes;
    T->n_nodes = 0;
    T->idx_local = malloc(N*sizeof(size_t));
    T->X_local = malloc(2*N*sizeof(size_t));
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
    double xmin, xmax, ymin, ymax = 0;
    get_min_max(X, N, &xmin, &xmax, &ymin, &ymax);
    // Expand the region for correct ball-within-region checking
    double dx = xmax-xmin;
    double dy = ymax-ymin;
    xmin -= 2*dx; xmax += 2*dx;
    ymin -= 2*dy; ymax += 2*dy;

    T->median_buffer = malloc(N*sizeof(double));

    kdtree_split(T, X, idx_list, N, xmin, xmax, ymin, ymax);

    free(T->median_buffer);
    T->median_buffer = NULL;
    T->n_nodes = T->next;
    free(idx_list);

    return T;
}

size_t kdtree_query_closest(kdtree_t * T, double * X)
{

    kdtree_node_t * N = T->root;
    while(N->node_left > 0)
    {
        int var = N->var;
        if(X[var] > N->pivot)
        {
            N = T->root+N->node_right;
        } else {
            N = T->root+N->node_left;
        }
    }
    assert(N->node_left == -1);
    assert(N->node_left == -1);
    double dmin = eudist(X, N->X);
    size_t imin = N->idx[0];
    for(size_t kk = 0; kk<N->n_points; kk++)
    {
        double d = eudist(X, N->X + 2*kk);
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
    if(x + r > node->xmax || x - r < node->xmin ||
       y + r > node->ymax || y - r < node->ymin)
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

#ifdef aq
    double rmax = aq_get_max_value(aq);
#else
    const double rmax = pqheap_get_max_value(T->pq);
#endif


    // find nearest x and nearest y
    double nx = 0;
    double ny = 0;
    const double x = Q[0];
    const double y = Q[1];
    if( x < node->xmin)
    {
        nx = node->xmin - x;
    } else if (x > node->xmax)
    {
        nx = x - node->xmax;
    }

    if( y < node->ymin)
    {
        nx = node->ymin - y;
    } else if (y > node->ymax)
    {
        ny = y - node->ymax;
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
    printf("Visiting node %6zu #=%7zu [% f, % f, % f, % f]\n", node->number, node->n_points,
           node->xmin, node->xmax, node->ymin, node->ymax);
    }

    #ifdef aq

    #else
    pqheap_t * pq = T->pq;
    #endif

    // If final node
    if(node->node_left == -1)
    {
        T->direct_path = 0;
        const double * X = node->X;
        const size_t * idx = node->idx;

        // Add all points
        for(size_t kk = 0; kk<node->n_points; kk++)
        {
            double d2 = eudist2(X+2*kk, Q);
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
    int var = node->var;
    int done = 0;
    if(Q[var] > node->pivot)
    {
        // correct direction
        if(T->direct_path || bounds_overlap_ball(T, T->root + node->node_right, Q))
        {
        done = kdtree_search(T, T->root + node->node_right, Q);
        if(done == 1)
        {
            return done;
        }
        }
        // "wrong direction"
        if(bounds_overlap_ball(T, T->root + node->node_left, Q))
        {
            done = kdtree_search(T, T->root + node->node_left, Q);
            if(done == 1)
            {
                return done;
            }
        }
    } else {
        // "correct" direction
        if(T->direct_path || bounds_overlap_ball(T, T->root + node->node_left, Q))
        {
        done = kdtree_search(T, T->root + node->node_left, Q);
        if(done)
        {
            return 1;
        }
        }

        // "wrong" direction
        if(bounds_overlap_ball(T, T->root + node->node_right, Q))
        {
          done = kdtree_search(T, T->root + node->node_right, Q);
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
    kdtree_search(T, T->root, Q);

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
