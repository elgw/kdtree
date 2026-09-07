#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

#include "kdtree.h"
#include "../src/pqheap.h"

typedef double fxx;
typedef int8_t i8;
typedef int32_t i32;
typedef uint32_t u32;
typedef uint64_t u64;

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


/* For doing parallel queries */

/* Query the nQ points in Q for the k nearest neighbors
 *
 * The returned matrix is kxN elements large and should be freed by
 * the caller.
 */
size_t * kdtree_query_knn_multi(kdtree_t * T,
                                const double * Q, size_t nQ,
                                int k, int ntheads);


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
    size_t * KNN = calloc(nQ*k, sizeof(size_t));
    assert(KNN != NULL);

    if(nthreads == 1)
    {
        for(size_t kk = 0; kk<nQ; kk++)
        {
            size_t * knn = kdtree_query_knn(T, Q+2*kk, k);
            memcpy(KNN + k*kk, knn, k*sizeof(size_t));
        }
        return KNN;
    } else {
        pthread_t * threads = calloc(nthreads, sizeof(pthread_t));
        assert(threads != NULL);
        _p_query_t * confs = calloc(nthreads, sizeof(_p_query_t));
        assert(confs != NULL);

        for(int kk = 0; kk<nthreads; kk++)
        {
            confs[kk].T = kdtree_copy_shallow(T);
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
            kdtree_free_shallow(confs[kk].T);
        }
        free(confs);
        free(threads);
    }

    return KNN;
}


/* Dynamic vector array */
struct dvarray {
    double * data;
    size_t n_used;
    size_t n_alloc;
    int ndim;
};

static void dvarray_n_more(struct dvarray * A, size_t nmore)
{
    if(A->n_used + nmore >= A->n_alloc)
    {
        size_t new_size = A->n_alloc + nmore;
        if(new_size < 1.2 *A->n_alloc)
        {
            new_size = 1.2*A->n_alloc;
        }
        double * t = realloc(A->data, A->ndim*new_size*sizeof(double));
        if(t == NULL) // This makes fanalyzer happy
        {
            free(A->data);
        }
        assert(t != NULL);
        A->data = t;


        assert(A->data != NULL);
        A->n_alloc = new_size;
    }
}

static void dvarray_insert_vector(struct dvarray * A, const double * X)
{
    dvarray_n_more(A, 1);
    memcpy(A->data + A->ndim*A->n_used,
           X,
           A->ndim*sizeof(double));
    A->n_used++;
    return;
}

static void dvarray_free(struct dvarray * A)
{
    free(A->data);
    free(A);
    return;
}

struct dvarray * dvarray_new(size_t n, int ndim)
{
    assert(n > 0);
    struct dvarray * A = calloc(1, sizeof(struct dvarray));
    assert(A != NULL);
    A->ndim = ndim;
    A->data = calloc(ndim*n, sizeof(double));
    assert(A->data != NULL);
    A->n_alloc = n;
    return A;
}

double * rand_points(size_t N, int ndim)
{
    double * X = calloc(ndim*N, sizeof(double));
    assert(X != NULL);
    for(size_t kk = 0; kk<ndim*N; kk++)
    {
        X[kk] = 1000 * (double) rand() / (double) RAND_MAX;
    }
    return X;
}

static double
eudist3_sq(const double * A, const double * B, int ndim)
{
    double d2 = 0;
    for(int kk = 0; kk < ndim; kk++){
        d2 += pow(A[kk]-B[kk], 2);
    }
    return d2;
}

static double
eudist3(const double * A, const double * B, int ndim)
{
    return sqrt(eudist3_sq(A, B, ndim));
}

static double timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

#ifdef __APPLE__
size_t get_peakMemoryKB(void)
{
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    return (size_t) round((double) r_usage.ru_maxrss/1024.0);
}
#endif

#ifndef __APPLE__
size_t get_peakMemoryKB(void)
{
    char * statfile = calloc(100, sizeof(char));
    assert(statfile != NULL);
    sprintf(statfile, "/proc/%d/status", getpid());
    FILE * sf = fopen(statfile, "r");
    if(sf == NULL)
    {
        fprintf(stderr, "Failed to open %s\n", statfile);
        free(statfile);
        return 0;
    }

    char * peakline = NULL;

    char * line = NULL;
    size_t len = 0;

    while( getline(&line, &len, sf) > 0)
    {
        if(strlen(line) > 6)
        {
            if(strncmp(line, "VmPeak", 6) == 0)
            {
                free(peakline);
                peakline = strdup(line);
                assert(peakline != NULL);
            }
        }
    }
    free(line);
    fclose(sf);
    free(statfile);

    // Parse the line starting with "VmPeak"
    // Seems like it is always in kB
    // (reference: fs/proc/task_mmu.c)
    // actually in kiB i.e., 1024 bytes
    // since the last three characters are ' kb' we can skip them and parse in between
    size_t peakMemoryKB = 0;
    //  printf("peakline: '%s'\n", peakline);
    if(peakline == NULL)
    {
        return 0;
    }
    if(strlen(peakline) > 11)
    {
        peakline[strlen(peakline) -4] = '\0';

        //    printf("peakline: '%s'\n", peakline+7);
        peakMemoryKB = (size_t) atol(peakline+7);
    }

    free(peakline);
    return peakMemoryKB;
}
#endif

void fprint_peakMemory(FILE * fout)
{
    size_t pm = get_peakMemoryKB();

    if(fout == NULL) fout = stdout;
    fprintf(fout, "peakMemory: %zu kiB\n", pm);

    return;
}


bool is_radially_sorted(const double * X, const double * Q, const size_t * N,
                        int ndim,
                        const int k)
{
    // Check that the distance to Q is increasing or constant
    double r0 = eudist3(X + ndim*N[0], Q, ndim);
    for(int kk = 1; kk < k; kk++)
    {
        double r1 = eudist3(X + ndim*N[kk], Q, ndim);
        if(r0 > r1)
        {
            return false;
        }
        r0 = r1;
    }
    return true;
}

bool is_unique(const size_t * N, const int k)
{
    for(int kk = 0; kk < k; kk++)
    {
        for(int ll = kk+1; ll < k; ll++)
        {
            if( N[kk] == N[ll] )
            {
                return false;
            }
        }
    }
    return true;
}

bool found_correct(double * X,
                   size_t N, int ndim, double * Q,
                   size_t * knn, int k)
{
    // if the last two distances are unique
    // i.e., not duplicates we should find
    // exactly k-1 points below r

    double r0 = eudist3(Q, X + ndim*knn[k-2], ndim);
    double r1 = eudist3(Q, X + ndim*knn[k-1], ndim);

    if(r0 == r1)
    {
        printf("last two points at equal distance, not checking if the result is correct\n");
        return true;
    }
    double r = 0.5*(r0+r1);
    int nfound = 0;

    for(size_t kk = 0; kk<N; kk++)
    {
        double rp = eudist3(Q, X + ndim*kk, ndim);

        if(rp < r)
        {
            nfound++;
        }
    }

    if(nfound == k-1)
    {
        return true;
    }
    printf("ERROR: Found %d point(s) < %f, expected %d\n", nfound, r, k-1);
    //assert(nfound == k-1);

    return false;
}

void test_threads(size_t N, int ndim, int k, int binsize)
{
    printf("\n--> test_threads(N=%zu, d=%d, binsize=%d)\n",
           N, k, binsize);
    double * X = rand_points(N, ndim);

    kdtree_t * T = kdtree_new(X, N, ndim, binsize);
    assert(T != NULL);
    // Timing with 1, ... 8 threads
    for(int nthreads = 1; nthreads < 9; nthreads++)
    {
        struct timespec tstart, tend;
        clock_gettime(CLOCK_REALTIME, &tstart);
        size_t * KNN = kdtree_query_knn_multi(T, X, N, k, nthreads);
        free(KNN);
        clock_gettime(CLOCK_REALTIME, &tend);
        double t_query = timespec_diff(&tend, &tstart);
        printf("Took %f s using %d threads\n", t_query, nthreads);
    }
    kdtree_free(T);
    free(X);
    return;
}

void basic_tests(size_t N, int ndim, int max_leaf_size)
{
    printf("\n--> basic_tests(N=%zu, ndim=%d, max_leaf_size=%d)\n",
           N, ndim, max_leaf_size);
    double * X = rand_points(N, ndim);

    printf("Create tree with zero points\n");
    kdtree_t * T = kdtree_new(NULL, 0, ndim, 10);
    assert(T == NULL);
    kdtree_free(T);
    T = NULL;

    printf("Create and free a Tree\n");
    T = kdtree_new(X, N, ndim, max_leaf_size);
    if(T == NULL)
    {
        printf("Could not construct a kd-tree\n");
        exit(EXIT_FAILURE);
    }
    kdtree_validate(T);

    kdtree_free(T); T = NULL;
    printf("done\n");

    free(X); X= NULL;

    printf("-- All points identical\n");
    X = calloc(N*ndim, sizeof(double));
    T = kdtree_new(X, N, ndim, max_leaf_size);
    size_t * idx = kdtree_query_knn(T, X, 5);
    printf("%zu, %zu, %zu, %zu, %zu",
           idx[0], idx[1], idx[2], idx[3], idx[4]);
    free(X); X = NULL;
    kdtree_free(T); T = NULL;
}

void kde_mean_ref(const double * X,
                  size_t N,
                  int ndim,
                  const double * Q,
                  double sigma,
                  double * mean_ref)
{
    double xmeank[3] = {0};
    double meank = 0;
    for(size_t kk = 0; kk < N; kk++)
    {
        double r = eudist3(X+3*kk, Q, ndim);
        double k = exp(-r*r/(2.0*sigma*sigma));
        //printf("r = %f, k = %f\n", r, k);
        meank += k;
        for(int ll = 0; ll < 3; ll++)
        {
            xmeank[ll] += k*X[3*kk + ll];
        }
    }
    //printf("meank = %f ", meank);
    for(int ll = 0; ll < 3; ll++)
    {
        mean_ref[ll] = xmeank[ll] / meank;
    }
}

void test_kdtree_kde_mean(size_t N, int ndim, int max_leaf_size)
{
    printf("\n--> test_kdtree_kde_mean(N=%zu, max_leaf_size=%d)\n",
           N, max_leaf_size);
    double * X = rand_points(N, ndim);


    kdtree_t * T = kdtree_new(X, N, ndim, max_leaf_size);
    if(T == NULL)
    {
        printf("Could not construct a kd-tree\n");
        exit(EXIT_FAILURE);
    }

    double sigma = 200.0;

    for(int kk = 0; kk < (int) N; kk++)
    {
        double mean[3] = {0};
        double * Q = X + kk*3;
        kdtree_kde_mean(T, Q, sigma, 0, mean);

        double mean_ref[3] = {0};
        kde_mean_ref(X,
                     N, ndim,
                     Q,
                     sigma, mean_ref);

        double err = fabs(eudist3(mean, mean_ref, ndim));
        if(err > 1e-4)
        {
            printf("Q = [%f, %f, %f] ", Q[0], Q[1], Q[2]);
            printf("mu = [%f, %f, %f] ", mean[0], mean[1], mean[2]);
            printf("ref = [%f, %f, %f]\n", mean_ref[0], mean_ref[1], mean_ref[2]);
            printf("abs err=%f\n", err);
            printf("Failure\n");
            exit(EXIT_FAILURE);
        }

    }



    kdtree_free(T); T = NULL;
    printf("done\n");

    free(X); X= NULL;
    kdtree_free(T); T = NULL;


}

static void
print_point(const double * P, int ndim)
{
    printf("(");
    for(int kk = 0; kk < ndim; kk++){
        printf("%f", P[kk]);
        if(kk + 1 < ndim) {
            printf(", ");
        } else {
            printf(")");
        }
    }
    return;
}

void print_query_and_result(const double * X,
                            const double * Q,
                            int ndim,
                            const size_t * idx,
                            size_t k)
{
    printf("Query point: ");
    print_point(Q, ndim);
    printf("\n");

    for(size_t kk = 0; kk< k; kk++)
    {
        printf("#%zu ", idx[kk]);
        print_point(X + ndim*idx[kk], ndim);
        printf("\n");
    }
    return;
}

void test_align_dots(size_t N, int ndim)
{
    printf("\n--> test_align_dots(%zu)\n", N);
    double sigma = 1;
    printf("    sigma=%f\n", sigma);

    double * A = rand_points(N, ndim);
    double * B = rand_points(N, ndim);

    double dx = 4.21;
    double dy = 2.67;
    double dz = -1.001;
    printf("    delta = (%.2f, %.2f, %.2f)\n", dx, dy, dz);
    size_t ncommon = 5;
    ncommon > N ? ncommon = N : 0;
    for(size_t kk = 0; kk < ncommon; kk++)
    {
        B[3*kk + 0] = A[3*kk + 0] + dx;
        B[3*kk + 1] = A[3*kk + 1] + dy;
        B[3*kk + 2] = A[3*kk + 2] + dz;
    }

    double pairing_radius = 10;
    printf("    Pairing radius: %f\n", pairing_radius);

    kdtree_t * TA = kdtree_new(A, N, ndim, 10);
    assert(TA != NULL);
    struct dvarray * arr = dvarray_new(N, ndim);
    size_t nfound_total = 0;
    for(size_t kk = 0; kk < N; kk++)
    {
        double * Q = B + 3*kk;
        size_t nfound = 0;
        size_t * A_idx = kdtree_query_radius(TA, Q, pairing_radius, &nfound);
        for(size_t ll = 0; ll < nfound; ll++)
        {
            double * P = A + 3*A_idx[ll];
            double D[3] = {0};
            D[0] = Q[0] - P[0];
            D[1] = Q[1] - P[1];
            D[2] = Q[2] - P[2];
            dvarray_insert_vector(arr, D);
        }
        nfound_total += nfound;
        free(A_idx);
    }
    printf("Found %zu pairs within the capture radius\n", nfound_total);
    kdtree_free(TA); TA = NULL;
    free(A); A = NULL;
    free(B); B = NULL;


    kdtree_t * TD = kdtree_new(arr->data, arr->n_used, ndim, 10);
    assert(TD != NULL);
    dvarray_free(arr);
    arr = NULL;

    double maxkde = 0;
    double maxpos[3] = {0};
    printf("Grid search\n");
    for(double x = -pairing_radius; x <= pairing_radius; x+=0.5) {
        for(double y = -pairing_radius; y <= pairing_radius; y+=0.5) {
            for(double z = -pairing_radius; z <= pairing_radius; z+=0.5) {
                double P[] = {x, y, z};
                double v = kdtree_kde(TD, P, sigma, 0);
                if(v > maxkde)
                {
                    maxkde = v;
                    memcpy(maxpos, P, 3*sizeof(double));
                }
            }
        }
    }

    printf("Max kde: %.1f, at (%.2f, %.2f, %.2f)\n",
           maxkde, maxpos[0], maxpos[1], maxpos[2]);

    /* Refinement over the grid search */
    double rs = 2*sigma; // Region size
    while(rs > 1e-3)
    {
        double center[3];
        memcpy(center, maxpos, 3*sizeof(double));
        for(double x = -rs; x <= rs; x+= rs/5.0) {
            for(double y = -rs; y <= rs; y+= rs/5.0) {
                for(double z = -rs; z <= rs; z+= rs/5.0) {
                    double P[] = {
                        x+center[0],
                        y+center[1],
                        z+center[2]};
                    double v = kdtree_kde(TD, P, sigma, 0);
                    //printf("%f, %f, %f -> %f\n", P[0], P[1], P[2], v);
                    if(v > maxkde)
                    {
                        maxkde = v;
                        memcpy(maxpos, P, 3*sizeof(double));
                    }
                }
            }
        }
        rs /= 2.0;
    }

    printf("Refined position: (%.2f, %.2f, %.2f) (kde=%.1f)\n",
           maxpos[0],maxpos[1], maxpos[2], maxkde);

    kdtree_free(TD);
    TD = NULL;

    return;
}

void test_query_radius(size_t N, int ndim, double radius)
{
    printf("\n--> test_query_radius(N=%zu, ndim=%d, r=%f)\n", N, ndim, radius);
    double * X = rand_points(N, ndim);
    kdtree_t * T = kdtree_new(X, N, ndim, 10);
    assert(T != NULL);
    size_t n = 0;

    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    size_t * idx = kdtree_query_radius(T, X, radius, &n);
    clock_gettime(CLOCK_REALTIME, &tend);
    printf("Found %zu points in %f s\n", n, timespec_diff(&tend, &tstart));
    size_t nshow = 10;
    if(n < nshow)
    {
        nshow = n;
    }
    for(size_t kk = 0; kk < nshow; kk++)
    {
        double * p = X + ndim*idx[kk];
        printf("(%f, %f, %f), r = %f\n", p[0], p[1], p[1], eudist3(p, X, ndim));
    }
    if(nshow < n)
    {
        printf("Showed %zu / %zu points\n", nshow, n);
    }
    clock_gettime(CLOCK_REALTIME, &tstart);
    size_t n_brute_force = 0;
    for(size_t kk = 0 ; kk<N; kk++)
    {
        if( eudist3(X, X+ndim*kk, ndim) < radius)
        {
            n_brute_force++;
        }
    }
    clock_gettime(CLOCK_REALTIME, &tend);
    printf("Brute force found %zu points in %f s\n", n_brute_force, timespec_diff(&tend, &tstart));
    if(n == n_brute_force)
    {
        printf("Agreement with brute force\n");
    } else {
        printf("Error: Brute force found %zu\n", n_brute_force);
        assert(0);
    }
    kdtree_free(T);
    free(idx);
    free(X);
    return;
}

// A number in [-1, 1]
static fxx inbox(void)
{
    return 2*(0.5 -((fxx) rand() / (fxx) RAND_MAX));
}

// n 3D points in [-1, 1]^3
static fxx * random_points(u32 n, u32 ndim)
{
    fxx * X = malloc(ndim*n*sizeof(fxx));
    assert(X != NULL);
    for(u32 kk = 0; kk < ndim*n; kk++) {
        X[kk] = inbox();
    }
    return X;
}

// See the correct collisions are found
// using kdtree_query_radius
static
void test_kdtree_query_radius_vs_bf(size_t N,
                                    int ndim,
                                    __attribute__((unused)) double vq)
{
    printf("test_collision(%zu, %f)\n", N, vq);
    fxx * X = random_points(N, ndim);
    assert(X != NULL);
    fxx radius = 2.0/cbrt(N);
    fxx radius2 = radius*radius;
    kdtree_t * T = kdtree_new(X, N, ndim, 5);
    if(T == NULL) {
        free(X);
        return;
    }
    assert(T != NULL);

    i8 * C = calloc(N, sizeof(i8));
    i8 * C2 = calloc(N, sizeof(i8));
    u64 n_found_total = 0;
    for(u32 kk = 0; kk < N; kk++) {

        for(u32 ll = 0; ll < N; ll++) {
            C[ll] = eudist3_sq(X+ndim*kk, X+ndim*ll, ndim) < radius2;
        }
        memset(C2, 0, N);
        size_t n_found;
        size_t * result = kdtree_query_radius(T, X+3*kk, radius, &n_found);
        for(size_t ll = 0; ll < n_found; ll++){
            C2[result[ll]] = 1;
            n_found_total++;
        }
        free(result);
        for(u32 ll = 0; ll < N; ll++) {
            if(C[ll] != C2[ll]) {
                printf("[%f, %f, %f] vs [%f, %f, %f]\n",
                       X[3*kk], X[3*kk+1], X[3*kk+2],
                       X[3*ll], X[3*ll+1], X[3*ll+2]);
                printf("distance = %f, bf=%d, kdtree=%d\n",
                       eudist3(X+ndim*kk, X+ndim*ll, ndim),
                       C[ll], C2[ll]);
                printf("test_collision failed\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    printf("Found %lu collisions, everything matches brute force\n", n_found_total);
    kdtree_free(T);
    free(C);
    free(C2);
    free(X);
}

typedef struct {
    i8 * C;
    size_t N;
} cb1_struct;

typedef struct {
    size_t N;
} cb2_struct;


void test_cb1(u32 u, u32 v, void* _data){
    cb1_struct * data = (cb1_struct*) _data;
    //printf("collisions? (%u, %u)\n", u, v);
    if(u > v){ u32 t = v; v = u; u = t;}
    data->C[u + data->N*v]++;
}

void test_cb2(__attribute__((unused)) u32 u,
              __attribute__((unused)) u32 v,
              void* _data){
    cb2_struct * data = (cb2_struct*) _data;
    data->N++;
}


void test_kdtree_collide(size_t N, int ndim)
{
    printf("test_kdtree_collide(N=%zu, ndim=%d)\n", N, ndim);
    fxx * X = random_points(N, ndim);
    assert(X != NULL);
    fxx radius = 2.0/cbrt(N);
    fxx radius2 = radius*radius;

    i8 * C1 = calloc(N*N, sizeof(i8));
    i8 * C2 = calloc(N*N, sizeof(i8));
    u64 nC1 = 0;
    for(u32 kk = 0; kk < N; kk++) {
        for(u32 ll = kk+1; ll < N; ll++) {
            C1[kk + N*ll] = eudist3_sq(X+3*kk, X+3*ll, ndim) < radius2;
            nC1 += C1[kk + N*ll];
        }
    }

    struct timespec t0, t1, t2, t3;;
    clock_gettime(CLOCK_REALTIME, &t0);
    kdtree_t * T = kdtree_new(X, N, ndim, 5);
    clock_gettime(CLOCK_REALTIME, &t1);
    if(T == NULL) {
        free(X);
        free(C1);
        free(C2);
        return;
    }

    assert(T != NULL);
    cb1_struct cb_data = {};
    cb_data.N = N;
    cb_data.C = C2;

    clock_gettime(CLOCK_REALTIME, &t2);
    kdtree_collide(T, radius, test_cb1, (void*) &cb_data);
    clock_gettime(CLOCK_REALTIME, &t3);

    u64 nC2 = 0;
    for(u32 kk = 0; kk < N; kk++) {
        for(u32 ll = kk+1; ll < N; ll++) {
            nC2 += C2[kk + N*ll];
        }
    }
    printf("nC2 = %lu\n", nC2);

    // validate
    for(u32 kk = 0; kk < N; kk++) {
        for(u32 ll = kk+1; ll < N; ll++) {
            if(C1[kk+N*ll] != C2[kk+N*ll]) {
                printf("[%f, %f, %f] vs [%f, %f, %f]\n",
                       X[3*kk], X[3*kk+1], X[3*kk+2],
                       X[3*ll], X[3*ll+1], X[3*ll+2]);
                printf("distance = %f, bf=%d, kdtree=%d\n",
                       eudist3(X+3*kk, X+3*ll, ndim),
                       C1[ll], C2[ll]);
                printf("test_collision failed\n");
                printf("query radius=%f\n", radius);
                exit(EXIT_FAILURE);
            }
        }
    }
    printf("Found %lu collisions, everything matches brute force\n", nC1);
    kdtree_free(T);
    free(C1);
    free(C2);
    free(X);
    printf("construct : %.3f [ms]\n", 1000.0*timespec_diff(&t1, &t0));
    printf("query     : %.3f [ms]\n", 1000.0*timespec_diff(&t3, &t2));
}

static void
benchmark_kdtree_collide(size_t N, int ndim)
{
    printf("test_kdtree_collide(%zu)\n", N);
    fxx * X = random_points(N, ndim);
    assert(X != NULL);
    fxx radius = 2.0/cbrt(N);

    struct timespec t0, t1, t2, t3;;
    clock_gettime(CLOCK_REALTIME, &t0);
    kdtree_t * T = kdtree_new(X, N, ndim, 5);
    clock_gettime(CLOCK_REALTIME, &t1);
    if(T == NULL) {
        free(X);
        return;
    }

    assert(T != NULL);
    cb2_struct cb_data = {};
    cb_data.N = 0;


    clock_gettime(CLOCK_REALTIME, &t2);
    kdtree_collide(T, radius, test_cb2, (void*) &cb_data);
    clock_gettime(CLOCK_REALTIME, &t3);

    kdtree_free(T);
    free(X);
    printf("construct : %.3f [ms]\n", 1000.0*timespec_diff(&t1, &t0));
    printf("query     : %.3f [ms]\n", 1000.0*timespec_diff(&t3, &t2));
    return;
}

void gen_benchmark_table(int ndim, int k, int binsize)
{
    printf("\n--> gen_benchmark_table(ndim=%d, k=%d, binsize=%d)\n",
           ndim, k, binsize);
    printf("| method | N    | t_construct [ms] | t_query [ms] | t_total [ms] |\n");
    printf("| ---    | ---: | ---:             | ---:         | --:          |\n");

        for(u32 N = 128; N < 2<<21; N*=2)
        {
            double * X = rand_points(N, ndim);

            struct timespec tstart, tend;
            clock_gettime(CLOCK_REALTIME, &tstart);
            kdtree_t * T = kdtree_new(X, N, ndim, binsize);
            if(T == NULL)
            {
                printf("Could not construct a kd-tree\n");
                exit(EXIT_FAILURE);
            }
            clock_gettime(CLOCK_REALTIME, &tend);
            double t_build_tree = timespec_diff(&tend, &tstart);


            clock_gettime(CLOCK_REALTIME, &tstart);
            size_t dummy = 0;
            for(size_t kk = 0; kk<N; kk++)
            {
                //printf("\n-> Q: %zu (%f, %f)\n", kk, X[2*kk], X[2*kk+1]);
                size_t * knn = kdtree_query_knn(T, X+ndim*kk, k);

                for(int kk = 0; kk<3; kk++)
                {
                    dummy += knn[kk];
                    //printf("%zu ", knn[kk]);
                }
                //printf("\n");
            }
            assert(dummy > 0);
            clock_gettime(CLOCK_REALTIME, &tend);
            double t_scan = timespec_diff(&tend, &tstart);

            printf("| kdtree | %u | %.3f | %.3f | %.3f |\n",
                   N,
                   1000.0*t_build_tree,
                   1000.0*t_scan,
                   1000.0*(t_build_tree + t_scan));
        }
    printf("\n");
}

void benchmark(size_t N, int ndim, int k, int binsize)
{
    printf("\n--> benchmark(N=%zu, k=%d, binsize=%d)\n",
           N, k, binsize);
    double * X = rand_points(N, ndim);

    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    kdtree_t * T = kdtree_new(X, N, ndim, binsize);
    if(T == NULL)
    {
        printf("Could not construct a kd-tree\n");
        exit(EXIT_FAILURE);
    }
    clock_gettime(CLOCK_REALTIME, &tend);
    double t_build_tree = timespec_diff(&tend, &tstart);
    printf("To build the kd-tree took %f s\n", t_build_tree);


    clock_gettime(CLOCK_REALTIME, &tstart);
    printf("-> %d-NN, all vs all\n", k);
    size_t dummy = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        //printf("\n-> Q: %zu (%f, %f)\n", kk, X[2*kk], X[2*kk+1]);
        size_t * knn = kdtree_query_knn(T, X+ndim*kk, k);

        for(int kk = 0; kk<3; kk++)
        {
            dummy += knn[kk];
            //printf("%zu ", knn[kk]);
        }
        //printf("\n");
    }
    assert(dummy > 0);

    clock_gettime(CLOCK_REALTIME, &tend);
    double t_all_knn = timespec_diff(&tend, &tstart);

    printf("To find %d-NN for all %zu points took %f s\n", k, N, t_all_knn);
    printf("Total time: %f s\n", t_build_tree + t_all_knn);


#ifndef NDEBUG
    printf("-> Validation\n");
    for(size_t kk = 0; kk<N; kk++)
    {
        if(kk % 1000 == 0)
        {
            printf("\r %zu / %zu", kk, N); fflush(stdout);
        }
        double * Q = X+ndim*kk;
        size_t * knn = kdtree_query_knn(T, Q, k);
        bool ok = true;
        if(knn != NULL)
        {
            ok = is_radially_sorted(X, Q, knn, ndim, k);
        }
        bool all_ok = true;
        if(!ok)
        {
            printf("\nERROR: Resulting match not ordered radially\n");
            all_ok = false;
        }

        ok = is_unique(knn, k);
        if(!ok)
        {
            printf("\nERROR: Resulting match has duplicates\n");
            all_ok = false;
        }
        ok = found_correct(X, N, ndim, Q, knn, k);
        if(!ok)
        {
            printf("\nERROR: Did not find the correct points\n");
            all_ok = false;
        }
        if(!all_ok)
        {
            print_query_and_result(X, Q, ndim, knn, k);

            printf("\n\nReturning\n\n");
            exit(EXIT_FAILURE);
        }
    }
    printf("\r %zu / %zu\n", N, N);
#endif
    kdtree_free(T);
    free(X);

    fprint_peakMemory(stdout);
}


int main(int argc, char ** argv)
{
    srand((unsigned) time(NULL));

    if(argc > 1) {
        if(strcmp(argv[1], "--table1") == 0) { // 3D
            int k = 5;
            int binsize = 35;
            gen_benchmark_table(2, k, binsize);
            exit(EXIT_SUCCESS);
        }
        if(strcmp(argv[1], "--table2") == 0) { // 2D
            int k = 5;
            int binsize = 35;
            gen_benchmark_table(3, k, binsize);
            exit(EXIT_SUCCESS);
        }
        if(strcmp(argv[1], "--table3") == 0) { // 2D
            int k = 5;
            int binsize = 35;
            gen_benchmark_table(7, k, binsize);
            exit(EXIT_SUCCESS);
        }
    }


    size_t N = 1000;
    int k = 5;
    int binsize = 20;
    if(argc > 1)
    {
        N = atol(argv[1]);
    }
    if(argc > 2)
    {
        k = atoi(argv[2]);
    }
    if(argc > 3)
    {
        binsize = atoi(argv[3]);
    }
    printf("N = %zu, k = %d, binsize = %d\n", N, k, binsize);

    benchmark(N, 3, k, binsize);
    benchmark_kdtree_collide(10000, 3);

    test_kdtree_query_radius_vs_bf(123, 3, 0.0);

    test_kdtree_query_radius_vs_bf(1000, 3, 0.0);
    test_kdtree_query_radius_vs_bf(1234, 3, 0.0);
    for(int ndim = 1; ndim < 6; ndim++)
    {
        basic_tests(N, ndim, binsize);
    }

    test_kdtree_kde_mean(N, 3, binsize);



    if(N > 100000 )
    {
        return EXIT_SUCCESS;
    }

    for(int ndim=1; ndim < 6; ndim++)
    {
        test_query_radius(N, ndim, 1);
        test_query_radius(N, ndim, 10);
        test_query_radius(N, ndim, 100);
        test_query_radius(N, ndim, 100);
    }
    test_align_dots(5000, 3);

    test_threads(N, 3, k, binsize);





    return EXIT_SUCCESS;
}
