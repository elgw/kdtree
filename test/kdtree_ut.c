#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>


#include "kdtree.h"

typedef double fxx;
typedef int8_t i8;
typedef int32_t i32;
typedef uint32_t u32;
typedef uint64_t u64;

typedef uint32_t kdtree_index;

/* Dynamic vector array */
struct dvarray {
    double * data;
    kdtree_index n_used;
    kdtree_index n_alloc;
    int ndim;
};

static void dvarray_n_more(struct dvarray * A, kdtree_index nmore)
{
    if(A->n_used + nmore >= A->n_alloc)
    {
        kdtree_index new_size = A->n_alloc + nmore;
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

struct dvarray * dvarray_new(kdtree_index n, int ndim)
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

double * rand_points(kdtree_index N, int ndim)
{
    double * X = calloc(ndim*N, sizeof(double));
    assert(X != NULL);
    for(kdtree_index kk = 0; kk<ndim*N; kk++)
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
kdtree_index get_peakMemoryKB(void)
{
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    return (kdtree_index) round((double) r_usage.ru_maxrss/1024.0);
}
#endif

#ifndef __APPLE__
kdtree_index get_peakMemoryKB(void)
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
    kdtree_index peakMemoryKB = 0;
    //  printf("peakline: '%s'\n", peakline);
    if(peakline == NULL)
    {
        return 0;
    }
    if(strlen(peakline) > 11)
    {
        peakline[strlen(peakline) -4] = '\0';

        //    printf("peakline: '%s'\n", peakline+7);
        peakMemoryKB = (kdtree_index) atol(peakline+7);
    }

    free(peakline);
    return peakMemoryKB;
}
#endif

void fprint_peakMemory(FILE * fout)
{
    kdtree_index pm = get_peakMemoryKB();

    if(fout == NULL) fout = stdout;
    fprintf(fout, "peakMemory: %u kiB\n", pm);

    return;
}

static void
basic_tests(kdtree_index N, int ndim, int max_leaf_size)
{
    printf("\n--> basic_tests(N=%u, ndim=%d, max_leaf_size=%d)\n",
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

    kdtree_free(T); T = NULL;
    printf("done\n");

    free(X); X= NULL;

    printf("-- All points identical\n");
    X = calloc(N*ndim, sizeof(double));
    T = kdtree_new(X, N, ndim, max_leaf_size);
    kdtree_index * idx = kdtree_query_knn(T, X, 5);
    printf("%u, %u, %u, %u, %u",
           idx[0], idx[1], idx[2], idx[3], idx[4]);
    free(idx);
    free(X); X = NULL;
    kdtree_free(T); T = NULL;
}

static void
kde_mean_ref(const double * X,
             kdtree_index N,
             int ndim,
             const double * Q,
             double sigma,
             double * mean_ref)
{
    double xmeank[3] = {0};
    double meank = 0;
    for(kdtree_index kk = 0; kk < N; kk++)
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

static void
test_kdtree_kde_mean(kdtree_index N, int ndim, int max_leaf_size)
{
    printf("\n--> test_kdtree_kde_mean(N=%u, max_leaf_size=%d)\n",
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

static void
test_align_dots(kdtree_index N, int ndim)
{
    printf("\n--> test_align_dots(%u)\n", N);
    double sigma = 1;
    printf("    sigma=%f\n", sigma);

    double * A = rand_points(N, ndim);
    double * B = rand_points(N, ndim);

    double dx = 4.21;
    double dy = 2.67;
    double dz = -1.001;
    printf("    delta = (%.2f, %.2f, %.2f)\n", dx, dy, dz);
    kdtree_index ncommon = 5;
    ncommon > N ? ncommon = N : 0;
    for(kdtree_index kk = 0; kk < ncommon; kk++)
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
    kdtree_index nfound_total = 0;
    u32 * A_idx = NULL;
    u32 A_idx_capacity = 0;
    for(kdtree_index kk = 0; kk < N; kk++)
    {
        double * Q = B + 3*kk;
        kdtree_index nfound = kdtree_query_radius(TA, Q, pairing_radius, &A_idx, &A_idx_capacity);
        for(kdtree_index ll = 0; ll < nfound; ll++)
        {
            double * P = A + 3*A_idx[ll];
            double D[3] = {0};
            D[0] = Q[0] - P[0];
            D[1] = Q[1] - P[1];
            D[2] = Q[2] - P[2];
            dvarray_insert_vector(arr, D);
        }
        nfound_total += nfound;
    }
    free(A_idx);

    printf("Found %u pairs within the capture radius\n", nfound_total);
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
void kdtree_query_radius__test(kdtree_index N,
                               int ndim,
                               __attribute__((unused)) double vq)
{
    fxx * X = random_points(N, ndim);
    assert(X != NULL);
    fxx radius = 2.0/cbrt(N);
    fxx radius2 = radius*radius;
    kdtree_t * T = kdtree_new(X, N, ndim, 5);
    if(T == NULL) {
        fprintf(stderr, "Failed to construct a tree\n");
        goto fail;
        return;
    }
    assert(T != NULL);

    i8 * C = calloc(N, sizeof(i8));
    assert(C != NULL);
    i8 * C2 = calloc(N, sizeof(i8));
    assert(C2 != NULL);
    u64 n_found_total = 0;
    u32 * result = NULL;
    u32 result_capacity = 0;
    // For each point, compare the result to brute force
    for(u32 kk = 0; kk < N; kk++) {
        // X + ndim*kk is checked against all other points
        for(u32 ll = 0; ll < N; ll++) {
            C[ll] = eudist3_sq(X+ndim*kk, X+ndim*ll, ndim) < radius2;
        }
        memset(C2, 0, N);
        kdtree_index n_found= kdtree_query_radius(T, X+ndim*kk, radius, &result, &result_capacity);

        for(kdtree_index ll = 0; ll < n_found; ll++){
            C2[result[ll]] = 1;
            n_found_total++;
        }
        for(u32 ll = 0; ll < N; ll++) {
            if(C[ll] != C2[ll]) {
                printf("X[%u]=", kk);
                print_point(X+ndim*kk, ndim);
                printf("X[%u]=", ll);
                print_point(X+ndim*ll, ndim);

                printf("query radius = %f, distance = %f, found by bf=%d, found by kdtree=%d\n",
                       radius,
                       eudist3(X+ndim*kk, X+ndim*ll, ndim),
                       C[ll], C2[ll]);
                goto fail;
            }
        }
    }
    free(result);
    kdtree_free(T);
    free(C);
    free(C2);
    free(X);
    return;
 fail:
    fprintf(stderr, "\nkdtree_query_radius__test FAILED, arguments: N=%u, ndim=%d\n", N, ndim);
    exit(EXIT_FAILURE);
}

typedef struct {
    i8 * C;
    kdtree_index N;
} cb1_struct;

typedef struct {
    kdtree_index N;
} cb2_struct;


static void test_cb1(u32 u, u32 v, void* _data){
    cb1_struct * data = (cb1_struct*) _data;
    //printf("collisions? (%u, %u)\n", u, v);
    if(u > v){ u32 t = v; v = u; u = t;}
    data->C[u + data->N*v]++;
}

static void test_cb2(__attribute__((unused)) u32 u,
                     __attribute__((unused)) u32 v,
                     void* _data){
    cb2_struct * data = (cb2_struct*) _data;
    data->N++;
}


static void test_kdtree_collide(kdtree_index N, int ndim)
{
    printf("test_kdtree_collide(N=%u, ndim=%d)\n", N, ndim);
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
    cb1_struct cb_data = {0};
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
benchmark_kdtree_collide(kdtree_index N, int ndim)
{
    printf("test_kdtree_collide(%u)\n", N);
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
    cb2_struct cb_data = {0};
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

static void
gen_benchmark_table_query_knn(int ndim, int k, int binsize)
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
        kdtree_index dummy = 0;
        for(kdtree_index kk = 0; kk<N; kk++)
        {
            //printf("\n-> Q: %u (%f, %f)\n", kk, X[2*kk], X[2*kk+1]);
            kdtree_index * knn = kdtree_query_knn(T, X+ndim*kk, k);

            for(int kk = 0; kk<3; kk++)
            {
                dummy += knn[kk];
                //printf("%u ", knn[kk]);
            }
            free(knn);
        }
        assert(dummy > 0);
        clock_gettime(CLOCK_REALTIME, &tend);
        double t_scan = timespec_diff(&tend, &tstart);

        printf("| kdtree | %u | %.3f | %.3f | %.3f |\n",
               N,
               1000.0*t_build_tree,
               1000.0*t_scan,
               1000.0*(t_build_tree + t_scan));
        free(X);
        kdtree_free(T);
    }
    printf("\n");
}

static void
gen_benchmark_table_query_distance(int ndim, int binsize)
{

    kdtree_index n_found_total = 0;
    printf("\n--> gen_benchmark_table_query_distance(ndim=%d, binsize=%d)\n",
           ndim, binsize);
    printf("| method | N    | t_construct [ms] | t_query [ms] | t_total [ms] |\n");
    printf("| ---    | ---: | ---:             | ---:         | --:          |\n");

    kdtree_index * knn = NULL;
    u32 knn_capacity = 0;
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
        double radius =  2.0/cbrt(N);

        for(kdtree_index kk = 0; kk<N; kk++)
        {
            //printf("\n-> Q: %u (%f, %f)\n", kk, X[2*kk], X[2*kk+1]);
            kdtree_index n_found = kdtree_query_radius(T, X+ndim*kk, radius, &knn, &knn_capacity);
            n_found_total += n_found;
        }

        clock_gettime(CLOCK_REALTIME, &tend);
        double t_scan = timespec_diff(&tend, &tstart);

        printf("| kdtree | %u | %.3f | %.3f | %.3f |\n",
               N,
               1000.0*t_build_tree,
               1000.0*t_scan,
               1000.0*(t_build_tree + t_scan));
        free(X);
        kdtree_free(T);
    }
    free(knn);
    printf("\n");
    printf("n_found=%u\n", n_found_total);
}


typedef struct {
    u32 idx;
    double d2;
} index_distance;


static int
cmp_index_distance(const void * _A, const void * _B)
{
    index_distance * A = (index_distance*) _A;
    index_distance * B = (index_distance*) _B;
    if(A->d2 < B->d2) {
        return -1;
    }
    if(A->d2 > B->d2) {
        return 1;
    }
    return 0;
}

static int
valid_knn_result(double * X,
                 u32 N, u32 ndim,
                 double * Q,
                 u32 * knn,
                 u32 k)
{
    index_distance * id = malloc(N*sizeof(index_distance));
    for(u32 kk = 0; kk < N; kk++){
        id[kk].idx = kk;
        id[kk].d2 = eudist3_sq(Q, X + kk*ndim, ndim);
    }

    qsort(id, N, sizeof(index_distance), cmp_index_distance);
    for(u32 kk = 0; kk < k; kk++)
    {
        if(id[kk].idx != knn[kk]){
            goto fail;
        }
    }
    return 1;

 fail:
    printf("\n");
    printf("Got:       ");
    for(u32 kk = 0; kk < k; kk++){
        double d = eudist3_sq(Q, X + knn[kk]*ndim, ndim);
        printf("%u(%f) ", knn[kk], sqrt(d));
    }
    printf("\n");
    printf("Should be: ");
    for(u32 kk = 0; kk < k; kk++){
        printf("%u(%f) ", id[kk].idx, sqrt(id[kk].d2));
    }
    printf("\n");

    return 0;
}

static void
kdtree_query_knn__test(kdtree_index N, int ndim, int k, int binsize)
{
    double * X = rand_points(N, ndim);

    kdtree_t * T = kdtree_new(X, N, ndim, binsize);
    if(T == NULL){
        printf("Could not construct a kd-tree\n");
        exit(EXIT_FAILURE);
    }

    for(kdtree_index kk = 0; kk<N; kk++) {
        double * Q = X + ndim*kk;
        kdtree_index * knn = kdtree_query_knn(T, Q, k);
        if(!valid_knn_result(X, N, ndim, Q, knn, k)){
            goto fail;
        }
        free(knn);
    }

    kdtree_free(T);
    free(X);
    return;
 fail:
    printf("\nFAILURE\nkdtree_query_knn (N=%u, ndim=%d, k=%d, binsize=%d)\n",
           N, ndim, k, binsize);
    exit(EXIT_FAILURE);
    return;
}


int main(int argc, char ** argv)
{
    srand((unsigned) time(NULL));

    if(argc > 1) {
        if(strcmp(argv[1], "--table1") == 0) { // 2D
            int k = 5;
            int binsize = 35;
            gen_benchmark_table_query_knn(2, k, binsize);
            exit(EXIT_SUCCESS);
        }
        if(strcmp(argv[1], "--table2") == 0) { // 2D
            int k = 5;
            int binsize = 35;
            gen_benchmark_table_query_knn(3, k, binsize);
            exit(EXIT_SUCCESS);
        }
        if(strcmp(argv[1], "--table3") == 0) { // 7D

            int k = 5;
            int binsize = 35;
            gen_benchmark_table_query_knn(7, k, binsize);
            exit(EXIT_SUCCESS);
        }
        if(strcmp(argv[1], "--table4") == 0) {
            int binsize = 35;
            gen_benchmark_table_query_distance(3, binsize);
            exit(EXIT_SUCCESS);
        }
    }

    printf("Testing 'kdtree_query_knn' with various inputs\n");
    for(int nd = 1; nd < 5; nd++){
        for(int n = 10; n < 300; n*=1.2){
            for(int bs = 1; bs < 6; bs++){
                for(int knn = 1; knn < 10; knn++){
                    printf("\rdimensions = %d, n = %d, bs = %d, knn = %d", nd, n, bs, knn);
                    fflush(stdout);
                    kdtree_query_knn__test(n, nd, knn, bs);
                }
            }
        }
    }
    printf("\r -- ok!                                                               \n");

    // kdtree_query_radius is tested against brute force
    // so keep the number of points low
    printf("Testing 'kdtree_query_radius' with various inputs\n");
    for(int nd = 1; nd < 5; nd++){
        for(int n = 10; n < 1000; n*=1.1){
            for(int bs = 1; bs < 6; bs++){
                printf("\rdimensions = %d, n = %d, bs = %d", nd, n, bs);
                fflush(stdout);
                kdtree_query_radius__test(n, nd, 0.0);
            }
        }
    }
    printf("\r -- ok!                                                               \n");

    return EXIT_SUCCESS;

    // Other things that were in progress at some point
    int N = 1000;
    int binsize = 5;

    benchmark_kdtree_collide(10000, 3);

    for(int ndim = 1; ndim < 6; ndim++)
    {
        basic_tests(N, ndim, binsize);
    }

    for(int bs = 1; bs < 10; bs++)
    {
        test_kdtree_kde_mean(N, 3, bs);

    }

    test_align_dots(5000, 3);

    return EXIT_SUCCESS;
}
