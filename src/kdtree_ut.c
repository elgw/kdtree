#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#include "../include/kdtree.h"
#include "pqheap.h"

#define DIM 3

/* Dynamic vector array */
struct dvarray {
    double * data;
    size_t n_used;
    size_t n_alloc;
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
        double * t = realloc(A->data, 3*new_size*sizeof(double));
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
    memcpy(A->data + 3*A->n_used,
           X,
           3*sizeof(double));
    A->n_used++;
    return;
}

static void dvarray_free(struct dvarray * A)
{
    free(A->data);
    free(A);
    return;
}

struct dvarray * dvarray_new(size_t n)
{
    assert(n > 0);
    struct dvarray * A = calloc(1, sizeof(struct dvarray));
    assert(A != NULL);
    A->data = calloc(3*n, sizeof(double));
    assert(A->data != NULL);
    A->n_alloc = n;
    return A;
}

double * rand_points(size_t N)
{
    double * X = calloc(DIM*N, sizeof(double));
    assert(X != NULL);
    for(size_t kk = 0; kk<DIM*N; kk++)
    {
        X[kk] = 1000 * (double) rand() / (double) RAND_MAX;
    }
    return X;
}

static double eudist3(const double * A, const double * B)
{
    return sqrt(pow(A[0]-B[0],2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2));
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


bool is_radially_sorted(const double * X, const double * Q, const size_t * N, const int k)
{
    // Check that the distance to Q is increasing or constant
    double r0 = eudist3(X + DIM*N[0], Q);
    for(int kk = 1; kk < k; kk++)
    {
        double r1 = eudist3(X + DIM*N[kk], Q);
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
                   size_t N, double * Q,
                   size_t * knn, int k)
{
    // if the last two distances are unique
    // i.e., not duplicates we should find
    // exactly k-1 points below r

    double r0 = eudist3(Q, X + DIM*knn[k-2]);
    double r1 = eudist3(Q, X + DIM*knn[k-1]);

    if(r0 == r1)
    {
        printf("last two points at equal distance, not checking if the result is correct\n");
        return true;
    }
    double r = 0.5*(r0+r1);
    int nfound = 0;

    for(size_t kk = 0; kk<N; kk++)
    {
        double rp = eudist3(Q, X + DIM*kk);

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

void threads(size_t N, int k, int binsize)
{
    double * X = rand_points(N);

    kdtree_t * T = kdtree_new(X, N, binsize);
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

void basic_tests(size_t N, int max_leaf_size)
{
    double * X = rand_points(N);

    printf("Create and free a Tree\n");
    kdtree_t * T = kdtree_new(X, N, max_leaf_size);
    if(T == NULL)
    {
        printf("Could not construct a kd-tree\n");
        exit(EXIT_FAILURE);
    }
    kdtree_validate(T);

    kdtree_free(T);
    printf("done\n");

    free(X);

}

void print_query_and_result(const double * X,
                       const double * Q,
                       const size_t * idx,
                       size_t k)
{
    printf("Query point: (%f, %f, %f)\n", Q[0], Q[1], Q[2]);
    for(size_t kk = 0; kk< k; kk++)
    {
        printf("#%zu (%f, %f, %f)\n", idx[kk],
               X[DIM*idx[kk]], X[DIM*idx[kk]+1], X[DIM*idx[kk]+2]);
    }
    return;
}

void test_align_dots(size_t N)
{
    printf("\n--> test_align_dots(%zu)\n", N);
    double * A = rand_points(N);
    double * B = calloc(3*N, sizeof(double));
    assert(B != NULL);

    double dx = 4.2;
    double dy = 2;
    double dz = -1;
    printf("delta = (%f, %f, %f)\n", dx, dy, dz);
    for(size_t kk = 0; kk < N; kk++)
    {
        B[3*kk + 0] = A[3*kk + 0] + dx;
        B[3*kk + 1] = A[3*kk + 1] + dy;
        B[3*kk + 2] = A[3*kk + 2] + dz;
    }

    double pairing_radius = 5;
    printf("Pairing radius: %f\n", pairing_radius);

    kdtree_t * TA = kdtree_new(A, N, 10);
    assert(TA != NULL);
    struct dvarray * arr = dvarray_new(N);
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
            //printf("%f %f %f, %f %f %f, %f %f %f\n", Q[0], Q[1], Q[2], P[0], P[1], P[2], D[0], D[1], D[2]);
            dvarray_insert_vector(arr, D);
        }
        nfound_total += nfound;
        free(A_idx);
    }
    printf("Found %zu pairs within the capture radius\n", nfound_total);
    kdtree_free(TA); TA = NULL;
    free(A); A = NULL;
    free(B); B = NULL;

    printf("\n");
    kdtree_t * TD = kdtree_new(arr->data, arr->n_used, 10);
    assert(TD != NULL);
    kdtree_print_info(TD);
    dvarray_free(arr);
    arr = NULL;

    double maxkde = 0;
    double maxpos[3] = {0};
    printf("Grid search\n");
    for(double x = -pairing_radius; x <= pairing_radius; x+=0.5) {
        for(double y = -pairing_radius; y <= pairing_radius; y+=0.5) {
            for(double z = -pairing_radius; z <= pairing_radius; z+=0.5) {
                double P[] = {x, y, z};
                double v = kdtree_kde(TD, P, 1);
                if(v > maxkde)
                {
                    maxkde = v;
                    memcpy(maxpos, P, 3*sizeof(double));
                }
            }
        }
    }
    // Write out image or numpy array https://developers.google.com/speed/webp/docs/api
    printf("Max kde: %f, at (%f, %f, %f)\n",
           maxkde, maxpos[0], maxpos[1], maxpos[2]);

    /* Refinement over the grid search */
    // GSL ... ?


    kdtree_free(TD);
    TD = NULL;

    /*
    size_t s = 1 + 2*pairing_radius;
    double * M = calloc(s*s, sizeof(double));
    for(size_t xx = 0; xx < s; xx++) {
        for(size_t yy = 0 ; yy < s; yy++) {
            double Q[3] = (double) xx -
            M[xx + s*yy] = kdtree_kde(TD, )
        }
    */

    return;
}

void test_query_radius(size_t N, double radius)
{
    printf("test_query_radius(N=%zu, r=%f)\n", N, radius);
    double * X = rand_points(N);
    kdtree_t * T = kdtree_new(X, N, 10);
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
        double * p = X + 3*idx[kk];
        printf("(%f, %f, %f), r = %f\n", p[0], p[1], p[1], eudist3(p, X));
    }
    if(nshow < n)
    {
        printf("Showed %zu / %zu points\n", nshow, n);
    }
    clock_gettime(CLOCK_REALTIME, &tstart);
    size_t n_brute_force = 0;
    for(size_t kk = 0 ; kk<N; kk++)
    {
        if( eudist3(X, X+3*kk) < radius)
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

void benchmark(size_t N, int k, int binsize)
{
    double * X = rand_points(N);

    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    kdtree_t * T = kdtree_new(X, N, binsize);
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
    for(size_t kk = 0; kk<N; kk++)
    {
        //printf("\n-> Q: %zu (%f, %f)\n", kk, X[2*kk], X[2*kk+1]);
        size_t * knn = kdtree_query_knn(T, X+DIM*kk, k);
        if(0){
            for(int kk = 0; kk<3; kk++)
            {
                printf("%zu ", knn[kk]);
            }
            printf("\n");
        }
        //getchar();
    }

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
        double * Q = X+DIM*kk;
        size_t * knn = kdtree_query_knn(T, Q, k);
        bool ok = true;
        if(knn != NULL)
        {
            ok = is_radially_sorted(X, Q, knn, k);
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
        ok = found_correct(X, N, Q, knn, k);
        if(!ok)
        {
            printf("\nERROR: Did not find the correct points\n");
            all_ok = false;
        }
        if(!all_ok)
        {
            print_query_and_result(X, Q, knn, k);

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

    benchmark(N, k, binsize);

    basic_tests(N, binsize);

    test_query_radius(N, 1);
    test_query_radius(N, 10);
    test_query_radius(N, 100);
    test_query_radius(N, 100);
    test_align_dots(5000);

    threads(N, k, binsize);



    return EXIT_SUCCESS;
}
