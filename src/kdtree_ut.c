#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>

#include "../include/kdtree.h"
#include "pqheap.h"

#define DIM 3

double * rand_points(size_t N)
{
    double * X = malloc(DIM*N*sizeof(double));
    for(size_t kk = 0; kk<DIM*N; kk++)
    {
        X[kk] = 1000 * (double) rand() / (double) RAND_MAX;
    }
    return X;
}

static double eudist3(const double * A, const double * B)
{
    return pow(A[0]-B[0],2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2);
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
  char * statfile = malloc(100*sizeof(char));
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
    kdtree_t * T = kdtree_new(X, N, 3, binsize);
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
    kdtree_free(&T);
    free(X);
    return;
}

void basic_tests(size_t N, int k, int binsize)
{
    double * X = rand_points(N);

    printf("Create and free a Tree\n");
    kdtree_t * T = kdtree_new(X, N, 3, binsize);
    if(T == NULL)
    {
        printf("Could not construct a kd-tree\n");
        exit(EXIT_FAILURE);
    }
    kdtree_validate(T);

    kdtree_free(&T);
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
               Q[DIM*idx[kk]],Q[DIM*idx[kk]+1], Q[DIM*idx[kk]+2]);
    }
    return;
}

void benchmark(size_t N, int k, int binsize)
{
    double * X = rand_points(N);

    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    kdtree_t * T = kdtree_new(X, N, DIM, binsize);
    if(T == NULL)
    {
        printf("Could not construct a kd-tree\n");
        exit(EXIT_FAILURE);
    }
    clock_gettime(CLOCK_REALTIME, &tend);
    double t_build_tree = timespec_diff(&tend, &tstart);
    printf("To build the kd-tree took %f s\n", t_build_tree);

    if(0){
        int next = 0;
        while(next >= 0)
        {
            kdtree_node_t node = T->nodes[next];
            printf("%d, %zu ", next, node.n_points);
            node_print_bbx(&node);
            next = node.node_right;
        }
    }

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



    printf("-> Validation\n");
    for(size_t kk = 0; kk<N; kk++)
    {
        if(kk % 1000 == 0)
        {
            printf("\r %zu / %zu", kk, N); fflush(stdout);
        }
        double * Q = X+DIM*kk;
        size_t * knn = kdtree_query_knn(T, Q, k);
        bool ok = is_radially_sorted(X, Q, knn, k);
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

            exit(EXIT_FAILURE);
        }
    }
    printf("\r %zu / %zu\n", N, N);
    kdtree_free(&T);
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


    basic_tests(N, k, binsize);


    return EXIT_SUCCESS;

    benchmark(N, k, binsize);


    threads(N, k, binsize);



    return EXIT_SUCCESS;
}
