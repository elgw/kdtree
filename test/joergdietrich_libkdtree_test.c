// git clone https://github.com/joergdietrich/libkdtree.git

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "kdtree.h"

static double timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}


float * rand_points(size_t N, int ndim)
{
    float * X = calloc(ndim*N, sizeof(float));
    assert(X != NULL);
    for(size_t kk = 0; kk<ndim*N; kk++)
    {
        X[kk] = 1000 * (float) rand() / (float) RAND_MAX;
    }
    return X;
}


int main(int argc, char ** argv)
{
    int ndim = 3;
    int k = 5;
    int nthreads = 1;
    printf("ndim=%d, kNN=%d\n", ndim, k);
    printf("threads=%d\n", nthreads);
    printf("| method | N    | t_construct [ms] | t_query [ms] | t_total [ms] |\n");
    printf("| ---    | ---: | ---:             | ---:         | --:          |\n");
    size_t n_found_total = 0;
    float min[3] = {0, 0, 0};
    float max[3] = {1000, 1000, 1000};
    for(int N = 128; N < 2<<21; N*=2){
        float * X = rand_points(N, ndim);

        int npoints = N;

        struct kd_point * pointlist = calloc(N, sizeof(kd_point));
        for(int kk = 0; kk < N; kk++){
            pointlist[kk].point = X+3*kk;
        }
        struct timespec t0, t1, t2, t3;

        clock_gettime(CLOCK_REALTIME, &t0); // t_create
        struct kdNode *kdTree;
        if ((kdTree = kd_buildTree(pointlist, npoints, NULL, NULL, min, max, ndim,
                                   nthreads)) == NULL) {
            fprintf(stderr, "Error building kd-tree\n");
            exit(EXIT_FAILURE);
        }
        clock_gettime(CLOCK_REALTIME, &t1);

        clock_gettime(CLOCK_REALTIME, &t2); // t_scan
        for(int kk = 0; kk < N; kk++)
        {
            float * point = X + 3*kk;
            float max_dist_sq = 1.0; // max_dist_sq
            struct pqueue *nearest = kd_qnearest(kdTree, point, &max_dist_sq, k, ndim);
            n_found_total += nearest->size;
            free(nearest);
        }
        clock_gettime(CLOCK_REALTIME, &t3);
        double t_create = 1000.0*timespec_diff(&t1, &t0);
        double t_scan = 1000.0*timespec_diff(&t3, &t2);
        double t_total = t_create + t_scan;
        printf("| joergdietrich | %d | %.3f | %.3f | %.3f |\n", N, t_create, t_scan, t_total);
        free(X);
        kd_destroyTree(kdTree, NULL);
    }

    printf("%zu\n", n_found_total); // so not everything is optimized away
    return EXIT_SUCCESS;
}
