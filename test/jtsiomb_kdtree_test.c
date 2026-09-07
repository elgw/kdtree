#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>


// Test https://github.com/jtsiomb/kdtree
#include "kdtree/kdtree.h"

static double timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
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

int main(int argc, char ** argv)
{
    printf("| method | N    | t_construct [ms] | t_query [ms] | t_total [ms] |\n");
    printf("| ---    | ---: | ---:             | ---:         | --:          |\n");
    size_t n_found_total = 0;
    for(int N = 128; N < 2<<21; N*=2){
        int ndim = 3;
        double * X = rand_points(N, ndim);
        struct timespec t0, t1, t2, t3;
        clock_gettime(CLOCK_REALTIME, &t0);
        void * kd = kd_create(3);
        for(int kk = 0; kk < N; kk++){
            double * x = X + kk*ndim;
            kd_insert3(kd, x[0], x[1], x[2], 0);
        }
        clock_gettime(CLOCK_REALTIME, &t1);

        double radius =  2.0/cbrt(N);
        clock_gettime(CLOCK_REALTIME, &t2);
        for(int kk = 0; kk < N; kk++)
        {
            const double * x = X + kk*ndim;
            void * set = kd_nearest_range3(kd, x[0], x[1], x[2], radius);
            n_found_total += kd_res_size(set);
            kd_res_free(set);
        }
        clock_gettime(CLOCK_REALTIME, &t3);
        double t_create = 1000.0*timespec_diff(&t1, &t0);
        double t_scan = 1000.0*timespec_diff(&t3, &t2);
        double t_total = t_create + t_scan;
        printf("| jtsiomb | %d | %.3f, | %.3f | %.3f |\n", N, t_create, t_scan, t_total);
        kd_free(kd);
        free(X);
    }

    printf("%zu\n", n_found_total); // so not everything is optimized away
    return EXIT_SUCCESS;
}
