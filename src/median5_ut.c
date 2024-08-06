#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_statistics_double.h>
#include "median5.h"


static double timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

void validate(size_t N)
{

    double * a = malloc(5*sizeof(double));
    for(size_t kk =0 ; kk<N; kk++)
    {
        for(int ll = 0; ll<5; ll++)
        {
            a[ll] = (double) rand() / (double) RAND_MAX;
        }
        double m5 = median5(a);
        int below = 0;
        int above = 0;
        for(int ll = 0; ll<5; ll++)
        {
            if(a[ll] < m5)
            {
                below++;
            }
            if(a[ll] > m5)
            {
                above++;
            }
        }
        if(above > 2 || below  > 2)
        {
            for(int ll = 0; ll<5; ll++)
            {
                printf("%f ", a[ll]);
            }
            printf(" -> %f\n", m5);
            assert(above < 3);
            assert(below < 3);
        }


    }
    return;
}

void bench(size_t N)
{
    double * X = malloc(N*sizeof(double));

    for(size_t kk = 0; kk<N; kk++)
    {
        X[kk] = (double) rand() / (double) RAND_MAX;
    }

    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);

    double * m5 = medians5(X, N, NULL);

    double mom =  gsl_stats_median(m5, 1, N/5);
    clock_gettime(CLOCK_REALTIME, &tend);
    double t0 = timespec_diff(&tend, &tstart);
    printf("median 5-medians of %zu elements took %f s\n", N, t0);
    printf("mom = %f\n", mom);
    clock_gettime(CLOCK_REALTIME, &tstart);
    double m =  gsl_stats_median(X, 1, N);
    clock_gettime(CLOCK_REALTIME, &tend);
    double t1 = timespec_diff(&tend, &tstart);
    printf("median of %zu elements took %f s\n", N, t1);
    printf("m = %f\n", m);

    free(m5);
    free(X);
}

int main(int argc, char ** argv)
{

    size_t N = 1000;
    if(argc > 1)
    {
        N = atol(argv[1]);
    }

    validate(N/10);
    bench(N);
    return EXIT_SUCCESS;
}
