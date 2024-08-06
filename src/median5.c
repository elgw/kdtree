#include "median5.h"

static inline double min(double x, double y)
{
    if(x<y)
    {
        return x;
    }
    return y;
}

double median5(double * a)
{
    // Median of 5 elements using 6 comparisons
    // Method by user "James Fingas"
    // https://www.ocf.berkeley.edu/~wwu/cgi-bin/yabb/YaBB.cgi?board=riddles_cs;action=display;num=1061827085
    // Visited 2021-06-14

    // 2. Make sure that a[1] < a[2], a[4]<a[5], a[1]<a[4]
    if(a[0] > a[1])
    {
        double t = a[0];
        a[0] = a[1];
        a[1] = t;
    }

    if (a[3] > a[4])
    {
        double t = a[3];
        a[3] = a[4];
        a[4] = t;
    }
    if(a[0] > a[3])
    {
        double t = a[0];
        a[0] = a[3];
        a[3] = t;
        t = a[1];
        a[1] = a[4];
        a[4] = t;
    }

    assert(a[0] < a[1]);
    assert(a[3] < a[4]);
    assert(a[0] < a[3]);

    // 3.
    if(a[2] > a[1])
    {
        if(a[1] < a[3])
        {
            return min(a[2], a[3]);
        } else {
            return min(a[1], a[4]);
        }
    } else {
        if(a[2] > a[3])
        {
            return min(a[2], a[4]);
        } else {
            return min(a[1], a[3]);
        }
    }
}

double * medians5(double * X, size_t N, double * _m5)
{
    double * m5 = _m5;
    if(m5 == NULL)
    {
        m5 = malloc(N/5*sizeof(double));
    }
    size_t w = 0;
    for(size_t kk = 0; kk<N; kk+=5)
    {
        m5[w++] = median5(X+kk);
    }
    return m5;
}
