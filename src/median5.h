#ifndef __median5_h__
#define __median5_h__

#include <assert.h>
#include <stdlib.h>

// median of five elements using 6 comparisons
double median5(double *);

// returns an array of size nX/5 containing the 5-medians
// the caller is responsible for freeing the returned array
// if M5 != NULL, the result is placed in M5 which has to below
// nX/5 elements.
double * medians5(double * X, size_t nX, double * M5);

#endif
