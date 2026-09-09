#pragma once
#include <stddef.h>

// quickselect -- find the kth smallest number in X
//
// Arguments:
// X the data points which will be rearranged by the function.
// nX number of elements of X
// k what element to select from sorted X
// n the value of element s of sorted X.

double quickselect(double * X, size_t nX, size_t k);
