#ifndef H_SCAN
#define H_SCAN

#include <stdlib.h>
#include <cilk/cilk.h>
#include <math.h>

typedef double (*bin_operator) (const double a, const double b);

double* scan(const double* a, const double* b, bin_operator plus, bin_operator cross, bin_operator companion);

#endif
