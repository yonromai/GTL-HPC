#ifndef H_SCAN
#define H_SCAN

#include <stdlib.h>
// #include <cilk/cilk.h>
#include <math.h>

//#define cilk_for for

typedef double (*bin_operator) (const double a, const double b);

int scan(double* out, const double* a, const double* b, const unsigned int size, bin_operator plus, bin_operator cross, bin_operator companion);

#endif
