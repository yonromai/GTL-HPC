#ifndef H_SCAN
#define H_SCAN

#include <stdlib.h>
#include <math.h>
// #include <cilk/cilk.h>
#define CILK_NWORKERS 2
// #define _Cilk_for for

typedef double (*bin_operator) (const double a, const double b);

int scan(double* out, const double* a, const double* b, const unsigned int size, bin_operator plus, bin_operator cross, bin_operator companion);

#endif
