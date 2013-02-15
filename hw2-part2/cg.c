/**
 *  \file cg.c
 *
 *  \brief Implements a simple serial conjugate gradient solver for
 *  linear systems.
 *
 *  \note To enable debugging at compile-time, #define the flag,
 *  "DEBUG_ME".
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cg.h"
#include "scan.h"

//#define DEBUG_ME

int calc_flags(int* out, const int* ptr, const unsigned int size);
int calc_values(int* out, const double* val, const double* x, const int* ind, const unsigned int size);

double original_dot (const double *x, const double *y, int n);
void par_axpy (double* dest, double alpha, const double* x, const double* y, int n);

int
cg (matvec_t matvec, const csr_t* Adata, const double* b,
    double *x, double rtol, int n, double* rhist, int maxiter)
{
  const int nbytes = n * sizeof(double);

  double bnorm2;              /* ||b||^2 */
  double rnorm2, rnorm2_old;  /* Residual norm squared */
  double rho;
  double alpha;

  double *s;                  /* Search direction */
  double *r;                  /* Residual         */
  double *z;                  /* Temporary vector */

  int i;                      /* Current iteration */

  s = (double *)malloc (nbytes); assert (s);
  r = (double *)malloc (nbytes); assert (r);
  z = (double *)malloc (nbytes); assert (z);

  bnorm2 = original_dot(b, b, n);
  memset (x, 0, nbytes);
  memcpy (r, b, nbytes);
  memcpy (s, b, nbytes);

  rnorm2 = original_dot (r, r, n);

  i = 0;
  do {

#if defined (DEBUG_ME)
    fprintf (stderr, "Iteration %i: ||r||_2 = %g\n", i, sqrt (rnorm2));
#endif

    rnorm2_old = rnorm2;

    matvec (z, Adata, s);
    alpha = rnorm2_old / original_dot(s, z, n);
    par_axpy (x, alpha, s, x, n);
    par_axpy (r, -alpha, z, r, n);
    rnorm2 = original_dot (r, r, n);
    par_axpy (s, rnorm2 / rnorm2_old, s, r, n);

    if (rhist != NULL)
      rhist[i] = sqrt(rnorm2 / bnorm2);
  }
  while (++i < maxiter && rnorm2 > bnorm2 * rtol * rtol);

  free (z);
  free (r);
  free (s);

  if (i >= maxiter)
    return -1;

  return i;
}

void
original_axpy (double* dest, double alpha, const double* x, const double* y, int n)
{
  int i;
  for (i = 0; i < n; ++i) {
    dest[i] = alpha * x[i] + y[i];
  }
}

void
par_axpy (double* dest, double alpha, const double* x, const double* y, int n)
{
  int i,j;
  
// #if defined (DEBUG_ME)
//   fprintf (stderr, "Calling axpy");
// #endif
  double beta = alpha;
  
  _Cilk_for(i = 0; i < n; ++i){
    dest[i] = alpha * x[i] + y[i];
  }

  // Peut etre que le fait qu'on y passe 2 fois plus de temps en faisant
  // deux axpy fait qu'on se mange une erreur?
  
  double * dest2 = malloc(n * sizeof(double));
  _Cilk_sync;
  original_axpy(dest2, beta, x, y, n);
  
  _Cilk_sync;
  for(j=0; j < n; ++j) {
    if(abs(dest2[i] - dest[i]) > 0.0001){
      fprintf (stderr, "ERROR - PAR_AXPY: expected[%d]: %f, actual[%d]: %f\n", i, dest2[i], i, dest[i]);
      fprintf (stderr, "beta=%lf, alpha=%lf, x[%d]: %f, y[%d]: %f\n", beta, alpha, i, x[i], i, y[i]);
    }
   }
   free(dest2);
}

// here are the operators needed by our generic scan implementation
double
add(const double a, const double b){
  return a + b;
}

double
mul(const double a, const double b){
  return a * b;
}

double 
original_dot(const double* x, const double* y, int n){
  int i;
  double sum = 0;
  for (i = 0; i < n; ++i)
    sum += x[i] * y[i];
  return sum;
}

double
par_dot (const double* x, const double* y, int n)
{
  double * prod = malloc(n * sizeof(double));
  double * out = malloc(n * sizeof(double));
  int i;
  double result;
  
  _Cilk_for(i = 0; i < n; ++i) {
    prod[i] = x[i] * y[i];
  }
  
  if(scan(out, prod, NULL, n, &add, NULL, NULL) == EXIT_FAILURE) {
    return EXIT_FAILURE;
  }
  result = out[n-1];
  free(out);
  free(prod);

  double expected = original_dot(x, y, n);
  if(abs(expected - result) > 0.0001)
    fprintf (stderr, "ERROR - PAR_DOT: expected: %f, actual: %f\n", expected, result);

  return expected;
}

/* eof */
