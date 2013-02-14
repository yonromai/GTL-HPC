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

  bnorm2 = dot (b, b, n);
  memset (x, 0, nbytes);
  memcpy (r, b, nbytes);
  memcpy (s, b, nbytes);

  rnorm2 = dot (r, r, n);

  i = 0;
  do {

#if defined (DEBUG_ME)
    fprintf (stderr, "Iteration %i: ||r||_2 = %g\n", i, sqrt (rnorm2));
#endif

    rnorm2_old = rnorm2;

    matvec (z, Adata, s);
    alpha = rnorm2_old / dot (s, z, n);
    axpy (x, alpha, s, x, n);
    axpy (r, -alpha, z, r, n);
    rnorm2 = dot (r, r, n);
    axpy (s, rnorm2 / rnorm2_old, s, r, n);

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
axpy (double* dest, double alpha, const double* x, const double* y, int n)
{
  int i;
  cilk_for(i = 0; i < n; ++i)
    dest[i] = alpha * x[i] + y[i];
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
dot (const double *x, const double *y, int n)
{
  double * prod = malloc(n * sizeof(double));
  int i;
  cilk_for(i = 0; i < n; ++i)
    prod[i] = x[i] * y[i];

  double * out = malloc(n * sizeof(double));
  if(scan(out, prod, NULL, n, &add, NULL, NULL) == 0)
    return EXIT_FAILURE;
  double sum = out[n-1];
  free(out);
  free(prod);
  return sum;
}

/* eof */
