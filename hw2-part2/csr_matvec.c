/**
 *  \file csr_matvec.c
 *
 *  \brief Implements a serial sparse matrix-vector multiply for a
 *  matrix using compressed sparse row (CSR) storage.
 */

#include <assert.h>
#include <strings.h>
#include "csr.h"

/* === Sequential implementation === */
void
csr_matvec__sequential (double* y, const csr_t* A, const double* x)
{
  assert (A);
  assert ((x && y) || !A->m);

  for (int i = 0; i < A->m; ++i) {
    int k;
    double y_i = 0;
    for (int k = A->ptr[i]; k < A->ptr[i+1]; ++k)
      y_i += A->val[k] * x[A->ind[k]];
    y[i] = y_i;
  }
}

/* === A Cilk Plus implementation based on par-for */
void
csr_matvec__parfor (double* y, const csr_t* A, const double* x)
{
  assert (A);
  assert ((x && y) || !A->m);

  _Cilk_for (int i = 0; i < A->m; ++i) {
    int k;
    double y_i = 0;
    for (k = A->ptr[i]; k < A->ptr[i+1]; ++k)
      y_i += A->val[k] * x[A->ind[k]];
    y[i] = y_i;
  }
}

/* === A Cilk Plus based on recursive spawns */
void
csr_matvec__recspawn (double* y, const csr_t* A, const double* x)
{
  assert (A);
  assert ((x && y) || !A->m);

  if (A->m <= 100) // base case; a "tuning" parameter
    csr_matvec__sequential (y, A, x);
  else {
    csr_t A1, A2;
    int m1 = A->m / 2;
    int m2 = A->m - m1;
    csr_instantiate (m1, A->ptr[m1], A->ptr, A->ind, A->val, &A1);
    csr_instantiate (m2, A->ptr[m2] - A->ptr[m1], A->ptr + m1, A->ind, A->val, &A2);
    _Cilk_spawn csr_matvec__recspawn (y, &A1, x);
    csr_matvec__recspawn (y + m1, &A2, x);
    _Cilk_sync;
  }
}

/* === Insert your segmented scan implementation here === */
void
csr_matvec__segscan (double* y, const csr_t* A, const double* x)
{
  csr_matvec__sequential (y, A, x);
}

/* eof */
