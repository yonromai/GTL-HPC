/**
 *  \file csr.c
 *
 *  \brief Implements support routines for CSR matrices.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Sparse matrix converter (SMC) utilities -- bebop.cs.berkeley.edu */
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>
#include <bebop/smc/csr_matrix.h>

#include "timer.h"
#include "csr.h"

static void make_well_conditioned (csr_t* A);
static void dump_vec__double (int n, const double* x, const char* tag);
static void dump_vec__int (int n, const int* x, const char* tag);

/* ====================================================================
 */
void
csr_instantiate (int m, int nnz, int* ptr, int* ind, double* val,
		 csr_t* A)
{
  assert (A);
  A->m = m;
  A->nnz = nnz;
  A->ptr = ptr;
  A->ind = ind;
  A->val = val;
}

/* ====================================================================
 */

void
csr_readHB (const char* file, csr_t* B)
{
  struct sparse_matrix_t* A = NULL;
  struct stopwatch_t* timer = NULL;
  int err;

  assert (file);

  timer = stopwatch_create ();
  assert (timer);

  fprintf (stderr, "Loading '%s'...", file); fflush (stderr);
  stopwatch_start (timer);
  A = load_sparse_matrix (HARWELL_BOEING, file); assert (A);
  stopwatch_stop (timer);
  fprintf (stderr, "done [%Lg secs].\n", stopwatch_elapsed (timer));

  fprintf (stderr, "Expanding to unsymmetric (if symmetric)..."); fflush (stderr);
  stopwatch_start (timer);
  err = sparse_matrix_expand_symmetric_storage (A); assert (!err);
  stopwatch_stop (timer);
  fprintf (stderr, "done [%Lg secs].\n", stopwatch_elapsed (timer));

  fprintf (stderr, "Converting to CSR..."); fflush (stderr);
  stopwatch_start (timer);
  sparse_matrix_convert (A, CSR);
  stopwatch_stop (timer);
  fprintf (stderr, "done [%Lg secs].\n", stopwatch_elapsed (timer));

  stopwatch_destroy (timer);

  assert (A->format == CSR);
  struct csr_matrix_t* A_csr = (struct csr_matrix_t *)A->repr;
  assert (A_csr);
  assert (A_csr->nnz == (A_csr->rowptr[A_csr->m] - A_csr->rowptr[0]));

  fprintf (stderr, "--> Matrix is %d x %d with %d non-zeros.\n",
	   A_csr->m, A_csr->n, A_csr->nnz);

  /* Copy and return A_csr */
  assert (B);
  B->m = A_csr->m;
  B->nnz = A_csr->nnz;

  B->ptr = (int *)malloc (A_csr->nnz * sizeof (int));
  assert (B->ptr);
  memcpy (B->ptr, A_csr->rowptr, (A_csr->n+1) * sizeof (int));

  B->ind = (int *)malloc (A_csr->nnz * sizeof (int));
  assert (B->ind);
  memcpy (B->ind, A_csr->colidx, A_csr->nnz * sizeof (int));

  B->val = (double *)malloc (A_csr->nnz * sizeof (double));
  assert (B->val);
  memcpy (B->val, A_csr->values, A_csr->nnz * sizeof (double));
  make_well_conditioned (B);

  destroy_sparse_matrix (A);
}

/* ====================================================================
 */

/**
 *  "Hacks" the values of A to make it appear
 *  well-conditioned. Assumes A has full diagonal entries.
 */
static
void
make_well_conditioned (csr_t* A)
{
  if (!A) return;
  int m = A->m;
  for (int i = 0; i < A->m; ++i) {
    double abs_sum = 0.0;
    for (int k = A->ptr[i]; k < A->ptr[i+1]; ++k)
      abs_sum += fabs (A->val[k]);
    for (int k = A->ptr[i]; k < A->ptr[i+1]; ++k) {
      if (A->ind[k] == i) /* diagonal entry */
	A->val[k] = -2.0; /* makes A diagonally dominant */
      else
	A->val[k] /= abs_sum;
    }
  }
}

/* ====================================================================
 */

void
csr_free (csr_t* A)
{
  if (A) {
    if (A->ptr) { free (A->ptr); A->ptr = 0; }
    if (A->ind) { free (A->ind); A->ind = 0; }
    if (A->val) { free (A->val); A->val = 0; }
  }
}

/* ====================================================================
 */

static
void
dump_vec__double (int n, const double* x, const char* tag)
{
  if (tag)
    fprintf (stderr, "=== Double-precision vector: %s [length %d] ===\n", tag, n);
  for (int i = 0; i < n; ++i) {
    fprintf (stderr, " %d:%g", i, x[i]);
  }
  fprintf (stderr, "\n");
}

static
void
dump_vec__int (int n, const int* x, const char* tag)
{
  if (tag)
    fprintf (stderr, "=== Integer vector: %s [length %d] ===\n", tag, n);
  for (int i = 0; i < n; ++i) {
    fprintf (stderr, " %d:%d", i, x[i]);
  }
  fprintf (stderr, "\n");
}

/* eof */
