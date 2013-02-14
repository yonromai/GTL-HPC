/**
 *  \file csr_matvec.c
 *
 *  \brief Implements a serial sparse matrix-vector multiply for a
 *  matrix using compressed sparse row (CSR) storage.
 */

#include <assert.h>
#include <strings.h>
#include "csr.h"
#include "scan.h"

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

// here are the operators needed to implement a matvec segscan
double
add(const double a, const double b){
  return a + b;
}

double
cross(const double a, const double b){
  return a * b;
}

double
companion(const double a, const double b){
  return a * b;
}

// f = [1 if i in ptr else 0 for i in range(len(val))]
int
calc_flags(double* out, const int* ptr, const unsigned int size_ptr, const unsigned int size_out) {
  int i = 0;
  
  _Cilk_for(i = 0; i <  size_ptr; ++i) {
    if (ptr[i] < size_out)
      out[ptr[i]] = 1;
  }
  
  return EXIT_SUCCESS;
}

// a = [val[i]*x[ind[i]] for i in range(len(val))]
int
calc_values(double* out, const double* val, const double* x, const int* ind, const unsigned int size) {
  int i = 0;
  
  _Cilk_for(i = 0; i <  size; ++i) {
    out[i] = val[i]*x[ind[i]];
  }
  
  return EXIT_SUCCESS;
}

/* x = [5,2,3]
mat = [[0, 3, 1],[2, 2, 0],[0, 0, 1]]
val = [3, 1, 2, 2, 1]
ind = [1, 2, 0, 1, 2]
ptr = [0, 2, 4, 5]

f = [1 if i in ptr else 0 for i in range(len(val))]

a = [val[i]*x[ind[i]] for i in range(len(val))]

print a
print f

s = Scan(f, a, op_plus, op_cross, op_companion)
c = s.scan()
[c[ptr[i]-1] for i in range(1,len(ptr))]
*/

/* === Insert your segmented scan implementation here === */
void
csr_matvec__segscan (double* y, const csr_t* A, const double* x)
{
  // csr_matvec__sequential (y, A, x);
  
  // We have to use flags of type double if we want to use our implementation of scan 
  double* flags = (double*) calloc(A->nnz, sizeof(double));
  double* values = (double*) malloc(A->nnz*sizeof(double));
  double* out = (double*) malloc(A->nnz*sizeof(double));
  
  calc_flags(flags, A->ptr, A->m, A->nnz);
  calc_values(values, A->val, x, A->ind, A->nnz);
  if(scan(out, flags, values, n, &add, &cross, &companion) == EXIT_SUCCESS) {
    int i = 0;
    _Cilk_for(i = 1; i < A->m; ++i) {
      y[i-1] = out[A->ptr[i]-1];
    }
  }
  
  free(flags);
  free(values);
  free(out);
}

/* eof */
