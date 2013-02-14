/**
 *  \file cg.h
 *
 *  \brief Interface for a simple serial conjugate gradient solver for
 *  linear systems.
 */

#ifndef INC_CG_H
#define INC_CG_H /*!< cg.h included. */

#include "csr.h"

/** CSR sparse matrix-vector multiply call-back routine. */
typedef void (*matvec_t) (double* y, const csr_t* Adata, const double* x);

/**
 *  \brief Compute alpha*x + y and store in dest (daxpy-like operation)
 *
 *  \param[out] dest - destination vector.  Can be the same as x or y.
 *  \param[in] alpha - scalar multiplier
 *  \param[in] x - vector input
 *  \param[in] y - vector input
 *  \param[in] n - vector size
 */
void axpy (double* dest, double alpha, const double* x, const double* y,
	   int n);

/**
 *  \brief Compute the dot product of two vectors x'*y
 *
 *  \param[in] x - vector input
 *  \param[in] y - vector input
 *  \param[in] n - vector size
 *
 *  \returns The dot product, x'*y
 */
double dot (const double *x, const double *y, int n);

/**
 *  \brief Solve Ax = b using the conjugate-gradient method.
 * 
 *  \param[in] matvec(y, Adata, x, n) - Multiplies A*x into y
 *  \param[in] Adata - Opaque pointer to data used by matvec
 *  \param[in] b - System right-hand side
 *  \param[out] x - Final result
 *  \param[in] rtol - Tolerance on the relative residual (||r||/||b||)
 *  \param[in] n - System size
 *  \param[out] rhist - Residual norm history.  Should be at least
 *  maxiter doubles if it is not NULL.
 *  \param[in] maxiter - Maximum number of iterations before returning
 *
 *  \returns Iteration count on success, -1 on failure
 */
int cg (const matvec_t, const csr_t* Adata, const double* b,
	double* x, double rtol, int n, double* rhist, int maxiter);

#endif                          /* CG_H */

/* eof */
