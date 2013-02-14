/**
 *  \file csr.h
 *
 *  \brief Implements a serial sparse matrix-vector multiply for a
 *  matrix using compressed sparse row (CSR) storage.
 */

#if !defined (INC_CSR_H)
#define INC_CSR_H /*!< csr.h included. */

typedef struct
{
  int m;       /*! no. of rows */
  int nnz;     /*!< no. of stored non-zeros */
  int* ptr;    /*! row pointers, 0-based indices */
  int* ind;    /*! column indices, 0-based */
  double* val; /*! stored values */
} csr_t;

/**
 *  \brief Instantiates a CSR matrix using the data given.
 *  \note This routine makes a __shallow-copy__ of the inputs.
 */
void csr_instantiate (int m, int nnz, int* ptr, int* ind, double* val,
		      csr_t* A);

/**
 *  \brief Reads the matrix from the specified file, which must be in
 *  Harwell-Boeing format, and returns it as a handle to a CSR
 *  matrix. If the matrix is stored symmetrically, it is expanded to
 *  full format, i.e., the non-zeros in both the lower and upper
 *  triangles are represented explicitly.
 */
void csr_readHB (const char* filename, csr_t* B);

/**
 *  \brief Performs a sparse matrix-vector multiply, y <- A*x, on a
 *  matrix A stored in CSR format.
 */
void csr_matvec__sequential (double* y, const csr_t* A, const double* x);
void csr_matvec__parfor (double* y, const csr_t* A, const double* x);
void csr_matvec__recspawn (double* y, const csr_t* A, const double* x);
void csr_matvec__segscan (double* y, const csr_t* A, const double* x);

/**
 *  \brief Free memory associated with a matrix structure.
 *  \note Does not free 'A' itself.
 */
void csr_free (csr_t* A);

#endif

/* eof */
