#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "scan.h"

#define A_SIZE 4

double add(const double d1, const double d2) {
	return d1+d2;
}

double diff(const double d1, const double d2) {
	return d1-d2;
}

double mul(const double d1, const double d2) {
	return d1*d2;
}

double cross(const double a, const double b){
  return (((int)b) == 1)?0:a;
}

double companion(const double a, const double b){
  return a || b;
}

// f = [1 if i in ptr else 0 for i in range(len(val))]
int calc_flags(double* out, const int* ptr, const unsigned int size_ptr, const unsigned int size_out) {
  int i = 0;
  
  _Cilk_for(i = 0; i <  size_ptr; ++i) {
    if (ptr[i] < size_out)
      out[ptr[i]] = 1;
  }
  
  return EXIT_SUCCESS;
}

// a = [val[i]*x[ind[i]] for i in range(len(val))]
int calc_values(double* out, const double* val, const double* x, const int* ind, const unsigned int size) {
  int i = 0;
  
  _Cilk_for(i = 0; i <  size; ++i) {
    out[i] = val[i]*x[ind[i]];
  }
  
  return EXIT_SUCCESS;
}

void
csr_matvec__segscan (double* y, const int nnz, const int m, const int* ptr, const double* val, const int* ind, const double* x)
{
  // csr_matvec__sequential (y, A, x);
  
  // We have to use flags of type double if we want to use our implementation of scan 
  double* flags = (double*) calloc(nnz, sizeof(double));
  double* values = (double*) malloc(nnz*sizeof(double));
  double* out = (double*) malloc(nnz*sizeof(double));
  
  calc_flags(flags, ptr, m, nnz);
  calc_values(values, val, x, ind, nnz);
  
  if(scan(out, flags, values, nnz, &add, &cross, &companion) == EXIT_SUCCESS) {
    int i;
    _Cilk_for(i = 1; i <= m; ++i) {
      y[i-1] = out[ptr[i]-1];
    }
  }
  
  free(flags);
  free(values);
  free(out);
}

//int scan(double* out, const double* a, const double* b, const unsigned int size, bin_operator plus, bin_operator cross, bin_operator companion)

int count_nnz (double** A) {
	int i, j;
	int nnz = 0;
	
	for (i = 0; i < A_SIZE; ++i) {
		for (j = 0; j < A_SIZE; ++j) {
			if (A[i][j] != 0) {
				nnz++;
			}
		}
	}
	
	return nnz;
}

void build(double** A, double* val, int* ptr, int* ind, const int m, const int nnz) {
	int i, j;
	int cpt = 0;
	
	for (i = 0; i < m; ++i) {
		ptr[i] = cpt;
		for (j = 0; j < m; ++j) {
			if (A[i][j] != 0) {
				val[cpt] = A[i][j];
				ind[cpt] = j;
				cpt++;
			}
		}
	}
	
	ptr[m] = nnz;
}

int main (void) {
	int i = 0, j = 0;
	double A_array[A_SIZE][A_SIZE] = {{2,0,0,1},{4,3,0,0},{0,0,7,0},{1,0,1,0}};
	double x_array[A_SIZE] = {3,2,-1,5};
	double** A_ptr = NULL;
	double* x_ptr = (double*) malloc (A_SIZE*sizeof(double));
	
	A_ptr = (double**) malloc(A_SIZE*sizeof(double*));
	_Cilk_for (i = 0; i < A_SIZE; ++i) {
		A_ptr[i] = (double*) malloc(A_SIZE*sizeof(double));
	}
	
	for (i = 0; i < A_SIZE; i++) {
		for (j = 0; j < A_SIZE; j++) {
			A_ptr[i][j] = A_array[i][j];
			printf("%lf, ", A_ptr[i][j]);
		}
		printf("\n");
	}
	
	printf("X:\n");
	for (i = 0; i < A_SIZE; i++) {
		x_ptr[i] = x_array[i];
		printf("%lf, ", x_ptr[i]);
	}
	printf("\n");
	
	int nnz = count_nnz(A_ptr);
	double* out = (double*) calloc (A_SIZE,sizeof(double));
	int* ptr = (int*) malloc ((A_SIZE+1)*sizeof(int));
	int* ind = (int*) malloc (nnz*sizeof(int));
	double* val = (double*) malloc (nnz*sizeof(double));
	
	build(A_ptr, val, ptr, ind, A_SIZE, nnz);
	
	/* _Cilk_for (i = 0; i < size; ++i) {
	  a[i] = i;
	  b[i] = i;
	  prod[i] = a[i] * b[i];
	} */

	csr_matvec__segscan (out, nnz, A_SIZE, ptr, val, ind, x_ptr);
	
	for (i = 0; i < A_SIZE; ++i) {
		printf("out[%d] = %lf\n", i, out[i]);
	}
	
	free(out);
	free(ptr);
	free(val);
	free(ind);
	for (i = 0; i < A_SIZE; i++) {
		free(A_ptr[i]);
	}
	free(A_ptr);
	free(x_ptr);
	
	/* double val = 0;
  for (i = 0; i < size; ++i) {
    val += a[i] * b[i];
    printf("%lf == out[%d] = %lf\n", val, i, out[i]);
  } */
	
	/* if (scan(out, prod, NULL, size, &add, NULL, NULL) == EXIT_SUCCESS) {
	  double val = 0;
	  for (i = 0; i < size; ++i) {
	    val += a[i] * b[i];
	    printf("%lf == out[%d] = %lf\n", val, i, out[i]);
	  }
	}
	
	free(out);
	free(a);
	free(b); */
	
	return 0;
}
