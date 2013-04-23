#ifdef MAGMA

#include <cuda_runtime_api.h>
#include <cublas.h>
#include "magma.h"
#include "magma_lapack.h"

#include "utils.h"
#include "utils_magma.h"

// magma_chol2inv_gpu: compute inverse of n by n pd symm mat A using Cholesky
// Note: A must be allocated on device
int magma_chol2inv_gpu(int n, double *A) {
	char uplo = 'U';
	int  info;

	// compute factorization
	magma_dpotrf_gpu(uplo, n, A, n, &info);
	if (info) {
		return(info);
	}

	// complete inverse
	magma_dpotri_gpu(uplo, n, A, n, &info);
	if (info) {
		return(info);
	}

	return(0);
}


// display triangle of a symmetric matrix from row/col lo to hi
void magma_disp_sym(double *mat, int n, int lo, int hi) {
	int i,j;

	for (i = lo; i < hi; i++) {
	  for (j = lo; j < hi; j++) {
	    if (j < i) {
	      MSG("     ");
	    } else {
	      MSG("%.2f ", mat[i + j*n]);
	    }
	  }

	  MSG("\n");
	}
}

#endif
