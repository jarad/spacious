#include <R.h>
#include <R_ext/Lapack.h>

#include "utils.h"

// chol2inv: compute inverse of n by n pd symm mat A using Cholesky
// assumes A is stored in packed upper triangular format
void chol2inv(int n, double *A) {
	char uplo = 'U';
	int  info;

	dpptrf_(&uplo, &n, A, &info);
	if (info) {
		MSG("Error with chol(A): info = %d\n", info);
	}

	// complete inverse
	dpptri_(&uplo, &n, A, &info);
	if (info) {
		MSG("Error with inv(chol(A)): info = %d\n", info);
	}
}

// display triangle of a symmetric matrix from row/col lo to hi
void disp_sym(double *mat, int lo, int hi) {
	int i,j;

	for (i = lo; i < hi; i++) { 
	  for (j = lo; j < hi; j++) {
	    if (j < i) {
	      MSG("     ");
	    } else {
	      MSG("%.2f ", mat[symi(i,j)]);
	    }
	  }

	  MSG("\n");
	}
}
