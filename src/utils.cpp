#include <R.h>
#include <R_ext/Lapack.h>

#include "utils.h"

// chol2inv: compute inverse of n by n pd symm mat A using Cholesky
// assumes A is stored in non-packed upper triangular format
int chol2inv(int n, double *A) {
	char uplo = 'U';
	int  info;

	// compute factorization
	dpotrf_(&uplo, &n, A, &n, &info);
	if (info) {
//		MSG("Error with chol(A): info = %d\n", info);
		return(info);
	}

	// complete inverse
	dpotri_(&uplo, &n, A, &n, &info);
	if (info) {
//		MSG("Error with inv(chol(A)): info = %d\n", info);
		return(info);
	}

	return(0);
}

// display triangle of a symmetric matrix from row/col lo to hi (packed format)
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

// display triangle of a symmetric matrix from row/col lo to hi (upper triangle stored)
void disp_usym(double *mat, int lo, int hi, int n) {
	int i,j;

	for (i = lo; i < hi; i++) { 
	  for (j = lo; j < hi; j++) {
	    if (j < i) {
	      MSG("     ");
	    } else {
	      MSG("%.2f ", mat[usymi(i,j,n)]);
	    }
	  }

	  MSG("\n");
	}
}
