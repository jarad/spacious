#include <R.h>
#include <R_ext/Lapack.h>

#include "utils.h"

// chol2inv: compute inverse of n by n pd symm mat A using Cholesky
// assumes A is stored in non-packed upper triangular format
int chol2inv(int n, double *A, bool do_log_det, double *log_det) {
	char uplo = 'U';
	int  info;

	// compute factorization
	dpotrf_(&uplo, &n, A, &n, &info);
	if (info) {
//		MSG("Error with chol(A): info = %d\n", info);
		return(info);
	}

	if (do_log_det) {
		// fill in log determinant
		*log_det = 0;

		for (int i = 0; i < n; i++)
			*log_det += log(A[i + i*n]);

		*log_det *= 2;
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

// compare doubles
int compare_double(const void *p1, const void *p2) {
	if ( *(double *)p1 > *(double *)p2 ) {
		return(1);
	} else if ( *(double *)p1 < *(double *)p2 ) {
		return(-1);
	} else {
		return(0);
	}
}

// compare ints
int compare_int(const void *p1, const void *p2) {
	return( ((st_int *)p1)->i - ((st_int *)p2)->i );
}

