#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

// MSG() takes the place of printf() and Rprintf().
// - Use printf when running from the command line (CLINE)
#ifdef CLINE
	#define MSG printf
#else
	#define MSG Rprintf
#endif

/*
 * NOTES:
 * - Symmetric matrices have packed storage (upper triangle): http://www.netlib.org/lapack/lug/node123.html
 */

// macro to determine location in packed storage of symmetric matrix
// - Note: this assumes a <= b
// - a=0, b=n can be used to identify # of elements needed to store an nxn matrix
//#define symi(a, b) ( (a) + (b)*((b)+1)/2)
// - We need to flexibly call this for a > b, so use inline instead
inline int symi(int a, int b) {
	if (a <= b) {
		return(a + b * (b+1)/2);
	} else {
		return(b + a * (a+1)/2);
	}
}

// full (non-packed) symmetric matrix indexing for upper triangle
inline int fsymi(int a, int b, int n) {
	if (a <= b) {
		return(a + b*n);
	} else {
		return(b + a*n);
	}
}

// function definitions
int chol2inv(int n, double *A);
void disp_sym(double *mat, int lo, int hi);

#endif
