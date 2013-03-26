// classes for working with various covariances matrices
#include <math.h>
#include "covs.h"

// set parameters to be fixed
void Cov::setFixed(int n, int *which, double *values) {
	mNfixed      = n;
	mFixed       = which;
	mFixedValues = values;
}

/*
 * Exponential covariance
 * Form: Sigma(i,j) = theta[0] * I(i == j) + theta[1] * exp(-d/theta[2])
 */

// fill in covariance
void CovExp::compute(double *Sigma, int n, double *theta, double *D, int offset) {
	int i,j;

	for (i = 0; i < n; i++) {
		Sigma[symi(i+offset,i+offset)] = theta[0] + theta[1];

		for (j = i+1; j < n; j++) {
			Sigma[symi(i+offset,j+offset)] = theta[1] * exp(-D[symi(i,j)]/theta[2]);
		}
	}
}

// fill in cross terms
void CovExp::cross(double *Sigma, int n, int m, double *theta, double *D) {
	int i,j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			Sigma[symi(i,j+n)] = theta[1] * exp(-D[i + j*n]/theta[2]);
		}

	}
}

// obtain partial derivatives for parameters
void CovExp::partials(double *theta) {
}

void CovExp::transformParams(double *theta) {
}
