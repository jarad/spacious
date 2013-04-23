// classes for working with various covariances matrices
#include <math.h>
#include "covs.h"

// compute covariance for a single block of locations
void Cov::compute(double *Sigma, double *theta, int n, double *D, bool packed) {
  compute(Sigma, theta, n, D, 0, NULL, NULL, packed);
}

// compute partials for a single block of locations
void Cov::partials(double *P, bool *diag, int param, double *theta, double *thetaT, int n, double *D) {
  partials(P, diag, param, theta, thetaT, n, D, 0, NULL, NULL);
}

/*
 * Exponential covariance
 * Form: Sigma(i,j) = theta[0] * I(i == j) + theta[1] * exp(-d/theta[2])
 * - On real scale: Sigma(i,j) = exp(thetaT[0]) * I(i == j) + exp(thetaT[1]) * exp(-d exp(theta[2]) )
 */
CovExp::CovExp() {
	mNparams = 3;
}

void CovExp::compute(double *Sigma, double *theta, int n1, double *D1, int n2, double *D2, double *Dc, bool packed) {
	int i,j;
	int index;
	int n = n1+n2;

	for (i = 0; i < n1; i++) {
		if (packed) { index = symi(i,i); } else { index = fsymi(i,i,n); }
		Sigma[index] = theta[0] + theta[1];

		for (j = i+1; j < n1; j++) {
			if (packed) { index = symi(i,j); } else { index = fsymi(i,j,n); }
			Sigma[index] = theta[1] * exp(-D1[symi(i,j)]/theta[2]);
		}
	}

	if (n2 > 0) {
		for (i = 0; i < n2; i++) {
			if (packed) { index = symi(i+n1,i+n1); } else { index = fsymi(i+n1,i+n1,n); }
			Sigma[index] = theta[0] + theta[1];

			for (j = i+1; j < n2; j++) {
				if (packed) { index = symi(i+n1,j+n1); } else { index = fsymi(i+n1,j+n1,n); }
				Sigma[index] = theta[1] * exp(-D2[symi(i,j)]/theta[2]);
			}
		}

		// cross terms
		for (i = 0; i < n1; i++) {
			for (j = 0; j < n2; j++) {
				if (packed) { index = symi(i,j+n1); } else { index = fsymi(i,j+n1,n); }
				Sigma[index] = theta[1] * exp(-Dc[i + j*n1]/theta[2]);
			}
		}
	}
}

void CovExp::partials(double *P, bool *diag, int param, double *theta, double *thetaT,
                      int n1, double *D1, int n2, double *D2, double *Dc) {
	int i,j;

	switch(param) {
		case 0:   // nugget
			for (i = 0; i < n1; i++) {
				P[symi(i,i)] = theta[0];

				for (j = i+1; j < n1; j++) {
					P[symi(i,j)] = 0;
				}
			}

			if (n2 > 0) {
				for (i = 0; i < n2; i++) {
					P[symi(i+n1,i+n1)] = theta[0];

					for (j = i+1; j < n2; j++) {
						P[symi(i+n1,j+n1)] = 0;
					}
				}

				// fill in cross terms
				for (i = 0; i < n1; i++) {
					for (j = 0; j < n2; j++) {
						P[symi(i,j+n1)] = 0;
					}
				}
			}

			*diag = true;
			break;

		case 1:   // partial sill
			for (i = 0; i < n1; i++) {
				P[symi(i,i)] = theta[1];

				for (j = i+1; j < n1; j++) {
					P[symi(i,j)] = theta[1]*exp(-D1[symi(i,j)]/theta[2]);
				}
			}

			if (n2 > 0) {
				for (i = 0; i < n2; i++) {
					P[symi(i+n1,i+n1)] = theta[1];

					for (j = i+1; j < n2; j++) {
						P[symi(i+n1,j+n1)] = theta[1]*exp(-D2[symi(i,j)]/theta[2]);
					}
				}

				if (Dc != NULL) {
					for (i = 0; i < n1; i++) {
						for (j = 0; j < n2; j++) {
							P[symi(i,j+n1)] = theta[1]*exp(-Dc[i+j*n1]/theta[2]);
						}
					}
				}
			}

			*diag = false;
			break;

		case 2:   // range
			for (i = 0; i < n1; i++) {
				P[symi(i,i)] = 0;

				for (j = i+1; j < n1; j++) {
					P[symi(i,j)] = -exp(thetaT[1]+thetaT[2]) * D1[symi(i,j)] * exp(-D1[symi(i,j)]/theta[2]);
				}
			}

			if (n2 > 0) {
				for (i = 0; i < n2; i++) {
					P[symi(i+n1,i+n1)] = 0;

					for (j = i+1; j < n2; j++) {
						P[symi(i+n1,j+n1)] = -exp(thetaT[1]+thetaT[2]) * D2[symi(i,j)] * exp(-D2[symi(i,j)]/theta[2]);
					}
				}

				if (Dc != NULL) {
					for (i = 0; i < n1; i++) {
						for (j = 0; j < n2; j++) {
							P[symi(i,j+n1)] = -exp(thetaT[1]+thetaT[2]) * Dc[i+j*n1] * exp(-Dc[i+j*n1]/theta[2]);
						}
					}
				}
			}

			*diag = false;
			break;
	}
}

void CovExp::transformToReal(double *theta) {
	theta[0] = log(theta[0]);
	theta[1] = log(theta[1]);
	theta[2] = -log(theta[2]);
}

void CovExp::transformFromReal(double *theta) {
	theta[0] = exp(theta[0]);
	theta[1] = exp(theta[1]);
	theta[2] = exp(-theta[2]);
}

