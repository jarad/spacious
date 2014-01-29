// classes for working with various covariances matrices
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include "covs.h"

// compute covariance for a single block of locations
void Cov::compute(double *Sigma, double *theta, int n, double *D, bool packed) {
  compute(Sigma, theta, n, D, 0, NULL, NULL, packed);
}

// compute partials for a single block of locations
void Cov::partials(double *P, bool *diag, int param, double *theta, double *thetaT, int n, double *D, bool packed) {
  partials(P, diag, param, theta, thetaT, n, D, 0, NULL, NULL, packed);
}

/*
 * Exponential covariance
 * Form: Sigma(i,j) = theta[0] * I(i == j) + theta[1] * exp(-d/theta[2])
 * - On real scale: Sigma(i,j) = exp(thetaT[0]) * I(i == j) + exp(thetaT[1]) * exp(-d exp(theta[2]) )
 */
CovExp::CovExp() {
	mNparams = 3;
}

void CovExp::compute(double *Sigma, double *theta, int n1, double *D1,
                     int n2, double *D2, double *Dc, bool packed, bool transpose) {
	int i,j;
	int index;
	int n = n1+n2;

	for (i = 0; i < n1; i++) {
		if (packed) { index = symi(i,i); } else { index = usymi(i,i,n); }
		Sigma[index] = theta[0] + theta[1];

		for (j = i+1; j < n1; j++) {
			if (packed) { index = symi(i,j); } else { index = usymi(i,j,n); }
			Sigma[index] = theta[1] * exp(-D1[symi(i,j)]/theta[2]);
		}
	}

	if (n2 > 0) {
		for (i = 0; i < n2; i++) {
			if (packed) { index = symi(i+n1,i+n1); } else { index = usymi(i+n1,i+n1,n); }
			Sigma[index] = theta[0] + theta[1];

			for (j = i+1; j < n2; j++) {
				if (packed) { index = symi(i+n1,j+n1); } else { index = usymi(i+n1,j+n1,n); }
				Sigma[index] = theta[1] * exp(-D2[symi(i,j)]/theta[2]);
			}
		}

		// cross terms
		computeCross(Sigma, theta, n1, n2, Dc, true, packed, transpose);
	}
}

void CovExp::computeCross(double *Sigma, double *theta, int n1, int n2, double *Dc, bool full, bool packed, bool transpose) {
	// fill Sigma with cross terms
	int i,j;
	int index;
	int n = n1+n2;

	if (full) {
		// we are filling in the (n1, n2) block of Sigma
		for (i = 0; i < n1; i++) {
			for (j = 0; j < n2; j++) {
				if (packed) { index = symi(i,j+n1); } else { index = usymi(i,j+n1,n); }

				if (transpose) Sigma[index] = theta[1] * exp(-Dc[j + i*n2]/theta[2]);
				else           Sigma[index] = theta[1] * exp(-Dc[i + j*n1]/theta[2]);
			}
		}
	} else {
		// we are filling in the n1 x n2 matrix Sigma that only has cross terms
		for (i = 0; i < n1; i++) {
			for (j = 0; j < n2; j++) {
				if (transpose) Sigma[i + j*n1] = theta[1] * exp(-Dc[j + i*n2]/theta[2]);
				else           Sigma[i + j*n1] = theta[1] * exp(-Dc[i + j*n1]/theta[2]);
			}
		}
	}
}

void CovExp::partials(double *P, bool *diag, int param, double *theta, double *thetaT,
                      int n1, double *D1, int n2, double *D2, double *Dc, bool packed) {
	int i,j;
	int index;
	int n = n1+n2;

	switch(param) {
		case 0:   // nugget
			for (i = 0; i < n1; i++) {
				if (packed) { index = symi(i,i); } else { index = usymi(i,i,n); }
				P[index] = theta[0];

				for (j = i+1; j < n1; j++) {
					if (packed) { index = symi(i,j); } else { index = usymi(i,j,n); }
					P[index] = 0;
				}
			}

			if (n2 > 0) {
				for (i = 0; i < n2; i++) {
					if (packed) { index = symi(i+n1,i+n1); } else { index = usymi(i+n1,i+n1,n); }
					P[index] = theta[0];

					for (j = i+1; j < n2; j++) {
						if (packed) { index = symi(i+n1,j+n1); } else { index = usymi(i+n1,j+n1,n); }
						P[index] = 0;
					}
				}

				// fill in cross terms
				for (i = 0; i < n1; i++) {
					for (j = 0; j < n2; j++) {
						if (packed) { index = symi(i,j+n1); } else { index = usymi(i,j+n1,n); }
						P[index] = 0;
					}
				}
			}

			*diag = true;
			break;

		case 1:   // partial sill
			for (i = 0; i < n1; i++) {
				if (packed) { index = symi(i,i); } else { index = usymi(i,i,n); }
				P[index] = theta[1];

				for (j = i+1; j < n1; j++) {
					if (packed) { index = symi(i,j); } else { index = usymi(i,j,n); }
					P[index] = theta[1]*exp(-D1[symi(i,j)]/theta[2]);
				}
			}

			if (n2 > 0) {
				for (i = 0; i < n2; i++) {
					if (packed) { index = symi(i+n1,i+n1); } else { index = usymi(i+n1,i+n1,n); }
					P[index] = theta[1];

					for (j = i+1; j < n2; j++) {
						if (packed) { index = symi(i+n1,j+n1); } else { index = usymi(i+n1,j+n1,n); }
						P[index] = theta[1]*exp(-D2[symi(i,j)]/theta[2]);
					}
				}

				for (i = 0; i < n1; i++) {
					for (j = 0; j < n2; j++) {
						if (packed) { index = symi(i,j+n1); } else { index = usymi(i,j+n1,n); }
						P[index] = theta[1]*exp(-Dc[i+j*n1]/theta[2]);
					}
				}
			}

			*diag = false;
			break;

		case 2:   // range
			for (i = 0; i < n1; i++) {
				if (packed) { index = symi(i,i); } else { index = usymi(i,i,n); }
				P[index] = 0;

				for (j = i+1; j < n1; j++) {
					if (packed) { index = symi(i,j); } else { index = usymi(i,j,n); }
					P[index] = -exp(thetaT[1]+thetaT[2]) * D1[symi(i,j)] * exp(-D1[symi(i,j)]/theta[2]);
				}
			}

			if (n2 > 0) {
				for (i = 0; i < n2; i++) {
					if (packed) { index = symi(i+n1,i+n1); } else { index = usymi(i+n1,i+n1,n); }
					P[index] = 0;

					for (j = i+1; j < n2; j++) {
						if (packed) { index = symi(i+n1,j+n1); } else { index = usymi(i+n1,j+n1,n); }
						P[index] = -exp(thetaT[1]+thetaT[2]) * D2[symi(i,j)] * exp(-D2[symi(i,j)]/theta[2]);
					}
				}

				for (i = 0; i < n1; i++) {
					for (j = 0; j < n2; j++) {
						if (packed) { index = symi(i,j+n1); } else { index = usymi(i,j+n1,n); }
						P[index] = -exp(thetaT[1]+thetaT[2]) * Dc[i+j*n1] * exp(-Dc[i+j*n1]/theta[2]);
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

/*
 * Matern covariance
 * Form: Sigma(i,j) = theta[0] * I(i == j) + theta[1] * rho(d; theta[2], theta[3])
 *   o rho(d; theta[2], theta[3]) = (2^(theta[3]-1) gammafn(theta[3]))^(-1) * (d/theta[2])^(theta[3]) bessel_k(d/theta[2], theta[3])
 * - On real scale: Sigma(i,j) = exp(thetaT[0]) * I(i == j) + exp(thetaT[1]) * exp(-d exp(theta[2]) )
 */
CovMatern::CovMatern() {
	mNparams = 4;
}

double CovMatern::rho(double d, double *theta, double *work, int partial) {
	// bessel_k isn't thread safe, so we need to provide our own work vector to bessel_k_ex
	double r = pow(d/theta[2], theta[3]) * bessel_k_ex(d/theta[2], theta[3]+partial, 1, work)/(pow(2.0, theta[3]-1.0) * gammafn(theta[3]));

	if (r != r) r = 1;  // handle if NaN

/*
	if (r < 0) {
		MSG("r < 0: %.4f / %.2f, %2.f, %.2f\n",
		    r, pow(d/theta[2], theta[3]), bessel_k_ex(d/theta[2], theta[3], 1, work), pow(2.0, theta[3]-1.0) * gammafn(theta[3]) );
	}
*/

	return(r);
}

void CovMatern::compute(double *Sigma, double *theta, int n1, double *D1,
                     int n2, double *D2, double *Dc, bool packed, bool transpose) {
	int i,j;
	int index;
	int n = n1+n2;
	double work[1+(long)floor(theta[3])];

	for (i = 0; i < n1; i++) {
		if (packed) { index = symi(i,i); } else { index = usymi(i,i,n); }
		Sigma[index] = theta[0] + theta[1];

		for (j = i+1; j < n1; j++) {
			if (packed) { index = symi(i,j); } else { index = usymi(i,j,n); }
			Sigma[index] = theta[1] * rho(D1[symi(i,j)], theta, work);
			if (Sigma[index] != Sigma[index]) {
				MSG("b1 nan: %.2f / %.2f, %.2f, %.2f, %.2f\n", D1[symi(i,j)], theta[0], theta[1], theta[2], theta[3]);
				Sigma[index] = theta[1];
			}
		}
	}

	if (n2 > 0) {
		for (i = 0; i < n2; i++) {
			if (packed) { index = symi(i+n1,i+n1); } else { index = usymi(i+n1,i+n1,n); }
			Sigma[index] = theta[0] + theta[1];

			for (j = i+1; j < n2; j++) {
				if (packed) { index = symi(i+n1,j+n1); } else { index = usymi(i+n1,j+n1,n); }
				Sigma[index] = theta[1] * rho(D2[symi(i,j)], theta, work);
				if (Sigma[index] != Sigma[index]) {
					MSG("(%d,%d) b2 nan: %.20f / %.2f, %.2f, %.2f, %.2f\n", i, j, D2[symi(i,j)], theta[0], theta[1], theta[2], theta[3]);
					Sigma[index] = theta[1];
				}
			}
		}

		// cross terms
		computeCross(Sigma, theta, n1, n2, Dc, true, packed, transpose);
	}
}

void CovMatern::computeCross(double *Sigma, double *theta, int n1, int n2, double *Dc, bool full, bool packed, bool transpose) {
	// fill Sigma with cross terms
	int i,j;
	int index;
	int n = n1+n2;
	double work[1+(long)floor(theta[3])];

	if (full) {
		// we are filling in the (n1, n2) block of Sigma
		for (i = 0; i < n1; i++) {
			for (j = 0; j < n2; j++) {
				if (packed) { index = symi(i,j+n1); } else { index = usymi(i,j+n1,n); }

				if (transpose) Sigma[index] = theta[1] * rho(Dc[j + i*n2], theta, work);
				else           Sigma[index] = theta[1] * rho(Dc[i + j*n1], theta, work);
				if (Sigma[index] != Sigma[index]) {
					MSG("cross nan: %.2f / %.2f, %.2f, %.2f, %.2f\n", Dc[j+i*n2], theta[0], theta[1], theta[2], theta[3]);
					Sigma[index] = theta[1];
				}
			}
		}
	} else {
		// we are filling in the n1 x n2 matrix Sigma that only has cross terms
		for (i = 0; i < n1; i++) {
			for (j = 0; j < n2; j++) {
				if (transpose) Sigma[i + j*n1] = theta[1] * rho(Dc[j + i*n2], theta, work);
				else           Sigma[i + j*n1] = theta[1] * rho(Dc[i + j*n1], theta, work);
			}
		}
	}
}

void CovMatern::partials(double *P, bool *diag, int param, double *theta, double *thetaT,
                         int n1, double *D1, int n2, double *D2, double *Dc, bool packed) {
	int i,j;
	int index;
	int n = n1+n2;
	double work[1+(long)floor(theta[3])];

	switch(param) {
		case 0:   // nugget
			for (i = 0; i < n1; i++) {
				if (packed) { index = symi(i,i); } else { index = usymi(i,i,n); }
				P[index] = theta[0];

				for (j = i+1; j < n1; j++) {
					if (packed) { index = symi(i,j); } else { index = usymi(i,j,n); }
					P[index] = 0;
				}
			}

			if (n2 > 0) {
				for (i = 0; i < n2; i++) {
					if (packed) { index = symi(i+n1,i+n1); } else { index = usymi(i+n1,i+n1,n); }
					P[index] = theta[0];

					for (j = i+1; j < n2; j++) {
						if (packed) { index = symi(i+n1,j+n1); } else { index = usymi(i+n1,j+n1,n); }
						P[index] = 0;
					}
				}

				// fill in cross terms
				for (i = 0; i < n1; i++) {
					for (j = 0; j < n2; j++) {
						if (packed) { index = symi(i,j+n1); } else { index = usymi(i,j+n1,n); }
						P[index] = 0;
					}
				}
			}

			*diag = true;
			break;

		case 1:   // partial sill
			for (i = 0; i < n1; i++) {
				if (packed) { index = symi(i,i); } else { index = usymi(i,i,n); }
				P[index] = theta[1];

				for (j = i+1; j < n1; j++) {
					if (packed) { index = symi(i,j); } else { index = usymi(i,j,n); }
					P[index] = theta[1]*rho(D1[symi(i,j)], theta, work);
				}
			}

			if (n2 > 0) {
				for (i = 0; i < n2; i++) {
					if (packed) { index = symi(i+n1,i+n1); } else { index = usymi(i+n1,i+n1,n); }
					P[index] = theta[1];

					for (j = i+1; j < n2; j++) {
						if (packed) { index = symi(i+n1,j+n1); } else { index = usymi(i+n1,j+n1,n); }
						P[index] = theta[1]*rho(D2[symi(i,j)], theta, work);
					}
				}

				for (i = 0; i < n1; i++) {
					for (j = 0; j < n2; j++) {
						if (packed) { index = symi(i,j+n1); } else { index = usymi(i,j+n1,n); }
						P[index] = theta[1]*rho(Dc[i+j*n1], theta, work);
					}
				}
			}

			*diag = false;
			break;

		case 2:   // range
			for (i = 0; i < n1; i++) {
				if (packed) { index = symi(i,i); } else { index = usymi(i,i,n); }
				P[index] = 0;

				for (j = i+1; j < n1; j++) {
					if (packed) { index = symi(i,j); } else { index = usymi(i,j,n); }
					P[index] = -theta[1]*(D1[symi(i,j)]/theta[2]) * rho(D1[symi(i,j)], theta, work, -1);
				}
			}

			if (n2 > 0) {
				for (i = 0; i < n2; i++) {
					if (packed) { index = symi(i+n1,i+n1); } else { index = usymi(i+n1,i+n1,n); }
					P[index] = 0;

					for (j = i+1; j < n2; j++) {
						if (packed) { index = symi(i+n1,j+n1); } else { index = usymi(i+n1,j+n1,n); }
						P[index] = -theta[1]*(D2[symi(i,j)]/theta[2]) * rho(D2[symi(i,j)], theta, work, -1);
					}
				}

				for (i = 0; i < n1; i++) {
					for (j = 0; j < n2; j++) {
						if (packed) { index = symi(i,j+n1); } else { index = usymi(i,j+n1,n); }
						P[index] = -theta[1]*(Dc[i+j*n1]/theta[2]) * rho(Dc[i+j*n1], theta, work, -1);
					}
				}
			}

			*diag = false;
			break;
	}
}

void CovMatern::transformToReal(double *theta) {
	theta[0] = log(theta[0]);
	theta[1] = log(theta[1]);
	theta[2] = -log(theta[2]);
	theta[3] = log(theta[3]);
}

void CovMatern::transformFromReal(double *theta) {
	theta[0] = exp(theta[0]);
	theta[1] = exp(theta[1]);
	theta[2] = exp(-theta[2]);
	theta[3] = exp(theta[3]);
}
