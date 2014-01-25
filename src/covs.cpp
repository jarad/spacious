// classes for working with various covariances matrices
#include <math.h>
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
				if (transpose) Sigma[j + i*n2] = theta[1] * exp(-Dc[j + i*n2]/theta[2]);
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

double CovMatern::rho(double d, double *theta) {
	return(
		pow(d*theta[2], theta[3]) * bessel_k(d*theta[2], theta[3], 1)/(pow(2.0, theta[3]-1.0) * gammafn(theta[3]))
	);
}

void CovMatern::compute(double *Sigma, double *theta, int n1, double *D1,
                     int n2, double *D2, double *Dc, bool packed, bool transpose) {
	int i,j;
	int index;
	int n = n1+n2;

	for (i = 0; i < n1; i++) {
		if (packed) { index = symi(i,i); } else { index = usymi(i,i,n); }
		Sigma[index] = theta[0] + theta[1];

		for (j = i+1; j < n1; j++) {
			if (packed) { index = symi(i,j); } else { index = usymi(i,j,n); }
			Sigma[index] = theta[1] * rho(D1[symi(i,j)], theta);
		}
	}

	if (n2 > 0) {
		for (i = 0; i < n2; i++) {
			if (packed) { index = symi(i+n1,i+n1); } else { index = usymi(i+n1,i+n1,n); }
			Sigma[index] = theta[0] + theta[1];

			for (j = i+1; j < n2; j++) {
				if (packed) { index = symi(i+n1,j+n1); } else { index = usymi(i+n1,j+n1,n); }
				Sigma[index] = theta[1] * rho(D2[symi(i,j)], theta);
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

	if (full) {
		// we are filling in the (n1, n2) block of Sigma
		for (i = 0; i < n1; i++) {
			for (j = 0; j < n2; j++) {
				if (packed) { index = symi(i,j+n1); } else { index = usymi(i,j+n1,n); }

				if (transpose) Sigma[index] = theta[1] * rho(Dc[j + i*n2], theta);
				else           Sigma[index] = theta[1] * rho(Dc[i + j*n1], theta);
			}
		}
	} else {
		// we are filling in the n1 x n2 matrix Sigma that only has cross terms
		for (i = 0; i < n1; i++) {
			for (j = 0; j < n2; j++) {
				if (transpose) Sigma[j + i*n2] = theta[1] * rho(Dc[j + i*n2], theta);
				else           Sigma[i + j*n1] = theta[1] * rho(Dc[i + j*n1], theta);
			}
		}
	}
}

void CovMatern::partials(double *P, bool *diag, int param, double *theta, double *thetaT,
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
/*
theta4 <- tsmooth(theta[4])
mid <- 2*D[in.pair,in.pair]*exp(theta[3])*sqrt(theta4)
rho <- mid^theta4 * besselK(mid, theta4)/(2^(theta4-1) * gamma(theta4))
rho[is.na(rho)] <- 1
exp(theta[2])*rho

		1.0/(pow(2.0, theta[3]-1.0) * gammafn(theta[3])) * pow(d/theta[2], theta[3]) * bessel_k(d/theta[2], theta[3], 1)
*/

			double mid;

			for (i = 0; i < n1; i++) {
				if (packed) { index = symi(i,i); } else { index = usymi(i,i,n); }
				P[index] = theta[1];

				for (j = i+1; j < n1; j++) {
					if (packed) { index = symi(i,j); } else { index = usymi(i,j,n); }
					P[index] = theta[1]*exp(-D1[symi(i,j)]/theta[2]);

					mid = 2*D1[symi(i,j)]*theta[2]*sqrt(theta[3]);
					P[index] = theta[1]*pow(mid, theta[3]) * bessel_k(mid, theta[3], 1)/( pow(2.0, theta[3]-1.0) * gammafn(theta[3]) );
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
