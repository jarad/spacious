// classes for working with various covariances matrices

#ifndef COVS_H
#define COVS_H

#include <math.h>
#include "utils.h"

// parent for all covariance types
class Cov {
public:
	// return number of params for this covariance type
	virtual int numParams() { return(mNparams); };

	// computes covariance matrix
	// - D1 is n1 by n1 holding distances between locations of interest
	// - Optionally, n2/D2 and nc/Dc can be used for filling in partials for a 2nd distance matrix and cross terms
	// - transpose = true uses Dc'
	virtual void compute(double *Sigma, double *theta, int n, double *D, bool packed = false);
	virtual void compute(double *Sigma, double *theta, int n1, double *D1,
	                     int n2, double *D2, double *Dc, bool packed = false, bool transpose = false) = 0;
	virtual void computeCross(double *Sigma, double *theta, int n1, int n2, double *Dc, bool full, bool packed = false, bool transpose = false) = 0;

	// obtain partial derivatives for parameters (in terms of on real line scale)
	// - Fills matrix of patrial derivatives P with respect to param
	// - theta is params on normal scale, thetaT on real line scale
	// - D1 is n1 by n1 holding distances between locations of interest
	// - Optionally, n2/D2 and nc/Dc can be used for filling in partials for a 2nd distance matrix and cross terms
	// - diag returns true if P is diagonal and false otherwise
	virtual void partials(double *P, bool *diag, int param, double *theta, double *thetaT, int n, double *D, bool packed = false);
	virtual void partials(double *P, bool *diag, int param, double *theta, double *thetaT,
	                      int n1, double *D1, int n2, double *D2, double *Dc, bool packed = false) = 0;

	// transform params to real line
	virtual void transformToReal(double *theta) = 0;
	// transform params from real line
	virtual void transformFromReal(double *theta) = 0;

protected:
	int     mNparams;
};

// exponential covariance
class CovExp : public Cov {
public:
	CovExp();

	virtual void compute(double *Sigma, double *theta, int n1, double *D1,
	                     int n2, double *D2, double *Dc, bool packed = false, bool transpose = false);
	virtual void computeCross(double *Sigma, double *theta, int n1, int n2, double *Dc, bool full, bool packed = false, bool transpose = false);
	virtual void partials(double *P, bool *diag, int param, double *theta, double *thetaT,
	                      int n1, double *D1, int n2, double *D2, double *Dc, bool packed = false);
	virtual void transformToReal(double *theta);
	virtual void transformFromReal(double *theta);

private:
};

// matern covariance
class CovMatern : public Cov {
public:
	CovMatern();

	virtual void compute(double *Sigma, double *theta, int n1, double *D1,
	                     int n2, double *D2, double *Dc, bool packed = false, bool transpose = false);
	virtual void computeCross(double *Sigma, double *theta, int n1, int n2, double *Dc, bool full, bool packed = false, bool transpose = false);
	virtual void partials(double *P, bool *diag, int param, double *theta, double *thetaT,
	                      int n1, double *D1, int n2, double *D2, double *Dc, bool packed = false);
	virtual void transformToReal(double *theta);
	virtual void transformFromReal(double *theta);

private:
	double rho(double d, double *theta, double *work, int partial=0);  // correlation based on distance d and parameters theta
};

#endif
