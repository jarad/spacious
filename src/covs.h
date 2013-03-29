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
	// - Sigma accessed with (i+offset, j+offset)
	// - D accessed with (i, j)
	virtual void compute(double *Sigma, int n, double *theta, double *D, int offset) = 0;

	// specialty function to compute cross covariance
	// - Assumes D holds the distances only between these cross terms
	// - D has dimensions n by m
	// - Sigma is assumed to be the full n+m covariance matrix, not just cross terms
	virtual void cross(double *Sigma, int n, int m, double *theta, double *D) = 0;

	// obtain partial derivatives for parameters (in terms of on real line scale)
	// - Fills matrix of patrial derivatives P with respect to param
	// - theta is params on normal scale, thetaT on real line scale
	// - D1 is n1 by n1 holding distances between locations of interest
	// - Optionally, n2/D2 and nc/Dc can be used for filling in partials for a 2nd distance matrix and cross terms
	// - diag returns true if P is diagonal and false otherwise
	virtual void partials(double *P, bool *diag, int param, double *theta, double *thetaT, int n, double *D);
	virtual void partials(double *P, bool *diag, int param, double *theta, double *thetaT,
	                      int n1, double *D1, int n2, double *D2, double *Dc) = 0;

	// transform params to real line
	virtual void transformToReal(double *theta) = 0;
	// transform params from real line
	virtual void transformFromReal(double *theta) = 0;

	// set parameters to be fixed
	void setFixed(int n, int *which, double *values);

protected:
	int     mNparams;
	int     mNfixed;
	int    *mFixed;
	double *mFixedValues;
};

// exponential covariance
class CovExp : public Cov {
public:
	CovExp();

	virtual void compute(double *Sigma, int n, double *theta, double *D, int offset);
	virtual void cross(double *Sigma, int n, int m, double *theta, double *D);
	virtual void partials(double *P, bool *diag, int param, double *theta, double *thetaT,
	                      int n1, double *D1, int n2, double *D2, double *Dc);
	virtual void transformToReal(double *theta);
	virtual void transformFromReal(double *theta);

private:
};

#endif
