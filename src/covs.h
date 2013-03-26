// classes for working with various covariances matrices

#ifndef COVS_H
#define COVS_H

#include <math.h>
#include "utils.h"

// parent for all covariance types
class Cov {
public:
	virtual int numParams() = 0;

	// computes covariance matrix
	// - Sigma accessed with (i+offset, j+offset)
	// - D accessed with (i, j)
	virtual void compute(double *Sigma, int n, double *theta, double *D, int offset) = 0;

	// specialty function to compute cross covariance
	// - Assumes D holds the distances only between these cross terms
	// - D has dimensions n by m
	// - Sigma is assumed to be the full n+m covariance matrix, not just cross terms
	virtual void cross(double *Sigma, int n, int m, double *theta, double *D) = 0;

	virtual void partials(double *theta) = 0;
	virtual void transformParams(double *theta) = 0;

	void setFixed(int n, int *which, double *values);

private:
	int     mNfixed;
	int    *mFixed;
	double *mFixedValues;
};

// exponential covariance
class CovExp : public Cov {
public:
	virtual int numParams() { return 3; }
	virtual void compute(double *Sigma, int n, double *theta, double *D, int offset);
	virtual void cross(double *Sigma, int n, int m, double *theta, double *D);
	virtual void partials(double *theta);
	virtual void transformParams(double *theta);

private:
};

#endif
