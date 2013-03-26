// block composite fitting class
#ifndef BLOCK_H
#define BLOCK_H

#include "covs.h"

class BlockComp {
public:
	BlockComp();
	~BlockComp();

	// types of likelihoods we can fit
	enum LikForm { Full, Pair, Block };
	enum CovType { Exp, Matern };

	void setConserveMemory(bool conserve) { mConsMem = conserve; }

	void initPointers();
	void setLikForm(LikForm form);
	void setCovType(CovType type) { mCovType = type; };
	void setData(int n, double *y, double *S, int nblocks, int *B, int p, double *X, int npairs, int *neighbors);
	void setInits(int ntheta, double *theta);

	// fit model with Fisher scoring
	bool fit(bool verbose);

	// methods to compute quantities of interest (requires fit)
	void computeStdErrs();
	void computeCLIC();

private:
	void cleanup();

	void updateBeta(Cov *cov);
	void updateTheta(Cov *cov);

	int     mN;          // number of observations
	double *mY;          // response

	double  *mS;          // spatial locations
	int      mNblocks;    // total number of blocks
	int     *mB;          // block membership
	int     *mNB;         // number of observations in each block
	int    **mWhichB;     // holds indicies for obs in each block
	double **mWithinD;    // within block distance matrices
	double **mBetweenD;   // between block distance matrices

	int     mP;          // number of covariates
	double *mX;          // model matrix

	int  mNpairs;        // number of block pairs
	int *mNeighbors;     // mNpairs by 2 matrix of neighbors
	int  mMaxPair;       // largest number of obs in each pair

	LikForm mLikForm;    // likelihood form
	CovType mCovType;    // covariance type

	bool   mConsMem;     // should memory be conserved?
	bool   mHasFit;      // do we have a fit?

	double mIterTol;     // control parameters
	int    mMaxIter;

	double *mThetaInits; // initial values
	double *mBeta;       // model parametes
	double *mTheta;
	double *mBetaT;      // transformed model parametes
	double *mThetaT;

	// update vars needed elsewhere
	double *mSigma;
	double *mBeta_A;
	double *mBeta_b;
};

#endif
