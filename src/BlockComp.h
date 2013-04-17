// block composite fitting class
#ifndef BLOCK_H
#define BLOCK_H

#include "covs.h"

#ifdef PTHREAD
typedef struct {
	int              id;
	class BlockComp *bc;
} pair_update_t;
#endif

class BlockComp {
public:
	BlockComp();
	BlockComp(int nthreads);
	~BlockComp();

	// types of likelihoods we can fit
	enum LikForm { Full, Pair, Block, IndBlock };
	enum CovType { Exp, Matern };

	// should we try to conserve memory during fit?
	void setConserveMemory(bool conserve) { mConsMem = conserve; }

	void init(int nthreads);
	void initPointers();
	void setLikForm(LikForm form);
	void setCovType(CovType type);
	void setData(int n, double *y, double *S, int nblocks, int *B, int p, double *X, int npairs, int *neighbors);
	void setInits(double *theta);
	void setFixed(bool *fixed, double *values);

	// fit model with Fisher scoring
	bool fit(bool verbose);

	// compute log likelihood at specified parameters
	double computeLogLik(double *beta, double *theta);

	// methods to compute quantities of interest (requires fit)
	void computeStdErrs();
	void computeCLIC();

	// extract items of interest
	void getBeta(double *beta);
	void getTheta(double *theta);

private:
	void cleanup();

	void setThreads(int nthreads);

	bool updateBeta();
	bool updateBetaPair(int pair, double *Sigma, double *A, double *b);
	bool updateTheta();

	int     mNthreads;   // number of processing threads

	int     mN;          // number of observations
	double *mY;          // response

	double  *mS;          // spatial locations
	int      mNblocks;    // total number of blocks
	int     *mB;          // block membership
	int     *mNB;         // number of observations in each block
	int    **mWhichB;     // holds indicies for obs in each block
	double **mWithinD;    // within block distance matrices
	double **mBetweenD;   // between block distance matrices

	int     mNbeta;      // number of covariates
	double *mX;          // model matrix

	int  mNpairs;        // number of block pairs
	int *mNeighbors;     // mNpairs by 2 matrix of neighbors
	int  mMaxPair;       // largest number of obs in each pair

	LikForm  mLikForm;    // likelihood form
	CovType  mCovType;    // covariance type
	Cov     *mCov;
	int      mNtheta;     // number of covariance parameters

	bool   mConsMem;     // should memory be conserved?
	bool   mHasFit;      // do we have a fit?
	bool   mConverged;   // did the fit converge?
	int    mIters;       // number of iters through fit

	double mIterTol;     // control parameters
	int    mMaxIter;

	double *mThetaInits;  // initial values for covariance params
	int     mNfixed;      // number fixed
	bool   *mFixed;       // which are fixed?
	double *mFixedVals;   // values of fixed params
	double *mBeta;        // model parametes
	double *mTheta;
	double *mThetaT;      // transformed model parametes

	// update vars
	double **mSigma;
	double  *mBeta_A;
	double  *mBeta_b;

	double **mTheta_W;
	double  *mTheta_H;
	double  *mTheta_P;

#ifdef PTHREAD
	static void *updateBetaThread(void *work);

	// variables for threading
	pthread_t      *mThreads;
	bool           *mThreadStatus;
	pair_update_t **mThreadWork;

	// thread specific update vars
	double **mBeta_A_t;
	double **mBeta_b_t;

	pthread_mutex_t mPairMutex;
	int             mPair_t;
#endif
};

#endif
