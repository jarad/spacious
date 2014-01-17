// block composite fitting class
#ifndef BLOCK_H
#define BLOCK_H

#include "covs.h"

#ifdef PTHREAD
typedef struct {
	int              id;
	double           log_lik;
	class BlockComp *bc;
} pair_update_t;
#endif

#ifdef CUDA
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

class BlockComp {
public:
	BlockComp();
	BlockComp(int nthreads, bool gpu);
	~BlockComp();

	// types of likelihoods we can fit
	enum LikForm { Full, Pair, Block, IndBlock };
	enum CovType { Exp, Matern };

	// should we try to conserve memory during fit?
	void setConserveMemory(bool conserve) { mConsMem = conserve; }

	void init(int nthreads, bool gpu);
	void initPointers();
	void setLikForm(LikForm form);
	void setCovType(CovType type);
	bool setData(int n, double *y, double *S, int nblocks, int *B, int p, double *X, int npairs, int *neighbors);
	void setInits(double *theta);
	void setFixed(bool *fixed, double *values);
	void computeWithinDistance(int n, const double *S, double *D);
	void computeBetweenDistance(int n1, const double *S1, int n2, const double *S2, double *D);

	// fit model with Fisher scoring
	bool fit(bool verbose);

	// methods to compute quantities of interest (requires fit)
	bool predict(int n_0, double *y_0, const double *newS, const int *newB, const double *newX,
	             bool do_sd, double *sd, bool local=false, int Nlocal=25);

	// extract items of interest
	void getBeta(double *beta);
	void getTheta(double *theta);
	bool getConverged() { return(mConverged); }
	int  getIters()     { return(mIters); }

	void getFitted(double *fitted);
	void getResiduals(double *resids);

private:
	void cleanup();

	void setThreads(int nthreads);

	bool computeLogLik(double *log_lik);
	bool computeLogLikPair(double *log_lik, int pair, double *Sigma, double *resids, double *q);

	void computeFitted();
	void computeFitted(int n, double *fitted, double *X);

	void computeResiduals();
	void computeStdErrs();
	void computeCLIC();

	bool invertFullCov(bool do_log_det=false, double *log_det=NULL);
	bool updateBeta();
	bool updateBetaPair(int pair, double *Sigma, double *A, double *b);
	bool updateTheta();
	bool updateThetaPair(int pair, double *Sigma, double **W, double *H, double *P, double *resids, double *q, double *u);

	bool blockPredict(int block, int n_0, double *y_0, const double *newS, const double *newX, bool do_sd, double *sd);

	int     mNthreads;   // number of processing threads

	int     mN;          // number of observations
	double *mY;          // response

	double  *mS;          // spatial locations
	int      mNblocks;    // total number of blocks
	int     *mB;          // block membership
	int     *mNB;         // number of observations in each block
	int    **mWhichB;     // holds indices for obs in each block
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

	bool   mGPU;         // use GPU?
	bool   mConsMem;     // should memory be conserved?
	bool   mConsMemGPU;  // should we conserve GPU memory?
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
	double *mIterBeta;    // ... at each iteration
	double *mIterTheta;
	double *mIterLogLik;

	double *mFitted;      // fitted values
	double *mResids;      // residuals

	// update vars
	double **mSigma;
	double  *mBeta_A;
	double  *mBeta_b;

	double **mTheta_W;
	double  *mTheta_H;
	double  *mTheta_P;

#ifdef PTHREAD
	static void *updateBetaThread(void *work);
	static void *updateThetaThread(void *work);
	static void *computeLogLikThread(void *work);

	// variables for threading
	pthread_t      *mThreads;
	bool           *mThreadStatus;
	pair_update_t **mThreadWork;

	// thread specific update vars
	double **mBeta_A_t;
	double **mBeta_b_t;

	double ***mTheta_W_t;
	double  **mTheta_H_t;
	double  **mTheta_P_t;
	double  **mTheta_u_t;

	pthread_mutex_t mPairMutex;
	int             mPair_t;
#endif

#ifdef CUDA
	cublasHandle_t   mCublasHandle;
	double          *mDevSigma;
// other candidates: mDevX; mDevBeta; mDevResids; mDevq;
#endif

};

#endif
