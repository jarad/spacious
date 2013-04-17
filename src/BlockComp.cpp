// Estimates block composite models with Fisher scoring
#include <stdio.h>
#include <unistd.h>

#ifdef PTHREAD
#include <pthread.h>
#endif

#include <R.h>
#include <R_ext/Lapack.h>

#include "BlockComp.h"
#include "covs.h"
#include "utils.h"

// constructor
BlockComp::BlockComp() {
	init(1);
}

BlockComp::BlockComp(int nthreads) {
#ifdef PTHREAD
	init(nthreads);
#else
	init(1);
#endif
}

BlockComp::~BlockComp() {
	cleanup();

	// only clear these to end
	delete mCov;
	free(mThetaInits);
	free(mFixed);
	free(mFixedVals);
}

void BlockComp::init(int nthreads) {
	initPointers();

	setThreads(nthreads);

	// only initalize these pointers to start
	mCov        = NULL;
	mThetaInits = NULL;
	mFixed      = NULL;
	mFixedVals  = NULL;

	// default booleans
	mConsMem   = false;
	mHasFit    = false;
	mConverged = false;

	// default setup
	setLikForm(Block);
	setCovType(Exp);

	// default control params
	mIterTol = 1e-3;
	mMaxIter = 100;
}

// initialize pointers
void BlockComp::initPointers() {
	mNB         = NULL;
	mWhichB     = NULL;
	mWithinD    = NULL;
	mBetweenD   = NULL;

	mBeta       = NULL;
	mTheta      = NULL;
	mThetaT     = NULL;

	mSigma      = NULL;
	mBeta_A     = NULL;
	mBeta_b     = NULL;

	mTheta_W    = NULL;
	mTheta_H    = NULL;
	mTheta_P    = NULL;

#ifdef PTHREAD
	mThreads      = NULL;
	mThreadStatus = NULL;
	mThreadWork   = NULL;

	mBeta_A_t     = NULL;
	mBeta_b_t     = NULL;

	mTheta_W_t    = NULL;
	mTheta_H_t    = NULL;
	mTheta_P_t    = NULL;
	mTheta_u_t    = NULL;
#endif
}

// free's allocated data
void BlockComp::cleanup() {
	int i;

	free(mNB);

	if (mWhichB != NULL) {
		for (i = 0; i < mNblocks; i++) {
			free(mWhichB[i]);
		}
		free(mWhichB);
	}

	if (mWithinD != NULL) {
		if (mLikForm == Block) {
			for (i = 0; i < mNblocks; i++) {
				free(mWithinD[i]);
			}
		} else if (mLikForm == Full || mLikForm == Pair) {
			free(mWithinD[0]);
		}

		free(mWithinD);
	}

	if (mBetweenD != NULL) {
		for (i = 0; i < mNpairs; i++) {
			free(mBetweenD[i]);
		}

		free(mBetweenD);
	}

	free(mBeta);
	free(mTheta);
	free(mThetaT);

	if (mSigma != NULL) {
		for (i = 0; i < mNthreads; i++) {
			free(mSigma[i]);
		}

		free(mSigma);
	}

	free(mBeta_A);
	free(mBeta_b);

	if (mTheta_W != NULL) {
		for (i = 0; i < mNtheta; i++) {
			free(mTheta_W[i]);
		}
		free(mTheta_W);
	}
	free(mTheta_H);
	free(mTheta_P);

#ifdef PTHREAD
	free(mThreads);
	free(mThreadStatus);

	if (mNthreads > 1) {
		if (mThreadWork != NULL) {
			for (i = 0; i < mNthreads; i++) {
				free(mThreadWork[i]);
			}
		}
		free(mThreadWork);

		if (mBeta_A_t != NULL && mBeta_b_t != NULL) {
			for (i = 0; i < mNthreads; i++) {
				free(mBeta_A_t[i]);
				free(mBeta_b_t[i]);
			}
		}

		free(mBeta_A_t);
		free(mBeta_b_t);

		if (mTheta_W_t != NULL && mTheta_H_t != NULL && mTheta_P_t && mTheta_u_t) {
			for (i = 0; i < mNthreads; i++) {
				for (int j = 0; j < mNtheta; j++) {
					free(mTheta_W_t[i][j]);
				}
				free(mTheta_W_t[i]);
				free(mTheta_H_t[i]);
				free(mTheta_P_t[i]);
				free(mTheta_u_t[i]);
			}
		}

		free(mTheta_W_t);
		free(mTheta_H_t);
		free(mTheta_P_t);
		free(mTheta_u_t);
	}
#endif

	initPointers();
}

void BlockComp::setThreads(int nthreads) {
	// number of processing threads
	mNthreads = nthreads;
}

void BlockComp::setLikForm(LikForm form) {
	// some cases will require a new setData() call,
	// so always cleanup and require setData() just in case
	cleanup();

	mLikForm = form;

	if (mLikForm == Full) {
		// make sure we only use one thread
		setThreads(1);
	}
}

void BlockComp::setCovType(CovType type) {
	mCovType = type;

	delete mCov;

	// obtain class for working with covariance type
	if (mCovType == Exp) {
		mCov = new CovExp();
	} else {
		MSG("Unknown covariance type\n");
	}

	// initial values need to be re-specified
	free(mThetaInits);
	mThetaInits = NULL;

	// default to no fixed parameters
	bool   fixed[mCov->numParams()];
	double values[mCov->numParams()];

	for (int i = 0; i < mCov->numParams(); i++) {
		fixed[i]  = false;
		values[i] = 0;
	}

	setFixed(fixed, values);
}

// specify data to fit model to
void BlockComp::setData(int n, double *y, double *S, int nblocks, int *B, int p, double *X, int npairs, int *neighbors) {
	int i,j,k;
	int blk1, blk2;
	int nIn1, nIn2;

#ifdef DEBUG
	printf("n=%d, y[0]=%.2f, nblocks=%d, B[0]=%d, p=%d, X[0]=%.2f, npairs=%d, neighbors[0]=%d, mHasFit=%d\n",
	       n, y[0], nblocks, B[0], p, X[0], npairs, neighbors[0], mHasFit);
#endif

	// require new fit
	mHasFit  = false;

	cleanup();

	// load the new data
	mN = n;
	mY = y;

	mS       = S;
	mNblocks = nblocks;
	mB       = B;

	mNbeta = p;
	mX     = X;

	mNpairs    = npairs;
	mNeighbors = neighbors;

	if (mLikForm == Block) {
		// compute number of observations in each block
		mNB = (int *)malloc(sizeof(int)*mNblocks);
		for (i = 0; i < mNblocks; i++) { mNB[i] = 0; }

		for (i = 0; i < mN; i++) {
			mNB[ mB[i] ]++;
		}

		// determine the largest number of obs in each pair
		mMaxPair = 0;
		for (i = 0; i < mNpairs; i++) {
			j = mNB[mNeighbors[i]] + mNB[mNeighbors[i+mNpairs]];
			if (j > mMaxPair) {
				mMaxPair = j;
			}
		}

		// save indicies for each block for easy access
		mWhichB = (int **)malloc(sizeof(int *)*mNblocks);
		for (i = 0; i < mNblocks; i++) {
			mWhichB[i] = (int *)malloc(sizeof(int)*mNB[i]);
		}

		// temporarily keep track of locations for mWhichB
		int locs[mNblocks];
		for (i = 0; i < mNblocks; i++) { locs[i] = 0; }

		// fill in the indicies
		for (i = 0; i < mN; i++) {
			mWhichB[ mB[i] ][ locs[ mB[i] ]++ ] = i;
		}
	} else if (mLikForm == Full) {
		mMaxPair = mN;
	} else if (mLikForm == Pair) {
		mMaxPair = 2;
	}

	// allocate largest Sigma we need
	mSigma = (double **)malloc(sizeof(double *)*mNthreads);
	for (i = 0; i < mNthreads; i++) {
		mSigma[i] = (double *)malloc(sizeof(double)*symi(0,mMaxPair));

		for (j = 0; j < symi(0,mMaxPair); j++) { mSigma[i][j] = 0; }
	}

#ifdef DEBUG
	MSG("max pair=%d with length %d\n", mMaxPair, symi(0, mMaxPair));
#endif

	if (!mConsMem) {
		// not trying to conserve memory, so store distances...

		if (mLikForm == Block) {
			// ... within blocks
			mWithinD = (double **)malloc(sizeof(double *)*mNblocks);
			for (i = 0; i < mNblocks; i++) {
				mWithinD[i] = (double *)malloc(sizeof(double)*symi(0, mNB[i]));
			}

			for (i = 0; i < mNblocks; i++) {
				for (j = 0; j < mNB[i]; j++) {
					for (k = j; k < mNB[i]; k++) {
						if (j == k) {
							mWithinD[i][ symi(j,j) ] = 0;
						} else {
							mWithinD[i][ symi(j,k) ] = sqrt( pow(mS[mWhichB[i][j]]-mS[mWhichB[i][k]], 2) + pow(mS[mWhichB[i][j]+mN]-mS[mWhichB[i][k]+mN], 2) );
						}
					}
				}
			}

			// ... between blocks
			mBetweenD = (double **)malloc(sizeof(double *)*mNpairs);
			for (i = 0; i < mNpairs; i++) {
				mBetweenD[i] = (double *)malloc(sizeof(double)* mNB[mNeighbors[i]]*mNB[mNeighbors[i+mNpairs]]);
			}

			for (i = 0; i < mNpairs; i++) {
				blk1 = mNeighbors[i];
				blk2 = mNeighbors[i+mNpairs];

				nIn1 = mNB[mNeighbors[i]];
				nIn2 = mNB[mNeighbors[i+mNpairs]];

				for (j = 0; j < nIn1; j++) {
					for (k = 0; k < nIn2; k++) {
						mBetweenD[i][j + k*nIn1] = sqrt(
							pow(mS[mWhichB[blk1][j]]-mS[mWhichB[blk2][k]], 2) + pow(mS[mWhichB[blk1][j]+mN]-mS[mWhichB[blk2][k]+mN], 2)
						);
					}
				}
			}

		} else if (mLikForm == Full || mLikForm == Pair) {
			// we need distances between all points for Full and Pair
			mWithinD    = (double **)malloc(sizeof(double *));
			mWithinD[0] = (double *)malloc(sizeof(double)*symi(0, mN));

			for (i = 0; i < mN; i++) {
				for (j = i; j < mN; j++) {
					if (i == j) {
						mWithinD[0][ symi(j,j) ] = 0;
					} else {
						mWithinD[0][ symi(i,j) ] = sqrt( pow(mS[i]-mS[j], 2) + pow(mS[i+mN] - mS[j+mN], 2) );
					}
				}
			}

		}

	}  // end memory conservation check
}

void BlockComp::setInits(double *theta) {
	free(mThetaInits);

	mThetaInits = (double *)malloc(sizeof(double)*mCov->numParams());
	for (int i = 0; i < mCov->numParams(); i++) {
		mThetaInits[i] = theta[i];
	}
}

// which covariance parameters are fixed?
void BlockComp::setFixed(bool *fixed, double *values) {
	// require a (new) fit
	mHasFit = false;

	free(mFixed);
	free(mFixedVals);

	mFixed     = (bool *)malloc(sizeof(bool)*mCov->numParams());
	mFixedVals = (double *)malloc(sizeof(double)*mCov->numParams());

	mNfixed = 0;
	for (int i = 0; i < mCov->numParams(); i++) {
		mFixed[i]     = fixed[i];
		mFixedVals[i] = values[i];

		if (mFixed[i]) {
			mNfixed++;
		}
	}
}

// fit model to data
bool BlockComp::fit(bool verbose) {
	int i;
	bool largeDiff;

#ifdef DEBUG
	MSG("Starting fit()...\n");
#endif

	// make sure we have initial values
	if (mThetaInits == NULL) {
		if (verbose) { MSG("Missing initial values for covariance parameters.\n"); }
		return(false);
	}

	// number of covariance parameters
	mNtheta = mCov->numParams();

	// initialize covariance params
	free(mTheta);
	free(mThetaT);
	mTheta  = (double *)malloc(sizeof(double) * mNtheta);
	mThetaT = (double *)malloc(sizeof(double) * mNtheta);

	// set initial theta and transform
	for (i = 0; i < mNtheta; i++) {
		if (mFixed[i]) {
			mTheta[i]  = mFixedVals[i];
		} else {
			mTheta[i]  = mThetaInits[i];
		}
		mThetaT[i] = mTheta[i];
	}
	mCov->transformToReal(mThetaT);

	// prepare variables for updating beta
	free(mBeta);
	free(mBeta_A);
	free(mBeta_b);
	mBeta   = (double *)malloc(sizeof(double)*mNbeta);
	mBeta_A = (double *)malloc(sizeof(double)*symi(0,mNbeta));
	mBeta_b = (double *)malloc(sizeof(double)*mNbeta);

	// prepare variables for updating theta
	if (mTheta_W != NULL) {
		for (i = 0; i < mNtheta; i++) {
			free(mTheta_W[i]);
		}
		free(mTheta_W);
	}
	free(mTheta_H);
	free(mTheta_P);
	mTheta_W = (double **)malloc(sizeof(double *)*mNtheta);
	mTheta_H = (double *)malloc(sizeof(double)*symi(0,mNtheta));
	mTheta_P = (double *)malloc(sizeof(double)*symi(0,mMaxPair));
	for (i = 0; i < mNtheta; i++) {
		mTheta_W[i] = (double *)malloc(sizeof(double)*pow(mMaxPair,2));
	}

#ifdef PTHREAD
	free(mThreads);
	free(mThreadStatus);

	if (mNthreads > 1) {
		MSG("TODO: order block pairs for threading!\n");

		// allocate space for threads
		mThreads      = (pthread_t *)malloc(sizeof(pthread_t)*mNthreads);
		mThreadStatus = (bool *)malloc(sizeof(bool)*mNthreads);

		// allocate thread work
		if (mThreadWork != NULL) {
			for (i = 0; i < mNthreads; i++) {
				free(mThreadWork[i]);
			}
		}
		free(mThreadWork);

		mThreadWork = (pair_update_t **)malloc(sizeof(pair_update_t *)*mNthreads);
		for (i = 0; i < mNthreads; i++) {
			mThreadWork[i] = (pair_update_t *)malloc(sizeof(pair_update_t)*mNthreads);
		}

		// allocate update variables specific to each thread
		if (mBeta_A_t != NULL && mBeta_b_t != NULL) {
			for (i = 0; i < mNthreads; i++) {
				free(mBeta_A_t[i]);
				free(mBeta_b_t[i]);
			}
		}
		free(mBeta_A_t);
		free(mBeta_b_t);

		mBeta_A_t = (double **)malloc(sizeof(double *)*mNthreads);
		mBeta_b_t = (double **)malloc(sizeof(double *)*mNthreads);
		for (i = 0; i < mNthreads; i++) {
			mBeta_A_t[i] = (double *)malloc(sizeof(double)*symi(0,mNbeta));
			mBeta_b_t[i] = (double *)malloc(sizeof(double)*mNbeta);
		}

		if (mTheta_W_t != NULL && mTheta_H_t != NULL && mTheta_P_t && mTheta_u_t) {
			for (i = 0; i < mNthreads; i++) {
				for (int j = 0; j < mNtheta; j++) {
					free(mTheta_W_t[i][j]);
				}
				free(mTheta_W_t[i]);
				free(mTheta_H_t[i]);
				free(mTheta_P_t[i]);
				free(mTheta_u_t[i]);
			}
		}

		free(mTheta_W_t);
		free(mTheta_H_t);
		free(mTheta_P_t);
		free(mTheta_u_t);

		mTheta_W_t = (double ***)malloc(sizeof(double **)*mNthreads);
		mTheta_H_t = (double **)malloc(sizeof(double *)*mNthreads);
		mTheta_P_t = (double **)malloc(sizeof(double *)*mNthreads);
		mTheta_u_t = (double **)malloc(sizeof(double *)*mNthreads);
		for (i = 0; i < mNthreads; i++) {
			mTheta_W_t[i] = (double **)malloc(sizeof(double)*mNtheta);
			for (int j = 0; j < mNtheta; j++) {
				mTheta_W_t[i][j] = (double *)malloc(sizeof(double)*pow(mMaxPair,2));
			}

			mTheta_H_t[i] = (double *)malloc(sizeof(double)*symi(0,mNtheta));
			mTheta_P_t[i] = (double *)malloc(sizeof(double)*symi(0,mMaxPair));
			mTheta_u_t[i] = (double *)malloc(sizeof(double)*mNtheta);
		}

	}

#endif

	// get initial beta
	if (!updateBeta()) {
		MSG("Unable to get initial values for beta\n");
		return(false);
	}

	// vectors to hold previous parameter values
	double *prevBeta   = (double *)malloc(sizeof(double)*mNbeta);
	double *prevThetaT = (double *)malloc(sizeof(double)*mNtheta);

	for (mIters = 0; mIters < mMaxIter; mIters++) {
		// store previous values
		for (i = 0; i < mNbeta; i++)  { prevBeta[i]   = mBeta[i]; }
		for (i = 0; i < mNtheta; i++) { prevThetaT[i] = mThetaT[i]; }

		// update covariance params
		if (!updateTheta()) {
			MSG("Unable to update theta at iteration %d\n", mIters+1);
			free(prevBeta); free(prevThetaT);
			return(false);
		}

		// update mean params
		if (!updateBeta()) {
			MSG("Unable to update beta at iteration %d\n", mIters+1);
			free(prevBeta); free(prevThetaT);
			return(false);
		}

		if (verbose) {
			// display iteration information
			MSG("iter=%d: ", mIters+1);
			MSG("beta: ");
			for (i = 0; i < mNbeta; i++) { MSG("%.2f ", mBeta[i]); }
			MSG("; theta: ");
			for (i = 0; i < mNtheta; i++) { MSG("%.2f ", mTheta[i]); }
			MSG("\n");
		}

		// detect convergence
		largeDiff = false;

		// check betas for convergence
		for (i = 0; i < mNbeta; i++) {
			if (fabs(mBeta[i]-prevBeta[i])/fabs(prevBeta[i]) > mIterTol) {
				largeDiff = true;
				break;
			}
		}

		if (!largeDiff) {
			// check thetas for convergence
			for (i = 0; i < mNtheta; i++) {
				if (fabs(mThetaT[i]-prevThetaT[i])/fabs(prevThetaT[i]) > mIterTol) {
					largeDiff = true;
					break;
				}
			}
		}

		if (!largeDiff) {
			// convergence!
			if (verbose) {
				MSG("Converged at iteration %d\n", mIters+1);
			}

			mConverged = true;

			break;
		}

	}

	// clean up
	free(prevBeta);
	free(prevThetaT);

	// we have a fit!
	mHasFit = true;

	return(true);
}

// update mean parameters
bool BlockComp::updateBeta() {
	int i,j,k,l;
	int selem,icol,jcol;
	int pair;

	// initialize A and b
	for (i = 0; i < symi(0,mNbeta); i++) { mBeta_A[i] = 0; }
	for (i = 0; i < mNbeta;         i++) { mBeta_b[i] = 0; }

	if (mLikForm == Block) {
		// update beta with block composite likelihood

#ifdef PTHREAD
		if (mNthreads <= 1) {
#endif
			// process each block pair in order
			for (pair = 0; pair < mNpairs; pair++) {
				if (!updateBetaPair(pair, mSigma[0], mBeta_A, mBeta_b)) {
					// error updating this pair
					return(false);
				}
			}
#ifdef PTHREAD
		} else {
			// use threads to process each block pair

			// setup mutex
			pthread_mutex_init(&mPairMutex, NULL);
			mPair_t    = 0;

			// create threads
			for (i = 0; i < mNthreads; i++) {
				mThreadStatus[i]      = true;
				mThreadWork[i]->id    = i;
				mThreadWork[i]->bc    = this;

				pthread_create(&mThreads[i], NULL,
					&BlockComp::updateBetaThread,
					(void *)mThreadWork[i]
				);

			}

			// wait for all threads to complete
			for (i = 0; i < mNthreads; i++) {
				pthread_join(mThreads[i], 0);
			}

			// destroy mutex
			pthread_mutex_destroy(&mPairMutex);

			// did we have any errors?
			for (i = 0; i < mNthreads; i++) {
				if (!mThreadStatus[i]) {
					MSG("updateBeta(): Error processing thread %d.\n", i);
					return(false);
				}
			}

			// combine results from each thread
			for (j = 0; j < mNthreads; j++) {
				for (i = 0; i < symi(0,mNbeta); i++) { mBeta_A[i] += mBeta_A_t[j][i]; }
				for (i = 0; i < mNbeta;         i++) { mBeta_b[i] += mBeta_b_t[j][i]; }
			}
		}
#endif

	} else if (mLikForm == IndBlock) {
		// update beta using independent blocks
		MSG("TODO: updateBeta(): mLikForm=IndBlock\n");
		return(false);
	} else if (mLikForm == Pair) {
		// update beta using pairwise composite likelihood
		MSG("TODO: updateBeta(): mLikForm=Pair\n");
		return(false);
	} else if (mLikForm == Full) {
		// update beta using full likelihood

		if (!mConsMem) {
			// fill in covariance matrix between these two blocks
			mCov->compute(mSigma[0], mTheta, mN, mWithinD[0]);
		} else {
			// we're conserving memory

			// fill in distance
			// fill in covariance
			MSG("TODO: updateBeta(): mConsMem=true, mLikForm=Full\n");
			return(false);
		}

		// invert Sigma
		if (chol2inv(mN, mSigma[0])) {
			MSG("updateBeta(): unable to invert Sigma\n");
			return(false);
		}

		// add contribution to A and b
		for (i = 0; i < mNbeta; i++) {
			icol = i*mN;

			for (j = i; j < mNbeta; j++) {
				selem = symi(i,j);
				jcol  = j*mN;

				for (k = 0; k < mN; k++) {
					for (l = 0; l < mN; l++) {
						mBeta_A[selem] += mX[l + icol] * mSigma[0][symi(l,k)] * mX[k + jcol];
						if (i == j) {
							mBeta_b[i] += mX[l + icol] * mSigma[0][symi(l,k)] * mY[k];
						}
					}
				}

			}
		}

	} // end different likelihood forms

	// compute beta
	if (chol2inv(mNbeta, mBeta_A)) {
		MSG("updateBeta(): unable to invert A\n");
		return(false);
	}

	for (i = 0; i < mNbeta; i++) {
		mBeta[i] = 0;

		for (j = 0; j < mNbeta; j++) {
			mBeta[i] += mBeta_A[symi(i,j)] * mBeta_b[j];
		}
	}

	return(true);
}

bool BlockComp::updateBetaPair(int pair, double *Sigma, double *A, double *b) {
	// update specified pair of blocks
	int i,j,k,l;
	int selem,icol,jcol;

	int blk1 = mNeighbors[pair];
	int blk2 = mNeighbors[pair+mNpairs];
	int N_in_pair = mNB[blk1] + mNB[blk2];

	if (!mConsMem) {
		// fill in covariance matrix between these two blocks
		mCov->compute(Sigma, mTheta, mNB[blk1], mWithinD[blk1], mNB[blk2], mWithinD[blk2], mBetweenD[pair]);
	} else {
		// we're conserving memory

		// fill in distance
		// fill in covariance
		MSG("TODO: updateBetaPair(): mConsMem=true, mLikForm=Block\n");
		return(false);
	}

	// invert Sigma
	if (chol2inv(N_in_pair, Sigma)) {
		MSG("updateBetaPair(): unable to invert Sigma\n");
		return(false);
	}

	// add contribution to A and b
	for (i = 0; i < mNbeta; i++) {
		icol = i*mN;

		for (j = i; j < mNbeta; j++) {
			selem = symi(i,j);
			jcol  = j*mN;

			for (k = 0; k < mNB[blk1]; k++) {     // block 1...
				for (l = 0; l < mNB[blk1]; l++) {   // with block 1
					A[selem] += mX[mWhichB[blk1][l] + icol] * Sigma[symi(l,k)] * mX[mWhichB[blk1][k] + jcol];
					if (i == j) {
						b[i] += mX[mWhichB[blk1][l] + icol] * Sigma[symi(l,k)] * mY[mWhichB[blk1][k]];
					}
				}

				for (l = 0; l < mNB[blk2]; l++) {   // with block 2
					A[selem] += mX[mWhichB[blk2][l] + icol] * Sigma[symi(l+mNB[blk1],k)] * mX[mWhichB[blk1][k] + jcol];
					if (i == j) {
						b[i] += mX[mWhichB[blk2][l] + icol] * Sigma[symi(l+mNB[blk1],k)] * mY[mWhichB[blk1][k]];
					}
				}
			}

			for (k = 0; k < mNB[blk2]; k++) {     // block 2...
				for (l = 0; l < mNB[blk1]; l++) {   // with block 1
					A[selem] += mX[mWhichB[blk1][l] + icol] * Sigma[symi(l,k+mNB[blk1])] * mX[mWhichB[blk2][k] + jcol];
					if (i == j) {
						b[i] += mX[mWhichB[blk1][l] + icol] * Sigma[symi(l,k+mNB[blk1])] * mY[mWhichB[blk2][k]];
					}
				}

				for (l = 0; l < mNB[blk2]; l++) {   // with block 2
					A[selem] += mX[mWhichB[blk2][l] + icol] * Sigma[symi(l+mNB[blk1],k+mNB[blk1])] * mX[mWhichB[blk2][k] + jcol];
					if (i == j) {
						b[i] += mX[mWhichB[blk2][l] + icol] * Sigma[symi(l+mNB[blk1],k+mNB[blk1])] * mY[mWhichB[blk2][k]];
					}
				}
			}

		}
	}

	return(true);
}

#ifdef PTHREAD
void *BlockComp::updateBetaThread(void *work) {
	pair_update_t *w = (pair_update_t *)work;
	int            id = w->id;
	BlockComp     *bc = w->bc;

	int i;
	int pair;

	// initialize A and b for this thread
	for (i = 0; i < symi(0,bc->mNbeta); i++) { bc->mBeta_A_t[id][i] = 0; }
	for (i = 0; i < bc->mNbeta;         i++) { bc->mBeta_b_t[id][i] = 0; }

	// process blocks
	while (1) {
		// get pair to process
		pthread_mutex_lock(&(bc->mPairMutex));
		pair = bc->mPair_t++;
		pthread_mutex_unlock(&(bc->mPairMutex));

		if (pair >= bc->mNpairs) {
			// no more to process
			break;
		}

		if (!bc->updateBetaPair(pair, bc->mSigma[id], bc->mBeta_A_t[id], bc->mBeta_b_t[id])) {
			// error updating this pair
			bc->mThreadStatus[id] = false;
			return(NULL);
		}
	}

	bc->mThreadStatus[id] = true;
	return(NULL);
}
#endif

bool BlockComp::updateTheta() {
	int i,j,k,l,c;
	int iH,jH;       // used to fill hessian properly when params fixed
	int pair;

	double resids[mMaxPair];
	double q[mMaxPair];
	double u[mNtheta];
	bool   diag;

	if (mNfixed == mNtheta) {
		// all are fixed -- no updates needed
		return(true);
	}

	// transform theta to real line
	for (i = 0; i < mNtheta; i++) { mThetaT[i] =  mTheta[i]; }
	mCov->transformToReal(mThetaT);

	// initialize u and H
	for (i = 0; i < mNtheta;         i++) { u[i] = 0; }
	for (i = 0; i < symi(0,mNtheta); i++) { mTheta_H[i] = 0; }

	if (mLikForm == Block) {
		// update theta with block composite likelihood

#ifdef PTHREAD
		if (mNthreads <= 1) {
#endif
			for (pair = 0; pair < mNpairs; pair++) {
				if (!updateBetaPair(pair, mSigma[0], mBeta_A, mBeta_b)) {
					// error updating this pair
					return(false);
				}
			}

			// process each block pair in order
			for (pair = 0; pair < mNpairs; pair++) {
				if (!updateThetaPair(pair, mSigma[0], mTheta_W, mTheta_H, mTheta_P, resids, q, u)) {
						// error updating this pair
						return(false);
				}
			}

#ifdef PTHREAD
		} else {
			// use threads to process each block pair

			// setup mutex
			pthread_mutex_init(&mPairMutex, NULL);
			mPair_t    = 0;

			// create threads
			for (i = 0; i < mNthreads; i++) {
				mThreadStatus[i]      = true;
				mThreadWork[i]->id    = i;
				mThreadWork[i]->bc    = this;

				pthread_create(&mThreads[i], NULL,
					&BlockComp::updateThetaThread,
					(void *)mThreadWork[i]
				);

			}

			// wait for all threads to complete
			for (i = 0; i < mNthreads; i++) {
				pthread_join(mThreads[i], 0);
			}

			// destroy mutex
			pthread_mutex_destroy(&mPairMutex);

			// did we have any errors?
			for (i = 0; i < mNthreads; i++) {
				if (!mThreadStatus[i]) {
					MSG("updateTheta(): Error processing thread %d.\n", i);
					return(false);
				}
			}

			// combine results from each thread
			for (j = 0; j < mNthreads; j++) {
				for (i = 0; i < symi(0, mNtheta-mNfixed); i++) { mTheta_H[i] += mTheta_H_t[j][i]; }
				for (i = 0; i < mNtheta;                  i++) { u[i]        += mTheta_u_t[j][i]; }
			}
		}
#endif

	} else if (mLikForm == IndBlock) {
		// update theta using independent blocks
		MSG("TODO: updateTheta(): mLikForm=IndBlock\n");
		return(false);
	} else if (mLikForm == Pair) {
		// update theta using pairwise composite likelihood
		MSG("TODO: updateTheta(): mLikForm=Pair\n");
		return(false);
	} else if (mLikForm == Full) {
		// update theta using full likelihood

		// note that updateBeta() already computes inv(Sigma) used here

		// initialize residuals and q
		for (i = 0; i < mN; i++) {
			resids[i] = 0;
			q[i]      = 0;
		}

		// compute resids = y - X'b
		for (i = 0; i < mN; i++) {
			for (j = 0; j < mNbeta; j++) {
				resids[i] += mX[i + j*mN] * mBeta[j];
			}

			resids[i] = mY[i] - resids[i];
		}

		// compute q = inv(Sigma) x resids
		for (i = 0; i < mN; i++) {
			for (j = 0; j < mN; j++) {
				q[i] += mSigma[0][symi(i,j)] * resids[j];
			}
		}

		// fill in W and u
		for (i = 0; i < mNtheta; i++) {
			if (mFixed[i]) {
				// this parameter is fixed, so skip these operations
				continue;
			}

			// get partial derivatives
			mCov->partials(mTheta_P, &diag, i, mTheta, mThetaT, mN, mWithinD[0]);

			// initialize W[i]
			for (j = 0; j < pow(mN, 2); j++) { mTheta_W[i][j] = 0; }

			// compute inv(Sigma) x P
			if (diag) {
				// take advantage of P being diagonal
				for (j = 0; j < mN; j++) {
					for (k = 0; k < mN; k++) {
						mTheta_W[i][j + k*mN] = mSigma[0][symi(j,k)] * mTheta_P[symi(k,k)];
					}
				}
			} else {
				for (j = 0; j < mN; j++) {
					for (k = 0; k < mN; k++) {
						for (l = 0; l < mN; l++) {
							mTheta_W[i][j + k*mN] += mSigma[0][symi(j,l)] * mTheta_P[symi(l,k)];
						}
					}
				}
			}

			// fill in u
			for (j = 0; j < mN; j++) {
				u[i] -= 0.5*mTheta_W[i][j+j*mN];
			}

			for (j = 0; j < mN; j++) {
				u[i] += 0.5*mTheta_P[symi(j,j)] * q[j] * q[j];
				for (k = j+1; k < mN; k++) {
					u[i] += mTheta_P[symi(j,k)] * q[j] * q[k];
				}
			}

		}

		// compute hessian
		iH = 0;
		for (i = 0; i < mNtheta; i++) {
			if (mFixed[i]) { continue; }

			jH = iH;
			for (j = i; j < mNtheta; j++) {
				if (mFixed[j]) { continue; }

				c = symi(iH,jH);

				// add in diagonal elements of W[i] x W[j]
				for (k = 0; k < mN; k++) {
					for (l = 0; l < mN; l++) {
						mTheta_H[c] += mTheta_W[i][k+l*mN] * mTheta_W[j][l+k*mN];
					}
				}

				jH++;
			}

			iH++;
		}

	}

	// hessian elements should be scaled by 1/2
	for (i = 0; i < symi(0, mNtheta-mNfixed); i++) {
		mTheta_H[i] *= 0.5;
	}

	// invert hessian
	if (chol2inv(mNtheta-mNfixed, mTheta_H)) {
		MSG("updateTheta(): Unable to invert H\n");
		return(false);
	}

	// update theta
	iH = 0;
	for (i = 0; i < mNtheta; i++) {
		if (mFixed[i]) { continue; }

		jH = 0;
		for (j = 0; j < mNtheta; j++) {
			if (mFixed[j]) { continue; }

			mThetaT[i] += mTheta_H[symi(iH,jH)] * u[j];

			jH++;
		}

		iH++;
	}

	// transform thetaT to original scale
	for (i = 0; i < mNtheta; i++) { mTheta[i] =  mThetaT[i]; }
	mCov->transformFromReal(mTheta);

	return(true);
}

bool BlockComp::updateThetaPair(int pair, double *Sigma, double **W, double *H, double *P, double *resids, double *q, double *u) {
	int blk1 = mNeighbors[pair];
	int blk2 = mNeighbors[pair+mNpairs];
	int N_in_pair = mNB[blk1] + mNB[blk2];

	int i,j,k,l;
	int iH,jH;       // used to fill hessian properly when params fixed
	int c;
	bool diag;

	if (!mConsMem) {
		// fill in covariance matrix between these two blocks
		mCov->compute(Sigma, mTheta, mNB[blk1], mWithinD[blk1], mNB[blk2], mWithinD[blk2], mBetweenD[pair]);
	} else {
		// we're conserving memory

		// fill in distance
		// fill in covariance
		MSG("TODO: updateThetaPair(): mConsMem=true, mLikForm=Block\n");
		return(false);
	}

	// invert Sigma
	if (chol2inv(N_in_pair, Sigma)) {
		MSG("updateThetaPair(): Unable to invert Sigma\n");
		return(false);
	}

	// initialize residuals and q
	for (i = 0; i < N_in_pair; i++) {
		resids[i] = 0;
		q[i]      = 0;
	}

	// compute resids = y - X'b
	for (i = 0; i < mNB[blk1]; i++) {   // block 1
		c = mWhichB[blk1][i];

		for (j = 0; j < mNbeta; j++) {
			resids[i] += mX[c + j*mN] * mBeta[j];
		}

		resids[i] = mY[c] - resids[i];
	}

	for (i = 0; i < mNB[blk2]; i++) {   // block 2
		c = mWhichB[blk2][i];

		for (j = 0; j < mNbeta; j++) {
			resids[i+mNB[blk1]] += mX[c + j*mN] * mBeta[j];
		}

		resids[i+mNB[blk1]] = mY[c] - resids[i+mNB[blk1]];
	}

	// compute q = inv(Sigma) x resids
	for (i = 0; i < N_in_pair; i++) {
		for (j = 0; j < N_in_pair; j++) {
			q[i] += Sigma[symi(i,j)] * resids[j];
		}
	}

	// fill in W and u
	for (i = 0; i < mNtheta; i++) {
		if (mFixed[i]) {
			// this parameter is fixed, so skip these operations
			continue;
		}

		// get partial derivatives
		mCov->partials(P, &diag, i, mTheta, mThetaT, mNB[blk1], mWithinD[blk1], mNB[blk2], mWithinD[blk2], mBetweenD[pair]);

		// initialize W[i]
		for (j = 0; j < pow(N_in_pair, 2); j++) { W[i][j] = 0; }

		// compute inv(Sigma) x P
		if (diag) {
			// take advantage of P being diagonal
			for (j = 0; j < N_in_pair; j++) {
				for (k = 0; k < N_in_pair; k++) {
					W[i][j + k*N_in_pair] = Sigma[symi(j,k)] * P[symi(k,k)];
				}
			}
		} else {
			for (j = 0; j < N_in_pair; j++) {
				for (k = 0; k < N_in_pair; k++) {
					for (l = 0; l < N_in_pair; l++) {
						W[i][j + k*N_in_pair] += Sigma[symi(j,l)] * P[symi(l,k)];
					}
				}
			}
		}

		// fill in u
		for (j = 0; j < N_in_pair; j++) {
			u[i] -= 0.5*W[i][j+j*N_in_pair];
		}

		for (j = 0; j < N_in_pair; j++) {
			u[i] += 0.5*P[symi(j,j)] * q[j] * q[j];
			for (k = j+1; k < N_in_pair; k++) {
				u[i] += P[symi(j,k)] * q[j] * q[k];
			}
		}

	}

	// compute hessian contributions
	iH = 0;
	for (i = 0; i < mNtheta; i++) {
		if (mFixed[i]) { continue; }

		jH = iH;
		for (j = i; j < mNtheta; j++) {
			if (mFixed[j]) { continue; }

			c = symi(iH,jH);

			// add in diagonal elements of W[i] x W[j]
			for (k = 0; k < N_in_pair; k++) {
				for (l = 0; l < N_in_pair; l++) {
					H[c] += W[i][k+l*N_in_pair] * W[j][l+k*N_in_pair];
				}
			}

			jH++;
		}

		iH++;
	}

	return(true);
}

#ifdef PTHREAD
void *BlockComp::updateThetaThread(void *work) {
	pair_update_t *w = (pair_update_t *)work;
	int            id = w->id;
	BlockComp     *bc = w->bc;

	int i;
	int pair;
	double resids[bc->mMaxPair];
	double q[bc->mMaxPair];

	// initialize u and H for this thread
	for (i = 0; i < bc->mNtheta;         i++) { bc->mTheta_u_t[id][i] = 0; }
	for (i = 0; i < symi(0,bc->mNtheta); i++) { bc->mTheta_H_t[id][i] = 0; }

	// process blocks
	while (1) {
		// get pair to process
		pthread_mutex_lock(&(bc->mPairMutex));
		pair = bc->mPair_t++;
		pthread_mutex_unlock(&(bc->mPairMutex));

		if (pair >= bc->mNpairs) {
			// no more to process
			break;
		}

		if (!bc->updateThetaPair(pair, bc->mSigma[id], bc->mTheta_W_t[id], bc->mTheta_H_t[id],
		                         bc->mTheta_P_t[id], resids, q, bc->mTheta_u_t[id])) {
			// error updating this pair
			bc->mThreadStatus[id] = false;
			return(NULL);
		}
	}

	bc->mThreadStatus[id] = true;
	return(NULL);
}
#endif

void BlockComp::getBeta(double *beta) {
	for (int i = 0; i < mNbeta; i++) {
		beta[i] = mBeta[i];
	}
}

void BlockComp::getTheta(double *theta) {
	for (int i = 0; i < mNtheta; i++) {
		theta[i] = mTheta[i];
	}
}

#ifndef BLK

// test the estimation on the data in test_data.h

#include "test_data.h"

int main(void) {
	BlockComp blk(4);

	// start blocks at 0
	int i;
	for (i = 0; i < test_n; i++) {
		test_B[i]--;
	}

	for (i = 0; i < test_npairs*2; i++) {
		test_neighbors[i]--;
	}

	blk.setLikForm(BlockComp::Block);
	//blk.setLikForm(BlockComp::Full);
	blk.setCovType(BlockComp::Exp);

	blk.setData(test_n, test_y, test_S, test_nblocks, test_B, test_p, test_X, test_npairs, test_neighbors);
	//blk.setData(test_n, test_y, test_S, test_nblocks, test_B, test_p, test_X, test_npairs, test_neighbors);

	//double inits[] = {0.25, 0.25, 0.5};
	//double inits[] = {0.5, 0.7243848, 0.2105263};
	bool fixed[] = {false, false, false};
	double vals[] = {0.5, 0.4, 0.15};

	blk.setInits(test_inits);
	blk.setFixed(fixed, vals);
	if (!blk.fit(true)) {
		MSG("Error with fit.\n");
	}

	return(0);
}

#endif
