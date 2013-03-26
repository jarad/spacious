// Estimates block composite models with Fisher scoring
#include <stdio.h>

#include <R.h>
#include <R_ext/Lapack.h>

#include "BlockComp.h"
#include "covs.h"
#include "utils.h"

// constructor
BlockComp::BlockComp() {
	mConsMem = false;
	mHasFit  = false;

	// default setup
	mLikForm = Block;
	mCovType = Exp;

	// default control params
	mIterTol = 1e-3;
	mMaxIter = 100;

	initPointers();
}

BlockComp::~BlockComp() {
	cleanup();
}

// initialize pointers
void BlockComp::initPointers() {
	mNB         = NULL;
	mWhichB     = NULL;
	mWithinD    = NULL;
	mBetweenD   = NULL;
	mThetaInits = NULL;

	mBeta       = NULL;
	mTheta      = NULL;
	mBetaT      = NULL;
	mThetaT     = NULL;

	mSigma      = NULL;
	mBeta_A     = NULL;
	mBeta_b     = NULL;
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

	free(mThetaInits);
	free(mBeta);
	free(mTheta);
	free(mBetaT);
	free(mThetaT);

	free(mSigma);
	free(mBeta_A);
	free(mBeta_b);

	initPointers();
}

void BlockComp::setLikForm(LikForm form) {
	// some cases will require a new setData() call,
	// so always cleanup and require setData() just in case
	cleanup();

	mLikForm = form;
}

// specify data to fit model to
void BlockComp::setData(int n, double *y, double *S, int nblocks, int *B, int p, double *X, int npairs, int *neighbors) {
	int i,j,k;
	int blk1, blk2;
	int nIn1, nIn2;

	// data is loaded... require new fit
	mHasFit  = false;

	cleanup();

	// load the new data
printf("n=%d, y[0]=%.2f, nblocks=%d, B[0]=%d, p=%d, X[0]=%.2f, npairs=%d, neighbors[0]=%d, mHasFit=%d\n",
	n, y[0], nblocks, B[0], p, X[0], npairs, neighbors[0], mHasFit);

	mN = n;
	mY = y;

	mS       = S;
	mNblocks = nblocks;
	mB       = B;

	mP = p;
	mX = X;

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
	mSigma = (double *)malloc(sizeof(double)*symi(0,mMaxPair));
	for (i = 0; i < symi(0,mMaxPair); i++) { mSigma[i] = 0; }
MSG("max pair=%d with length %d\n", mMaxPair, symi(0, mMaxPair));

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

void BlockComp::setInits(int ntheta, double *theta) {
	free(mThetaInits);

	mThetaInits = (double *)malloc(sizeof(double)*ntheta);
	for (int i = 0; i < ntheta; i++) {
		mThetaInits[i] = theta[i];
	}
}

// fit model to data
bool BlockComp::fit(bool verbose) {
	int i;
	Cov *cov;
	int nParams;

	if (verbose) { MSG("Starting fit()...\n"); }

	// make sure we have initial values
	if (mThetaInits == NULL) {
		if (verbose) { MSG("Missing initial values for spatial parameters.\n"); }
		return(false);
	}

	// obtain class for working with covariance type
	if (mCovType == Exp) {
		cov = new CovExp();
	} else {
		return(false);
	}

	// number of sptial parameters
	nParams = cov->numParams();

	// initialize spatial params
	free(mTheta);
	mTheta = (double *)malloc(sizeof(double) * nParams);
	for (i = 0; i < nParams; i++) {
		mTheta[i] = mThetaInits[i];
	}

	// get initial beta based on initial theta
	free(mBeta);
	free(mBeta_A);
	free(mBeta_b);
	mBeta   = (double *)malloc(sizeof(double)*mP);
	mBeta_A = (double *)malloc(sizeof(double)*symi(0,mP));
	mBeta_b = (double *)malloc(sizeof(double)*mP);

	updateBeta(cov);

	delete cov;

	return(true);
}

// update mean parameters
void BlockComp::updateBeta(Cov *cov) {
	int i,j,k,c,e;
	int pair;
	int blk1,blk2;
	int N_in_pair,elem,kelem;

	// initialize A and b
	for (i = 0; i < symi(0,mP); i++) { mBeta_A[i] = 0; }
	for (i = 0; i < mP;         i++) { mBeta_b[i] = 0; }

	for (pair = 0; pair < mNpairs; pair++) {
		blk1 = mNeighbors[pair];
		blk2 = mNeighbors[pair+mNpairs];
		N_in_pair = mNB[blk1] + mNB[blk2];

		if (!mConsMem) {
			// fill in covariance matrix between these two blocks
			cov->compute(mSigma, mNB[blk1], mTheta, mWithinD[blk1], 0);
			cov->compute(mSigma, mNB[blk2], mTheta, mWithinD[blk2], mNB[blk1]);
			cov->cross(mSigma, mNB[blk1], mNB[blk2], mTheta, mBetweenD[pair]);

			// invert Sigma
			chol2inv(N_in_pair, mSigma);

			// add contribution to A and b
			for (i = 0; i < mP; i++) {
				for (j = i; j < mP; j++) {
					// A_ij...
					c = symi(i,j);

					for (k = 0; k < N_in_pair; k++) {
						// mSigma_.k
						if (k < mNB[blk1]) {
							kelem = mWhichB[blk1][k];
						} else {
							kelem = mWhichB[blk2][k-mNB[blk1]];
						}

						for (e = 0; e < N_in_pair; e++) {
							if (e < mNB[blk1]) {
								elem = mWhichB[blk1][e];
							} else {
								elem = mWhichB[blk2][e-mNB[blk1]];
							}

							mBeta_A[c] += mX[elem + i*mN] * mSigma[symi(e,k)] * mX[kelem + j*mN];

							if (i == j) {
								mBeta_b[i] += mX[elem + i*mN] * mSigma[symi(e,k)] * mY[kelem];
							}
						}

					}

				}
			}
		} else {
MSG("TODO\n");
		}
	}

	disp_sym(mBeta_A, 0, mP);
	for (i = 0; i < mP; i++) { MSG("%.2f ", mBeta_b[i]); } MSG("\n");

	// compute beta
	chol2inv(mP, mBeta_A);
	for (i = 0; i < mP; i++) {
		mBeta[i] = 0;

		for (j = 0; j < mP; j++) {
			mBeta[i] += mBeta_A[i + j*mP] * mBeta_b[j];
		}
	}

	for (i = 0; i < mP; i++) { MSG("%.2f ", mBeta[i]); } MSG("\n");
}

#ifndef BLK

// test the estimation on the data in test_data.h

#include "test_data.h"

int main(void) {
	BlockComp blk;

/*
	CovExp    cov;
	int    i,j;
	double theta[] = { log(1), log(1), log(1/0.2) };

	double test_D[ symi(0, test_n) ];
	double test_Sigma[ symi(0,test_n) ];

	// compute distances
	for (i = 0; i < test_n-1; i++) {
		test_D[ symi(i,i) ] = 0;

		for (j = i+1; j < test_n; j++) {
			test_D[ symi(i,j) ] = sqrt( pow(test_S[i]-test_S[j],2) + pow(test_S[i+test_n]-test_S[j+test_n],2) );
		}
	}

	// fill test_Sigma
	cov.compute(test_Sigma, test_n, theta, test_D);

MSG("%.2f, %.2f, %.2f\n", test_D[symi(0,0)], test_D[symi(0,1)], test_D[symi(1,1)]);
MSG("%.2f, %.2f, %.2f\n", test_Sigma[symi(0,0)], test_Sigma[symi(0,1)], test_Sigma[symi(1,1)]);

	// compute Chol(test_Sigma)
	char uplo = 'U';
	int info;

	dpptrf_(&uplo, &test_n, test_Sigma, &info);
	if (info) {
		MSG("Error with chol(): info = %d\n", info);
	}

MSG("info = %d\n", info);
MSG("%.2f, %.2f, %.2f\n", test_Sigma[symi(0,0)], test_Sigma[symi(0,1)], test_Sigma[symi(1,1)]);

	// invert Chol(test_Sigma)
	dpptri_(&uplo, &test_n, test_Sigma, &info);
	if (info) {
		MSG("Error with chol invert(): info = %d\n", info);
	}

MSG("info = %d\n", info);
MSG("%.2f, %.2f, %.2f\n", test_Sigma[symi(0,0)], test_Sigma[symi(0,1)], test_Sigma[symi(1,1)]);
return(0);
*/

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
	double inits[] = {0.5, 0.7243848, 0.2105263};
	blk.setInits(3, inits);
	blk.fit(true);

	return(0);
}
#endif
