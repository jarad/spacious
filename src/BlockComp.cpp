// Estimates block composite models with Fisher scoring
#include <stdio.h>
#include <unistd.h>
#include <time.h>

#include <R.h>

#include "BlockComp.h"
#include "covs.h"
#include "utils.h"

#ifdef PTHREAD
#include <pthread.h>
#endif

#ifdef CUDA
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "utils_cuda.h"
#endif

#include <R_ext/Lapack.h>

// constructor
BlockComp::BlockComp() {
	init(1, false);
}

BlockComp::BlockComp(int nthreads, bool gpu) {
#ifdef PTHREAD
	init(nthreads, gpu);
#else
	init(1, gpu);
#endif
}

BlockComp::~BlockComp() {
	cleanup();

	// only clear these to end
	delete mCov;
	free(mThetaInits);
	free(mFixed);
	free(mFixedVals);

#ifdef CUDA
	if (mGPU) {
		cublasDestroy(mCublasHandle);
		cudaDeviceReset();
	}
#endif
}

void BlockComp::init(int nthreads, bool gpu) {
	initPointers();

	if (nthreads < 1) nthreads = 1;   // make sure we have a least one thread

	setThreads(nthreads);

	// only initalize these pointers to start
	mCov        = NULL;
	mThetaInits = NULL;
	mFixed      = NULL;
	mFixedVals  = NULL;

	// default booleans
	mGPU        = gpu;
	mConsMem    = false;
	mConsMemGPU = false;
	mHasFit     = false;
	mConverged  = false;

	// default setup
	setLikForm(Block);
	setCovType(Exp);

	// default control params
	mIterTol = 1e-5;
	mMaxIter = 20;

#ifdef CUDA
	if (mGPU) {

		if (cublasCreate(&mCublasHandle) != CUBLAS_STATUS_SUCCESS) {
#ifdef CLINE
			MSG("init(): Unable to initialize CUBLAS.\n");
			exit(-1);
#else
			error("init(): Unable to initialize CUBLAS.\n");
#endif
		}

#ifdef DEBUG
		cuda_devices();
#endif
	}
#endif
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
	mIterBeta   = NULL;
	mIterTheta  = NULL;
	mIterLogLik = NULL;

	mFitted     = NULL;
	mResids     = NULL;

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

#ifdef CUDA
	mDevSigma     = NULL;
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
	free(mIterBeta);
	free(mIterTheta);
	free(mIterLogLik);

	free(mFitted);
	free(mResids);

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

#ifdef CUDA
	if (mGPU) {
		if (mDevSigma != NULL) {
			if (cudaFree(mDevSigma) != cudaSuccess) {
				MSG("cleanup(): Unable to free device memory for mDevSigma: %s\n", cudaGetErrorString(cudaGetLastError()));
			}
		}
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
bool BlockComp::setData(int n, double *y, double *S, int nblocks, int *B, int p, double *X, int npairs, int *neighbors) {
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
		mSigma[i] = (double *)malloc(sizeof(double)*mMaxPair*mMaxPair);

		for (j = 0; j < mMaxPair*mMaxPair; j++) { mSigma[i][j] = 0; }

#ifdef CUDA
		if (mGPU) {
			if (cudaMalloc( (void **) &mDevSigma, mMaxPair*mMaxPair*sizeof(double) ) != cudaSuccess) {
				MSG("setData(): Unable to allocate device memory for mDevSigma: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}
		}
#endif

#ifdef DEBUG
	MSG("max pair=%d with length %d\n", mMaxPair, mMaxPair*mMaxPair);
#endif
	}

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

			computeWithinDistance(mN, mS, mWithinD[0]);
		}

	}  // end memory conservation check

	return(true);
}

void BlockComp::computeWithinDistance(int n, const double *S, double *D) {
	int i,j;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			if (i == j) {
				D[ symi(j,j) ] = 0;
			} else {
				D[ symi(i,j) ] = sqrt( pow(S[i]-S[j], 2) + pow(S[i+n] - S[j+n], 2) );
			}
		}
	}
}

void BlockComp::computeBetweenDistance(int n1, const double *S1, int n2, const double *S2, double *D) {
	int i,j;

	for (i = 0; i < n1; i++) {
		for (j = 0; j < n2; j++) {
			D[i + j*n1] = sqrt( pow(S1[i]-S2[j], 2) + pow(S1[i+n1]-S2[j+n2], 2) );
		}
	}
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

// compute log-likelihood
bool BlockComp::computeLogLik(double *log_lik) {
	int i;
	double log_det;
	double q[mMaxPair];

	char   cN = 'N';
	double p1 = 1.0;
	int    i1 = 1;

	// initialization
	log_det = 0;
	*log_lik = 0;

	if (mLikForm == Block) {
		int pair;
		double resids[mMaxPair];

#ifdef PTHREAD
		if (mNthreads <= 1) {
#endif
			// process each block pair in order
			for (pair = 0; pair < mNpairs; pair++) {
				if (!computeLogLikPair(log_lik, pair, mSigma[0], resids, q)) {
						// error updating this pair
						MSG("Unable to compute log-likelihood for pair %d\n", pair);
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
				mThreadStatus[i]        = true;
				mThreadWork[i]->id      = i;
				mThreadWork[i]->log_lik = 0;
				mThreadWork[i]->bc      = this;

				pthread_create(&mThreads[i], NULL,
					&BlockComp::computeLogLikThread,
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
					MSG("computeLogLik(): Error processing thread %d.\n", i);
					return(false);
				}
			}

			// combine results from each thread
			for (i = 0; i < mNthreads; i++) {
				*log_lik += mThreadWork[i]->log_lik;
			}
		}
#endif

		*log_lik /= -2;
	} else if (mLikForm == Full) {
		for (i = 0; i < mN; i++) q[i] = 0;

		// invert covariance
		if (!invertFullCov(true, &log_det)) return(false);

		// add determinant piece to log-likelihood
		*log_lik += log_det;
#ifdef CUDA
		if (!mGPU) {
#endif
			computeResiduals();

			// compute q = inv(Sigma) x resids
			dgemv_(&cN, &mN, &mN, &p1, mSigma[0], &mN, mResids, &i1, &p1, q, &i1);
#ifdef CUDA
		} else {
#endif
			// use GPU
			double p1 = 1.0;
			double n1 = -1.0;
			int    i1 = 1;
			cublasStatus_t status;

			double *devX;
			double *devBeta;
			double *devResids;
			double *devq;

			// copy host Sigma to device
			if (cudaMemcpy(mDevSigma, mSigma[0], mN*mN*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("computeLogLik(): unable to copy Sigma to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// allocate space on device
			if (cudaMalloc((void **)&devX, mN*mNbeta*sizeof(double)) != cudaSuccess) {
				MSG("computeLogLik(): unable to allocate space for X on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devBeta, mNbeta*sizeof(double)) != cudaSuccess) {
				MSG("computeLogLik(): unable to allocate space for beta on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devResids, mN*sizeof(double)) != cudaSuccess) {
				MSG("computeLogLik(): unable to allocate space for resids on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devq, mN*sizeof(double)) != cudaSuccess) {
				MSG("computeLogLik(): unable to allocate space for q on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// copy host data to device
			if (cudaMemcpy(devX, mX, mN*mNbeta*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("computeLogLik(): unable to copy X to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMemcpy(devBeta, mBeta, mNbeta*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("computeLogLik(): unable to copy beta to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMemcpy(devResids, mY, mN*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("computeLogLik(): unable to copy resids to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMemcpy(devq, q, mN*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("computeLogLik(): unable to copy q to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// compute resids = y - Xb
			status = cublasDgemv(mCublasHandle, CUBLAS_OP_N, mN, mNbeta,
			                     &n1, devX, mN, devBeta, i1, &p1, devResids, i1);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("computeLogLik(): unable to call cublasDgemv(): %d\n", status);
				return(false);
			}

			// compute q = inv(Sigma) x resids
			status = cublasDgemv(mCublasHandle, CUBLAS_OP_N, mN, mN,
			                     &p1, mDevSigma, mN, devResids, i1, &p1, devq, i1);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("computeLogLik(): unable to call cublasDgemv(): %d\n", status);
				return(false);
			}

			// get back resids
			if (cudaMemcpy(mResids, devResids, mN*sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) {
				MSG("computeLogLik(): unable to copy resids to host: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// get back q
			if (cudaMemcpy(q, devq, mN*sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) {
				MSG("computeLogLik(): unable to copy q to host: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			cudaFree(devX);
			cudaFree(devBeta);
			cudaFree(devResids);
			cudaFree(devq);

#ifdef CUDA
		}
#endif

		// add in resids x inv(Sigma) x resids piece to log-likelihood
		for (i = 0; i < mN; i++)
			*log_lik += mResids[i] * q[i];

		*log_lik /= -2;
	}

	return(true);
}

bool BlockComp::computeLogLikPair(double *log_lik, int pair, double *Sigma, double *resids, double *q) {
	int blk1 = mNeighbors[pair];
	int blk2 = mNeighbors[pair+mNpairs];
	int N_in_pair = mNB[blk1] + mNB[blk2];

	int i,j;
	int c;

	if (!mConsMem) {
		// fill in covariance matrix between these two blocks
		mCov->compute(Sigma, mTheta, mNB[blk1], mWithinD[blk1], mNB[blk2], mWithinD[blk2], mBetweenD[pair]);
	} else {
		// we're conserving memory

		// fill in distance
		// fill in covariance
		MSG("TODO: computeLogLikPair(): mConsMem=true, mLikForm=Block\n");
		return(false);
	}

	// invert Sigma
	if (chol2inv(N_in_pair, Sigma)) {
		MSG("computeLogLikPair(): Unable to invert Sigma\n");
		return(false);
	}

	// fill in lower triangle
	for (i = 0; i < N_in_pair; i++) {
		for (j = 0; j < N_in_pair; j++) {
			Sigma[lsymi(i,j,N_in_pair)] = Sigma[usymi(i,j,N_in_pair)];
		}
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

	char   cN = 'N';
	double p1 = 1.0;
	double z = 0;
	int    i1 = 1;

	// compute q = inv(Sigma) x resids
	dgemv_(&cN, &N_in_pair, &N_in_pair, &p1, Sigma, &N_in_pair, resids, &i1, &z, q, &i1);

	// add in resids x inv(Sigma) x resids piece to log-likelihood
	for (i = 0; i < N_in_pair; i++)
		*log_lik += resids[i] * q[i];

	return(true);
}

#ifdef PTHREAD
void *BlockComp::computeLogLikThread(void *work) {
	pair_update_t *w = (pair_update_t *)work;
	int            id = w->id;
	BlockComp     *bc = w->bc;

	int pair;
	double resids[bc->mMaxPair];
	double q[bc->mMaxPair];

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

		if (!bc->computeLogLikPair(&w->log_lik, pair, bc->mSigma[id], resids, q)) {
			// error updating this pair
			bc->mThreadStatus[id] = false;
			return(NULL);
		}
	}

	bc->mThreadStatus[id] = true;
	return(NULL);
}
#endif


// fit model to data
bool BlockComp::fit(bool verbose) {
	int i;

#ifdef DEBUG
	MSG("Starting fit() (# of threads=%d, use gpu=%d)...\n", mNthreads, mGPU);
#endif

	// allocate space for residuals and fitted
	free(mResids);
	free(mFitted);
	mResids = (double *)malloc(sizeof(double)*mN);
	mFitted = (double *)malloc(sizeof(double)*mN);

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
	mBeta_A = (double *)malloc(sizeof(double)*mNbeta*mNbeta);
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
	mTheta_H = (double *)malloc(sizeof(double)*mNtheta*mNtheta);
	mTheta_P = (double *)malloc(sizeof(double)*mMaxPair*mMaxPair);
	for (i = 0; i < mNtheta; i++) {
		mTheta_W[i] = (double *)malloc(sizeof(double)*mMaxPair*mMaxPair);
	}

	// prepare where to save at each iteration
	free(mIterBeta);
	free(mIterTheta);
	free(mIterLogLik);
	mIterBeta   = (double *)malloc(sizeof(double)*mNbeta*(mMaxIter+1));
	mIterTheta  = (double *)malloc(sizeof(double)*mNtheta*(mMaxIter+1));
	mIterLogLik = (double *)malloc(sizeof(double)*(mMaxIter+1));

#ifdef PTHREAD
	free(mThreads);
	free(mThreadStatus);

	if (mNthreads > 1) {
#ifdef DEBUG
		MSG("TODO: order block pairs for threading!\n");
#endif

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
			mThreadWork[i] = (pair_update_t *)malloc(sizeof(pair_update_t));
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
			mBeta_A_t[i] = (double *)malloc(sizeof(double)*mNbeta*mNbeta);
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
				mTheta_W_t[i][j] = (double *)malloc(sizeof(double)*mMaxPair*mMaxPair);
			}

			mTheta_H_t[i] = (double *)malloc(sizeof(double)*mNtheta*mNtheta);
			mTheta_P_t[i] = (double *)malloc(sizeof(double)*mMaxPair*mMaxPair);
			mTheta_u_t[i] = (double *)malloc(sizeof(double)*mNtheta);
		}

	}

#endif

	if (mLikForm == Full) {
		// invert full covariance once to start
		if (!invertFullCov()) {
			MSG("Unable to invert full covariance using initial values.\n");
			return(false);
		}
	}

	// get initial beta
	if (!updateBeta()) {
		MSG("Unable to get initial values for beta\n");
		return(false);
	}

	// compute log likelihood
	double log_lik;
	if (!computeLogLik(&log_lik)) return(false);

	// save initial values
	for (i = 0; i < mNbeta; i++)  { mIterBeta[i]   = mBeta[i];  }
	for (i = 0; i < mNtheta; i++) { mIterTheta[i]  = mTheta[i]; }
	mIterLogLik[0] = log_lik;

	for (mIters = 0; mIters < mMaxIter; mIters++) {
		// update covariance params
		if (!updateTheta()) {
			MSG("Unable to update theta at iteration %d\n", mIters+1);
			return(false);
		}

		// update mean params
		if (!updateBeta()) {
			MSG("Unable to update beta at iteration %d\n", mIters+1);
			return(false);
		}

		// compute log-likelihood for updated params
		computeLogLik(&log_lik);

		// save values at this iteration
		for (i = 0; i < mNbeta; i++)  { mIterBeta[i + mNbeta*(mIters+1)]  = mBeta[i];  }
		for (i = 0; i < mNtheta; i++) { mIterTheta[i + mNbeta*(mIters+1)] = mTheta[i]; }
		mIterLogLik[mIters+1] = log_lik;

		if (verbose) {
			// display iteration information
			MSG("iter=%d: ll: %.2f; ", mIters+1, log_lik);
			MSG("beta: ");
			for (i = 0; i < mNbeta; i++) { MSG("%.2f ", mBeta[i]); }
			MSG("; theta: ");
			for (i = 0; i < mNtheta; i++) { MSG("%.2f ", mTheta[i]); }
			MSG("\n");
		}

		// check log-likelihood for convergence
		if (fabs(mIterLogLik[mIters+1] - mIterLogLik[mIters])/fabs(0.1 + mIterLogLik[mIters+1]) <= mIterTol) {
			// convergence!
			if (verbose) {
				MSG("Converged at iteration %d\n", mIters+1);
			}

			mConverged = true;
			mIters++;

			break;
		}

	}

	// we have a fit!
	mHasFit = true;

	return(true);
}

// invert full covariance matrix
bool BlockComp::invertFullCov(bool do_log_det, double *log_det) {
	int i,j;

	// fill in covariance matrix
	mCov->compute(mSigma[0], mTheta, mN, mWithinD[0]);

#ifdef CUDA
	if (!mGPU) {  // compiled with cuda, but don't use GPU...
#endif
		if (chol2inv(mN, mSigma[0], do_log_det, log_det)) {
			MSG("invertFullCov(): unable to invert Sigma\n");
			return(false);
		}
#ifdef CUDA
	} else {
		// use GPU

		// copy host Sigma to device
		if (cudaMemcpy(mDevSigma, mSigma[0], mN*mN*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
			MSG("invertFullCov(): unable to copy Sigma to device: %s\n", cudaGetErrorString(cudaGetLastError()));
			return(false);
		}

		// invert Sigma with GPU
		if (cuda_chol2inv(mCublasHandle, mN, mDevSigma, do_log_det, log_det)) {
			MSG("invertFullCov(): unable to invert Sigma\n");
			return(false);
		}

		// copy device Sigma to host
		if (cudaMemcpy(mSigma[0], mDevSigma, mN*mN*sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) {
			MSG("invertFullCov(): unable to copy Sigma to host: %s\n", cudaGetErrorString(cudaGetLastError()));
			return(false);
		}

	}
#endif

	// fill in lower triangle of Sigma
	for (i = 0; i < mN; i++) {
		for (j = i+1; j < mN; j++) {
			mSigma[0][lsymi(i,j,mN)] = mSigma[0][usymi(i,j,mN)];
		}
	}

	return(true);
}

// update mean parameters
bool BlockComp::updateBeta() {
	int i,j;
	int pair;
#ifdef CUDA
	cublasStatus_t status;
#endif

	// initialize A and b
	for (i = 0; i < mNbeta*mNbeta; i++) { mBeta_A[i] = 0; }
	for (i = 0; i < mNbeta;        i++) { mBeta_b[i] = 0; }

	if (mLikForm == Block) {
		// update beta with block composite likelihood

#ifdef PTHREAD
		if (mNthreads <= 1) {
#endif
			// process each block pair in order
			for (pair = 0; pair < mNpairs; pair++) {
				if (!updateBetaPair(pair, mSigma[0], mBeta_A, mBeta_b)) {
					// error updating this pair
					MSG("Unable to update beta for pair %d\n", pair);
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
				for (i = 0; i < mNbeta*mNbeta; i++) { mBeta_A[i] += mBeta_A_t[j][i]; }
				for (i = 0; i < mNbeta;        i++) { mBeta_b[i] += mBeta_b_t[j][i]; }
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

		// note that computeLogLik() already computes inv(Sigma) used here

		// compute A = X'inv(Sigma)X and b = X'inv(Sigma)y
#ifdef CUDA
		if (!mGPU) {  // compiled with cuda, but don't use GPU...
#endif

			double *dXtS = (double *)malloc(sizeof(double)*mNbeta*mN);
			for (i = 0; i < mNbeta*mN; i++) { dXtS[i] = 0; }
			char   cT = 'T';
			char   cN = 'N';
			double p1 = 1.0;
			int    i1 = 1;

			// compute X'inv(Sigma)
			dgemm_(&cT, &cN, &mNbeta, &mN, &mN, &p1, mX, &mN, mSigma[0], &mN, &p1, dXtS, &mNbeta);

			// compute A = X'inv(Sigma)X
			dgemm_(&cN, &cN, &mNbeta, &mNbeta, &mN, &p1, dXtS, &mNbeta, mX, &mN, &p1, mBeta_A, &mNbeta);

			// compute b = X'inv(Sigma)y
			dgemv_(&cN, &mNbeta, &mN, &p1, dXtS, &mNbeta, mY, &i1, &p1, mBeta_b, &i1);

			free(dXtS);

#ifdef CUDA
		} else {
			// use GPU

			double *devX;
			double *devXtS;
			double *devy;
			double *devA;
			double *devb;
			double p1 = 1.0;
			int    i1 = 1;

			// allocate space on device
			if (cudaMalloc((void **)&devX, mN*mNbeta*sizeof(double)) != cudaSuccess) {
				MSG("updateBeta(): unable to allocate space for X on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devXtS, mNbeta*mN*sizeof(double)) != cudaSuccess) {
				MSG("updateBeta(): unable to allocate space for X'inv(Sigma) on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devy, mN*sizeof(double)) != cudaSuccess) {
				MSG("updateBeta(): unable to allocate space for Y on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devA, mNbeta*mNbeta*sizeof(double)) != cudaSuccess) {
				MSG("updateBeta(): unable to allocate space for A on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devb, mNbeta*sizeof(double)) != cudaSuccess) {
				MSG("updateBeta(): unable to allocate space for b on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// copy host Sigma to device
			if (cudaMemcpy(mDevSigma, mSigma[0], mN*mN*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("updateBeta(): unable to copy Sigma to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// copy host X to device
			if (cudaMemcpy(devX, mX, mN*mNbeta*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("updateBeta(): unable to copy X to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// copy host Y to device
			if (cudaMemcpy(devy, mY, mN*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("updateBeta(): unable to copy Y to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// initialize device data
			if (cudaMemset(devXtS, 0, mNbeta*mN*sizeof(double)) != cudaSuccess) {
				MSG("updateBeta(): unable to initialize X'inv(Sigma) on device\n");
				return(false);
			}

			if (cudaMemset(devA, 0, mNbeta*mNbeta*sizeof(double)) != cudaSuccess) {
				MSG("updateBeta(): unable to initialize A on device\n");
				return(false);
			}

			if (cudaMemset(devb, 0, mNbeta*sizeof(double)) != cudaSuccess) {
				MSG("updateBeta(): unable to initialize b on device\n");
				return(false);
			}

			// compute X'inv(Sigma)
			status = cublasDgemm(mCublasHandle, CUBLAS_OP_T, CUBLAS_OP_N, mNbeta, mN, mN,
                           &p1, devX, mN, mDevSigma, mN, &p1, devXtS, mNbeta);
      if (status != CUBLAS_STATUS_SUCCESS) {
        MSG("updateBeta(): unable to call cublasDgemm(): %d\n", status);
        return(false);
      }

			// compute A = X'inv(Sigma)X
			status = cublasDgemm(mCublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, mNbeta, mNbeta, mN,
                           &p1, devXtS, mNbeta, devX, mN, &p1, devA, mNbeta);
      if (status != CUBLAS_STATUS_SUCCESS) {
        MSG("updateBeta(): unable to call cublasDgemm(): %d\n", status);
        return(false);
      }

			// copy A from device
			if (cudaMemcpy(mBeta_A, devA, mNbeta*mNbeta*sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) {
				MSG("updateBeta(): unable to copy A from device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// compute b = X'inv(Sigma)y
			status = cublasDgemv(mCublasHandle, CUBLAS_OP_N, mNbeta, mN,
			                     &p1, devXtS, mNbeta, devy, i1, &p1, devb, i1);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("updateBeta(): unable to call cublasDgemv(): %d\n", status);
				return(false);
			}

			// copy b from device
			if (cudaMemcpy(mBeta_b, devb, mNbeta*sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) {
				MSG("updateBeta(): unable to copy b from device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			cudaFree(devX);
			cudaFree(devXtS);
			cudaFree(devy);
			cudaFree(devA);
			cudaFree(devb);
		}
#endif

	} // end different likelihood forms

	// compute beta
	if (chol2inv(mNbeta, mBeta_A)) {
		MSG("updateBeta(): unable to invert A\n");
		return(false);
	}

	for (i = 0; i < mNbeta; i++) {
		mBeta[i] = 0;

		for (j = 0; j < mNbeta; j++) {
			mBeta[i] += mBeta_A[usymi(i,j,mNbeta)] * mBeta_b[j];
		}
	}

	return(true);
}

bool BlockComp::updateBetaPair(int pair, double *Sigma, double *A, double *b) {
	// update specified pair of blocks
	int blk1 = mNeighbors[pair];
	int blk2 = mNeighbors[pair+mNpairs];
	int N_in_pair = mNB[blk1] + mNB[blk2];
	int i,j,icol,jcol,k,l,selem;

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
			selem = usymi(i,j,mNbeta);
			jcol  = j*mN;

			for (k = 0; k < mNB[blk1]; k++) {     // block 1...
				for (l = 0; l < mNB[blk1]; l++) {   // with block 1
					A[selem] += mX[mWhichB[blk1][l] + icol] * Sigma[usymi(l,k,N_in_pair)] * mX[mWhichB[blk1][k] + jcol];
					if (i == j) {
						b[i] += mX[mWhichB[blk1][l] + icol] * Sigma[usymi(l,k,N_in_pair)] * mY[mWhichB[blk1][k]];
					}
				}

				for (l = 0; l < mNB[blk2]; l++) {   // with block 2
					A[selem] += mX[mWhichB[blk2][l] + icol] * Sigma[usymi(l+mNB[blk1],k,N_in_pair)] * mX[mWhichB[blk1][k] + jcol];
					if (i == j) {
						b[i] += mX[mWhichB[blk2][l] + icol] * Sigma[usymi(l+mNB[blk1],k,N_in_pair)] * mY[mWhichB[blk1][k]];
					}
				}
			}

			for (k = 0; k < mNB[blk2]; k++) {     // block 2...
				for (l = 0; l < mNB[blk1]; l++) {   // with block 1
					A[selem] += mX[mWhichB[blk1][l] + icol] * Sigma[usymi(l,k+mNB[blk1],N_in_pair)] * mX[mWhichB[blk2][k] + jcol];
					if (i == j) {
						b[i] += mX[mWhichB[blk1][l] + icol] * Sigma[usymi(l,k+mNB[blk1],N_in_pair)] * mY[mWhichB[blk2][k]];
					}
				}

				for (l = 0; l < mNB[blk2]; l++) {   // with block 2
					A[selem] += mX[mWhichB[blk2][l] + icol] * Sigma[usymi(l+mNB[blk1],k+mNB[blk1],N_in_pair)] * mX[mWhichB[blk2][k] + jcol];
					if (i == j) {
						b[i] += mX[mWhichB[blk2][l] + icol] * Sigma[usymi(l+mNB[blk1],k+mNB[blk1],N_in_pair)] * mY[mWhichB[blk2][k]];
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
	for (i = 0; i < bc->mNbeta*bc->mNbeta; i++) { bc->mBeta_A_t[id][i] = 0; }
	for (i = 0; i < bc->mNbeta;            i++) { bc->mBeta_b_t[id][i] = 0; }

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
	for (i = 0; i < mNtheta*mNtheta; i++) { mTheta_H[i] = 0; }

	if (mLikForm == Block) {
		// update theta with block composite likelihood

#ifdef PTHREAD
		if (mNthreads <= 1) {
#endif

			// process each block pair in order
			for (pair = 0; pair < mNpairs; pair++) {
				if (!updateThetaPair(pair, mSigma[0], mTheta_W, mTheta_H, mTheta_P, resids, q, u)) {
						// error updating this pair
						MSG("Unable to update theta for pair %d\n", pair);
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
				for (i = 0; i < mNtheta*mNtheta; i++) { mTheta_H[i] += mTheta_H_t[j][i]; }
				for (i = 0; i < mNtheta;         i++) { u[i]        += mTheta_u_t[j][i]; }
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

		// note that logLik() already computes inv(Sigma) and residuals used here
		for (i = 0; i < mN; i++) q[i] = 0;

#ifdef CUDA
		if (!mGPU) {
#endif
			char   cN = 'N';
			double p1 = 1.0;
			double z = 0;
			int    i1 = 1;

			// compute q = inv(Sigma) x resids
			dgemv_(&cN, &mN, &mN, &p1, mSigma[0], &mN, mResids, &i1, &p1, q, &i1);

			// fill in W and u
			// possible parallel: each of these is independent
			for (i = 0; i < mNtheta; i++) {
				if (mFixed[i]) {
					// this parameter is fixed, so skip these operations
					continue;
				}

				// get partial derivatives
				mCov->partials(mTheta_P, &diag, i, mTheta, mThetaT, mN, mWithinD[0]);

				// compute inv(Sigma) x P
				if (diag) {
					// take advantage of P being diagonal

					// initialize W[i]
					for (j = 0; j < mN*mN; j++) mTheta_W[i][j] = 0;

					for (j = 0; j < mN; j++) {
						for (k = 0; k < mN; k++) {
							mTheta_W[i][j + k*mN] = mSigma[0][usymi(j,k,mN)] * mTheta_P[usymi(k,k,mN)];
						}
					}
				} else {
					// fill in lower triangle
					for (j = 0; j < mN; j++) {
						for (k = j+1; k < mN; k++) {
							mTheta_P[lsymi(j,k,mN)] = mTheta_P[usymi(j,k,mN)];
						}
					}

					dgemm_(&cN, &cN, &mN, &mN, &mN, &p1, mSigma[0], &mN, mTheta_P, &mN, &z, mTheta_W[i], &mN);
				}

				// fill in u
				for (j = 0; j < mN; j++) {
					u[i] -= 0.5*mTheta_W[i][j+j*mN];
				}

				for (j = 0; j < mN; j++) {
					u[i] += 0.5*mTheta_P[usymi(j,j,mN)] * q[j] * q[j];
					for (k = j+1; k < mN; k++) {
						u[i] += mTheta_P[usymi(j,k,mN)] * q[j] * q[k];
					}
				}
			}
#ifdef CUDA
		} else {
			// use GPU
			double p1 = 1.0;
			double n1 = -1.0;
			int    i1 = 1;
			cublasStatus_t status;

			double *devX;
			double *devBeta;
			double *devResids;
			double *devq;
			double *devP;
			double *devW;

			// copy host Sigma to device
			if (cudaMemcpy(mDevSigma, mSigma[0], mN*mN*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("updateTheta(): unable to copy Sigma to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// allocate space on device
			if (cudaMalloc((void **)&devX, mN*mNbeta*sizeof(double)) != cudaSuccess) {
				MSG("updateTheta(): unable to allocate space for X on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devBeta, mNbeta*sizeof(double)) != cudaSuccess) {
				MSG("updateTheta(): unable to allocate space for beta on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devResids, mN*sizeof(double)) != cudaSuccess) {
				MSG("updateTheta(): unable to allocate space for resids on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devq, mN*sizeof(double)) != cudaSuccess) {
				MSG("updateTheta(): unable to allocate space for q on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devP, mN*mN*sizeof(double)) != cudaSuccess) {
				MSG("updateTheta(): unable to allocate space for P on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMalloc((void **)&devW, mN*mN*sizeof(double)) != cudaSuccess) {
				MSG("updateTheta(): unable to allocate space for P on device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// copy host data to device
			if (cudaMemcpy(devX, mX, mN*mNbeta*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("updateTheta(): unable to copy X to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMemcpy(devBeta, mBeta, mNbeta*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("updateTheta(): unable to copy beta to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMemcpy(devResids, mY, mN*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("updateTheta(): unable to copy resids to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			if (cudaMemcpy(devq, q, mN*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				MSG("updateTheta(): unable to copy q to device: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// compute resids = y - Xb
			status = cublasDgemv(mCublasHandle, CUBLAS_OP_N, mN, mNbeta,
			                     &n1, devX, mN, devBeta, i1, &p1, devResids, i1);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("updateTheta(): unable to call cublasDgemv(): %d\n", status);
				return(false);
			}

			// compute q = inv(Sigma) x resids
			status = cublasDgemv(mCublasHandle, CUBLAS_OP_N, mN, mN,
			                     &p1, mDevSigma, mN, devResids, i1, &p1, devq, i1);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("updateTheta(): unable to call cublasDgemv(): %d\n", status);
				return(false);
			}

			// get back q
			if (cudaMemcpy(q, devq, mN*sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) {
				MSG("updateTheta(): unable to copy q to host: %s\n", cudaGetErrorString(cudaGetLastError()));
				return(false);
			}

			// fill in W and u
			// possible parallel: each of these is independent
			for (i = 0; i < mNtheta; i++) {
				if (mFixed[i]) {
					// this parameter is fixed, so skip these operations
					continue;
				}

				// get partial derivatives
				mCov->partials(mTheta_P, &diag, i, mTheta, mThetaT, mN, mWithinD[0]);

				// compute inv(Sigma) x P
				if (diag) {
					// take advantage of P being diagonal
					for (j = 0; j < mN; j++) {
						for (k = 0; k < mN; k++) {
							mTheta_W[i][j + k*mN] = mSigma[0][usymi(j,k,mN)] * mTheta_P[usymi(k,k,mN)];
						}
					}
				} else {
					// fill in lower triangle
					for (j = 0; j < mN; j++) {
						for (k = j+1; k < mN; k++) {
							mTheta_P[lsymi(j,k,mN)] = mTheta_P[usymi(j,k,mN)];
						}
					}

					// transfer P to device
					if (cudaMemcpy(devP, mTheta_P, mN*mN*sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
						MSG("updateTheta(): unable to copy P to device: %s\n", cudaGetErrorString(cudaGetLastError()));
						return(false);
					}

					// initialize W[i]
					if (cudaMemset(devW, 0, mN*mN*sizeof(double)) != cudaSuccess) {
						MSG("updateTheta(): unable to initialize W on device\n");
						return(false);
					}

					// perform the multiply
					status = cublasDgemm(mCublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, mN, mN, mN,
		                           &p1, mDevSigma, mN, devP, mN, &p1, devW, mN);
		      if (status != CUBLAS_STATUS_SUCCESS) {
		        MSG("updateTheta(): unable to call cublasDgemm(): %d\n", status);
		        return(false);
		      }

					if (cudaMemcpy(mTheta_W[i], devW, mN*mN*sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) {
						MSG("updateTheta(): unable to copy W to host: %s\n", cudaGetErrorString(cudaGetLastError()));
						return(false);
					}
				}

				// fill in u
				for (j = 0; j < mN; j++) {
					u[i] -= 0.5*mTheta_W[i][j+j*mN];
				}

				for (j = 0; j < mN; j++) {
					u[i] += 0.5*mTheta_P[usymi(j,j,mN)] * q[j] * q[j];
					for (k = j+1; k < mN; k++) {
						u[i] += mTheta_P[usymi(j,k,mN)] * q[j] * q[k];
					}
				}
			}

			cudaFree(devX);
			cudaFree(devBeta);
			cudaFree(devResids);
			cudaFree(devq);
			cudaFree(devP);
			cudaFree(devW);

		}
#endif

		// compute hessian
		iH = 0;
		for (i = 0; i < mNtheta; i++) {
			if (mFixed[i]) { continue; }

			jH = iH;
			for (j = i; j < mNtheta; j++) {
				if (mFixed[j]) { continue; }

				c = usymi(iH,jH,mNtheta-mNfixed);

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

	} // end different likelihood forms

	// hessian elements should be scaled by 1/2
	for (i = 0; i < mNtheta*mNtheta; i++) {
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

			mThetaT[i] += mTheta_H[usymi(iH,jH,mNtheta-mNfixed)] * u[j];

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

	// fill in lower triangle
	for (i = 0; i < N_in_pair; i++) {
		for (j = 0; j < N_in_pair; j++) {
			Sigma[lsymi(i,j,N_in_pair)] = Sigma[usymi(i,j,N_in_pair)];
		}
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

	char   cN = 'N';
	double p1 = 1.0;
	double z = 0;
	int    i1 = 1;

	// compute q = inv(Sigma) x resids
	dgemv_(&cN, &N_in_pair, &N_in_pair, &p1, Sigma, &N_in_pair, resids, &i1, &z, q, &i1);

	// fill in W and u
	for (i = 0; i < mNtheta; i++) {
		if (mFixed[i]) {
			// this parameter is fixed, so skip these operations
			continue;
		}

		// get partial derivatives
		mCov->partials(P, &diag, i, mTheta, mThetaT, mNB[blk1], mWithinD[blk1], mNB[blk2], mWithinD[blk2], mBetweenD[pair]);

		// compute inv(Sigma) x P
		if (diag) {
			// take advantage of P being diagonal

			// initialize W[i]
			for (j = 0; j < N_in_pair*N_in_pair; j++) W[i][j] = 0;

			for (j = 0; j < N_in_pair; j++) {
				for (k = 0; k < N_in_pair; k++) {
					W[i][j + k*N_in_pair] = Sigma[usymi(j,k,N_in_pair)] * P[usymi(k,k,N_in_pair)];
				}
			}
		} else {
			// fill in lower triangle
			for (j = 0; j < N_in_pair; j++) {
				for (k = j+1; k < N_in_pair; k++) {
					P[lsymi(j,k,N_in_pair)] = P[usymi(j,k,N_in_pair)];
				}
			}

			dgemm_(&cN, &cN, &N_in_pair, &N_in_pair, &N_in_pair, &p1, Sigma, &N_in_pair, P, &N_in_pair, &z, W[i], &N_in_pair);
		}

		// fill in u
		for (j = 0; j < N_in_pair; j++) {
			u[i] -= 0.5*W[i][j+j*N_in_pair];
		}

		for (j = 0; j < N_in_pair; j++) {
			u[i] += 0.5*P[usymi(j,j,N_in_pair)] * q[j] * q[j];
			for (k = j+1; k < N_in_pair; k++) {
				u[i] += P[usymi(j,k,N_in_pair)] * q[j] * q[k];
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

			c = usymi(iH,jH,mNtheta-mNfixed);

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
	for (i = 0; i < bc->mNtheta;             i++) { bc->mTheta_u_t[id][i] = 0; }
	for (i = 0; i < bc->mNtheta*bc->mNtheta; i++) { bc->mTheta_H_t[id][i] = 0; }

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

void BlockComp::computeFitted() {
	// computed fitted values to model data

	computeFitted(mN, mFitted, mX);
}

void BlockComp::computeFitted(int n, double *fitted, double *X) {
	// compute fitted values for covariates in X
	char   cN = 'N';
	double p1 = 1.0;
	double z  = 0.0;
	int    i1 = 1;

	// compute fitted = Xb
	dgemv_(&cN, &n, &mNbeta, &p1, X, &n, mBeta, &i1, &z, fitted, &i1);
}

void BlockComp::computeResiduals() {
	// compute residuals for model data

	// make sure we have fitted values
	computeFitted();

	for (int i = 0; i < mN; i++) {
		mResids[i] = mY[i]-mFitted[i];
	}
}

void BlockComp::getFitted(double *fitted) {
	computeFitted();

	if (mFitted == NULL) return;

	for (int i = 0; i < mN; i++) fitted[i] = mFitted[i];
}

void BlockComp::getResiduals(double *resids) {
	computeResiduals();

	if (mResids == NULL) return;

	for (int i = 0; i < mN; i++) resids[i] = mResids[i];
}

bool BlockComp::predict(int n_0, double *y_0, const double *newS, const int *newB, const double *newX,
                        bool do_sd, double *sd, bool local, int Nlocal) {

	if (!mHasFit) {
		MSG("predict() requires a fitted model.\n");
		return(false);
	}

	int i,j;

	// computed fitted and residuals
	computeResiduals();

	if (local) {
		// predict with local kriging
		MSG("TODO: predict(): implement local kriging\n");
		return(false);
	} else if (mLikForm == Block) {
		// predict with block composite likelihood

		int ix,jx;
		int nNewB[mNblocks];
		int nMaxInB = 0;

		// compute how many prediction sites are in each block
		for (i = 0; i < mNblocks; i++) nNewB[i] = 0;
		for (i = 0; i < n_0; i++) nNewB[newB[i]]++;

		// compute largest number of prediction sites in a block
		for (i = 0; i < mNblocks; i++) if (nNewB[i] > nMaxInB) nMaxInB = nNewB[i];

		// store locations in each block
		int **newWhichB = (int **)malloc(sizeof(int *)*mNblocks);
		for (i = 0; i < mNblocks; i++) {
			if (nNewB[i] == 0) continue;
			newWhichB[i] = (int *)malloc(sizeof(int)*nNewB[i]);
		}

		int locs[mNblocks];
		for (i = 0; i < mNblocks; i++) { locs[i] = 0; }
		for (i = 0; i < n_0; i++) newWhichB[ newB[i] ][ locs[newB[i]]++ ] = i;

		// variables re-used for prediction in each block
		double y_block[nMaxInB];
		double X_block[nMaxInB*mNbeta];
		double S_block[2*nMaxInB];
		double sd_block[nMaxInB];

		for (i = 0; i < nMaxInB; i++)        y_block[i] = sd_block[i] = 0;
		for (i = 0; i < nMaxInB*mNbeta; i++) X_block[i] = 0;
		for (i = 0; i < 2*nMaxInB; i++)      S_block[i] = 0;

		for (i = 0; i < mNblocks; i++) {
			if (nNewB[i] <= 0) continue;

#ifdef DEBUG
			MSG("%d predictions at block %d\n", nNewB[i], i);
#endif

			// fill in S
			for (j = 0; j < nNewB[i]; j++) {
				S_block[j]          = newS[newWhichB[i][j]];
				S_block[j+nNewB[i]] = newS[newWhichB[i][j]+n_0];
			}

			// fill in X
			for (jx = 0; jx < mNbeta; jx++) {
				for (ix = 0; ix < nNewB[i]; ix++) {
					X_block[ix + jx*nNewB[i]] = newX[newWhichB[i][ix] + jx*n_0];
				}
			}

			if (!blockPredict(i, nNewB[i], y_block, S_block, X_block, do_sd, sd_block)) {
				MSG("Unable to predict sites in block %d\n", i);

				for (j = 0; j < mNblocks; j++) {
					if (nNewB[j] == 0) continue;
					free(newWhichB[j]);
				}
				free(newWhichB);

				return(false);
			}

			// fill in results
			for (ix = 0; ix < nNewB[i]; ix++) y_0[ newWhichB[i][ix] ] = y_block[ix];

			if (do_sd) {
				for (ix = 0; ix < nNewB[i]; ix++) sd[ newWhichB[i][ix] ] = sd_block[ix];
			}
		}

		for (i = 0; i < mNblocks; i++) {
			if (nNewB[i] == 0) continue;
			free(newWhichB[i]);
		}
		free(newWhichB);
	} else if (mLikForm == Full) {
		// predict using full likelihood

		// make room for fitted + predicted locations
		int combN = mN+n_0;
		double *combSigma = (double *)malloc(combN*combN*sizeof(double));
		double *combBetweenD = (double *)malloc(mN*n_0*sizeof(double));
		double *newD = (double *)malloc(symi(0,n_0)*sizeof(double));

		// fill in distance
		computeWithinDistance(n_0, newS, newD);
		computeBetweenDistance(mN, mS, n_0, newS, combBetweenD);

		// fill in covariance matrix
		mCov->compute(combSigma, mTheta, mN, mWithinD[0], n_0, newD, combBetweenD);

		// invert fitted Sigma
		for (i = 0; i < mN; i++) {
			for (j = i; j < mN; j++) {
				mSigma[0][usymi(i,j,mN)] = combSigma[usymi(i,j,combN)];
			}
		}

		if (chol2inv(mN, mSigma[0])) {
			MSG("Unable to invert fitted Sigma\n");
			free(combSigma); free(combBetweenD); free(newD);
			return(false);
		}

		// fill in lower tringular piece after inversion
		for (i = 0; i < mN; i++) {
			for (j = i+1; j < mN; j++) {
				mSigma[0][lsymi(i,j,mN)] = mSigma[0][usymi(i,j,mN)];
			}
		}

		char   cN = 'N';
		double p1 = 1.0;
		double z  = 0.0;
		int    i1 = 1;
		double corrvec[mN];

		// compute new y_0 = Xb
		dgemv_(&cN, &n_0, &mNbeta, &p1, newX, &n_0, mBeta, &i1, &z, y_0, &i1);

		// compute inv(Sigma) x resids
		dgemv_(&cN, &mN, &mN, &p1, mSigma[0], &mN, mResids, &i1, &z, corrvec, &i1);

		// finish y_0 = Xb + Sigma[new,fit] inv(Sigma[fit,fit]) resids
		double covBetween[n_0*mN];
		for (j = 0; j < mN; j++) {
			for (i = 0; i < n_0; i++) {
				covBetween[i + j*n_0] = combSigma[usymi(i+mN,j,combN)];
			}
		}
		dgemv_(&cN, &n_0, &mN, &p1, covBetween, &n_0, corrvec, &i1, &p1, y_0, &i1);

		if (do_sd && sd != NULL) {
			// fill sd with prediction variances

			// invert combined Sigma
			if (chol2inv(combN, combSigma)) {
				MSG("Unable to invert combined Sigma\n");
				free(combSigma); free(combBetweenD); free(newD);
				return(false);
			}

			// invert the new data block of inv(combSigma)
			for (i = 0; i < n_0; i++) {
				for (j = i; j < n_0; j++) {
					mSigma[0][usymi(i,j,n_0)] = combSigma[usymi(i+mN,j+mN,combN)];
				}
			}

			if (chol2inv(n_0, mSigma[0])) {
				MSG("Unable to invert new data block of inv(combSigma)\n");
				free(combSigma); free(combBetweenD); free(newD);
				return(false);
			}

			// fill in std devs
			for (i = 0; i < n_0; i++) {
				sd[i] = sqrt(mSigma[0][usymi(i,i,n_0)]);
			}
		}

		free(combSigma);
		free(combBetweenD);
		free(newD);
	}

	return(true);
}

// predict at sites in block
bool BlockComp::blockPredict(int block, int n_0, double *y_0,
	const double *newS, const double *newX, bool do_sd, double *sd) {

	char   cN = 'N';
	char   cT = 'T';
	double p1 = 1.0;
	int    i1 = 1;
	double zero = 0.0;

	int i,j,k;
	int ip,ip2;
	int neighbor,neighbor2;
	int N_in_pair,N_in_pair2;
	bool trans,trans2;

	int    N_block_new = n_0+mNB[block];
	double A_0[n_0*n_0];
	double b_0[n_0];
	double J_0[n_0*n_0];
	double B_0[n_0*N_block_new];

	double pairD[symi(0, mMaxPair+n_0)];
	double pairSigma[(mMaxPair+n_0)*(mMaxPair+n_0)];
	double pairD2[symi(0, mMaxPair+n_0)];
	double pairSigma2[(mMaxPair+n_0)*(mMaxPair+n_0)];
	double neighborSigma[mMaxPair*mMaxPair];

	double predMat[n_0*mMaxPair];
	double predMat2[n_0*mMaxPair];
	double pairResids[mMaxPair];

	// initialize
	for (i = 0; i < n_0*n_0; i++)         A_0[i] = J_0[i] = 0;
	for (i = 0; i < n_0; i++)             b_0[i] = 0;
	for (i = 0; i < n_0*N_block_new; i++) B_0[i] = 0;

	for (i = 0; i < symi(0, mMaxPair+n_0); i++)         pairD[i] = pairD2[i] = 0;
	for (i = 0; i < (mMaxPair+n_0)*(mMaxPair+n_0); i++) pairSigma[i] = pairSigma2[i] = 0;
	for (i = 0; i < mMaxPair*mMaxPair; i++)             neighborSigma[i] = 0;

	for (i = 0; i < n_0*mMaxPair; i++) predMat[i] = predMat2[i] = 0;
	for (i = 0; i < mMaxPair; i++)     pairResids[i] = 0;

	// fill in distance matrix for prediction sites and this block
	// ... within prediction
	for (i = 0; i < n_0; i++)
		for (j = i; j < n_0; j++)
			pairD[symi(i,j)] = sqrt( pow(newS[i]-newS[j], 2) + pow(newS[i+n_0]-newS[j+n_0], 2) );
	// ... between prediction and block
	for (i = 0; i < n_0; i++)
		for (j = 0; j < mNB[block]; j++)
			pairD[symi(i,n_0+j)] = sqrt( pow(newS[i]-mS[mWhichB[block][j]], 2) + pow(newS[i+n_0]-mS[mWhichB[block][j]+mN], 2) );
	// ... within block
	for (i = 0; i < mNB[block]; i++)
		for (j = 0; j < mNB[block]; j++)
			pairD[symi(i+n_0,j+n_0)] = mWithinD[block][symi(i,j)];

	for (i = 0; i < symi(0, mMaxPair+n_0); i++) pairD2[i] = pairD[i];

	for (i = 0; i < mNB[block]; i++) pairResids[i] = mResids[mWhichB[block][i]];

	// for each neighboring block...
	for (ip = 0; ip < mNpairs; ip++) {
		// get the neighbor
		if (mNeighbors[ip] == block) {
			neighbor = mNeighbors[ip + mNpairs];
		} else if (mNeighbors[ip + mNpairs] == block) {
			neighbor = mNeighbors[ip];
		} else {
			continue;   // not a neighbor
		}

		N_in_pair = N_block_new+mNB[neighbor];
		if (block > neighbor) trans = true;
		else                  trans = false;

		// create distance matrix for this pair...
		// ... between prediction and neighbor
		for (i = 0; i < n_0; i++)
			for (j = 0; j < mNB[neighbor]; j++)
				pairD[symi(i,j+N_block_new)] = sqrt( pow(newS[i]-mS[mWhichB[neighbor][j]], 2) + pow(newS[i+n_0]-mS[mWhichB[neighbor][j]+mN], 2) );
		// ... within neighbor
		for (i = 0; i < mNB[neighbor]; i++)
			for (j = i; j < mNB[neighbor]; j++)
				pairD[symi(i+N_block_new,j+N_block_new)] = mWithinD[neighbor][symi(i,j)];
		// ... between block and neighbor
		for (i = 0; i < mNB[block]; i++)
			for (j = 0; j < mNB[neighbor]; j++)
				pairD[symi(i+n_0,j+N_block_new)] = trans ? mBetweenD[ip][j + i*mNB[neighbor]] : mBetweenD[ip][i + j*mNB[block]];

		// compute covariance between pair
		mCov->compute(pairSigma, mTheta, N_in_pair, pairD);

		// invert pairSigma
		if (chol2inv(N_in_pair, pairSigma)) {
			MSG("blockPredict(): Unable to invert Sigma for pair (%d,%d)\n", block, neighbor);
			return(false);
		}

		// add prediction sites of inv(pairSigma) to A_0
		for (i = 0; i < n_0; i++)
			for (j =  0; j < n_0; j++)
				A_0[usymi(i,j,n_0)] += pairSigma[usymi(i,j,N_in_pair)];

		// add to b_0: inv(Sigma_pair)[pred,c(block,neighbor)] x resids[c(block,neighbor)]
		for (i = 0; i < n_0; i++)
			for (j = 0; j < N_in_pair; j++)
				predMat[i+j*n_0] = pairSigma[usymi(i,j,N_in_pair)];

		for (i = 0; i < mNB[neighbor]; i++)
			pairResids[i+mNB[block]] = mResids[mWhichB[neighbor][i]];

		dgemv_(&cN, &n_0, &N_in_pair, &p1, predMat, &n_0, pairResids, &i1, &p1, b_0, &i1);

		if (do_sd) {

			for (i = 0; i < n_0; i++) {
				for (j = 0; j < n_0; j++)        B_0[i+j*n_0]       += pairSigma[usymi(i,j,N_in_pair)];
				for (j = 0; j < mNB[block]; j++) B_0[i+(j+n_0)*n_0] += pairSigma[usymi(i,j+n_0,N_in_pair)];
			}

			// for each neighboring block...
			for (ip2 = 0; ip2 < mNpairs; ip2++) {
				// get the neighbor
				if (mNeighbors[ip2] == block) {
					neighbor2 = mNeighbors[ip2 + mNpairs];
				} else if (mNeighbors[ip2 + mNpairs] == block) {
					neighbor2 = mNeighbors[ip2];
				} else {
					continue;   // not a neighbor
				}

				N_in_pair2 = N_block_new+mNB[neighbor2];
				if (block > neighbor2) trans2 = true;
				else                   trans2 = false;

				// create distance matrix for this pair...
				// ... between prediction and neighbor
				for (i = 0; i < n_0; i++)
					for (j = 0; j < mNB[neighbor2]; j++)
						pairD2[symi(i,j+N_block_new)] = sqrt( pow(newS[i]-mS[mWhichB[neighbor2][j]], 2) + pow(newS[i+n_0]-mS[mWhichB[neighbor2][j]+mN], 2) );
				// ... within neighbor
				for (i = 0; i < mNB[neighbor2]; i++)
					for (j = i; j < mNB[neighbor2]; j++)
						pairD2[symi(i+N_block_new,j+N_block_new)] = mWithinD[neighbor2][symi(i,j)];
				// ... between block and neighbor
				for (i = 0; i < mNB[block]; i++)
					for (j = 0; j < mNB[neighbor2]; j++)
						pairD2[symi(i+n_0,j+N_block_new)] = trans2 ? mBetweenD[ip2][j + i*mNB[neighbor2]] : mBetweenD[ip2][i + j*mNB[block]];

				// compute covariance between pair
				mCov->compute(pairSigma2, mTheta, N_in_pair2, pairD2);

				// invert pairSigma2
				if (chol2inv(N_in_pair2, pairSigma2)) {
					MSG("blockPredict(): Unable to invert Sigma2 for pair (%d,%d)\n", block, neighbor2);
					return(false);
				}

				// compute covariance between neighbor and neighbor2

				// are neighbor and neighbor2 neighbors?
				int neighbor_pair = -1;
				double *neighborBetweenD;
				for (i = 0; i < mNpairs; i++) {
					if ( (mNeighbors[i] == neighbor && mNeighbors[i + mNpairs] == neighbor2) ||
					     (mNeighbors[i] == neighbor2 && mNeighbors[i + mNpairs] == neighbor) ) {
						neighbor_pair = i;
						break;
					}
				}

				if (neighbor_pair < 0) {
					// we need to compute distance between these blocks
					neighborBetweenD = (double *)malloc(sizeof(double)*mNB[neighbor]*mNB[neighbor2]);

					for (i = 0; i < mNB[neighbor]; i++) { 
						for (j = 0; j < mNB[neighbor2]; j++) {
							neighborBetweenD[i + j*mNB[neighbor]] = sqrt(
								pow(mS[mWhichB[neighbor][i]]-mS[mWhichB[neighbor2][j]], 2) + pow(mS[mWhichB[neighbor][i]+mN]-mS[mWhichB[neighbor2][j]+mN], 2)
							);
						}
					}

					mCov->computeCross(neighborSigma, mTheta, mNB[neighbor], mNB[neighbor2], neighborBetweenD, false);

					free(neighborBetweenD);
				} else {
					// we already have distance between these blocks
					if (neighbor > neighbor2) trans2 = true;
					else                      trans2 = false;

					mCov->computeCross(neighborSigma, mTheta, mNB[neighbor], mNB[neighbor2], mBetweenD[neighbor_pair], false, false, trans2);
				}

				// add to J_0: inv(Sigma_pair)[pred,neighbor] x Sigma_n1n2[neighbor,neighbor2] x inv(Sigma_pair2)[neighbor2,pred]
				dgemm_(&cN, &cN, &n_0, &mNB[neighbor2], &mNB[neighbor], &p1, predMat, &n_0, neighborSigma, &mNB[neighbor], &zero, predMat2, &n_0);

				for (i = 0; i < n_0; i++)
					for (j = 0; j < N_in_pair2; j++)
						predMat[i+j*n_0] = pairSigma2[usymi(i,j,N_in_pair2)];

				dgemm_(&cN, &cN, &n_0, &n_0, &mNB[neighbor2], &p1, predMat2, &n_0, predMat, &mNB[neighbor2], &p1, J_0, &n_0);

			}

		}
	}  // end neighbor

	// invert A_0
	if (chol2inv(n_0, A_0)) {
		MSG("blockPredict(): unable to invert A_0\n");
		return(false);
	}

	// fill in lower part of A_0
	for (i = 0; i < n_0; i++)
		for (j = i+1; j < n_0; j++)
			A_0[lsymi(i,j,n_0)] = A_0[usymi(i,j,n_0)];

	// negate b_0
	for (i = 0; i < n_0; i++) b_0[i] = -b_0[i];

	// set y_0[pred] = newX[pred] x beta + inv(A_0) x b_0
	dgemv_(&cN, &n_0, &mNbeta, &p1, newX, &n_0, mBeta, &i1, &zero, y_0, &i1);
	dgemv_(&cN, &n_0, &n_0, &p1, A_0, &n_0, b_0, &i1, &p1, y_0, &i1);

	if (do_sd) {
		// negate A_0
		for (i = 0; i < n_0*n_0; i++) A_0[i] = -A_0[i];

		// add to J_0: B_0 x Sigma_block x t(B_0)...
		// ...compute covariance for this block
		mCov->compute(pairSigma, mTheta, N_block_new, pairD);

		// fill in lower triangular piece
		for (i = 0; i < N_block_new; i++)
			for (j = i+1; j < N_block_new; j++)
				pairSigma[lsymi(i,j,N_block_new)] = pairSigma[usymi(i,j,N_block_new)];

		// ... multiplty B_0 x pairSigma x t(B_0)
		dgemm_(&cN, &cN, &n_0, &mNB[block], &N_block_new, &p1, B_0, &n_0, pairSigma, &N_block_new, &zero, predMat, &n_0);
		dgemm_(&cN, &cT, &n_0, &n_0, &N_block_new, &p1, predMat, &n_0, B_0, &n_0, &p1, J_0, &n_0);

		// for each neighboring block...
		for (ip = 0; ip < mNpairs; ip++) {
			// get the neighbor
			if (mNeighbors[ip] == block) {
				neighbor = mNeighbors[ip + mNpairs];
			} else if (mNeighbors[ip + mNpairs] == block) {
				neighbor = mNeighbors[ip];
			} else {
				continue;   // not a neighbor
			}

			N_in_pair = N_block_new+mNB[neighbor];
			if (block > neighbor) trans = true;
			else                  trans = false;

			// create distance matrix for this pair...
			// ... between prediction and neighbor
			for (i = 0; i < n_0; i++)
				for (j = 0; j < mNB[neighbor]; j++)
					pairD[symi(i,j+N_block_new)] = sqrt( pow(newS[i]-mS[mWhichB[neighbor][j]], 2) + pow(newS[i+n_0]-mS[mWhichB[neighbor][j]+mN], 2) );
			// ... within neighbor
			for (i = 0; i < mNB[neighbor]; i++)
				for (j = i; j < mNB[neighbor]; j++)
					pairD[symi(i+N_block_new,j+N_block_new)] = mWithinD[neighbor][symi(i,j)];
			// ... between block and neighbor
			for (i = 0; i < mNB[block]; i++)
				for (j = 0; j < mNB[neighbor]; j++)
					pairD[symi(i+n_0,j+N_block_new)] = trans ? mBetweenD[ip][j + i*mNB[neighbor]] : mBetweenD[ip][i + j*mNB[block]];

			// compute covariance between pair
			mCov->compute(pairSigma, mTheta, N_in_pair, pairD);

//			for (i = 0; i < N_in_pair;

			// store pairSigma in pairSigma2 for inversion
			for (i = 0; i < N_in_pair; i++)
				for (j = i; j < N_in_pair; j++)
					pairSigma2[usymi(i,j,N_in_pair)] = pairSigma[usymi(i,j,N_in_pair)];

			// invert pairSigma = inv(pairSigma2)
			if (chol2inv(N_in_pair, pairSigma2)) {
				MSG("blockPredict(): Unable to invert Sigma for pair (%d,%d)\n", block, neighbor);
				return(false);
			}

			// fill in lower parts of pairSigma and pairSigma2
			for (i = 0; i < N_in_pair; i++) {
				for (j = i+1; j < N_in_pair; j++) {
					pairSigma[lsymi(i,j,N_in_pair)]  = pairSigma[usymi(i,j,N_in_pair)];
					pairSigma2[lsymi(i,j,N_in_pair)] = pairSigma2[usymi(i,j,N_in_pair)];
				}
			}

			// add 2 x B_0 x sigmaPair[block,neighbor] x invSigmaPair[neighbor,pred] to J_0
			for (i = 0; i < N_block_new; i++) {
				for (j = 0; j < n_0; j++) {
					predMat[i + j*N_block_new] = 0;
					for (k = 0; k < mNB[neighbor]; k++) {
						predMat[i + j*N_block_new] += pairSigma[usymi(i,k,N_in_pair)]*pairSigma2[usymi(k,j,N_in_pair)];
					}
				}
			}

			for (i = 0; i < n_0; i++) {
				for (j = 0; j < n_0; j++) {
					for (k = 0; k < N_block_new; k++) {
						J_0[i + j*n_0] += 2 * B_0[i + k*n_0]*predMat[k + j*N_block_new];
					}
				}
			}

		}

		// invert J_0
		if (chol2inv(n_0, J_0)) {
			MSG("blockPredict(): Unable to invert J_0\n");
			return(false);
		}
		// fill in lower part
		for (i = 0; i < n_0; i++) for (j = i+1; j < n_0; j++) J_0[lsymi(i,j,n_0)] = J_0[usymi(i,j,n_0)];

		// compute A_0 x inv(J_0) x A_0
		double S_0[n_0*n_0];

		for (i = 0; i < n_0; i++) for (j = 0; j < n_0; j++) S_0[i + j*n_0] = 0;
		for (i = 0; i < n_0; i++)
			for (j = 0; j < n_0; j++)
				for (k = 0; k < n_0; k++)
					S_0[i + j*n_0] += J_0[i + k*n_0]*A_0[k + j*n_0];

		for (i = 0; i < n_0; i++) for (j = 0; j < n_0; j++) J_0[i + j*n_0] = 0;
		for (i = 0; i < n_0; i++)
			for (j = 0; j < n_0; j++)
				for (k = 0; k < n_0; k++)
					J_0[i + j*n_0] += A_0[i + k*n_0]*S_0[k + j*n_0];

		// invert J_0
		if (chol2inv(n_0, J_0)) {
			MSG("blockPredict(): Unable to invert J_0\n");
			return(false);
		}

		for (i = 0; i < n_0; i++) {
			sd[i] = sqrt(J_0[i + i*n_0]);
		}
	} // end do_sd

/*
	- For each block we have prediction sites in...
		o We must construct A_0 (n_0 by n_0) and b_0 (n_0 by 1)
		o For prediction SDs, we must construct J_0 (n_0 by n_0) and B_0 (n_0 by n_0+n_block)
		o For each neighboring block (neighbor):
			- Invert Cov(Y_pair) to get inv(Sigma_pair)
			- Add the n_0 block of inv(Sigma_pair) to A_0
			- Add to b_0: inv(Sigma_pair)[pred,block] x resids_block + inv(Sigma_pair)[pred,neighbor] x resids_neighbor
				o Can probably write this as: inv(Sigma_pair)[pred,c(block,neighbor)] x resids[c(block,neighbor)]
			- For prediction SDs:
				o Add to B_0: inv(Sigma_pair)[pred,c(pred,block)]
				o For each neighboring block (neighbor2):
					- Invert Cov(Y_pair2) to get inv(Sigma_pair2)
					- Compute Cov(neighbor,neighbor2)
					- Add to J_0: inv(Sigma_pair)[pred,neighbor] x Sigma_n1n2[neighbor,neighbor2] x inv(Sigma_pair2)[n2,pred]
							J_0 <- J_0 + invSigma.n1[1:n.in.new,n.in.new+n.in.obs+1:n.in.n1] %*%
								Sigma.n1n2[1:n.in.n1,n.in.n1+1:n.in.n2] %*% invSigma.n2[n.in.new+n.in.obs+1:n.in.n2,1:n.in.new]
		o Invert A_0
		o Negate b_0
		o Take y_0[pred] = newX[pred] x beta + inv(A_0) x b_0
		o For prediction SDs:
			- Take H_0 = -A_0
			- Add to J_0: B_0 x Sigma_block x t(B_0)
					Sigma <- compute_cov(object$cov, theta, D[in.b,in.b])
					J_0   <- J_0 + B_0 %*% Sigma %*% t(B_0)
			- For each neighboring block (neighbor):
				o Compute inv(Sigma[c(new,block,neighbor)])
				o Add to J_0: 2 x B_0 x Sigma[block,neighbor] x inv(Sigma)[neighbor,pred]
			- Compute prediction SDs: sqrt(diag( inv(H_0 x inv(J_0) x H_0) ))

*/

	return(true);
}

void BlockComp::getBeta(double *beta) {
	for (int i = 0; i < mNbeta; i++) beta[i] = mBeta[i];
}

void BlockComp::getTheta(double *theta) {
	for (int i = 0; i < mNtheta; i++) theta[i] = mTheta[i];
}

#ifdef CLINE

// test the estimation on the data in test_data.h

#include "test_data.h"

void test_bc(BlockComp::LikForm lf, int nthreads, bool gpu) {
	BlockComp blk(nthreads, gpu);

	//blk.setLikForm(BlockComp::Block);
	//blk.setLikForm(BlockComp::Full);
	blk.setLikForm(lf);
	blk.setCovType(BlockComp::Exp);

	if (!blk.setData(test_n, test_y, test_S, test_nblocks, test_B, test_p, test_X, test_npairs, test_neighbors)) {
		MSG("Error allocating data for fit.\n");
		return;
	}
	//blk.setData(test_n, test_y, test_S, test_nblocks, test_B, test_p, test_X, test_npairs, test_neighbors);

	//double inits[] = {0.25, 0.25, 0.5};
	//double inits[] = {0.5, 0.7243848, 0.2105263};
	bool fixed[] = {false, false, false};
	double vals[] = {0.5, 0.4, 0.15};

	blk.setInits(test_inits);
	blk.setFixed(fixed, vals);

	if (!blk.fit(true)) {
		MSG("Error with fit.\n");
		return;
	}

	//double newS[] = { 0.75, 0.80 };
/*
	double newS[] = { 0.6, 0.85 };
	double newX[] = { 1.0, 0.0 };
	double newY[] = { 0.0 };
	double newSD[] = { 0.0 };
	blk.predict(1, newY, newS, newX, true, newSD);
MSG("predicted y_0: %.2f (%.3f)\n", newY[0], newSD[0]);
*/
	double newS[] = { 0.75, 0.6, 0.80, 0.85 };
	int newB[] = { 1, 1 };
	double newX[] = { 1.0, 1.0, 0.0, 0.0 };
	double newY[] = { 0, 0 };
	double newSD[] = { 0, 0 };

blk.predict(2, newY, newS, newB, newX, false, newSD);
MSG("predicted y_0: %.2f (%.3f); %.2f (%.3f)\n", newY[0], newSD[0], newY[1], newSD[1]);
newY[0]=0;newY[1]=0;
newSD[0]=0;newSD[1]=0;
blk.predict(2, newY, newS, newB, newX, true, newSD);
MSG("predicted y_0: %.2f (%.3f); %.2f (%.3f)\n", newY[0], newSD[0], newY[1], newSD[1]);

}

int main(void) {
	// start blocks at 0
	int i;
	for (i = 0; i < test_n; i++) {
		test_B[i]--;
	}

	for (i = 0; i < test_npairs*2; i++) {
		test_neighbors[i]--;
	}

	clock_t t1;
/*
	MSG("GPU\n");
	t1 = clock();
	test_bc(4, false);
	MSG("--> Done (%.2fsec)\n", (double)(clock() - t1)/CLOCKS_PER_SEC);
	MSG("=========================================================================\n");
	MSG("CPU\n");
*/

/*
	MSG("GPU full\n");
	t1 = clock();
	test_bc(BlockComp::Full, 4, true);
	MSG("--> Done (%.2fsec)\n", (double)(clock() - t1)/CLOCKS_PER_SEC);

	MSG("CPU full\n");
	t1 = clock();
	test_bc(BlockComp::Full, 4, false);
	MSG("--> Done (%.2fsec)\n", (double)(clock() - t1)/CLOCKS_PER_SEC);
*/

	MSG("CPU blocks, no threads\n");
	t1 = clock();
	test_bc(BlockComp::Block, 1, false);
	MSG("--> Done (%.2fsec)\n", (double)(clock() - t1)/CLOCKS_PER_SEC);

	MSG("CPU blocks, 4 threads\n");
	t1 = clock();
	test_bc(BlockComp::Block, 4, false);
	MSG("--> Done (%.2fsec)\n", (double)(clock() - t1)/CLOCKS_PER_SEC);
/*
*/

	return(0);
}

#endif // end main()
