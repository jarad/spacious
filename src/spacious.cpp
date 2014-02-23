#include <stdio.h>
#include <string.h>
#include <R.h>
#include "BlockComp.h"
#include "utils.h"

extern "C" {

// spacious_fit function to be called from R
void spacious_fit(
	// input data
	double *y, double *X, double *S, int *B, int *neighbors,
	int *n, int *p, int *nblocks, int *npairs,
	// type of model to fit
	char **lik_form, char **cov,
	// parameter estimates and convergence info
	double *theta, int *theta_fixed, double *beta,
	bool *convergence, int *nIter,
	// standard errors
	double *se_beta, double *se_theta,
	double *vcov_beta, double *vcov_theta,
	// fitted values and residuals
	double *fitted, double *resids,
	// values of theta and log likelihood at each iteration
	double *iters_theta, double *iters_ll,
	// fitting control parameters
	bool *verbose, double *tol, int *max_iter, bool *compute_se,
	// parallelization options
	int *nthreads, bool *gpu
) {

	BlockComp *blk;

	blk = new BlockComp(nthreads[0], gpu[0]);

	// set likelihood form
	if (strcmp(lik_form[0], "block") == 0) {
		blk->setLikForm(BlockComp::Block);
	} else if (strcmp(lik_form[0], "full") == 0) {
		blk->setLikForm(BlockComp::Full);
	} else if (strcmp(lik_form[0], "pair") == 0) {
		blk->setLikForm(BlockComp::Pair);
	} else {
		error("Unknown likelihood form: %s\n", lik_form[0]);
		delete blk;
		return;
	}

	// set covariance model
	if (strcmp(cov[0], "exp") == 0) {
		blk->setCovType(BlockComp::Exp);
	} else if (strcmp(cov[0], "matern") == 0) {
		blk->setCovType(BlockComp::Matern);
	} else {
		error("Unknown covariance type: %s\n", cov[0]);
		delete blk;
		return;
	}

	// setup data for fit
	if (!blk->setData(*n, y, S, *nblocks, B, *p, X, *npairs, neighbors)) {
		error("Error with fit: unable to set data\n");
		delete blk;
		return;
	}

	// set initial values
	blk->setInits(theta);

	// set fixed parameters
	bool isFixed[blk->getNumTheta()];
	for (int i = 0; i < blk->getNumTheta(); i++) {
		if (theta_fixed[i] == 1)
			isFixed[i] = true;
		else
			isFixed[i] = false;
	}
	blk->setFixed(isFixed, theta);

	// set tolerance
	blk->setTol(*tol);

	// set max iters
	blk->setMaxIter(*max_iter);

	// fit!
	if (!blk->fit(verbose[0])) {
		error("Error with fit.\n");
		delete blk;
		return;
	}

	// get estimates
	blk->getBeta(beta);
	blk->getTheta(theta);

	//blk->getBetaIter(iterBeta);
	blk->getThetaIter(iters_theta);
	blk->getLogLikIter(iters_ll);

	// get convergence info
	convergence[0] = blk->getConverged();
	nIter[0]       = blk->getIters();

	// get fitted and residuals
	blk->getFitted(fitted);
	blk->getResiduals(resids);

	// get standard errors
	blk->getStdErrs(se_beta, vcov_beta, se_theta, vcov_theta);

	delete blk;
	return;
}

// spacious_predict function to be called from R
void spacious_predict(
	// new site information
	int *n_0, double *y_0, double *newS, int *newB, double *newX,
	bool *do_sd, double *sd, bool *local, int *Nlocal,
	// data used in fit
	double *y, double *X, double *S, int *B, int *neighbors,
	int *n, int *p, int *nblocks, int *npairs,
	// covariance type
	char **lik_form, char **cov,
	// fitted values
	double *beta, double *theta,
	// parallelization options
	int *nthreads, bool *gpu
) {

	BlockComp *blk;

	blk = new BlockComp(nthreads[0], gpu[0]);

	// set likelihood form
	if (strcmp(lik_form[0], "block") == 0) {
		blk->setLikForm(BlockComp::Block);
	} else if (strcmp(lik_form[0], "full") == 0) {
		blk->setLikForm(BlockComp::Full);
	} else if (strcmp(lik_form[0], "pair") == 0) {
		blk->setLikForm(BlockComp::Pair);
	} else {
		error("Unknown likelihood form: %s\n", lik_form[0]);

		delete blk;
		return;
	}

	// set covariance model
	if (strcmp(cov[0], "exp") == 0) {
		blk->setCovType(BlockComp::Exp);
	} else if (strcmp(cov[0], "matern") == 0) {
		blk->setCovType(BlockComp::Matern);
	} else {
		error("Unknown covariance type: %s\n", cov[0]);
		delete blk;
		return;
	}

	// set data
	if (!blk->setData(*n, y, S, *nblocks, B, *p, X, *npairs, neighbors)) {
		error("Error with predict: unable to set data\n");
		delete blk;
		return;
	}

	// set fit
	if (!blk->setFit(beta, theta)) {
		error("Error with predict: unable to set fit\n");
		delete blk;
		return;
	}

	// get predictions
	if (!blk->predict(*n_0, y_0, newS, newB, newX, *do_sd, sd, *local, *Nlocal)) {
		error("Error with predict.\n");
		delete blk;
		return;
	}

	delete blk;
	return;
}

} // end extern "C"
