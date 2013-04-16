#include <stdio.h>
#include <string.h>
#include <R.h>
#include "BlockComp.h"
#include "utils.h"

extern "C" {

// spacious_fit function to be called from R
void spacious_fit(double *y, double *X, double *S, int *B, int *neighbors,
                  int *n, int *p, int *nblocks, int *npairs,
                  char **lik_form, char **cov,
                  double *theta, bool *theta_fixed, double *beta,
                  bool *verbose, double *tol, int *max_iter) {

	BlockComp blk;

	// set likelihood form
	if (strcmp(lik_form[0], "block") == 0) {
		blk.setLikForm(BlockComp::Block);
	} else if (strcmp(lik_form[0], "full") == 0) {
		blk.setLikForm(BlockComp::Full);
	} else if (strcmp(lik_form[0], "pair") == 0) {
		blk.setLikForm(BlockComp::Pair);
	} else {
		error("Unknown likelihood form: %s\n", lik_form[0]);
		return;
	}

	// set covariance model
	if (strcmp(cov[0], "exp") == 0) {
		blk.setCovType(BlockComp::Exp);
	} else if (strcmp(cov[0], "matern") == 0) {
		blk.setCovType(BlockComp::Matern);
	} else {
		error("Unknown covariance type: %s\n", cov[0]);
		return;
	}

	// setup data for fit
	blk.setData(n[0], y, S, nblocks[0], B, p[0], X, npairs[0], neighbors);

	// set initial values
	blk.setInits(theta);

	// set fixed parameters
	blk.setFixed(theta_fixed, theta);

	// fit!
	if (!blk.fit(verbose[0])) {
		error("Error with fit.\n");
		return;
	}

	// get estimates
	blk.getBeta(beta);
	blk.getTheta(theta);

	// TODO: set tolerance and max iters

	return;

}

} // end extern "C"
