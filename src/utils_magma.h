#ifndef UTILS_MAGMA_H
#define UTILS_MAGMA_H

// function definitions
int magma_chol2inv_gpu(int n, double *A);
void magma_disp_sym(double *mat, int n, int lo, int hi);

#endif
