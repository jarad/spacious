#ifdef CUDA
// Code to accelerate matrix operations on GPUs.
// Ideas taken from Matt Wheeler, netlib.org, and MAGMA.

#include <R.h>
#include <R_ext/Lapack.h>

#include <cuda_runtime_api.h>
#include <cublas_v2.h>

#include "utils.h"
#include "utils_cuda.h"

void cuda_devices() {
	int ndevices;
	int i;
	cudaDeviceProp prop;

	// get number of devices
	cudaGetDeviceCount(&ndevices);

	for (i = 0; i < ndevices; i++) {
		cudaGetDeviceProperties(&prop, i);

		MSG("%d: %s, %.1f MHz clock, %.1f MB memory, compute capability %d.%d\n",
		    i, prop.name, prop.clockRate/1000.0, prop.totalGlobalMem/(1024.0*1024.0), prop.major, prop.minor);
	}
}

// cuda_chol2inv: compute inverse of n by n pd symm mat A using Cholesky
// Note: A must be allocated on device
int cuda_chol2inv(cublasHandle_t handle, int n, double *A) {
	int  info;

	// compute factorization
	cuda_dpotrf(handle, n, A, &info);
	if (info) {
		return(info);
	}

	// complete inverse
	cuda_dpotri(handle, n, A, &info);
	if (info) {
		return(info);
	}

	return(0);
}

// compute Cholesky factorization (upper triangle) of A
void cuda_dpotrf(cublasHandle_t handle, int n, double *A, int *info) {
	char           uplo = 'U';
	double         n1 = -1.0, p1=1.0;
	double        *work;
	int            i,ibs;
	int            bs = cuda_block_size();
	cudaStream_t   streams[2];
	cublasStatus_t status;

	if (n < bs) {
		// small n; use CPU
		MSG("TODO: cuda_dpotrf(): handle small n case\n");
		*info = -1;
		return;
	}

	// create streams
	for (i = 0; i < 2; i++) {
		if (cudaStreamCreate(&streams[i]) != cudaSuccess) {
			MSG("cuda_dpotrf(): unable to create cuda stream %d: %s\n", i, cudaGetErrorString(cudaGetLastError()));
			*info = -1;
			return;
		}
	}

	// allocate workspace
	if (cudaMallocHost((void **)&work, bs*bs*sizeof(double)) != cudaSuccess) {
		MSG("cuda_dpotrf(): unable to allocate workspace: %s\n", cudaGetErrorString(cudaGetLastError()));
		*info = -1;
		return;
	}

	for (i = 0; i < n; i += bs) {
		ibs = imin(bs, n-i);

		// perform update for this block
		status = cublasDsyrk(handle, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_T, ibs, i,
		                     &n1, &A[fsymi(0,i,n)], n, &p1, &A[fsymi(i,i,n)], n);
		if (status != CUBLAS_STATUS_SUCCESS) {
			MSG("cuda_dpotrf(): unable to call cublasDsyrk(): %d\n", status);
			*info = -1;
			break;
		}

		status = cublasGetMatrixAsync(ibs, ibs, sizeof(double), &A[fsymi(i,i,n)], n, work, ibs, streams[1]);
		if (status != CUBLAS_STATUS_SUCCESS) {
			MSG("cuda_dpotrf(): unable to call cublasGetMatrixAsync(): %d\n", status);
			*info = -1;
			break;
		}

		if ((i+ibs) < n) {
			// operate on row for this block
			status = cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, ibs, n-i-ibs, i,
			                     &n1, &A[fsymi(0,i,n)], n, &A[fsymi(0,i+ibs,n)], n, &p1, &A[fsymi(i,i+ibs,n)], n);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("cuda_dpotrf(): unable to call cublasDgemm(): %d\n", status);
				*info = -1;
				break;
			}
		}

		if (cudaStreamSynchronize(streams[1]) != cudaSuccess) {
			MSG("cuda_dpotrf(): unable to synchronize cuda stream: %s\n", cudaGetErrorString(cudaGetLastError()));
			*info = -1;
			break;
		}

		// factorize workspace
		dpotrf_(&uplo, &ibs, work, &ibs, info);
		if (*info != 0) {
			// unable to factorize workspace
			MSG("cuda_dpotrf(): unable to factorize workspace: %d, %d\n", *info, i);
			*info = *info + i;
			break;
		}

		status = cublasSetMatrixAsync(ibs, ibs, sizeof(double), work, ibs, &A[fsymi(i,i,n)], n, streams[0]);
		if (status != CUBLAS_STATUS_SUCCESS) {
			MSG("cuda_dpotrf(): unable to call cublasSetMatrixAsync(): %d\n", status);
			*info = -1;
			break;
		}

		if ((i+ibs) < n) {
			// run the triangular solve
			status = cublasDtrsm(handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT,
				                     ibs, n-i-ibs, &p1, &A[fsymi(i,i,n)], n, &A[fsymi(i,i+ibs,n)], n);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("cuda_dpotrf(): unable to call cublasDtrsm(): %d\n", status);
				*info = -1;
				break;
			}

		}

	}

	// destroy streams
	for (i = 0; i < 2; i++) {
		cudaStreamDestroy(streams[i]);
	}

	// free workspace
	cudaFreeHost(work);
}

// compute inverse from Cholesky factorization (upper triangle) of A
void cuda_dpotri(cublasHandle_t handle, int n, double *A, int *info) {
	// invert upper triangle
	cuda_dtrtri(handle, CUBLAS_DIAG_NON_UNIT, n, A, n, info);

	if (*info == 0) {
		// complete inverse with inv(U)inv(U)'
		cuda_dlauum(handle, n, A, n, info);
	}
}

// invert upper triangle of A
void cuda_dtrtri(cublasHandle_t handle, cublasDiagType_t diag, int n, double *A, int lda, int *info) {
	char           uplo = 'U';
	char           c_diag = 'N';
	double         n1 = -1.0, p1=1.0;
	double        *work;
	int            i,ibs;
	int            bs = cuda_block_size();
	cudaStream_t   streams[2];
	cublasStatus_t status;

	if (n < bs) {
		// small n; use CPU
		MSG("TODO: cuda_dtrtri(): handle small n case\n");
		*info = -1;
		return;
	}

	// create streams
	for (i = 0; i < 2; i++) {
		if (cudaStreamCreate(&streams[i]) != cudaSuccess) {
			MSG("cuda_dpotrf(): unable to create cuda stream %d: %s\n", i, cudaGetErrorString(cudaGetLastError()));
			*info = -1;
			return;
		}
	}

	// allocate workspace
	if (cudaMallocHost((void **)&work, bs*bs*sizeof(double)) != cudaSuccess) {
		MSG("cuda_dtrtri(): unable to allocate workspace: %s\n", cudaGetErrorString(cudaGetLastError()));
		*info = -1;
		return;
	}

	for (i = 0; i < n; i += bs) {
		ibs = imin(bs, n-i);

		status = cublasDtrmm(handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT,
		                     i, ibs, &p1, &A[0], lda, &A[fsymi(0,i,n)], lda, &A[fsymi(0,i,n)], lda);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("cuda_dtrtri(): unable to call cublasDtrmm(): %d\n", status);
				*info = -1;
				break;
			}

		status = cublasDtrsm(handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT,
		                     i, ibs, &n1, &A[fsymi(i,i,n)], lda, &A[fsymi(0,i,n)], lda);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("cuda_dtrtri(): unable to call cublasDtrsm(): %d\n", status);
				*info = -1;
				break;
			}

		status = cublasGetMatrixAsync(ibs, ibs, sizeof(double), &A[fsymi(i,i,n)], lda, work, ibs, streams[1]);
		if (status != CUBLAS_STATUS_SUCCESS) {
			MSG("cuda_dtrtri(): unable to call cublasGetMatrixAsync(): %d\n", status);
			*info = -1;
			break;
		}

		if (cudaStreamSynchronize(streams[1]) != cudaSuccess) {
			MSG("cuda_dtrtri(): unable to synchronize cuda stream: %s\n", cudaGetErrorString(cudaGetLastError()));
			*info = -1;
			break;
		}

		// invert workspace
		dtrtri_(&uplo, &c_diag, &ibs, work, &ibs, info);
		if (*info != 0) {
			// unable to factorize workspace
			MSG("cuda_dtrtri(): unable to invert workspace: %d, %d\n", *info, i);
			break;
		}

		status = cublasSetMatrixAsync(ibs, ibs, sizeof(double), work, ibs, &A[fsymi(i,i,n)], n, streams[0]);
		if (status != CUBLAS_STATUS_SUCCESS) {
			MSG("cuda_dtrtri(): unable to call cublasSetMatrixAsync(): %d\n", status);
			*info = -1;
			break;
		}

	}

	// destroy streams
	for (i = 0; i < 2; i++) {
		cudaStreamDestroy(streams[i]);
	}

	// free workspace
	cudaFreeHost(work);
}

// compute product of upper triangle of A with it's transpose
void cuda_dlauum(cublasHandle_t handle, int n, double *A, int lda, int *info) {
	char           uplo = 'U';
	double         p1=1.0;
	double        *work;
	int            i,ibs;
	int            bs = cuda_block_size();
	cublasStatus_t status;

	if (n < bs) {
		// small n; use CPU
		MSG("TODO: cuda_dlauum(): handle small n case\n");
		*info = -1;
		return;
	}

	// allocate workspace
	if (cudaMallocHost((void **)&work, bs*bs*sizeof(double)) != cudaSuccess) {
		MSG("cuda_dlauum(): unable to allocate workspace: %s\n", cudaGetErrorString(cudaGetLastError()));
		*info = -1;
		return;
	}

	for (i = 0; i < n; i += bs) {
		ibs = imin(bs, n-i);

		status = cublasDtrmm(handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT,
		                     i, ibs, &p1, &A[fsymi(i,i,n)], lda, &A[fsymi(0,i,n)], lda, &A[fsymi(0,i,n)], lda);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("cuda_dlauum(): unable to call cublasDtrmm(): %d\n", status);
				*info = -1;
				break;
			}

		status = cublasGetMatrix(ibs, ibs, sizeof(double), &A[fsymi(i,i,n)], lda, work, ibs);
		if (status != CUBLAS_STATUS_SUCCESS) {
			MSG("cuda_dlauum(): unable to call cublasGetMatrix(): %d\n", status);
			*info = -1;
			break;
		}

		// invert workspace
		dlauum_(&uplo, &ibs, work, &ibs, info);
		if (*info != 0) {
			// unable to factorize workspace
			MSG("cuda_dlauum(): unable to invert workspace: %d, %d\n", *info, i);
			break;
		}

		status = cublasSetMatrix(ibs, ibs, sizeof(double), work, ibs, &A[fsymi(i,i,n)], lda);
		if (status != CUBLAS_STATUS_SUCCESS) {
			MSG("cuda_dlauum(): unable to call cublasSetMatrix(): %d\n", status);
			*info = -1;
			break;
		}

		if ((i+ibs) < n) {
			status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, i, ibs, n-i-ibs,
			                     &p1, &A[fsymi(0,i+ibs,n)], lda, &A[fsymi(i,i+ibs,n)], lda, &p1, &A[fsymi(0,i,n)], lda);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("cuda_dlauum(): unable to call cublasDgemm(): %d\n", status);
				*info = -1;
				break;
			}

			status = cublasDsyrk(handle, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, ibs, n-i-ibs,
		                     &p1, &A[fsymi(i,i+ibs,n)], lda, &p1, &A[fsymi(i,i,n)], lda);
			if (status != CUBLAS_STATUS_SUCCESS) {
				MSG("cuda_dlauum(): unable to call cublasDsyrk(): %d\n", status);
				*info = -1;
				break;
			}

		}
	}

	// free workspace
	cudaFreeHost(work);
}

#endif
