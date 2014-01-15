#ifndef UTILS_CUDA_H
#define UTILS_CUDA_H

#ifndef CUDA_BS
#define CUDA_BS 256
#endif

// function definitions
void cuda_devices();
int  cuda_chol2inv(cublasHandle_t handle, int n, double *A, bool do_log_det=false, double *log_det=NULL);
void cuda_dpotrf(cublasHandle_t handle, int n, double *A, int *info);
void cuda_dpotri(cublasHandle_t handle, int n, double *A, int *info);
void cuda_dtrtri(cublasHandle_t handle, cublasDiagType_t diag, int n, double *A, int lda, int *info);
void cuda_dlauum(cublasHandle_t handle, int n, double *A, int lda, int *info);

inline int cuda_block_size() {
	return(CUDA_BS);
}

// kernels
void cuda_log_det(double *v, double *A, int n);
__global__ void cuda_log_det_k(double *v, double *A, int n);

#endif
