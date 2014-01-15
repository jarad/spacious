#include <cuda_runtime_api.h>
#include <cublas_v2.h>

#include <R.h>

#include "utils.h"
#include "utils_cuda.h"

void cuda_log_det(double *v, double *A, int n) {
	int b = ceil((double)n/CUDA_BS);
//printf("b=%d\n", b);
	cuda_log_det_k<<<ceil((double)n/CUDA_BS), CUDA_BS>>>(v, A, n);
}

__global__ void cuda_log_det_k(double *v, double *A, int n) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n) v[i] = log(A[i + i*n]);
}
