/*
This Code is provied to be freely used, distributed, or modified.
However it comes without warranty of any kind.
Matt Wheeler 2011 

Functions to calculate the covariance and derivatives of the covariance
for spatially correlated data goverend by Matern spatial correlation.
*/

#ifndef KERNELS
#define KERNELS

#define CUDA_MAX_X 32
//Change depending on compute capacity:
// 1.x: 16
// 2.x: 32
#define CUDA_MAX_Y 32


//Calculates the covariance matrix.  Overwrites _Cx with the calculated values
__global__ void cuda_make_cov(DATA_TYPE *_Cx, DATA_TYPE *_Cy, DATA_TYPE sig2, 
                              DATA_TYPE phi, DATA_TYPE tau2, unsigned int ld) {

  int x = blockDim.x*blockIdx.x + threadIdx.x;
  int y = blockDim.y*blockIdx.y + threadIdx.y;

  DATA_TYPE covVal = 0.;

  __shared__ volatile DATA_TYPE Cx[CUDA_MAX_X][CUDA_MAX_Y];
  __shared__ volatile DATA_TYPE Cy[CUDA_MAX_X][CUDA_MAX_Y];


  Cx[threadIdx.x][threadIdx.y] = _Cx[y*ld + x];
  Cy[threadIdx.x][threadIdx.y] = _Cy[y*ld + x];

  covVal = sqrt(pow(Cx[threadIdx.x][threadIdx.y],2) + 
                pow(Cy[threadIdx.x][threadIdx.y],2))*phi;
  covVal = sig2*(1+covVal) * exp(-covVal);
  if(x == y) covVal += tau2;

  _Cx[y*ld + x] = covVal;
}

//Calcuates the derivative with respect to the nugget parameter.
//Stores value in _der. 
__global__ void cuda_make_ddtau2(DATA_TYPE *_Cx, DATA_TYPE *_Cy, DATA_TYPE *_der, DATA_TYPE sig2,
                                      DATA_TYPE phi, DATA_TYPE tau2, unsigned int ld) {

  int x = blockDim.x*blockIdx.x + threadIdx.x;
  int y = blockDim.y*blockIdx.y + threadIdx.y;

  DATA_TYPE covVal;

  __shared__ DATA_TYPE Cx[CUDA_MAX_X][CUDA_MAX_Y];
  __shared__ DATA_TYPE Cy[CUDA_MAX_X][CUDA_MAX_Y];

  Cx[threadIdx.x][threadIdx.y] = _Cx[y*ld + x];
  Cy[threadIdx.x][threadIdx.y] = _Cy[y*ld + x];

  covVal = sqrt(pow(Cx[threadIdx.x][threadIdx.y],2) + 
                pow(Cy[threadIdx.x][threadIdx.y],2))*phi;

  covVal = sig2*(1+covVal) * exp(-covVal);
 
  _der[y*ld + x] = -covVal; 
}

//Calcuates the derivative with respect to the phi parameter.
//Stores value in _der. 
__global__ void cuda_make_ddphi(DATA_TYPE *_Cx, DATA_TYPE *_Cy, DATA_TYPE *_der, DATA_TYPE sig2,
                                      DATA_TYPE phi, DATA_TYPE tau2, unsigned int ld) {

  int x = blockDim.x*blockIdx.x + threadIdx.x;
  int y = blockDim.y*blockIdx.y + threadIdx.y;

  DATA_TYPE covVal;

  __shared__ DATA_TYPE Cx[CUDA_MAX_X][CUDA_MAX_Y];
  __shared__ DATA_TYPE Cy[CUDA_MAX_X][CUDA_MAX_Y];

  Cx[threadIdx.x][threadIdx.y] = _Cx[y*ld + x];
  Cy[threadIdx.x][threadIdx.y] = _Cy[y*ld + x];

  covVal = sqrt(pow(Cx[threadIdx.x][threadIdx.y],2) + 
                pow(Cy[threadIdx.x][threadIdx.y],2))*phi;

  _der[y*ld + x] = -sig2*(exp(-covVal)*covVal*covVal);
}

//Calcuates the derivative with respect to the sigma^2 parameter.
//Stores value in _der. 
__global__ void cuda_make_ddsig2(DATA_TYPE *_der, DATA_TYPE tau2, unsigned int ld) {

  int x = blockDim.x*blockIdx.x + threadIdx.x;
  int y = blockDim.y*blockIdx.y + threadIdx.y;

  if(x == y) {
    _der[y*ld + x] = -tau2;
  }
}

#endif
