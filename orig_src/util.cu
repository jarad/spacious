/*
This Code is provied to be freely used, distributed, or modified.
However it comes without warranty of any kind.
Matt Wheeler 2011 

Convenience functions to solve a linear system Ax = B and calculate
the trace of a matrix.  All functions use the GPU.
*/

#include <cublas.h>
#include <cutil_inline.h>

#include "util.h"
#include "Block.h"

void cholesky_cuda(DATA_TYPE*,int,int );

//Uses Cholesky Decomposition to solve d_C{x} = d_B (or x = d_C{^-1}d_B)
//Overwrites both d_C and d_B so watch out!
void leftLinCholSovle(DATA_TYPE *d_C, DATA_TYPE *d_B, const int& size, const int& cols,
                      const bool& needsChol) {

  if(needsChol) {
    DATA_TYPE *h_C = new DATA_TYPE[size*size];
    cutilSafeCall( cudaMemcpy(h_C, d_C, size*size*sizeof(DATA_TYPE), cudaMemcpyDeviceToHost) );
    cholesky_cuda(h_C, size, BLOCK_N);
    cutilSafeCall( cudaMemcpy(d_C, h_C, size*size*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
    delete [] h_C;
  }

#ifdef DOUBLE_PRECISION
  cublasDtrsm('L','L','N','N',size,cols,1.,d_C,size,d_B,size);
  cublasDtrsm('L','L','T','N',size,cols,1.,d_C,size,d_B,size);
#else
  cublasStrsm('L','L','N','N',size,cols,1.,d_C,size,d_B,size);
  cublasStrsm('L','L','T','N',size,cols,1.,d_C,size,d_B,size);
#endif

}

//Calculates the trace of size-by-sze matrix
DATA_TYPE trace(DATA_TYPE *mat, const int& size) {
  DATA_TYPE *h_ones = new DATA_TYPE[size], *d_ones,
        ans;
  int i;
  cutilSafeCall( cudaMalloc(&d_ones, size*sizeof(DATA_TYPE)) );

  for(i = 0; i < size; i++)
    h_ones[i] = 1.; 

  cutilSafeCall( cudaMemcpy(d_ones, h_ones, size*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );

#ifdef DOUBLE_PRECISION
  ans = cublasDdot(size,mat,size+1,d_ones,1);
#else
  ans = cublasSdot(size,mat,size+1,d_ones,1);
#endif

  cutilSafeCall( cudaFree(d_ones) );
  delete [] h_ones;

  return ans;
}

//Calculates the trace of the product of two size-by-size matrices
DATA_TYPE trace(DATA_TYPE *matA, DATA_TYPE *matB, const int& size) {
  DATA_TYPE  *mat,
        ans = 0.;

  //Which method is faster depends on the number of cols/rows
  if(size <= 2048) { 
    cutilSafeCall( cudaMalloc(&mat, size*size*sizeof(DATA_TYPE)) );
#ifdef DOUBLE_PRECISION
    cublasDgemm('N','N',size,size,size,1.,matA,size,matB,size,0,mat,size);
#else
    cublasSgemm('N','N',size,size,size,1.,matA,size,matB,size,0,mat,size);
#endif

    ans = trace(mat,size);
  
    cutilSafeCall( cudaFree(mat) );
  } else {
    int i;
    for(i = 0; i < size; i++)
#ifdef DOUBLE_PRECISION
      ans += cublasDdot(size,&matA[i],size,&matB[i*size],1);
#else
      ans += cublasSdot(size,&matA[i],size,&matB[i*size],1);
#endif
  }

  return ans;
}

