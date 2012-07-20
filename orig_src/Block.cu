/*
This Code is provied to be freely used, distributed, or modified.
However it comes without warranty of any kind.
Matt Wheeler 2011 

Block class definiation.  Defines a "block" of spatial data consisting
of locations and measured data.  Uses NVIDIA GPUs for most of the
calculations.
*/

#include "Block.h"
#include "util.h"

#include <cutil_inline.h>
#include <cublas.h>

#include <assert.h>

#include "kernels.cu"
#include "cholesky_kernel.cu"

using namespace std;

Block::Block() {
  allocatedOnDevice = false;
  allocatedOnHost = false;
  allocatedDataOnDevice = false;
  allocatedDataOnHost = false;

  size = 0;
}

Block::~Block(){
  if(allocatedOnHost) {
    delete [] h_x;
    delete [] h_y;
  }

  if(allocatedDataOnHost)
    delete [] h_data;

  if(allocatedOnDevice) {
    cutilSafeCall( cudaFree(d_x) );
    cutilSafeCall( cudaFree(d_y) );
  }
}

void Block::setData(DATA_TYPE* _data, const unsigned int& sz) {
  if(allocatedDataOnHost) 
    delete [] h_data;

  size = sz;
  h_data = new DATA_TYPE[size];
  allocatedDataOnHost = true;

  memcpy(h_data,_data,size*sizeof(DATA_TYPE));
}
 
void Block::setXYs(DATA_TYPE *xs, DATA_TYPE *ys, const unsigned int& sz) {
  if(allocatedOnHost) {
    delete [] h_x;
    delete [] h_y;
  }

  size = sz;
  h_x = new DATA_TYPE[size];
  h_y = new DATA_TYPE[size];
  allocatedOnHost = true;

  memcpy(h_x, xs, sz*sizeof(DATA_TYPE));
  memcpy(h_y, ys, sz*sizeof(DATA_TYPE));
}

void Block::randomFill(const unsigned int& s, const DATA_TYPE& xLower, const DATA_TYPE& yLower,
                       const DATA_TYPE& del) {
  mt = MTRand(MT_SEED + (int)xLower*s + (int)yLower);
  size = s;

  h_x = new DATA_TYPE[size];
  h_y = new DATA_TYPE[size];
  allocatedOnHost = true;

  int i;
  for(i = 0; i < size; i++) {
    h_x[i] = xLower*del + del*mt.rand();
    h_y[i] = yLower*del + del*mt.rand();
  }
}

void Block::generateCov(DATA_TYPE *lt, DATA_TYPE* cov) {
  DATA_TYPE sig2 = exp(-lt[0]),
        phi = exp(lt[1]),
        tau2 = exp(-lt[2]),
        *h_ones = new DATA_TYPE[size], *d_ones,
        *d_Cx,
        *d_Cy;


  int i;
  for(i = 0; i < size; i++) h_ones[i] = 1.;
  cutilSafeCall( cudaMalloc((void**)&d_ones,size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemcpy(d_ones, h_ones, size*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );

  cutilSafeCall( cudaMalloc((void**)&d_Cx,size*size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemset(d_Cx, 0, size*size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMalloc((void**)&d_x,size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemcpy(d_x, h_x, size*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );

#ifdef DOUBLE_PRECISION
  cublasDger(size,size,1.,d_x,1,d_ones,1,d_Cx,size);
  cublasDger(size,size,-1.,d_ones,1,d_x,1,d_Cx,size);
#else
  cublasSger(size,size,1.,d_x,1,d_ones,1,d_Cx,size);
  cublasSger(size,size,-1.,d_ones,1,d_x,1,d_Cx,size);
#endif

  cutilSafeCall( cudaMalloc((void**)&d_Cy,size*size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemset(d_Cy, 0, size*size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMalloc((void**)&d_y,size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemcpy(d_y, h_y, size*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
#ifdef DOUBLE_PRECISION
  cublasDger(size,size,1.,d_y,1,d_ones,1,d_Cy,size);
  cublasDger(size,size,-1.,d_ones,1,d_y,1,d_Cy,size);
#else
  cublasSger(size,size,1.,d_y,1,d_ones,1,d_Cy,size);
  cublasSger(size,size,-1.,d_ones,1,d_y,1,d_Cy,size);
#endif
  unsigned int threadx = 32,
               thready = 32;
  dim3 threads(threadx,thready),
       blocks(size/threadx,size/thready);

  cuda_make_cov<<<blocks,threads>>>(d_Cx,d_Cy,sig2,phi,tau2,size);
  cudaThreadSynchronize();
  cutilCheckMsg("cuda_make_cov failed\n");

  cutilSafeCall( cudaMemcpy(cov, d_Cx, size*size*sizeof(DATA_TYPE), cudaMemcpyDeviceToHost) );

  cutilSafeCall( cudaFree(d_ones) );
  cutilSafeCall( cudaFree(d_Cx) );
  cutilSafeCall( cudaFree(d_Cy) );
  cutilSafeCall( cudaFree(d_x) );
  cutilSafeCall( cudaFree(d_y) );
  delete [] h_ones;
}

void Block::generateCovAndDerivatives(DATA_TYPE *lt, DATA_TYPE *cov, DATA_TYPE *der) {
  DATA_TYPE sig2 = exp(-lt[0]),
        phi = exp(lt[1]),
        tau2 = exp(-lt[2]),
        *h_ones = new DATA_TYPE[size], *d_ones,
        *d_Cy,
        *d_der;
  int i;
  for(i = 0; i < size; i++) h_ones[i] = 1.;
  cutilSafeCall( cudaMalloc((void**)&d_ones,size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemcpy(d_ones, h_ones, size*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );

  cutilSafeCall( cudaMemset(cov, 0, size*size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMalloc((void**)&d_x,size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemcpy(d_x, h_x, size*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
#ifdef DOUBLE_PRECISION
  cublasDger(size,size,1.,d_x,1,d_ones,1,cov,size);
  cublasDger(size,size,-1.,d_ones,1,d_x,1,cov,size);
#else
  cublasSger(size,size,1.,d_x,1,d_ones,1,cov,size);
  cublasSger(size,size,-1.,d_ones,1,d_x,1,cov,size);
#endif

  cutilSafeCall( cudaMalloc((void**)&d_Cy,size*size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemset(d_Cy, 0, size*size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMalloc((void**)&d_y,size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemcpy(d_y, h_y, size*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
#ifdef DOUBLE_PRECISION
  cublasDger(size,size,1.,d_y,1,d_ones,1,d_Cy,size);
  cublasDger(size,size,-1.,d_ones,1,d_y,1,d_Cy,size);
#else
  cublasSger(size,size,1.,d_y,1,d_ones,1,d_Cy,size);
  cublasSger(size,size,-1.,d_ones,1,d_y,1,d_Cy,size);
#endif

  cutilSafeCall( cudaMalloc((void**)&d_der,size*size*sizeof(DATA_TYPE)) );
  unsigned int threadx = 32,
               thready = 32;
  dim3 threads(threadx,thready),
       blocks(size/threadx,size/thready);

  cutilSafeCall( cudaMemset(d_der, 0., size*size*sizeof(DATA_TYPE)) );
  cuda_make_ddtau2<<<blocks,threads>>>(cov,d_Cy,d_der,sig2,phi,tau2,size);
  cudaThreadSynchronize();
  cutilCheckMsg("cuda_make_ddtau2 failed\n");
  cutilSafeCall( cudaMemcpy(der, d_der, size*size*sizeof(DATA_TYPE), cudaMemcpyDeviceToHost) );

  cutilSafeCall( cudaMemset(d_der, 0., size*size*sizeof(DATA_TYPE)) );
  cuda_make_ddphi<<<blocks,threads>>>(cov,d_Cy,d_der,sig2,phi,tau2,size);
  cudaThreadSynchronize();
  cutilCheckMsg("cuda_make_ddphi failed\n");
  cutilSafeCall( cudaMemcpy(&der[size*size], d_der, size*size*sizeof(DATA_TYPE), cudaMemcpyDeviceToHost) );

  cutilSafeCall( cudaMemset(d_der, 0., size*size*sizeof(DATA_TYPE)) );
  cuda_make_ddsig2<<<blocks,threads>>>(d_der,tau2,size);
  cudaThreadSynchronize();
  cutilCheckMsg("cuda_make_ddsig2 failed\n");
  cutilSafeCall( cudaMemcpy(&der[2*size*size], d_der, size*size*sizeof(DATA_TYPE), cudaMemcpyDeviceToHost) );

  cuda_make_cov<<<blocks,threads>>>(cov,d_Cy,sig2,phi,tau2,size);
  cudaThreadSynchronize();
  cutilCheckMsg("cuda_make_cov failed\n");

  cudaThreadSynchronize();

  cutilSafeCall( cudaFree(d_ones) );
  cutilSafeCall( cudaFree(d_Cy) );
  cutilSafeCall( cudaFree(d_x) );
  cutilSafeCall( cudaFree(d_y) );
  cutilSafeCall( cudaFree(d_der) );
  delete [] h_ones;
}

DATA_TYPE Block::generateLik(DATA_TYPE *lt) {
  assert(size > 0 && lt && h_data); 
  DATA_TYPE lik = 0.;
  DATA_TYPE *C = new DATA_TYPE[size*size], *d_C,
        *h_xtCx = new DATA_TYPE[1], *d_xtCx,
        *d_Cx;
  cutilSafeCall( cudaMalloc(&d_C, size*size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMalloc(&d_data, size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMalloc(&d_Cx, size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMalloc(&d_xtCx, sizeof(DATA_TYPE)) );
  generateCov(lt,C);

  cholesky_cuda(C, size, BLOCK_N);

  DATA_TYPE logdet = 0;
  int j;
  for(j = 0; j < size; j++)
    logdet += log(C[j*size + j]);
  logdet *= 2;

  cutilSafeCall( cudaMemcpy(d_C, C, size*size*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
  assert(allocatedDataOnHost);
  cutilSafeCall( cudaMemcpy(d_data, h_data, size*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaMemcpy(d_Cx, h_data, size*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
#ifdef DOUBLE_PRECISION
  cublasDtrsm('L','L','N','N',size,1,1.,d_C,size,d_Cx,size);
  cublasDtrsm('L','L','T','N',size,1,1.,d_C,size,d_Cx,size);
  cublasDgemv('T',size,1,1.,d_data,size,d_Cx,1,0,d_xtCx,1);
#else
  cublasStrsm('L','L','N','N',size,1,1.,d_C,size,d_Cx,size);
  cublasStrsm('L','L','T','N',size,1,1.,d_C,size,d_Cx,size);
  cublasSgemv('T',size,1,1.,d_data,size,d_Cx,1,0,d_xtCx,1);
#endif

  cutilSafeCall( cudaMemcpy(h_xtCx, d_xtCx, sizeof(DATA_TYPE), cudaMemcpyDeviceToHost) );

  lik = -.5*(h_xtCx[0] + abs(logdet));

  cutilSafeCall( cudaFree(d_data) );
  cutilSafeCall( cudaFree(d_C) );
  cutilSafeCall( cudaFree(d_Cx) );
  cutilSafeCall( cudaFree(d_xtCx) );
  delete [] C;
  delete [] h_xtCx;
  return lik;
}

DATA_TYPE Block::generateCompositeLik(DATA_TYPE *lt, Block& neighbor) {

  assert(neighbor.getSize() == getSize() && lt);

  DATA_TYPE *xs   = new DATA_TYPE[2*size],
        *ys   = new DATA_TYPE[2*size],
        *data = new DATA_TYPE[2*size];

  memcpy(xs, h_x, size*sizeof(DATA_TYPE));
  memcpy(xs + size, neighbor.getXs(), size*sizeof(DATA_TYPE));

  memcpy(ys, h_y, size*sizeof(DATA_TYPE));
  memcpy(ys + size, neighbor.getYs(), size*sizeof(DATA_TYPE));

  memcpy(data, h_data, size*sizeof(DATA_TYPE));
  memcpy(data + size, neighbor.getData(), size*sizeof(DATA_TYPE));

  Block compositeBlock;
  compositeBlock.setXYs(xs,ys,2*size);
  compositeBlock.setData(data,2*size);

  delete [] xs;
  delete [] ys;
  delete [] data;

  return compositeBlock.generateLik(lt);
}

void Block::addCompositeScoreHessContrib(DATA_TYPE *lt, Block& neighbor,
                                         DATA_TYPE *derL, DATA_TYPE **derE2) {

  assert(neighbor.getSize() == getSize() && lt);

  DATA_TYPE *xs    = new DATA_TYPE[2*size],
        *ys    = new DATA_TYPE[2*size],
        *data  = new DATA_TYPE[2*size],
        *h_der_and_uu_vv = new DATA_TYPE[4*size*size*SP_PARAM_SIZE],
        *d_curr_der_uu_vv,
        *d_C,
        *d_u1,
        *d_y;
  unsigned int compSize;

  memcpy(xs, h_x, size*sizeof(DATA_TYPE));
  memcpy(xs + size, neighbor.getXs(), size*sizeof(DATA_TYPE));

  memcpy(ys, h_y, size*sizeof(DATA_TYPE));
  memcpy(ys + size, neighbor.getYs(), size*sizeof(DATA_TYPE));

  memcpy(data, h_data, size*sizeof(DATA_TYPE));
  memcpy(data + size, neighbor.getData(), size*sizeof(DATA_TYPE));

  compSize = 2*size;
  Block compositeBlock;
  compositeBlock.setXYs(xs,ys,compSize);
  compositeBlock.setData(data,compSize);

  cutilSafeCall( cudaMalloc(&d_C, compSize*compSize*sizeof(DATA_TYPE)) );
  compositeBlock.generateCovAndDerivatives(lt,d_C,h_der_and_uu_vv);

  cutilSafeCall( cudaMalloc(&d_u1, compSize*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemcpy(d_u1, data, compSize*sizeof(DATA_TYPE), 
                            cudaMemcpyHostToDevice) );

  leftLinCholSovle(d_C,d_u1,compSize,1,/*needsChol*/true);

  cutilSafeCall( cudaMalloc(&d_y, compSize*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMalloc(&d_curr_der_uu_vv, 2*compSize*compSize*sizeof(DATA_TYPE)) );
  int kl, lm;
  for(kl = 0; kl < SP_PARAM_SIZE; kl++) {
    cutilSafeCall( cudaMemcpy(&d_curr_der_uu_vv[0], &h_der_and_uu_vv[kl*compSize*compSize], 
                              compSize*compSize*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
#ifdef DOUBLE_PRECISION
    cublasDgemv('T',compSize,compSize,1.,&d_curr_der_uu_vv[0],compSize,d_u1,1,0.,d_y,1);
    derL[kl] += .5*cublasDdot(compSize,d_y,1,d_u1,1);
#else
    cublasSgemv('T',compSize,compSize,1.,&d_curr_der_uu_vv[0],compSize,d_u1,1,0.,d_y,1);
    derL[kl] += .5*cublasSdot(compSize,d_y,1,d_u1,1);
#endif
  }
 
  for(kl = 0; kl < SP_PARAM_SIZE; kl++) {
    cutilSafeCall( cudaMemcpy(&d_curr_der_uu_vv[0], &h_der_and_uu_vv[kl*compSize*compSize], 
                              compSize*compSize*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
    leftLinCholSovle(d_C,&d_curr_der_uu_vv[0],compSize,compSize,/*needsChol*/false);
    cutilSafeCall( cudaMemcpy(&h_der_and_uu_vv[kl*compSize*compSize], &d_curr_der_uu_vv[0], 
                              compSize*compSize*sizeof(DATA_TYPE), cudaMemcpyDeviceToHost) );
  }
  cutilSafeCall( cudaFree(d_C) );

  for(kl = 0; kl < SP_PARAM_SIZE; kl++) {
    cutilSafeCall( cudaMemcpy(d_curr_der_uu_vv, &h_der_and_uu_vv[kl*compSize*compSize], 
                              compSize*compSize*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
    derL[kl] -= .5*trace(d_curr_der_uu_vv,compSize);
    derE2[kl][kl] += .5*trace(d_curr_der_uu_vv,d_curr_der_uu_vv,compSize);

    for(lm = kl + 1; lm < SP_PARAM_SIZE; lm++) {
      cutilSafeCall( cudaMemcpy(&d_curr_der_uu_vv[compSize*compSize], 
                                &h_der_and_uu_vv[lm*compSize*compSize], 
                                compSize*compSize*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
      derE2[kl][lm] += .5*trace(d_curr_der_uu_vv,&d_curr_der_uu_vv[compSize*compSize],compSize);
    }
  }

  cutilSafeCall( cudaFree(d_curr_der_uu_vv) );
  cutilSafeCall( cudaFree(d_u1) );
  cutilSafeCall( cudaFree(d_y) );
  delete [] xs;
  delete [] ys;
  delete [] data;
  delete [] h_der_and_uu_vv;
}

void Block::makeScoreHess(DATA_TYPE *lt, DATA_TYPE *derL, DATA_TYPE **derE2) {

  DATA_TYPE 
        *h_der_and_uu_vv = new DATA_TYPE[size*size*SP_PARAM_SIZE],
        *d_curr_der_uu_vv,
        *d_C,
        *d_u1,
        *d_y;
  unsigned int compSize;

  compSize = size;

  cutilSafeCall( cudaMalloc(&d_C, compSize*compSize*sizeof(DATA_TYPE)) );
  generateCovAndDerivatives(lt,d_C,h_der_and_uu_vv);

  cutilSafeCall( cudaMalloc(&d_u1, compSize*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemcpy(d_u1, h_data, compSize*sizeof(DATA_TYPE), 
                            cudaMemcpyHostToDevice) );

  leftLinCholSovle(d_C,d_u1,compSize,1,/*needsChol*/true);

  cutilSafeCall( cudaMalloc(&d_y, compSize*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMalloc(&d_curr_der_uu_vv, 2*compSize*compSize*sizeof(DATA_TYPE)) );
  int kl, lm;
  for(kl = 0; kl < SP_PARAM_SIZE; kl++) {
    cutilSafeCall( cudaMemcpy(&d_curr_der_uu_vv[0], &h_der_and_uu_vv[kl*compSize*compSize], 
                              compSize*compSize*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
#ifdef DOUBLE_PRECISION
    cublasDgemv('T',compSize,compSize,1.,&d_curr_der_uu_vv[0],compSize,d_u1,1,0.,d_y,1);
    derL[kl] += .5*cublasDdot(compSize,d_y,1,d_u1,1);
#else
    cublasSgemv('T',compSize,compSize,1.,&d_curr_der_uu_vv[0],compSize,d_u1,1,0.,d_y,1);
    derL[kl] += .5*cublasSdot(compSize,d_y,1,d_u1,1);
#endif
  }
 
  for(kl = 0; kl < SP_PARAM_SIZE; kl++) {
    cutilSafeCall( cudaMemcpy(&d_curr_der_uu_vv[0], &h_der_and_uu_vv[kl*compSize*compSize], 
                              compSize*compSize*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
    leftLinCholSovle(d_C,&d_curr_der_uu_vv[0],compSize,compSize,/*needsChol*/false);
    cutilSafeCall( cudaMemcpy(&h_der_and_uu_vv[kl*compSize*compSize], &d_curr_der_uu_vv[0], 
                              compSize*compSize*sizeof(DATA_TYPE), cudaMemcpyDeviceToHost) );
  }
  cutilSafeCall( cudaFree(d_C) );

  for(kl = 0; kl < SP_PARAM_SIZE; kl++) {
    cutilSafeCall( cudaMemcpy(d_curr_der_uu_vv, &h_der_and_uu_vv[kl*compSize*compSize], 
                              compSize*compSize*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
    derL[kl] -= .5*trace(d_curr_der_uu_vv,compSize);
    derE2[kl][kl] += .5*trace(d_curr_der_uu_vv,d_curr_der_uu_vv,compSize);

    for(lm = kl + 1; lm < SP_PARAM_SIZE; lm++) {
      cutilSafeCall( cudaMemcpy(&d_curr_der_uu_vv[compSize*compSize], 
                                &h_der_and_uu_vv[lm*compSize*compSize], 
                                compSize*compSize*sizeof(DATA_TYPE), cudaMemcpyHostToDevice) );
      derE2[kl][lm] += .5*trace(d_curr_der_uu_vv,&d_curr_der_uu_vv[compSize*compSize],compSize);
    }
  }

  cutilSafeCall( cudaFree(d_curr_der_uu_vv) );
  cutilSafeCall( cudaFree(d_u1) );
  cutilSafeCall( cudaFree(d_y) );
  delete [] h_der_and_uu_vv;
}

void Block::generateData(DATA_TYPE *lt, const DATA_TYPE& xLower, const DATA_TYPE& yLower) {
  assert(size > 0 && h_x && h_y);

  mt = MTRand(MT_DB_SEED + (int)xLower*size + (int)yLower);
  DATA_TYPE *h_cov = new DATA_TYPE[size*size], *d_cov,
        *h_randVec = new DATA_TYPE[size];
  h_data = new DATA_TYPE[size];
  allocatedDataOnHost = true;
  generateCov(lt,h_cov);

  cholesky_cuda(h_cov,size,BLOCK_N);

  int i;
  DATA_TYPE x1, x2, w;
  for(i = 0; i < size; i+=2) {
    do {
      x1 = 2.*mt.rand() - 1.;
      x2 = 2.*mt.rand() - 1.;
      w = x1*x1 + x2*x2;
    } while (w >= 1.);
    w = sqrt( ( -2.*log(w))/w);
    
    h_randVec[i] = x1*w;
    h_randVec[i+1] = x2*w;
  } 

  cutilSafeCall( cudaMalloc((void**)&d_data,size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemcpy(d_data, h_randVec, size*sizeof(DATA_TYPE), 
                 cudaMemcpyHostToDevice) );

  cutilSafeCall( cudaMalloc((void**)&d_cov,size*size*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMemcpy(d_cov, h_cov, size*size*sizeof(DATA_TYPE), 
                 cudaMemcpyHostToDevice) );

#ifdef DOUBLE_PRECISION
  cublasDtrmv('L','N','N',size,d_cov,size,d_data,1);
#else
  cublasStrmv('L','N','N',size,d_cov,size,d_data,1);
#endif
  cutilSafeCall( cudaMemcpy(h_data, d_data, size*sizeof(DATA_TYPE), 
                 cudaMemcpyDeviceToHost) );

  cutilSafeCall( cudaFree(d_cov) );
  cutilSafeCall( cudaFree(d_data) );
  delete [] h_cov;
  delete [] h_randVec;
}


