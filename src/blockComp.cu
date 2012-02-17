/*
This Code is provied to be freely used, distributed, or modified.
However it comes without warranty of any kind.
Matt Wheeler 2011 

This is CUDA code used to evaluate the time required to compute the block
composite likelihood and mle of a 2D grid of spatially correlated points
governed by a Matern covariance structure.  The paper describing this is:
  Estimation and prediction in spatial models with block 
  composite likelihoods using parallel computing
  Jo Eidsvik, Benjamin A. Shaby, Brian J. Reich, Matthew Wheeler, Jarad Niemi
*/

#include <cublas.h>
#include <cutil_inline.h>
#include <iostream>

#include <iostream>
#include <iomanip>

#include "MersenneTwister.h"

#include "Block.h"
#include "util.h"

using namespace std;

//Convienence function that determines and iterates through the block pairs
//during a Fisher iteration.  If a non-square structure is desired with
//your blocks, this will need to be heavily edited.
void calcScoreAndHessian(Block** blocks, unsigned int& M1, DATA_TYPE *lt,
                         DATA_TYPE *derL, DATA_TYPE **derE2);

int main(/*int argc, char **argv*/)
{

  cublasInit();
  MTRand mt(MT_SEED);

  //Note that this code example does not set any convergence criterion for
  //the Fisher iterations.  That will have to be created if desired.
  unsigned int numFisherIters = 1;

  //////////////////////////////////////////////////////////////////////////////
  //Setting up the grid structure:
  //c is the number of points per block. Currently it must be a multiple of 32.
  //M1 is the dimension of the grid, e.g. M1 = 5 means a 5x5 grid of blocks, 
  //each of size c
  //n is total number of points equal to c*M1^2
  //It should not be necessary to force a square grid structure, however
  //this has not been tested and you will have to manage your sets of blocks explicitly.
  //////////////////////////////////////////////////////////////////////////////
  unsigned int c = 1024,
               M1 = 2;
  unsigned int n = c*M1*M1;


  //////////////////////////////////////////////////////////////////////////////
  //Here is where the data is generated or loaded.  Currently we only generate
  //random data that is not correlated between blocks.  If it is desired to
  //load data, use:
  // setXYs(...) instead of randomFill(...)
  // setData(...) instead of generateData(...)
  //with the desired data set. 
  //////////////////////////////////////////////////////////////////////////////

  //Matern covariance parameters
  DATA_TYPE lt[SP_PARAM_SIZE];
  lt[0] = 0;
  lt[1] = lt[2] = 2;

  //Create the grid
  Block **blocks = new Block*[M1];
  blocks[0] = new Block[M1*M1];

  //Fill the grid with location and measured data
  int i, j;
  for(j = 0; j < M1; j++) {
    blocks[0][j].randomFill(c,0.,(DATA_TYPE)j,1./(DATA_TYPE)M1);
    blocks[0][j].generateData(lt,0.,(DATA_TYPE)j);
  }
  for(i = 1; i < M1; i++) {
    blocks[i] = blocks[i-1] + M1;
    for(j = 0; j < M1; j++) {
      blocks[i][j].randomFill(c,(DATA_TYPE)i,(DATA_TYPE)j,1./(DATA_TYPE)M1);
      blocks[i][j].generateData(lt,(DATA_TYPE)i,(DATA_TYPE)j);
    }
  }
  //////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////
  //Eval Block Composite Likelihood
  DATA_TYPE blockCompLik = 0.;
  for(i = 0; i < M1 - 1; i++) {
    for(j = 0; j < M1; j++) {
      //lower right neighbor, if applicable
      if(j > 0) {
        blockCompLik += blocks[i][j].generateCompositeLik(lt,blocks[i+1][j-1]);
      }
      //upper and upper right neighbor, if applicable
      if((j + 1) % M1 != 0) {
        blockCompLik += blocks[i][j].generateCompositeLik(lt,blocks[i+1][j+1]);
        blockCompLik += blocks[i][j].generateCompositeLik(lt,blocks[i][j+1]);
      }
      //right neighbor
      blockCompLik += blocks[i][j].generateCompositeLik(lt,blocks[i+1][j]);
    }
  }
  std::cout << blockCompLik << "\n";
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  //Iterate to the MLE

  //Instantiate score and hessian arrays.
  DATA_TYPE *derL = new DATA_TYPE[SP_PARAM_SIZE], *d_derL,
        **derE2 = new DATA_TYPE*[SP_PARAM_SIZE], *d_derE2;
  derE2[0] = new DATA_TYPE[SP_PARAM_SIZE*SP_PARAM_SIZE];
  for(i = 1; i < SP_PARAM_SIZE; i++) 
    derE2[i] = derE2[i-1] + SP_PARAM_SIZE;

  //Instantiate paramater array to be estimated.
  DATA_TYPE ltold[SP_PARAM_SIZE],
        ltnew[SP_PARAM_SIZE], *d_ltnew;
  for(i = 0; i < SP_PARAM_SIZE; i++) {
    ltold[i] = lt[i];
  }

  //Allocate space on the GPU device for above quantities
  cutilSafeCall( cudaMalloc(&d_derL, SP_PARAM_SIZE*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMalloc(&d_ltnew, SP_PARAM_SIZE*sizeof(DATA_TYPE)) );
  cutilSafeCall( cudaMalloc(&d_derE2, SP_PARAM_SIZE*SP_PARAM_SIZE*sizeof(DATA_TYPE)) );

  unsigned int hTimer;
  DATA_TYPE gpuTime;
  cutilCheckError( cutCreateTimer(&hTimer) );
  cutilCheckError( cutResetTimer(hTimer) );
  cutilCheckError( cutStartTimer(hTimer) );

  //Fisher iterations
  for(i = 0; i < numFisherIters; i++) {
    memset(derL,0,SP_PARAM_SIZE*sizeof(DATA_TYPE));
    memset(derE2[0],0,SP_PARAM_SIZE*SP_PARAM_SIZE*sizeof(DATA_TYPE));
    calcScoreAndHessian(blocks,M1,ltold,derL,derE2);

    cutilSafeCall( cudaMemcpy(d_derL, derL, SP_PARAM_SIZE*sizeof(DATA_TYPE),
                              cudaMemcpyHostToDevice) );
    cutilSafeCall( cudaMemcpy(d_derE2, derE2[0], SP_PARAM_SIZE*SP_PARAM_SIZE*sizeof(DATA_TYPE),
                              cudaMemcpyHostToDevice) );

    leftLinCholSovle(d_derE2,d_derL,SP_PARAM_SIZE,1,/*needsChol*/true);

    cutilSafeCall( cudaMemcpy(ltnew, d_derL, SP_PARAM_SIZE*sizeof(DATA_TYPE),
                              cudaMemcpyDeviceToHost) );
    for(j = 0; j < SP_PARAM_SIZE; j++)
      ltnew[j] += ltold[j];

    for(j = 0; j < SP_PARAM_SIZE; j++)
      ltold[j] = ltnew[j];
  }
  cutilCheckError( cutStopTimer(hTimer) );
  gpuTime = 1.0e-3 * cutGetTimerValue(hTimer);
  std::cout << "n = " << n << ", c = " << c << ", time = " << gpuTime << std::endl;
  //////////////////////////////////////////////////////////////////////////////

  //clean up
  cublasShutdown();

  cutilSafeCall( cudaFree(d_derL) );
  cutilSafeCall( cudaFree(d_ltnew) );
  cutilSafeCall( cudaFree(d_derE2) );
  delete [] derE2[0];
  delete [] derE2;
  delete [] derL;
  delete [] blocks[0];
  delete [] blocks;
}

void calcScoreAndHessian(Block** blocks, unsigned int& M1, DATA_TYPE *lt,
                         DATA_TYPE *derL, DATA_TYPE **derE2) {
  int i, j;
  for(i = 0; i < M1 - 1; i++) {
    for(j = 0; j < M1; j++) {
      //lower right neighbor, if applicable
      if(j > 0) {
        blocks[i][j].addCompositeScoreHessContrib(lt,blocks[i+1][j-1],
                                                  derL,derE2);
      }
      //upper and upper right neighbor, if applicable
      if((j + 1) % M1 != 0) {
        blocks[i][j].addCompositeScoreHessContrib(lt,blocks[i+1][j+1],
                                                  derL,derE2);
        blocks[i][j].addCompositeScoreHessContrib(lt,blocks[i][j+1],
                                                  derL,derE2);
      }
      //right neighbor
      blocks[i][j].addCompositeScoreHessContrib(lt,blocks[i+1][j],
                                                derL,derE2);
    }
  }
}
