/*
This Code is provied to be freely used, distributed, or modified.
However it comes without warranty of any kind.
Matt Wheeler 2011 

Cholesky decomposition implemented on a GPU
*/


#ifndef CHOLESKY_KERNEL
#define CHOLESKY_KERNEL

//Needs to change based on compute capability of the GPU device:
//2.x: 32
//1.x: 22
#define CUDA_MAX_BLOCK_SIZE 32

//Evaluate cholesky decomposition iteratively on a block.  Used
//here to calculate decomposition on the upper-leftmost block,
//along the block diagonal, that has yet to be evaluated.
__global__ void cuda_chol_iter(DATA_TYPE *m, int n , int boffset ) {
  int k ;
  int x = threadIdx.x ;
  int y = threadIdx.y ;
  int bsize = blockDim.y ;

  __shared__ DATA_TYPE b[CUDA_MAX_BLOCK_SIZE][CUDA_MAX_BLOCK_SIZE] ;
  b[x][y] = 0;
  if(x >= y) {
    b[x][y] = m[ (y+boffset)*n + boffset + x ] ;
  }

  for (k = 0; k < bsize; k++) {
    __syncthreads( );
    if ( y == k ) {
      DATA_TYPE fac = rsqrt(b[y][y]);
      if ( x >= y ) {
        b[x][y] *= fac;
      }
    }
    __syncthreads ( ) ;
    if ( y > k && x >= y ) {
      b[x][y] -= b[y][k] * b[x][k];
    }
  }

  __syncthreads ( ) ;
  m[(boffset+y)*n + boffset + x] = b[x][y] ;
}

//Once the above function is run on the appropriate block, all values directly
//below that block must be updated according to the Cholesky algorithm.
//This function does that.
__global__ void cuda_solve_strip(DATA_TYPE *m, int n , int boffset) {
  int k ;
  int x = threadIdx.x;
  int y = threadIdx.y;
  int bsize = blockDim.y;
  int by = (blockIdx.y+1)*bsize + boffset;

  __shared__ volatile DATA_TYPE b[CUDA_MAX_BLOCK_SIZE][CUDA_MAX_BLOCK_SIZE];
  __shared__ DATA_TYPE ttl[CUDA_MAX_BLOCK_SIZE][CUDA_MAX_BLOCK_SIZE];

  b[x][y] = m[by + x + (boffset+y)*n];
  ttl[x][y] = m[(boffset + y)*n + boffset + x];
  __syncthreads( );

  for(k = 0; k < bsize; k++) {
    if(y == k) {
      b[x][y] /= ttl[k][k];
    }
    __syncthreads();

    if(y > k) {
      b[x][y] -= b[x][k]*ttl[y][k];
    }
    __syncthreads();
  }

  m[by + x + (boffset + y)*n] = b[x][y];
}

//Primary cholesky decomposition function.  To be called to find the
//decomposition of matrix m of size nxn.  The argument bsize is
//used to define the block sizes used on the GPU.  Overwrites
//m with the decomposition.
void cholesky_cuda(DATA_TYPE *m, int n , int bsize) {
  DATA_TYPE *dm;
  int size = n*n*sizeof(DATA_TYPE);
  int iter = (n + bsize - 1) / bsize;
  int i ;
  dim3 threads(bsize, bsize);
  cudaMalloc(&dm, size);
  cudaMemcpy(dm, m, size, cudaMemcpyHostToDevice);
  for (i = 0; i < iter; i++) {
    // top-most block of current vertical strip
    cuda_chol_iter<<<1,threads>>>(dm, n, i*bsize);
    cudaThreadSynchronize();
    cutilCheckMsg("cuda_chol_iter failed\n");
    cudaThreadSynchronize ( );

    if(iter - i - 1 > 0) {
      dim3 strip(1, iter - i - 1, 1);
      //Solve rest of current strip.
      cuda_solve_strip<<<strip,threads>>>(dm, n, i*bsize);
      cudaThreadSynchronize();
      cutilCheckMsg("cuda_solve_strip failed\n");
      cudaThreadSynchronize ( ) ;
      //Contribution of current vertical strip to next strip over
#ifdef DOUBLE_PRECISION
      cublasDsyrk( 'L', 'N', bsize, bsize, -1, 
  		  	&dm[(i+1)*bsize+(i*bsize)*n], n, 1, 
  		  	&dm[(i+1)*bsize+((i+1)*bsize)*n], n );
#else
      cublasSsyrk( 'L', 'N', bsize, bsize, -1, 
  		  	&dm[(i+1)*bsize+(i*bsize)*n], n, 1, 
  		  	&dm[(i+1)*bsize+((i+1)*bsize)*n], n );
#endif
      cudaThreadSynchronize();
      cutilCheckMsg("chol ssyrk 1 failed\n");

      if(n-((i+2)*bsize) > 0) {  
#ifdef DOUBLE_PRECISION
        cublasDgemm( 'N', 'T', n-((i+2)*bsize), bsize, bsize, -1, 
    	  	  	&dm[(i+2)*bsize+(i*bsize)*n], n, 
    	  	  	&dm[(i+1)*bsize+(i*bsize)*n], n, 1, 
    	  	  	&dm[(i+2)*bsize+(i+1)*bsize*n], n );
    
        //Contribution of current vertical strip to lower
        //right symm matrix
        cublasDsyrk( 'L', 'N', n-((i+2)*bsize), bsize, -1, 
    	      		&dm[(i+2)*bsize+(i*bsize)*n], n, 1, 
       	      		&dm[(i+2)*bsize*(n+1)], n );
#else
        cublasSgemm( 'N', 'T', n-((i+2)*bsize), bsize, bsize, -1, 
    	  	  	&dm[(i+2)*bsize+(i*bsize)*n], n, 
    	  	  	&dm[(i+1)*bsize+(i*bsize)*n], n, 1, 
    	  	  	&dm[(i+2)*bsize+(i+1)*bsize*n], n );
    
        //Contribution of current vertical strip to lower
        //right symm matrix
        cublasSsyrk( 'L', 'N', n-((i+2)*bsize), bsize, -1, 
    	      		&dm[(i+2)*bsize+(i*bsize)*n], n, 1, 
       	      		&dm[(i+2)*bsize*(n+1)], n );
#endif
        cudaThreadSynchronize();
        cutilCheckMsg("chol ssyrk 2/gemm failed\n");
      }
    }
  }

  cudaMemcpy(m, dm, size, cudaMemcpyDeviceToHost);
  cudaFree(dm);
}

#endif
