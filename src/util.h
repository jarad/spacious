/*
This Code is provied to be freely used, distributed, or modified.
However it comes without warranty of any kind.
Matt Wheeler 2011 

Convenience functions to solve a linear system Ax = B and calculate
the trace of a matrix.  All functions use the GPU.
*/

#ifndef UTIL_H
#define UTIL_H

//#define DOUBLE_PRECISION 1
#ifdef DOUBLE_PRECISION
#define DATA_TYPE double
#else
#define DATA_TYPE float
#endif

//Uses Cholesky Decomposition to solve d_C{x} = d_B (or x = d_C{^-1}d_B)
//Overwrites both d_C and d_B so watch out!
void leftLinCholSovle(DATA_TYPE *d_C, DATA_TYPE *d_B, const int& size, const int& cols,
                      const bool& needsChol = true);

//Calculates the trace of size-by-size matrix
DATA_TYPE trace(DATA_TYPE *mat, const int& size);

//Calculates the trace of the product of two size-by-size matrices
DATA_TYPE trace(DATA_TYPE *matA, DATA_TYPE *matB, const int& size);

#endif
