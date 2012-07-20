/*
This Code is provied to be freely used, distributed, or modified.
However it comes without warranty of any kind.
Matt Wheeler 2011 

Block class definiation.  Defines a "block" of spatial data consisting
of locations and measured data.  Uses NVIDIA GPUs for most of the
calculations.
*/

#ifndef BLOCK_H
#define BLOCK_H

#define MT_SEED 1009
#define MT_DB_SEED 811
#define BLOCK_N 32
#define SP_PARAM_SIZE 3

#include "MersenneTwister.h"

#include "util.h"

class Block {
  public:
    Block();
    virtual ~Block();

    DATA_TYPE* getXs() { return h_x; }
    DATA_TYPE* getYs() { return h_y; }

    DATA_TYPE* getData() { return h_data; }
    void setData(DATA_TYPE*,const unsigned int&);

    unsigned int getSize() { return size; }
   
    void setXYs(DATA_TYPE *Xs, DATA_TYPE *Ys, const unsigned int& size);

    //Generates the Covariance matrix (and derivatives if applicable) 
    //based on the paramaters in lt and stores the result in cov
    //(and der if applicable) based on the locations and data stored
    //in the calling block.
    void generateCov(DATA_TYPE *lt, DATA_TYPE *cov);
    void generateCovAndDerivatives(DATA_TYPE *lt,DATA_TYPE*,DATA_TYPE*);

    //Generates and returns the likelihood based on the data 
    //and locations in the calling block
    DATA_TYPE generateLik(DATA_TYPE*);

    //Generates and returns the composite likelihood based on the data 
    //and locations in the calling block as well as the data and locations
    //in the block "neighbor."
    DATA_TYPE generateCompositeLik(DATA_TYPE *lt, Block& neighbor);

    //Adds the composite contribution to the score and hessian of the
    //block pair defined by the calling block and the block passed
    //in the "neighbor" argument.
    void addCompositeScoreHessContrib(DATA_TYPE *lt, Block& neighbor,
                                      DATA_TYPE *derL, DATA_TYPE **derE2);

    //Generates the score (derL) and hessian (derE2) based on the data
    //and locations in the calling block.
    void makeScoreHess(DATA_TYPE *lt, DATA_TYPE *derL, DATA_TYPE **derE2);
   
    //Convenience functions used to generate data randomly.  Used in
    //testing and timing calculations. 
    void generateData(DATA_TYPE*, const DATA_TYPE&, const DATA_TYPE&);
    void randomFill(const unsigned int&, const DATA_TYPE&, const DATA_TYPE&,
                    const DATA_TYPE&);

  protected:
    DATA_TYPE *h_x, *d_x,
          *h_y, *d_y,
          *h_data, *d_data;
    unsigned int size;
    MTRand mt;
    bool allocatedOnDevice,
         allocatedOnHost,
         allocatedDataOnDevice,
         allocatedDataOnHost;
};

#endif
