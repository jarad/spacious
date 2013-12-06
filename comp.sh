#!/bin/bash

INC="-I/usr/local/lib64/R/include -I/usr/local/cuda-5.0/include -I/home/rjp/include"
#LIB="-L/usr/local/lib64/R/lib -L/usr/local/cuda-5.0/lib64 -L/home/rjp/lib -lR -lRblas -lRlapack -lm -pthread -lcublas -lcudart -lmagma -lmagmablas -lmagma"
LIB="-L/usr/local/lib64/R/lib -L/usr/local/cuda-5.0/lib64 -L/home/rjp/lib -lR -lRblas -lRlapack -lm -pthread -lcublas -lcudart"
#FLAGS="-DDEBUG -DCLINE -DPTHREAD -DADD_ -DCUBLAS_GFORTRAN -DHAVE_CUBLAS -DGPUSHMEM=200 -DCUDA -Wall -ggdb -O3"
FLAGS="-DDEBUG -DCLINE -DPTHREAD -DCUDA -Wall -ggdb -O3 -m64"
#FLAGS="-DDEBUG -DCLINE -DPTHREAD -DCUDA"
#FLAGS="-DDEBUG -DCLINE -Wall -ggdb";

#g++ -DCLINE -Wall -I/usr/local/lib64/R/include -L/usr/local/lib64/R/lib -c src/BlockComp.cpp
#g++ -DCLINE -Wall -I/usr/local/lib64/R/include -L/usr/local/lib64/R/lib -c src/covs.cpp

#echo "nvcc -c -o utils_cuda.o src/utils_cuda.cu"
#nvcc -c -o utils_cuda.o src/utils_cuda.cu

echo "g++ $FLAGS $INC -o block_test src/BlockComp.cpp src/covs.cpp src/utils.cpp src/utils_cuda.cpp utils_cuda.o $LIB"
g++ $FLAGS $INC -o block_test src/BlockComp.cpp src/covs.cpp src/utils.cpp src/utils_cuda.cpp utils_cuda.o $LIB

#echo "nvcc -Xcompiler $FLAGS $INC -Xlinker $LIB -o block_test src/BlockComp.cpp src/covs.cpp src/utils.cpp src/utils_cuda.cpp"
#nvcc -Xcompiler $FLAGS $INC -Xlinker $LIB src/BlockComp.cpp src/covs.cpp src/utils.cpp src/utils_cuda.cpp -o blocktest
#/usr/local/cuda/bin/nvcc -gencode arch=compute_10,code=sm_10 -gencode arch=compute_11,code=sm_11 -gencode arch=compute_12,code=sm_12 -gencode arch=compute_13,code=sm_13 -gencode arch=compute_20,code=sm_20 -shared -Xlinker  -L/usr/local/lib64/R/lib -lR -L/usr/local/cuda/lib64 -lcublas kendall.o classification.o   rinterface.o mi.o sort.o granger.o qrdecomp.o correlation.o hcluster.o distance.o matmult.o lsfit.o cuseful.o -o gputools.so

