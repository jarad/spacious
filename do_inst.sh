#!/bin/bash

#R CMD INSTALL -l ~/Rlib .
R CMD INSTALL --configure-args="--with-cuda --with-CUDA_BS=256 --with-debug" -l ~/Rlib .
