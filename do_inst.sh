#!/bin/bash

#R CMD INSTALL -l ~/Rlib .
R CMD INSTALL --configure-args="--with-cuda" -l ~/Rlib .
