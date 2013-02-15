#!/bin/sh -x

git pull
make clean
make
./cg ~rvuduc3/matrices/bfwb62/bfwb62.rb 2> err.out