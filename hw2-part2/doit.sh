#!/bin/sh -x

git pull
make clean
make
./cg ~rvuduc3/matrices/s3dkt3m2/s3dkt3m2.rb 2> err.out