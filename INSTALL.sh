#!/usr/bin/bash
#This files compiles the FORTRAN programs in snake_src
cd ./snake_src
make
make -f Makefile_SHA 
make -f Makefile_VG
mkdir ../bin
mv snakeSHAS ../bin/
mv snakeSHA ../bin/
mv snakeVGRACE ../bin/
cd ../bin
#Perform test
./snakeSHAS et ../snake_test/grdata.nc ../snake_test/test_flm.nc ../snake_test/test_check.nc 60
