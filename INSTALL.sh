#!/bin/bash
#This files compiles the FORTRAN programs in snake_src
$FORTLIBS
cd ./snake_src
make
make -f Makefile_SHA 
make -f Makefile_VG
make -f Makefile_SHS
mkdir ../bin
mv snakeSHAS ../bin/
mv snakeSHA ../bin/
mv snakeVGRACE ../bin/
mv snakeSHS ../bin/
cd ../bin
#Perform test
./snakeSHAS et ../snake_test/grdata.nc ../snake_test/test_flm.nc ../snake_test/test_check.nc 60
