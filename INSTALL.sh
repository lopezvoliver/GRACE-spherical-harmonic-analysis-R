#!/usr/bin/bash
#This files compiles the FORTRAN programs in snake_src
cd ./snake_src
make
make -f Makefile_SHA 
make -f Makefile_VG
mv snakeSHAS ../
mv snakeSHA ../
mv snakeVGRACE ../
