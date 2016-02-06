#!/bin/bash
cd ./CSR
for file in *.gz;
do gunzip $file;
done;

