#!/bin/bash
#cd ./RAW
for file in *.gz;
do gunzip $file;
done;

