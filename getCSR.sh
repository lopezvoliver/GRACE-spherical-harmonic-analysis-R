#! /bin/bash
mkdir ./CSR
cd ./CSR
HOST=podaac.jpl.nasa.gov
USER=anonymous
PASS=anonymous
ftp -inv $HOST << EOF
user $USER $PASS
cd /allData/grace/L2/CSR/RL05/
mget *0060_0005.gz*
