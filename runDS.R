#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2){
   cat('Usage: ./runDS.R input.nc output.nc\n')
   quit()
}
fnamei = args[1]
fnameo = args[2]

source("./ncDSfilter.R") 
#fnamei = "deltaGSM_UTCSR.nc"
#fnameo = "deltaGSMDS.nc"
ncdfDSfilter(fnamei,fnameo,c('clm','slm'),c('clm','slm'))
command<-sprintf('ncdump -h %s',fnameo)
system(command)



