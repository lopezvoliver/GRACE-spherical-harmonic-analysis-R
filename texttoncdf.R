#!/usr/bin/env Rscript
args<- commandArgs(trailingOnly=TRUE)
if(length(args)<2){
cat('Usage: ./textoncdf.R output.nc path_to_coef\n')
quit()
}
fnameo=args[1]
local_path = args[2]
library(ncdf4)

fortran2double<-function(v){
                #Given a character vector v (e.g. v[1:3]>
                #[1] "  0.5837063244D-03" " -0.3435992171D-03" " -0.6436518226D-03"
                nc=nchar(v[1])
                exps<-as.numeric(substr(v,nc-2,nc))
                charf<-substr(v,1,nc-4)
                charF<-gsub('\\s','',charf)
                values<-as.numeric(charF) 
                val<-values*10^(exps)
                return(val)
}

# Read and convert time periods from file names
GRACE_csr_periods<-dir(local_path,pattern='GSM-2_')
GCSR1<-as.Date(substr(GRACE_csr_periods,7,13),format='%Y%j')
GCSR2<-as.Date(substr(GRACE_csr_periods,15,21),format='%Y%j')
gracemid=as.Date(colMeans(t(cbind(GCSR1,GCSR2))),origin='1970/01/01')#GRACE midpoint

day1=GCSR1
day2=GCSR2

#File names to read
fils=list.files(local_path,pattern='*GSM-2_*')

#we will store the coefficients it as a small netcdf file with dimensions l,m,t
#Define dimensions	
dimt<-ncdim_def('time','GRACE month',1:length(day1))
diml<-ncdim_def('degree','l',0:60)
dimm<-ncdim_def('order','m',0:60)
	
#Define variables 
ncday1=ncvar_def('day1','days since 2000-01-01',list(dimt),-9999,longname='GRACE day start',prec='integer')
ncday2=ncvar_def('day2','days since 2000-01-01',list(dimt),-9999,longname='GRACE day end',prec='integer')	
ncoef1=ncvar_def('clm','geoid height',list(diml,dimm,dimt),-9999,longname='Clm',prec='float')
ncoef2=ncvar_def('slm','geoid height',list(diml,dimm,dimt),-9999,longname='Slm',prec='float')
	
#Define nc file
ncout<-nc_create(fnameo,list(ncday1,ncday2,ncoef1,ncoef2),force_v4=TRUE)

cat(sprintf('\n Reading GRACE coefficients and writing to file %s\n',fnameo))

#save day periods 
ncvar_put(ncout,'day1',as.numeric(day1-as.Date('2000-01-01')))
ncvar_put(ncout,'day2',as.numeric(day2-as.Date('2000-01-01')))

#loop through files and save coefficients
for (k in 1:length(fils)){
	cat(sprintf('\n Reading file %s \n',fils[k]))
	fname=sprintf('%s%s',local_path,fils[k])
	allfile<-readLines(fname) #reads all file line per line
	#identify how many lines to skip:
	skip.lines=0; cond=TRUE; j=1;
	while(cond==TRUE){
		ff<-substr(allfile[j],1,6) #first five characters of line
		if (!(ff=='GRCOF2')){ 
			skip.lines=skip.lines+1; 		
		} else{
			cond=FALSE
		}
		j=j+1
	}
	nl=length(allfile)
	#Extract coefficients from file
	deg<-as.numeric(substr(allfile[(skip.lines+1):nl],10,11))
	ord<-as.numeric(substr(allfile[(skip.lines+1):nl],15,16))
	coef1<-substr(allfile[(skip.lines+1):nl],18,35)
	coef2<-substr(allfile[(skip.lines+1):nl],37,54)
	C<-fortran2double(coef1)
	S<-fortran2double(coef2)
	CLM=array(NA,dim=c(max(deg)+1,max(ord)+1,1))
	SLM=array(NA,dim=c(max(deg)+1,max(ord)+1,1))
	#convert to matrix form
	for (j in 1:length(deg)){
		CLM[deg[j]+1,ord[j]+1,1]=C[j]
		SLM[deg[j]+1,ord[j]+1,1]=S[j]
	}
	
	#save coefficients 
	CLM[is.na(CLM)]=-9999
	SLM[is.na(SLM)]=-9999
	cat(sprintf('\n Writing to %s\n',fnameo))
	ncvar_put(ncout,'clm',CLM,start=c(1,1,k),count=c(61,61,1))
	ncvar_put(ncout,'slm',SLM,start=c(1,1,k),count=c(61,61,1))
}
#Close file
nc_close(ncout)
command<-sprintf('ncdump -h %s',fnameo)
system(command)
