#!/usr/bin/env Rscript 							 	
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2){
      cat('Usage: ./C11.R output.txt path_to_coef\n')
   quit()
}
fnameo = args[1]
#local_path = '~/Dropbox/WORK/SH/GRACECOEF/'
local_path = args[2]
# We need the path where the raw coefficients are stored, to read the date information:

# Read and convert time periods from file names
GRACE_csr_periods<-dir(local_path,pattern='GSM-2_')
GCSR1<-as.Date(substr(GRACE_csr_periods,7,13),format='%Y%j')
GCSR2<-as.Date(substr(GRACE_csr_periods,15,21),format='%Y%j')
gracemid=as.Date(colMeans(t(cbind(GCSR1,GCSR2))),origin='1970/01/01')#GRACE midpoint

day1=GCSR1
day2=GCSR2


#################################################################################
#										#
# GRACE GSM files represent the full gravity field as Stokes coefficients	# 
#							   Clm(t), Slm(t)	#
# Here we apply the following post-processing needed to obtain the signal in 	#
# equivalent water height:							#
#										#
# 4) Replace Degree 1 (C11, S11, C10) coefficients (geocenter estimates)        #
#       with those estimated by Swenson, Chambers, and Wahr (2008)              #
#  Source: ftp://podaac.jpl.nasa.gov/allData/tellus/L2/degree_1/deg1_coef.txt   # 
deg1_path='ftp://podaac.jpl.nasa.gov/allData/tellus/L2/degree_1/deg1_coef.txt'
maxdiff=17
cat(sprintf('Reading %s\n',deg1_path))
deg1file<-readLines(deg1_path)
#identify how many lines to skip (header):
f6=as.numeric(substr(deg1file,1,6))
line1=which(!is.na(f6))[1]
deg1file=deg1file[line1:length(deg1file)]
month<-unique(as.numeric(substr(deg1file,1,6)))
#add day '15' to 'month'
midpt<-as.Date(sprintf('%s15',month),format='%Y%m%d')
#find the corresponding GRACE month (same strategy as step 4)
corrind=array(NA,dim=length(midpt))
for (k in 1:length(midpt)){
           #get distance to each gracemidpoint
           dist=sqrt((as.numeric(gracemid-midpt[k]))^2)
        if (min(dist)<maxdiff){
                           corrind[k]=which(dist==min(dist))
                }
}
validind=which(!is.na(corrind))
#Get coefficient C10 and C11 from file:
C1110=as.numeric(substr(deg1file,17,17+19-1))
C10=C1110[as.logical(1:length(C1110) %%2)]
C11=C1110[as.logical(!(1:length(C1110) %%2))]
#Coefficient S11:
S1110=as.numeric(substr(deg1file,17+19,17+19+19-1))
S11=S1110[as.logical(!(1:length(S1110) %%2))]

#remove temporal mean:
C10=C10-mean(C10)
C11=C11-mean(C11)
S11=S11-mean(S11)

out<-cbind(corrind,C10)
out<-cbind(out,C11)
out<-cbind(out,S11)

write.table(x=out,file=fnameo,row.names=FALSE,col.names=FALSE)
command<-sprintf('Saved to %s \n',fnameo)
cat(command)
command<-sprintf('cat %s',fnameo)
cat(command)
cat('\n')
system(command)
