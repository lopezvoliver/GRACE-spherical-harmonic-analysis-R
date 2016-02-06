#!/usr/bin/env Rscript 							 	
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2){
      cat('Usage: ./C20_slr.R output.txt path_to_coef\n')
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
# 3) Replace Degree 2, order 0 coefficients by those from SLR analysis 		#
#						(Cheng and Tapley 2004)		#
#  Source: ftp://ftp.csr.utexas.edu/pub/slr/degree_2/C20_RL05.txt          	#
#  They are stored as a text file
slr_path='ftp://ftp.csr.utexas.edu/pub/slr/degree_2/C20_RL05.txt'

#						(Cheng and Tapley 2004)		
#  Source: ftp://ftp.csr.utexas.edu/pub/slr/degree_2/C20_RL05.txt          
cat(sprintf('Readind %s',slr_path))
slrfile<-readLines(slr_path)
#identify how many lines to skip (header):
skip.lines=0; cond=TRUE; j=1;
while(cond==TRUE){
	fc<-substr(slrfile[j],1,1) #first character
	if ((fc=='#')){ 
		skip.lines=skip.lines+1; 		
	} else{
		cond=FALSE
	}
	j=j+1
}
nl=length(slrfile)


#Extract coefficients from ftp file 
C20<-as.numeric(substr(slrfile[(skip.lines+1):nl],13,29))
#substract mean
C20<-C20-mean(C20)
#subset C20 to the corresponding dates in GRACE:
C20d1<-substr(slrfile[(skip.lines+1):nl],57,64)  #SLR calendar day 1
C20d2<-substr(slrfile[(skip.lines+1):nl],73,80)	 #SLR calendar day 2
C20d1<-as.Date(C20d1,format='%Y%m%d')		 
C20d2<-as.Date(C20d2,format='%Y%m%d')
C20mid<-as.Date(colMeans(t(cbind(C20d1,C20d2))),origin='1970/01/01') #SLR midpoint


#For each SLR midpoint, find the closest one in GRACE midpoint (within a certain range)
maxdiff=17 #if the closest is more than 15 days offset, then it is set to NA
corrind=array(NA,dim=length(C20mid))
for (k in 1:length(C20mid)){
	#get distance to each gracemidpoint
	dist=sqrt((as.numeric(gracemid-C20mid[k]))^2)
	if (min(dist)<maxdiff){
		corrind[k]=which(dist==min(dist))
	}
}
validslrj=which(!is.na(corrind))
C20mid=C20mid[validslrj]
corrind=corrind[validslrj]
C20=C20[validslrj]  

out<-cbind(corrind,C20)  #corrind will give the GRACE month number
write.table(x=out,file=fnameo,row.names=FALSE,col.names=FALSE)
command<-sprintf('Saved to %s \n',fnameo)
cat(command)
command<-sprintf('cat %s',fnameo)
cat(command)
cat('\n')
system(command)
