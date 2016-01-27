#################################################################################
#										#
#			      R  netcdf Spherical Harmonics 			#
#			Correlated-error filter (de-striping)			#
#				   Oliver Lopez					#
#					2014					#
#										#
#			       Output: netcdf file 				#
#							 			#
#										#
#The Coefficients of the gravity field changes of the same order and same	#
#parity in degrees are correlated with each other. 				#
#The correlated error can be reduced by high-pass filtering the correlated	#
#coefficients as a function of degree. 						#
#Smooth the Stokes coefficients for a particular order (m) with a quadratic	#
#polynomial in a moving window of width w centered about degree l. The		#
#polynomial coefficients are obtained by least-squares.				#
#										#
#References: Swenson and Wahr (2006) and Duan et al. (2009) 			#
#################################################################################
#requires: ncdf4
#Input: name of SH coefficients file e.g. AMAZONAS.CSIRO.MCOEF.nc
#which contains dimensions l,m and time
#and two variables named: Clm (l,m,time) and Slm(l,m,time)
#Output: netcdf file (e.g. AMAZONAS.CSIRO.MCOEFDS.nc)
#-----------------Function to get batches of coefficients for a given m
#m must be in 5<=m<=(ltrunc-4*alpha) where alpha=(window-1)/2
library(ncdf4)
ncdfDSfilter<-function (fnamei,fnameo,vari=c('clm','slm'),varo=c('Clm','Slm'),start=NA,count=NA) {
    cat(sprintf('NCDF Correlated error filter (Swenson and Wahr 2006, Duan et al. 2009)\n'))
    cat(sprintf('R script written by Oliver Lopez, 2014\n'))
    

    #open file and retrieve information
	nci<-nc_open(fnamei)

        if (is.na(start)){
           tfoo<-ncvar_get(nci,'time')
           start=c(1,1,1)
           count=c(-1,-1,length(tfoo))
        }else{
           start=c(1,1,start)
           count=c(-1,-1,count)
        }


        time<-ncvar_get(nci,'time',start[3],count[3])
	Clm<-ncvar_get(nci,vari[1],start,count)
	Slm<-ncvar_get(nci,vari[2],start,count)
	nc_close(nci)



    cat(sprintf('Succesfully loaded coefficients file %s\n',fnamei))
	cat(sprintf('Preparing time-independent components\n'))

	#Pre-allocating space
	Clm_ds=Clm
	Slm_ds=Slm
   	#Time independent components
	ltrunc=dim(Clm)[1]-1
	p = 2 #Maximum order of polynomial fit

	#---------Define window widths---------------------
	K = 10 #Empirical constant for window width
	A = 30 #Empirical constant for window width

	mvec = seq(0,ltrunc,1) #Order(m) array
	
	wemp<-array(0,dim=c(length(mvec),2))
	wemp[,1]=A*exp(-mvec/K)
	wemp[,2]=5
	wemp<-apply(wemp,1,max)
	#Empirical window width as a function of order
	#Swenson and Wahr (personal communication, 2008, cited in Duan et al. 2009)

	w = 2*floor(wemp/2)+1 #Round to nearest integer
	
	#---------------------------------------------------
	indx <- function (l,m) {
        #Errors
        if (l<m) stop("Degree cannot be lower than order")
        if (l>ltrunc) stop("'degree' must be <=ltrunc")
        
        window <-w[m+1]
        alpha <-(window-1)/2
        if (m>(ltrunc-4*alpha)) stop("m must be <ltrunc-2*(w-1)")
        
        index<-seq((l-2*alpha),(l+2*alpha),by=2)
        
        #Near-tesseral
        if (l<m+2*alpha) {
            #Check parity of m
            parm=(m%%2)
            #Check parity of l
            parl=(l%%2)
            
            a=m
            b=m+2*window-2
            index<-seq(a,b,by=2)
            
            if(parl!=parm) index<-index+1
            
            
        }
        
        #Near-maximum degree
        if (l>ltrunc-2*alpha) {
            #Check parity of l
            if (l%%2) {
                index<-seq(ltrunc-4*alpha,ltrunc,by=2)-1  #odd
            } else {index<-seq(ltrunc-4*alpha,ltrunc,by=2)}
            
        } 
        index
    }

	#Decorrelation is done for the SCs in order m=5 and above 	
	#Swenson and Wahr (personal communication, 2008, cited in Duan et al. 2009)
	
	#For each m, set window width and then: {
		#For each l, get batch of coefficients to be considered for the polynomial fit, considering that: 
		#for near-tesseral (l<m+(w-1)) or near-maximum degrees (l>lmax-(w-1), 'window' coefficients of the same parity as l with lowest or highest degrees are used
		#for the rest, 'window' coefficients (centered at l) are used.	
		#Also, coefficients l<m do not exist, neither do coefficients m>lmax
	#Once the coefficients to be fit are obtained, perform a polynomial fit (i.e. get the coefficients Q_{lm}^i)
	#The correlated part to be removed for each l,m combination is sum_{i=0}^p Q_{lm}^i l^i
	
    #Time cycle
    nts=dim(Clm)[3]
    cat(sprintf("Applying filter to coefficients\n"))
    cat(sprintf("Number of time steps: %d\n",nts))
    cat(sprintf('Current time step: '))

    #Temporal variables
    QC<-array(0,dim=(p+1))
    QS<-array(0,dim=(p+1))
	
    #######################Main cycle#######################
    for (ts in 1:nts){
        cat(sprintf('t=%d ',ts))
        #Main cycle  (m=5 to 40)  #as in GRACE-documentation "Converting Release-04 Gravity Coefficients into Maps of Equivalent Water Thickness"
	for (m in 5:40) {
                if (nts==1){
         	Cbatchall=Clm[(m+1):(ltrunc+1),m+1]
		Sbatchall=Slm[(m+1):(ltrunc+1),m+1]

                }else{
		Cbatchall=Clm[(m+1):(ltrunc+1),m+1,ts]
		Sbatchall=Slm[(m+1):(ltrunc+1),m+1,ts]
                }	
		for (l in m:ltrunc){
			index<-indx(l,m)	
			if (!all.equal(w[m+1],length(index))) stop("error! length(index) is not the same as window")
			if (min(index)<m) stop("error! attempting index for l less than m")
			#get coefficients to fit
                        if (nts==1){
                        cbatch = Clm[index+1,m+1]
                        sbatch = Slm[index+1,m+1]
                        }else{
			cbatch=Clm[index+1,m+1,ts]
			sbatch=Slm[index+1,m+1,ts]
                        }
			#Get fit 
			fitC<-lm(cbatch~poly(index,p,raw=TRUE))$coefficients
			QC[1]=fitC[[1]]; QC[2]=fitC[[2]] ;QC[3]=fitC[[3]]
            
			fitS<-lm(sbatch~poly(index,p,raw=TRUE))$coefficients
			QS[1]=fitS[[1]]; QS[2]=fitS[[2]]; QS[3]=fitS[[3]]
			
			#Remove correlated part from coefficient
                        if (nts==1){
                        Clm_ds[l+1,m+1]=Clm[l+1,m+1]-(QC[1]+QC[2]*l+QC[3]*l^2)
                        Slm_ds[l+1,m+1]=Slm[l+1,m+1]-(QS[1]+QS[2]*l+QS[3]*l^2)
                        }else{

			Clm_ds[l+1,m+1,ts]=Clm[l+1,m+1,ts]-(QC[1]+QC[2]*l+QC[3]*l^2)
			Slm_ds[l+1,m+1,ts]=Slm[l+1,m+1,ts]-(QS[1]+QS[2]*l+QS[3]*l^2)
                        }
		}
	}

	for (m in 41:ltrunc) {				#Every coefficient above n=40, m=40 is set to zero
		for (l in m:ltrunc){	
                        if (nts==1){
                           Clm_ds[l+1,m+1]=0
                           Slm_ds[l+1,m+1]=0
                        }else{
			Clm_ds[l+1,m+1,ts]=0
			Slm_ds[l+1,m+1,ts]=0
                        }
		}
	}

    }
    #######################End Main cycle#######################
    
#send output to new netcdf file
cat(sprintf('Succesfully applied de-striping filter to coefficients from file %s\n',fnamei))

#list(C=Clm_ds,S=Slm_ds)

#Coefficients are now stored in Clm_ds and Slm_ds
#and have dimensions l,m,t (degree, order, time)

#we will store them as a small netcdf file with dimensions l,m,t
#Define dimensions
if (nts==1) time=1
dimt<-ncdim_def('time','GRACE month',time)
diml<-ncdim_def('l','degree',0:ltrunc)
dimm<-ncdim_def('m','order',0:ltrunc)

#Define variables
ncoef1=ncvar_def(varo[1],'mm/day',list(diml,dimm,dimt),-9999,longname='coefficient 1',prec='float')
ncoef2=ncvar_def(varo[2],'mm/day',list(diml,dimm,dimt),-9999,longname='coeffiicent 2',prec='float')

#Define nc file
ncout<-nc_create(fnameo,list(ncoef1,ncoef2),force_v4=TRUE)

cat(sprintf('Creating output file %s\n',fnameo))

#save coefficients
Clm_ds[is.na(Clm_ds)]=-9999
Slm_ds[is.na(Slm_ds)]=-9999
nt=length(time)
for (i in 1:nt){
    if(nts==1){
    ncvar_put(ncout,varo[1],Clm_ds,start=c(1,1,1),count=c(-1,-1,1))
    ncvar_put(ncout,varo[2],Slm_ds,start=c(1,1,1),,count=c(-1,-1,1))

    }else{
    ncvar_put(ncout,varo[1],Clm_ds[,,i],start=c(1,1,i),count=c(-1,-1,1))
    ncvar_put(ncout,varo[2],Slm_ds[,,i],start=c(1,1,i),count=c(-1,-1,1))
    }
}

#Close file
nc_close(ncout)

cat(sprintf('Succesfully created file  %s\n',fnameo))
}
