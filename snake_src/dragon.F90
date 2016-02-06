!====================================================================
! Source code by Wang et al. (2006)
! Modified by Oliver Lopez (2015) to read input from netcdf files
! and to perform analysis/synthesis on the same grid on different
! time steps
!
! Module name: dragon
! Purpose: Spherical harmonic analysis and synthesis using dragon 
!           by Wang et al. (2006)
! Called by:    snake_driver.F90
!====================================================================
      MODULE dragon_module
         USE netcdf
         USE snake_common_module
         USE snake_data_module, only: nc_abort
         implicit none
         real, allocatable :: v(:), sti, stj   
         
         integer :: ncid, varid, status
         CONTAINS
!====================================================================
!
! Name: DRGSHA
! Purpose: Spherical harmonic analysis using dragon (Wang et al. 2006)
!     Called from: snake_driver
!     Input: Geographic data (currently common)
!            reads netcdf step by step
!     Output: coefficients flmi, flmr, etc.. (currently flm as complex common)
!     Output: variable is stored in common VAR 

      SUBROUTINE DRGSHA
      implicit double precision(a-h,o-z)
      integer ::  imax, jmax, l, m, I, J
      real :: dcta, dfai, thetaa, thetab, phai, phaib
      double complex f1,f2
!      double complex, allocatable :: flm(:,:) !moved to common
      dimension coew(4),fw(4)
      double precision, dimension(:,:,:),allocatable::  fa,fb
      double precision, dimension(:,:),allocatable:: f
      double precision, dimension(:,:,:),allocatable:: coe
      real, dimension(:,:), allocatable :: TEMPVAR, LOCALVAR 
      character fout*100, dum1*100, format_string*100
         print *, 'This is dragon analysis'
 
      imax = size(lat)-1
      jmax = size(lon)-1
      thetaa = 90 - maxval(lat)
      thetab = 90 - minval(lat)
      phaia = minval(lon)
      phaib = maxval(lon)
      dcta = (thetab-thetaa) /imax
      dfai = (phaib - phaia) / jmax
      if (.not. allocated(v)) ALLOCATE (v(0:lmax))
      ALLOCATE (coe(imax,jmax,4))
      ALLOCATE (fa(0:(lmax+1),0:(lmax+1),imax))
      ALLOCATE (fb(0:(lmax+1),0:(lmax+1),imax))

      if (allocated(flm)) then
         print *, 'flm was already allocated'
         deallocate(flm)
         if(allocated(VAR)) then
            print *, 'VAR was already allocated'
            deallocate(VAR)
         endif
      endif

      ALLOCATE (flm(0:lmax,0:lmax,numt))
!     ------------Physical parameters-------------------
      pi=4.d0*datan(1.d0)
      arc=pi/180.d0

      tta=thetaa*arc
      pha=phaia*arc
      dc=dcta*arc
      df=dfai*arc
!     ------------Square root in equation (10)-----------
!     Compute and store in memory (variable 'v')
!      print *, 'square root in eq 10'
      do l=1,lmax
      v1=0.d0
        do i=1,l
          v1=v1+0.5d0*dlog(dble(2*i-1))-0.5d0*dlog(dble(2*i))
        end do
      v(l)=dexp(v1)
      end do
      v(0)=1
!       -----------Three-looped cycle computation of------------
!      print *, 'Three-looped cycle computation of ST'
!       -----------ST(l,m,x) and SXT(l,m,x) (equation 19)-------
      do i=1,imax               !i is for each latitude
        cta1=tta+(i-1)*dc
        cta2=tta+i*dc
        x1=dcos(cta1)           !cos(90-theta) at node below
        x2=dcos(cta2)           !cos(90-theta) at node above
        do m=0,lmax             !m is for each order m
          ! Initial values (equation 20), where sti is STm,m and
          ! stj is STm+1,m
          sti=stmm(m,x1,x2)             !calls stmm
          stj=-v(m)*dsqrt(dble(2*m+1))/dble(m+2)*((1.d0-x2*x2)**&
          (dble(m+2)/dble(2))-(1.d0-x1*x1)**(dble(m+2)/dble(2)))
          do l=m,lmax           !l is for each order l
            !Calculate st and sxt by calling
            !subroutine stsxt (equations 19.1 and 19.2)
            !Here, st and sxt will have only one value (l,m,x)
            call stsxt(st,sxt,l,m,x1,x2)
            !Store st and sxt
            fa(l,m,i)=st
            fb(l,m,i)=sxt
          end do
        end do
      end do

      ALLOCATE(LOCALVAR(0:imax,0:jmax))
      ALLOCATE(TEMPVAR(0:imax,0:jmax))
      ALLOCATE(VAR(numlat,numlon,numt))
      !initialize VAR to zeros
      do i=0,imax
         do j=0,jmax
            LOCALVAR(i,j)=0
         enddo
      enddo
!      print *, 'numlon is ', numlon
!      print*, 'numlat is ', numlat
!      print *, 'numt is ', numt
!      print *, 'flipvar is', flip_var
!      print *, 'checklon is ', check_lon
!     Open netcdf file
         status = nf90_open(path = FNAMEI, mode = nf90_nowrite, &
                  ncid = ncid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error opening input netcdf file'//&
                     ' (SUBROUTINE DRGSHA)')
         status = nf90_inq_varid(ncid, VARNAME, varid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding '//VARNAME//&
                     ' (SUBROUTINE DRGSHA)')
          status = nf90_get_att(ncid, varid, "_FillValue",miss_value)
               if (status /= NF90_NOERR) CALL nc_abort (status, &
                     'Error getting missing value (SUBROUTINE DRGSHA)')
      print *, 'Missing value is ', miss_value
      print *, "Latord is ", latord

      do k=1,numt
     print *, "Processing step number ",k, " out of ",numt
      print *, "Reading slice"
      if (flip_var .eq. 1) then
         !flip_Var =1 then variable needs to be flipped (from lon, lat, time to lat, lon, time)
      status = nf90_get_var(ncid, varid, LOCALVAR , start = (/1,1,k/), &
               count = (/numlon, numlat, 1/), map = (/numlat, 1/))
               if (status /= NF90_NOERR) CALL nc_abort &
                (status, 'Error reading variable (SUBROUTINE DRGSHA)')
               
      else
      status = nf90_get_var(ncid, varid, LOCALVAR, start = (/1,1,k/), &
               count = (/numlat, numlon, 1/))
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading variable '//&
                      ' (SUBROUTINE DRGSHA)')
 
      endif
       if (check_lon .eq. 1) then
      print *, "Shifting longitude-wise"
        TEMPVAR = LOCALVAR(:,0:(pos_split-2))   ! Store left half
        LOCALVAR(:,0:(pos_split-2))=LOCALVAR(:,(pos_split-1):(numlon-1))!Right-half to left half
        LOCALVAR(:,(pos_split-1):(numlon-1)) = TEMPVAR !stored left half to right-half
      endif

      if (latord .ge. 0) then
!     Flip VAR upside down in lat dimension
         print *, "Flipping latitude-wise"
      LOCALVAR = LOCALVAR(ubound(LOCALVAR,1):lbound(LOCALVAR,1):-1,:)
      lat = lat(ubound(lat,1):lbound(lat,1):-1)
      endif
      ! Set all missing values to zero:
      forall(J=0:jmax,I=0:imax,LOCALVAR(I,J) .eq. miss_value)&
            LOCALVAR(I,J)=0.0

       ! Save LOCALVAR in global VAR:
       VAR(:,:,k) = LOCALVAR

      print *, "computing coefficients of bilinear functions"
!       ----------Coefficients of bilinear functions
!                 stored in coe(i,j,1..4)  (A=1,B=2,C=3,D=4)
      do i=1,imax
        cta1=tta+(i-1)*dc
        cta2=tta+i*dc
        x1=dcos(cta1)
        x2=dcos(cta2)

        do j=1,jmax
         fai1=pha+(j-1)*df
         fai2=pha+j*df

         fw(1)=LOCALVAR(i-1,j-1)
         fw(2)=LOCALVAR(i,j-1)
         fw(3)=LOCALVAR(i,j)
         fw(4)=LOCALVAR(i-1,j)
         !A,B,C,D are given below respectively
         coe(i,j,1)=(fai1*(fw(4)-fw(3))+fai2*(fw(2)-fw(1)))&
           /(x2-x1)/(fai2-fai1)
         coe(i,j,2)=(x1*(fw(2)-fw(3))+x2*(fw(4)-fw(1)))&
           /(x2-x1)/(fai2-fai1)
         coe(i,j,3)=(fw(1)+fw(3)-fw(2)-fw(4))/(x2-x1)/(fai2-fai1)
         coe(i,j,4)=(x2*fai2*fw(1)-x1*fai2*fw(2)&
            +x1*fai1*fw(3)-x2*fai1*fw(4))/(x2-x1)/(fai2-fai1)
       end do
      end do

!      print *, "ALL GOOD"
!       ------------ Main Computation of f(l,m)------------
      print *, "Computing flm"
!  Initialize variable flm (set to zero):
      do m=0,lmax
        do l=m,lmax
           flm(l,m,k)=(0.d0,0.d0)
        end do
      end do
      do i=1,imax
        do m=0,lmax
          do l=m,lmax
            !Get st and sxt from memory:
            st=fa(l,m,i)
            sxt=fb(l,m,i)
            do j=1,jmax
              fai1=pha+(j-1)*df
              fai2=pha+j*df
              beta1=dmod(dble(m)*fai1,2.d0*pi)
              beta2=dmod(dble(m)*fai2,2.d0*pi)
              coew(1)=coe(i,j,1)
              coew(2)=coe(i,j,2)
              coew(3)=coe(i,j,3)
              coew(4)=coe(i,j,4)
              !calls subroutine SE:
              call SE(m,coew,fai1,fai2,beta1,beta2,f1,f2)
              !Calculation of Stokes coefficients flm (add to sum):
          !                     (i.e. equation 14)
              flm(l,m,k)=flm(l,m,k)+(-1.d0)**(m+1)*&
                 dsqrt(dble(1+2*l)/4.d0/pi)*(f1*sxt+f2*st)
            end do
          end do
        end do
      end do

       enddo  !Time cycle
!       Close the netcdf files
      status = nf90_close(ncid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error closing file '//&
                      ' (SUBROUTINE DRGSHA)')
 
      print *, "Dragon analysis done"
      ENDSUBROUTINE DRGSHA
!====================================================================
!
! Name: DRGSHS
! Purpose: Spherical harmonic synthesis using dragon (Wang et al. 2006)
!     Called from: snake_driver
!     Input: Coefficients must be in complex form flm (common)
!     Output: variable is stored in common VARBACK 
       
      SUBROUTINE DRGSHS
      implicit double precision(a-h,o-z)
      integer imax, jmax, l, m
      real :: dcta, dfai, thetaa, thetab, phai, phaib

      real :: res_lat, res_lon
      
      real, dimension(:,:), allocatable :: fbackr
      double complex, allocatable :: fback(:,:)

      print *, 'This is dragon synthesis'

!!!!!!!!   Check if lat or lon vector have been allocated.
         !! Default vector if not..
         if (.not. allocated(lon)) then
            print *, "Preparing a default global space"
          res_lat = 0.5
          res_lon = 0.5
          numlon = 360*2
          numlat = 180*2
          print*, "numlat set to ", numlat
          print*, "numlon set to ", numlon
          print*, "resolution set to ", res_lon
           lat = (/(90-res_lat/2-i*res_lat,i=0,numlat-1,1)/)
           lon = (/(0+res_lon/2+i*res_lon,i=0,numlon-1,1)/)
         endif
!!!!!!!!!!!



      imax = size(lat)-1
      jmax = size(lon)-1
      thetaa = 90 - maxval(lat)
      thetab = 90 - minval(lat)
      phaia = minval(lon)
      phaib = maxval(lon)
      dcta = (thetab-thetaa) /imax
      dfai = (phaib - phaia) / jmax
!     ------------Physical parameters-------------------
      pi=4.d0*datan(1.d0)
      arc=pi/180.d0

      tta=thetaa*arc
      pha=phaia*arc
      dc=dcta*arc
      df=dfai*arc
      if (.not. allocated(v)) ALLOCATE (v(0:lmax))

!     ------------Square root in equation (10)-----------
!     Compute and store in memory (variable 'v')
!      print *, 'square root in eq 10'
      do l=1,lmax
      v1=0.d0
        do i=1,l
          v1=v1+0.5d0*dlog(dble(2*i-1))-0.5d0*dlog(dble(2*i))
        end do
      v(l)=dexp(v1)
      end do
      v(0)=1

      ALLOCATE(VARBACK(numlat,numlon,numt))
      ALLOCATE(fback(0:imax,0:jmax))
      ALLOCATE(fbackr(0:imax,0:jmax))

      do k=1,numt
      write(*,'(f8.2,a25)')k/dble(numt)*100,' % finished for synthesis'
       do i=0,imax
       do j=0,jmax
        fback(i,j)=(0.d0,0.d0)
        fbackr(i,j)=0.d0
        VARBACK(i+1,j+1,k)=0.d0 !Initialize VARBACK..
       end do
       end do

      do i=0,imax
!      write(*,'(f8.2,a25)')i/dble(imax)*100,' % finished for synthesis'
      ctai=tta+i*dc
      xi=dcos(ctai)
      do m=0,lmax
      do l=m,lmax
      a=tlm(l,m,xi)*(-1)**m*dsqrt((2*l+1.d0)/4/pi)
      do j=0,jmax     
      faij=pha+j*df
      beta=dmod(dble(m)*faij,2.d0*pi)
      fback(i,j)=fback(i,j)+a*flm(l,m,k)*cdexp((0.d0,1.d0)*beta)
      if(m.ne.0)then
      fback(i,j)=fback(i,j)+a*dconjg(flm(l,m,k))*cdexp((0.d0,-1.d0)*&
                  beta)
      endif
      fbackr(i,j)=dreal(fback(i,j))
!      VARBACK(ubound(VARBACK,1)-i,j+1,k)=fbackr(i,j)
      VARBACK(i+1,j+1,k)=fbackr(i,j)
      end do !j cycle
      end do !l cycle
      end do !m cycle
      end do !i cycle
      enddo !Time cycle

     

      END SUBROUTINE DRGSHS

      subroutine stsxt(st,sxt,l,m,x1,x2)
!     this subroutine calculates Eq.(19) of the paper
!     st and sxt are the output variables
         implicit double precision(a-h,o-z)
         integer l, m
         dimension a(3),b(3),c(3)

         if(l.eq.m)then

         if(m.eq.0)then
         a(1)=v(m)
         b(1)=v(m)
         else
         a(1)=v(m)*(1.d0-x1*x1)**(dble(m)/2.d0)
         b(1)=v(m)*(1.d0-x2*x2)**(dble(m)/2.d0)
         endif

         c(1)=sti
         st=c(1)

         sxt=v(m)*(-0.5d0/(1.d0+dble(m)/2.d0))&
            *((1.d0-x2*x2)**(dble(m+2)/2.d0)-(1.d0-x1*x1)&
            **(dble(m+2)/2.d0))

         elseif(l.eq.(m+1))then

           if(m.eq.0)then
         a(2)=x1*dsqrt(dble(2*m+1))*v(m)
         b(2)=x2*dsqrt(dble(2*m+1))*v(m)
         else
         a(2)=x1*(1.d0-x1*x1)**(dble(m)/2.d0)&
                   *dsqrt(dble(2*m+1))*v(m)
         b(2)=x2*(1.d0-x2*x2)**(dble(m)/2.d0)&
                   *dsqrt(dble(2*m+1))*v(m)
         endif

         c(2)=stj
         st=c(2)
         sxt=-1.d0/(l+2)*((1-x2*x2)*b(2)-(1.d0-x1*x1)*a(2))&
         +dsqrt(dble(l*l-m*m))/(l+2)*sti
         else
         a(3)=-dsqrt(dble((l-m-1)*(l+m-1))/dble((l+m)*(l-m)))*a(1)&
          +x1*dble(2*l-1)/dsqrt(dble((l-m)*(l+m)))*a(2)
         b(3)=-dsqrt(dble((l-m-1)*(l+m-1))/dble((l+m)*(l-m)))*b(1)&
          +x2*dble(2*l-1)/dsqrt(dble((l-m)*(l+m)))*b(2)
         a(1)=a(2)
         a(2)=a(3)
         b(1)=b(2)
         b(2)=b(3)

         c(3)=-dble(2*l-1)/dble(l+1)/dsqrt(dble(l*l-m*m))&
          *((1.d0-x2*x2)*b(1)-(1.d0-x1*x1)*a(1))&
          +dble(l-2)/dble(l+1)*dsqrt(dble((l+m-1)*(l-m-1))&
          /dble((l-m)*(l+m)))*c(1)
         c(1)=c(2)
         c(2)=c(3)
         st=c(3)

         sxt=-(1.d0-x2*x2)/dble(l+2)*b(2)+(1.d0-x1*x1)/dble(l+2)*a(2)&
         +dsqrt(dble(l*l-m*m))/dble(l+2)*c(1)
         endif
         return
         end
!     -----------stmm (equation 20.1)-------
      double precision function stmm(m,x1,x2)
!     this subroutine calculates first equation of Eqs.(20) of the paper
!     stmm is the output variable
         implicit double precision(a-h,o-z)
         integer m
         dimension a(2),b(2)

         if(m.eq.0)then
         a(1)=x2-x1
         w=a(1)
         endif

         if(m.eq.1)then
         b(1)=dsqrt(2.d0)/4.d0*(-(dacos(x2)-dacos(x1))+&
               dsqrt(1.d0-x2*x2)*x2-dsqrt(1.d0-x1*x1)*x1)
         w=b(1)
         endif

         if(m.ge.2.and.mod(m,2).eq.0)then
         a(2)=v(m)/dble(m+1)*((1.d0-x2*x2)**(dble(m)/2.d0)*x2&
                         -(1.d0-x1*x1)**(dble(m)/2.d0)*x1)&
                         +dble(m)/dble(m+1)*dsqrt(dble((2*m-1)*(2*m-3))&
                         /dble(2*m*(2*m-2)))*a(1)
         w=a(2)
         a(1)=a(2)
         endif
         if(m.ge.3.and.mod(m,2).ne.0)then
         b(2)=v(m)/dble(m+1)*((1.d0-x2*x2)**(dble(m)/2.d0)*x2&
                         -(1.d0-x1*x1)**(dble(m)/2.d0)*x1)&
                         +dble(m)/dble(m+1)*dsqrt(dble((2*m-1)*(2*m-3))&
                         /dble(2*m*(2*m-2)))*b(1)
         w=b(2)
         b(1)=b(2)
         endif
         stmm=w
         return
         end

!     -----------SE (equation 15)---------------
      subroutine SE(m,coe,fai1,fai2,beta1,beta2,xi1,xi2)
!this subroutine calculates Eq.(15) of the paper
!coe is the inputed coeficients A,B,C,D of bilinear function
!xi1,xi2 are output of Eq.(15) when A,C and D,B used respectively
        implicit double precision(a-h,o-z)
        dimension coe(4)
        integer m
        double complex xi1,xi2,xi
        a=coe(1)
        b=coe(2)
        c=coe(3)
        d=coe(4)
        xi=(0.0d0,1.0d0)
        if(m.eq.0)then
          xi1=a*(fai2-fai1)+0.5d0*c*(fai2**2-fai1**2)
          xi2=d*(fai2-fai1)+0.5d0*b*(fai2**2-fai1**2)
        else
          xi1=((a+c*fai2)/dble(m)*xi+c/dble(m*m))*cdexp(-xi*beta2)&
          -((a+c*fai1)/dble(m)*xi+c/dble(m*m))*cdexp(-xi*beta1)
          xi2=((d+b*fai2)/dble(m)*xi+b/dble(m*m))*cdexp(-xi*beta2)&
          -((d+b*fai1)/dble(m)*xi+b/dble(m*m))*cdexp(-xi*beta1)
        endif
        return
        end

 
      double precision function tlm(l,m,x)
!     this function calculates Eq.(7) of the paper
!     tlm is the output variable
        implicit double precision(a-h,o-z)
        integer l,m
        dimension a(3)
        if(m.lt.0.or.m.gt.l.or.dabs(x).gt.1.d0)&
            write(*,*)' bad arguments!'
          pi=4.d0*datan(1.d0)

         if(l.eq.m)then

         if(m.eq.0)then
         a(1)=1.d0
         else
         a(1)=(1.d0-x*x)**(dble(m)/2.d0)
         endif

         w=a(1)

         elseif(l.eq.(m+1))then

         if(m.eq.0)a(2)=x*dsqrt(dble(2*m+1))
         if(m.ne.0)a(2)=x*(1.d0-x*x)**(dble(m)/2.d0)*dsqrt(dble(2*m+1))

         w=a(2)

         else

         a(3)=-dsqrt(dble((l-m-1)*(l+m-1))/dble((l+m)*(l-m)))*a(1)&
             +x*dble(2*l-1)/dsqrt(dble((l-m)*(l+m)))*a(2)
         a(1)=a(2)
         a(2)=a(3)
         w=a(3)

         endif

         tlm=v(m)*w
         return
         end

      END MODULE dragon_module
