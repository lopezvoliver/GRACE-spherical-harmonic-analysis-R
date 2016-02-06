!====================================================================
! Source code by Oliver Lopez
!
! Module name: snake_flmout
! Purpose: Output results from the analysis
! Called by:    snake_driver.F90
!====================================================================
      MODULE snake_flmout_module
         USE netcdf
         USE snake_common_module
         implicit none
         character fout*100,dum1*100,format_string*100
         integer, allocatable :: lmvalues(:)
         integer ncido, ldimid, mdimid, todimid, &
                 lvarid, mvarid, tovarid, &
                 frvarid, fivarid, avarid, bvarid, &
                 cvarid, svarid, status, l, m
         double precision :: fac, pi
         double precision, allocatable :: flmr(:,:), &
                                          flmi(:,:), &
                                          alm(:,:), &
                                          blm(:,:), &
                                          clm(:,:), &
                                          slm(:,:)
          
         CONTAINS

         SUBROUTINE FLMOUT
! Output coefficients to netcdf AND .dat files (for faster comparison with dragon)
      print *, "Preparing output netcdf file"
!    ---------------------PREPARE OUTPUT NETCDF FILE---------
      ALLOCATE (lmvalues(0:lmax))
      lmvalues = (/(i, i=0, lmax, 1)/)
      ALLOCATE (flmr(0:lmax,0:lmax))
      ALLOCATE (flmi(0:lmax,0:lmax))
      ALLOCATE (alm(0:lmax, 0:lmax))
      ALLOCATE (blm(0:lmax, 0:lmax))
      ALLOCATE (clm(0:lmax, 0:lmax))
      ALLOCATE (slm(0:lmax, 0:lmax))
   
      pi = 4.d0*datan(1.d0)

     ! Initialize all coefficients to NA (-9999)
 
      status = nf90_create(FNAMEO, nf90_hdf5, ncido)
      ! Define dimensions: l, m, t
      status = nf90_def_dim(ncido, "l", lmax+1, ldimid)
      status = nf90_def_dim(ncido, "m", lmax+1, mdimid)
      status = nf90_def_dim(ncido, "time", numt, todimid)
      ! Define variables: l, m, t, flmr, flmi
      status = nf90_def_var(ncido, "l", nf90_int, (/ldimid/), lvarid)
      status = nf90_def_var(ncido, "m", nf90_int, (/mdimid/), mvarid)
      status = nf90_def_var(ncido, "time", nf90_int, (/todimid/), &
       tovarid)
      status = nf90_def_var(ncido, "flmr", nf90_double, (/ldimid, & 
       mdimid, todimid/), frvarid)
      status = nf90_def_var(ncido, "flmi", nf90_double, (/ldimid, & 
       mdimid, todimid/), fivarid)
      status = nf90_def_var(ncido, "alm", nf90_double, (/ldimid, & 
       mdimid, todimid/), avarid)
      status = nf90_def_var(ncido, "blm", nf90_double, (/ldimid, & 
       mdimid, todimid/), bvarid)
      status = nf90_def_var(ncido, "clm", nf90_double, (/ldimid, & 
       mdimid, todimid/), cvarid)
      status = nf90_def_var(ncido, "slm", nf90_double, (/ldimid, & 
       mdimid, todimid/), svarid) 
      ! Define and assign attributes
      status = nf90_put_att(ncido, frvarid, "_FillValue", &
               -9999.0d0)
      status = nf90_put_att(ncido, fivarid, "_FillValue", &
               -9999.0d0)
      status = nf90_put_att(ncido, avarid, "_FillValue", &
               -9999.0d0)
      status = nf90_put_att(ncido, bvarid, "_FillValue", &
               -9999.0d0)
      status = nf90_put_att(ncido, cvarid, "_FillValue", &
               -9999.0d0)
      status = nf90_put_att(ncido, svarid, "_FillValue", &
               -9999.0d0)
      ! End define mode
      status = nf90_enddef(ncido)
!      print *, "Status of end define mode ",status
      ! Write dimension variables (l, m, time)
!      print *, lmvalues
      status = nf90_put_var(ncido, lvarid,lmvalues)
!      print *, "Status of lmvalues", status
      status = nf90_put_var(ncido, mvarid,lmvalues)
      status = nf90_put_var(ncido, tovarid, time)


     !Initialize coefficients to NA:
      do l =0,lmax
       do m=l+1,lmax 
        flmr(l,m) = -9999.0
        flmi(l,m) = -9999.0
        alm(l,m) = -9999.0
        blm(l,m) = -9999.0
        clm(l,m) = -9999.0
        slm(l,m) = -9999.0
      end do
      end do


! Time cycle for outputting coefficients into netcdf file:
      do k=1,numt
     print *, "time step ", k, "out of ", numt
!       ----Output of coefficients-----
!       (copied from ncdragon)
!            -(3 types)-:
!     flm - defined in (1) of the (Wang) paper
!     alm, blm - defined in (A1) of the paper
!     clm,slm - defined in (A2) of the paper
      format_string="(I0.3)"
      write(dum1,format_string)k
      fout = trim("tmp")//trim(dum1)//".dat"
      open(20,file=fout)
      do l=0,lmax
      do m=0,l

!     Factor (l-m)!/(l+m)!
!     (needed for equations A.4 and A.5)
      fac=1.d0
      do i=l-m+1,l+m
        fac=fac/dble(i)
      end do

      !This needs to be checked according to dragon:
!      flmr(l,m) = aout(l+1,m+1,k)  !aout, bout, declared in common module
!      flmi(l,m) = bout(l+1,m+1,k)  !they are the output of spk

!     flm(l,m,k) contains the coefficients as a double complex variable
      flmr(l,m) = dreal(flm(l,m,k))
      flmi(l,m) = dimag(flm(l,m,k))

!     Equations A.4 and A.5 for the case m=0 and otherwise:
      if(m.eq.0)then
      alm(l,m)=dsqrt(dble(2*l+1)/4.d0/pi)*flm(l,m,k)!A.4
      blm(l,m)=0.d0
      clm(l,m)=flm(l,m,k)/dsqrt(4.d0*pi)         !A.5
      slm(l,m)=0.d0
      else
      alm(l,m)=(-1)**m*dsqrt(dble(2*l+1)/pi*fac)*flmr(l,m) !A.4
      blm(l,m)=(-1)**(m+1)*dsqrt(dble(2*l+1)/pi*fac)*flmi(l,m) !A.4
      clm(l,m)=(-1)**m*flmr(l,m)/dsqrt(2.d0*pi)     !A.5
      slm(l,m)=(-1)**(m+1)*flmi(l,m)/dsqrt(2.d0*pi) !A.5
      endif
!      Write 1 value to txt file:
      write(20,'(1x,2i6,6d18.10)')l,m,flmr(l,m),flmi(l,m),&
                                  alm(l,m),blm(l,m),clm(l,m),slm(l,m)
      end do
      end do
      close(20)

!     Write all values to netcdf file (in this time step):
      print *, "Writing to output file"
      status = nf90_put_var(ncido, frvarid, flmr, & 
               (/1,1,k/),(/lmax+1, lmax+1, 1/))
      status = nf90_put_var(ncido, fivarid, flmi, &
               (/1,1,k/),(/lmax+1, lmax+1, 1/))
      status = nf90_put_var(ncido, avarid, alm, & 
               (/1,1,k/),(/lmax+1, lmax+1, 1/))
      status = nf90_put_var(ncido, bvarid, blm, & 
               (/1,1,k/),(/lmax+1, lmax+1, 1/))
      status = nf90_put_var(ncido, cvarid, clm, & 
               (/1,1,k/),(/lmax+1, lmax+1, 1/))
      status = nf90_put_var(ncido, svarid, slm, & 
               (/1,1,k/),(/lmax+1, lmax+1, 1/))

      end do
!       Close the netcdf file
      status = nf90_close(ncido)
         END SUBROUTINE FLMOUT


      END MODULE snake_flmout_module

