!====================================================================
! Source code by Oliver Lopez
!
! Module name: snake_data
! Purpose: Read input data from a netcdf file for the main program
!          snake_driver.F90
!====================================================================
      MODULE snake_data_module
         USE netcdf
         USE snake_common_module
      implicit none
      !Switches
      integer ::  cell_c, &
                  status
      !Integers related to reading netcdf (Ids):
      integer ::  ncid, &
!                  numlon, &
!                  numlat, &
!                  numt, &
                  varid, &
                  londimid, &
                  latdimid, &
                  tdimid, &
                  lonvarid, &
                  latvarid, &
                  tvarid, &
                  nlat, &
                  nlon
                  
      integer, dimension(nf90_max_var_dims) :: dimIDs, dids 
      real, allocatable :: lons(:), &
                           lats(:)
!                              lon(:), &
!                              lat(:), &
!                              time(:), &
      real ::  res_lat, &
               res_lon

      CONTAINS
!====================================================================
!
! Name: nc_abort
! Purpose: Stop the program due to netcdf-related error
!          Display error from netcdf library and the caller
!     Called from: get_geo_data
!     Input: status of netcdf and Error from caller
!     Output: None, it stops the program

      SUBROUTINE nc_abort(status, message)
         CHARACTER (LEN=*), INTENT(IN) :: message
         INTEGER, INTENT(IN) :: status
         PRINT *, message !error from caller
         PRINT *, NF90_STRERROR(status) !netcdf error
         stop
      END SUBROUTINE nc_abort

!====================================================================
!
! Name: get_flm_data
! Purpose: Opens netcdf file containing the flm input data
!     Called from: snake_driver
!     Calls: nc_abort
!
!     Input: name of input file
!     Output: data is stored in common flm complex variable
      SUBROUTINE get_flm_data(FNAMEI)
         CHARACTER (LEN=*), INTENT(IN) :: FNAMEI
         integer ::  flmivarid, &
                     flmrvarid, &
                     timeid, &
                     lid, &
                     mid, &
                     numl, &
                     numm

         double precision, allocatable :: flmr(:,:,:), &
                                          flmi(:,:,:) 
                
         status = NF90_OPEN(path = FNAMEI, mode = nf90_nowrite, &
                  ncid = ncid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error opening input netcdf file'//&
                     ' (SUBROUTINE get_flm_data)')
         status = NF90_INQ_VARID(ncid, 'flmr', flmrvarid)
               if (status /= NF90_NOERR) then
                  print *, 'Warning: looking for clm variable'
                  status = NF90_INQ_VARID(ncid, "clm",flmrvarid)
                     if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding flmr or clm variable'//&
                     ' (SUBROUTINE get_flm_data)')
               endif
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding flmr variable'//&
                     ' (SUBROUTINE get_flm_data)')                  

         status = NF90_INQ_VARID(ncid, 'flmi', flmivarid)
               if (status /= NF90_NOERR) then
                  print *, 'Warning: looking for slm variable'
                  status = NF90_INQ_VARID(ncid, "slm",flmivarid)
                     if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding flmi or slm variable'//&
                     ' (SUBROUTINE get_flm_data)')
               endif
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding flmi variable'//&
                     ' (SUBROUTINE get_flm_data)')                  

        status = NF90_INQ_DIMID(ncid, "l", lid)
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_DIMID(ncid, "degree",lid)
                     if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding l dimension'//&
                     ' (SUBROUTINE get_flm_data)')
               endif
           status = NF90_INQUIRE_DIMENSION(ncid, lid, len = numl)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining l dimension '//&
                      ' (SUBROUTINE get_flm_data)') 
         status = NF90_INQ_DIMID(ncid, "m", mid)
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_DIMID(ncid, "order", mid)
                     if (status /=NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding m dimension'//&
                     ' (SUBROUTINE get_flm_data)')
               endif
            status = NF90_INQUIRE_DIMENSION(ncid, mid, len = numm)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining m dimension '//&
                      ' (SUBROUTINE get_flm_data)') 
         status = NF90_INQ_DIMID(ncid, "time", timeid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding time variable'//&
                     ' (SUBROUTINE get_flm_data)')
          status = NF90_INQUIRE_DIMENSION(ncid, timeid, len = numt)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining time dimension '//&
                      ' (SUBROUTINE get_flm_data)') 
               if (numl .ne. numm) then
                  print *, 'l dimension is not equal to m dimension'
                  stop
                endif

         lmax = numl - 1
         print *, 'lmax is ', lmax
         print *, 'numt is ', numt

         
         ALLOCATE(flmr(0:lmax,0:lmax,numt))
         ALLOCATE(flmi(0:lmax,0:lmax,numt))
         if (allocated(flm)) then
           print *, 'Warning: flm was already allocated'
           deallocate(flm)
         endif
         ALLOCATE (flm(0:lmax,0:lmax,numt))


           status = NF90_GET_VAR(ncid, flmrvarid, flmr )
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading flmr variable'//&
                      ' (SUBROUTINE get_flm_data)') 
           status = NF90_GET_VAR(ncid, flmivarid, flmi )
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading flmi variable'//&
                      ' (SUBROUTINE get_flm_data)')                 
 
          do k=1,numt
               flm(:,:,k) = CMPLX(flmr(:,:,k),flmi(:,:,k))
          enddo

          print *, flm(0,0,1)

         ALLOCATE (time(numt))
          status = NF90_INQ_VARID(ncid, "time", tvarid)
              if (status /= NF90_NOERR) CALL nc_abort &
                    (status, 'Error reading time dimension '//&
                     ' (SUBROUTINE get_flm_data)') 
          status = NF90_GET_VAR(ncid, tvarid, time)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading time variable'//&
                      ' (SUBROUTINE get_flm_data)')
 
      END SUBROUTINE get_flm_data

!====================================================================
!
! Name: get_clm_data
! Purpose: Opens netcdf file containing clm input data (e.g. GRACE)
!          and convert to flm data 
!     Called from: snake_driver
!     Calls: nc_abort
!
!     Input: name of input file
!     Output: data is stored in common flm complex variable
!     Update: addition of a gaussian filter for GRACE coefficients

      SUBROUTINE get_clm_data(FNAMEI)
         CHARACTER (LEN=*), INTENT(IN) :: FNAMEI
         integer ::  clmvarid, &
                     slmvarid, &
                     timeid, &
                     lid, &
                     mid, &
                     numl, &
                     numm

         double precision, allocatable :: clm(:,:,:), &
                                          slm(:,:,:), &
                                          flmr(:,:), &
                                          flmi(:,:)
         real :: wl, ae
         integer :: l,m
         double precision :: fac, pi
         ae = 6371! mean radius of Earth in km
         pi = 4.d0*datan(1.d0)

         status = NF90_OPEN(path = FNAMEI, mode = nf90_nowrite, &
                  ncid = ncid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error opening input netcdf file'//&
                     ' (SUBROUTINE get_clm_data)')
          status = NF90_INQ_VARID(ncid, "clm", clmvarid)
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_VARID(ncid, "Clm",clmvarid)
                     if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding clm variable'//&
                     ' (SUBROUTINE get_clm_data)')
               endif
           status = NF90_INQ_VARID(ncid, "slm", slmvarid)
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_VARID(ncid, "Slm",slmvarid)
                     if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding slm variable'//&
                     ' (SUBROUTINE get_clm_data)')
               endif
                         
       status = NF90_INQ_DIMID(ncid, "l", lid)
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_DIMID(ncid, "degree",lid)
                     if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding l dimension'//&
                     ' (SUBROUTINE get_clm_data)')
               endif
 
          status = NF90_INQUIRE_DIMENSION(ncid, lid, len = numl)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining l dimension '//&
                      ' (SUBROUTINE get_clm_data)') 
              status = NF90_INQ_DIMID(ncid, "m", mid)
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_DIMID(ncid, "order", mid)
                     if (status /=NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding m dimension'//&
                     ' (SUBROUTINE get_clm_data)')
               endif
 
           status = NF90_INQUIRE_DIMENSION(ncid, mid, len = numm)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining m dimension '//&
                      ' (SUBROUTINE get_clm_data)') 
         status = NF90_INQ_DIMID(ncid, "time", timeid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding time variable'//&
                     ' (SUBROUTINE get_clm_data)')
          status = NF90_INQUIRE_DIMENSION(ncid, timeid, len = numt)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining time dimension '//&
                      ' (SUBROUTINE get_clm_data)') 
               if (numl .ne. numm) then
                  print *, 'l dimension is not equal to m dimension'
                  stop
                endif

         lmax = numl - 1
         print *, 'lmax is ', lmax
         print *, 'numt is ', numt
 
         ALLOCATE(clm(0:lmax,0:lmax,numt))
         ALLOCATE(slm(0:lmax,0:lmax,numt))
         ALLOCATE(flmr(0:lmax,0:lmax))
         ALLOCATE(flmi(0:lmax,0:lmax))
        
         if (allocated(flm)) then
           print *, 'Warning: flm was already allocated'
           deallocate(flm)
         endif
         ALLOCATE (flm(0:lmax,0:lmax,numt))


           status = NF90_GET_VAR(ncid, clmvarid, clm )
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading Clm variable'//&
                      ' (SUBROUTINE get_clm_data)') 
           status = NF90_GET_VAR(ncid, slmvarid, flmi )
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading Slm variable'//&
                      ' (SUBROUTINE get_clm_data)')                 
 
          do k=1,numt
               !Convert clm and slm to flmr and flmi:
               do l=0,lmax
                ! Gaussian filter:
                ! Wl = exp(-(nr/a)^2/4ln2) 
                !wl = exp(-(l*gausr/ae)**2/(4*log(2.0)))
               do m=l+1,lmax
                  flmr(l,m)=-9999.0
                  flmi(l,m)=-9999.0
               enddo
               do m=0,l

               fac=1.d0
               do i=l-m+1,l+m
                  fac=fac/dble(i)
               end do

                  if(m.eq.0)then
                     flmr(l,m)=clm(l,m,k)*dsqrt(4.d0*pi)
                     flmi(l,m)=0.d0
                  else
                     flmr(l,m)=clm(l,m,k)*(-1)**(-m)*dsqrt(2.d0*pi)
                     flmi(l,m)=slm(l,m,k)*(-1)**(-m-1)*dsqrt(2.d0*pi)
                  endif
               enddo
               enddo
               !send to global variable:
               flm(:,:,k) = CMPLX(flmr(:,:),flmi(:,:))
          enddo

          print *, flm(0,0,1)

         ALLOCATE (time(numt))
          status = NF90_INQ_VARID(ncid, "time", tvarid)
              if (status /= NF90_NOERR) CALL nc_abort &
                    (status, 'Error reading time dimension '//&
                     ' (SUBROUTINE get_clm_data)') 
          status = NF90_GET_VAR(ncid, tvarid, time)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading time variable'//&
                      ' (SUBROUTINE get_clm_data)')
 
      END SUBROUTINE get_clm_data 

!====================================================================
!
! Name: get_grace_data
! Purpose: Opens netcdf file containing GRACE coefficients
!          and convert to flm data 
      ! assuming the name is (1) clm, slm, (2) Clm, Slm
! IMPORTANT: These are not the same 'clm,slm' defined as in dragon      
!     Called from: snake_driver
!     Calls: nc_abort
!
!     Input: name of input file
!     Output: data is stored in common flm complex variable
!     Update: addition of a gaussian filter for GRACE coefficients

      SUBROUTINE get_grace_data(FNAMEI)
         CHARACTER (LEN=*), INTENT(IN) :: FNAMEI
         integer ::  almvarid, &
                     blmvarid, &
                     timeid, &
                     lid, &
                     mid, &
                     numl, &
                     numm

         double precision, allocatable :: alm(:,:,:), &
                                          blm(:,:,:), &
                                          flmr(:,:), &
                                          flmi(:,:)
         integer :: l,m
         double precision :: fac, pi
         pi = 4.d0*datan(1.d0)

         status = NF90_OPEN(path = FNAMEI, mode = nf90_nowrite, &
                  ncid = ncid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error opening input netcdf file'//&
                     ' (SUBROUTINE get_grace_data)')
          status = NF90_INQ_VARID(ncid, "clm", almvarid)
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_VARID(ncid, "Clm",almvarid)
                     if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding alm variable'//&
                     ' (SUBROUTINE get_grace_data)')
               endif
           status = NF90_INQ_VARID(ncid, "slm", blmvarid)
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_VARID(ncid, "Slm",blmvarid)
                     if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding blm variable'//&
                     ' (SUBROUTINE get_grace_data)')
               endif
                         
       status = NF90_INQ_DIMID(ncid, "l", lid)
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_DIMID(ncid, "degree",lid)
                     if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding l dimension'//&
                     ' (SUBROUTINE get_grace_data)')
               endif
 
          status = NF90_INQUIRE_DIMENSION(ncid, lid, len = numl)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining l dimension '//&
                      ' (SUBROUTINE get_grace_data)') 
              status = NF90_INQ_DIMID(ncid, "m", mid)
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_DIMID(ncid, "order", mid)
                     if (status /=NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding m dimension'//&
                     ' (SUBROUTINE get_grace_data)')
               endif
 
           status = NF90_INQUIRE_DIMENSION(ncid, mid, len = numm)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining m dimension '//&
                      ' (SUBROUTINE get_grace_data)') 
         status = NF90_INQ_DIMID(ncid, "time", timeid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding time variable'//&
                     ' (SUBROUTINE get_grace_data)')
          status = NF90_INQUIRE_DIMENSION(ncid, timeid, len = numt)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining time dimension '//&
                      ' (SUBROUTINE get_grace_data)') 
               if (numl .ne. numm) then
                  print *, 'l dimension is not equal to m dimension'
                  stop
                endif

         lmax = numl - 1
         print *, 'lmax is ', lmax
         print *, 'numt is ', numt
 
         ALLOCATE(alm(0:lmax,0:lmax,numt))
         ALLOCATE(blm(0:lmax,0:lmax,numt))
         ALLOCATE(flmr(0:lmax,0:lmax))
         ALLOCATE(flmi(0:lmax,0:lmax))
        
         if (allocated(flm)) then
           print *, 'Warning: flm was already allocated'
           deallocate(flm)
         endif
         ALLOCATE (flm(0:lmax,0:lmax,numt))


           status = NF90_GET_VAR(ncid, almvarid, alm )
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading alm variable'//&
                      ' (SUBROUTINE get_grace_data)') 
           status = NF90_GET_VAR(ncid, blmvarid, blm )
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading blm variable'//&
                      ' (SUBROUTINE get_grace_data)')                 
 
          do k=1,numt
               do l=0,lmax
                do m=l+1,lmax
                  flmr(l,m)=-9999.0
                  flmi(l,m)=-9999.0
               enddo
               do m=0,l

               fac=1.d0
               do i=l-m+1,l+m
                  fac=fac/dble(i)
               end do

                  if(m.eq.0)then
                      flmr(l,m)=alm(l,m,k)
                      !Inversion of A.4
                      flmi(l,m)=0.d0
                  else
                     flmr(l,m)=alm(l,m,k)*(-1)**(-m)
                     flmi(l,m)= blm(l,m,k)*(-1)**(-m-1)
                  endif
               enddo
               enddo
               !send to global variable:
               flm(:,:,k) = CMPLX(flmr(:,:),flmi(:,:))
          enddo

          print *, flm(0,0,1)

         ALLOCATE (time(numt))
          status = NF90_INQ_VARID(ncid, "time", tvarid)
              if (status /= NF90_NOERR) CALL nc_abort &
                    (status, 'Error reading time dimension '//&
                     ' (SUBROUTINE get_grace_data)') 
          status = NF90_GET_VAR(ncid, tvarid, time)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading time variable'//&
                      ' (SUBROUTINE get_grace_data)')
 
      END SUBROUTINE get_grace_data 
      
!====================================================================
!
! Name: get_geo_data
! Purpose: Opens netcdf file containing the input data and determines:
!     1. Spatiotemporal details - number of grid cells, latitudes, 
!                                longitudes, number of steps
!     2. Transforms data if necessary (rotation and/or dimension
!        permutation)
!
!     Called from: snake_driver
!     Calls: nc_abort
!
!     Input: name of input file
!     Output: ordered data with dimensions lat, lon, time

      SUBROUTINE get_geo_data(FNAMEI)
         CHARACTER (LEN=*), INTENT(IN) :: FNAMEI
!         print *, "get_geo_data(",FNAMEI,")"

         status = NF90_OPEN(path = FNAMEI, mode = nf90_nowrite, &
                  ncid = ncid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error opening input netcdf file'//&
                     ' (SUBROUTINE get_geo_data)')
         status = NF90_INQ_VARID(ncid, VARNAME, varid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding '//VARNAME//&
                     ' (SUBROUTINE get_geo_data)')
!         status = NF90_INQ_DIMID(ncid, "lon", dimIDs(1))
!               if (status /= NF90_NOERR) CALL nc_abort &
!                     (status, 'Error finding lon variable'//&
!                     ' (SUBROUTINE get_geo_data)')
 
         status = NF90_INQ_DIMID(ncid, "lon", dimIDs(1))
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_VARID(ncid, "Longitude",dimIDs(1))
                     if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding lon variable'//&
                     ' (SUBROUTINE get_geo_data)')
               endif
 
!         status = NF90_INQ_DIMID(ncid, "lat", dimIDs(2))
!               if (status /= NF90_NOERR) CALL nc_abort &
!                     (status, 'Error finding lat variable'//&
!                     ' (SUBROUTINE get_geo_data)')
 
         status = NF90_INQ_DIMID(ncid, "lat", dimIDs(2))
               if (status /= NF90_NOERR) then
                  status = NF90_INQ_VARID(ncid, "Latitude",dimIDs(2))
                     if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding lat variable'//&
                     ' (SUBROUTINE get_geo_data)')
               endif
 
         status = NF90_INQ_DIMID(ncid, "time", dimIDs(3))
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error finding time variable'//&
                     ' (SUBROUTINE get_geo_data)')
 
         status = NF90_INQUIRE_VARIABLE(ncid, varid, dimids = dids)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining dimensions in '//&
                     VARNAME//' (SUBROUTINE get_geo_data)')
 
!     Check if first dimension is lat or lon
         print *, dimIDs(1:3)
         print *, dids(1:3)
         if (dids(1) .eq. dimIDs(1))then  
          flip_var = 1
         else
          flip_var = 0
         endif
!     Determine number of longitude, latitude and time points
         status = NF90_INQUIRE_DIMENSION(ncid, dimIDs(1), len = numlon)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining lon dimension '//&
                      ' (SUBROUTINE get_geo_data)')
         status = NF90_INQUIRE_DIMENSION(ncid, dimIDs(2), len = numlat)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining lat dimension '//&
                      ' (SUBROUTINE get_geo_data)') 
         status = NF90_INQUIRE_DIMENSION(ncid, dimIDs(3), len = numt)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error determining time dimension '//&
                      ' (SUBROUTINE get_geo_data)') 

! Delete following line..it was added Only for testing!!! : 
!                   numt = 1

         print *, "Number of longitude points: ", numlon
         print *, "Number of latitude points: ", numlat
         print *, "Number of time points: ", numt
!     Allocate longitude, latitude and time vectors
         if (allocated(lon)) then
            print *, 'Warning: lon was already allocated'
            deallocate(lon)
         endif

         if(allocated(lat)) then
            print *, 'Warning: lat was already allocated'
            deallocate(lat)
         endif

         if(allocated(time)) then
            print * , 'Warning: time was already allocated'
            deallocate(time)
         endif
         allocate(lon(numlon))
         allocate(lat(numlat))
         allocate(time(numt))
!     Retrieve longitude, latitude and time vectors         
         status = NF90_INQ_VARID(ncid, "lon", lonvarid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading lon dimension '//&
                      ' (SUBROUTINE get_geo_data)') 
         status = NF90_INQ_VARID(ncid, "lat", latvarid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading lat dimension '//&
                      ' (SUBROUTINE get_geo_data)') 
         status = NF90_INQ_VARID(ncid, "time", tvarid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading time dimension '//&
                      ' (SUBROUTINE get_geo_data)') 
         status = NF90_GET_VAR(ncid, lonvarid, lon)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading lon variable'//&
                      ' (SUBROUTINE get_geo_data)') 
         status = NF90_GET_VAR(ncid, latvarid, lat)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading lat variable '//&
                      ' (SUBROUTINE get_geo_data)') 
         status = NF90_GET_VAR(ncid, tvarid, time)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading time variable'//&
                      ' (SUBROUTINE get_geo_data)')
         status = NF90_CLOSE(ncid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error closing netcdf file '//&
                      ' (SUBROUTINE get_all_data)') 


!     Check longitude resolution                      
         res_lon = lon(2) - lon(1)
         if (res_lon .le. 0) then
             print *, "Error in lon dimension: lon must be in &
                      increasing order"
            stop
         endif
!     Check latitude resolution
         res_lat = lat(2) - lat(1)
         latord = res_lat
         res_lon = abs(res_lon)
         res_lat = abs(res_lat)
!     Check if lat is cell-centered or node
!     And determine number of points assuming a global data set
         if (mod(90-res_lat/2-lat(1),res_lat) .eq. 0)then
           print *, "Latitude is cell-centered"
           nlat = nint(180/res_lat) 
           cell_c = 1
         elseif (mod(90-lat(1),res_lat).eq. 0)then
           print *, "Latitude is node-based"
           nlat = nint(180/res_lat +1)
           cell_c = 0
         else
           print *, "Warning: resolution in lat and lon are different"
           print *, "Assuming the data is GLOBAL, node-based,&
                    starting at 90"
           nlat = numlat
           nlon = numlon
           cell_c = 2
!          stop
         endif  !If latitude is cell-centered

         if(cell_c .ne. 2)then
         ! Check if we need to reshape the longitudes
         if (minval(lon,1) .lt. 0) then
        ! If lon does not passes by 0 (now 180), then no need to do this
            if (maxval(lon,1) .le. 180)then
              check_lon = 0
            else
         check_lon = 1
         print *, "Lon range detected to be in -180:180 format. &
                   This will be changed to 0:360 format"
         lon = lon + 180  !This only changes the lon vector, but when SUBVAR is retrieved,
                          ! then we will need to rearrange it
 
              !Calculate where the variable is 'split' (where it passes through 0)
              ! this depends if it is cell-centered or node-based
              if (cell_c .eq. 1)then
                 pos_split = 180/res_lon + res_lon/2 - lon(1) +1
              else
                 pos_split = 180/res_lon - lon(1) + 1
              endif
             
             print *, "Position to split: ",pos_split
             
             endif !If lon passes by 0
             else
               check_lon = 0
          endif  !If minimum of lon is less than 0


         nlon = 360/res_lon
         endif ! if cell_ce is not 2 

         print *, "The resolution in lon is  ", res_lon
         print *, "The resolution in lat is ", res_lat
         print *, "Cell_c is ", cell_c
         print *, "For a global data set,"
         print *, "the number of lon points should be ", nlon
         print *, "the number of lat points should be ", nlat
 
!        At this point we have the geo data. Nothing else to be done.
      END SUBROUTINE get_geo_data



!====================================================================
!
! Name: get_all_data
! Purpose: Retrieve the whole variable
! This is used in case of using SPHEREPACK
!     Called from: snake_driver, if using SPHEREPACK ONLY
!     Calls: nc_abort
!
!     Input: name of input file
!     Output: Data is stored in variable VAR
         SUBROUTINE get_all_data(FNAMEI)
           CHARACTER (LEN=*), INTENT(IN) :: FNAMEI

           !Local variables
           integer :: xpos, ypos, xposf, yposf
           real, dimension(:,:,:), allocatable :: TEMPVAR, SUBVAR
           
 
           if(.not. allocated(lons)) ALLOCATE(lons(nlon))
           if(.not. allocated(lats)) ALLOCATE(lats(nlat))

!     Case where it is node-based.. to do: add the cell-centered case
           lats = (/(90-i*res_lat,i=0,nlat-1,1)/)
           lons = (/(0+i*res_lon,i=0,nlon-1,1)/)
!     Determine position of lon in lons and lat in lats
           xpos = lon(1)/res_lon +1      ! assuming the first position in vector will be '1'
           if (latord .gt. 0)then
              ypos = (90-lat(numlat))/res_lat+1 !lat originally in increasing order
           else
              ypos = (90-lat(1))/res_lat +1 ! lat in decreasing order
           endif
           xposf = numlon - xpos + 1
           yposf = numlat - ypos + 1
      !Copy lons and lats to lon lat vectors
           DEALLOCATE(lon)
           DEALLOCATE(lat)
      
           ALLOCATE(lon(nlon))
           ALLOCATE(lat(nlat))

           if (cell_c .ne. 2)then          
            lon=lons              
            lat=lats                 ! assuming the data set is as should be
          else
            xpos = 1            !this if is to bypass the modification of the variable 
            xposf = nlon        ! in case the resolutions are different
            ypos = 1            !assuming the data set is in correct order
            yposf = nlat
          endif

          ALLOCATE(SUBVAR(numlat,numlon,numt))
          if(allocated(VAR)) DEALLOCATE(VAR)
          ALLOCATE(VAR(nlat,nlon,numt))

!     Initialize VAR to zeros
          do i=1,nlat
             do j=1,nlon
                do k=1,numt
                   VAR(i,j,k)=0      
                enddo
             enddo
          enddo
     
         status = NF90_OPEN(path = FNAMEI, mode = nf90_nowrite, &
                  ncid = ncid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error opening input netcdf file'//&
                     ' (SUBROUTINE get_all_data)')
 
          print *, "Reading variable (full) into memory."
          print *, "This should take a while if the size is large. "
          status = NF90_INQ_VARID(ncid, VARNAME, varid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error closing netcdf file '//&
                      ' (SUBROUTINE get_all_data)') 


          status = NF90_GET_ATT(ncid, varid, "_FillValue",miss_value)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error closing netcdf file '//&
                      ' (SUBROUTINE get_all_data)') 


          if (flip_var .eq. 1) then
             !flip_Var =1 then variable nees to be flipped (from lon, lat, time to lat, lon, time)
          print *, "Variable is being read and transposed"
          status = NF90_GET_VAR(ncid, varid, SUBVAR , &
                   start = (/1,1,1/), count = (/numlon, numlat, &
                   numt/), map = (/numlat, 1, numlon*numlat/)) 
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading variable '//&
                      ' (SUBROUTINE get_all_data)') 
          ! by using map = (/numlat, 1, numlon*numlat/). If we wanted to read it as is, it would be
          !          map = (/1, numlat, numlon*numlat/), but in that case map wouldn't have to be specified
          else
          print *, "Variable is being read"
          status = NF90_GET_VAR(ncid, varid, SUBVAR)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error reading variable '//&
                      ' (SUBROUTINE get_all_data)') 
          endif
      
          print *, "Missing values are: ", miss_value
          print *, "Read netCDF file correctly"
          status = NF90_CLOSE(ncid)
               if (status /= NF90_NOERR) CALL nc_abort &
                     (status, 'Error closing netcdf file '//&
                      ' (SUBROUTINE get_all_data)') 

          ALLOCATE(TEMPVAR(numlat,pos_split-1,numt))
          if (check_lon .eq. 1) then
           print *, "Splitting SUBVAR at 180, at position ", pos_split
           TEMPVAR = SUBVAR(:,1:(pos_split-1),:)   ! Store left half
           SUBVAR(:,1:(pos_split-1),:)=SUBVAR(:,pos_split:numlon,:)  !Right-half to left half
           SUBVAR(:,pos_split:numlon,:) = TEMPVAR !stored left half to right-half
         endif

         !Pass SUBVAR to VAR in correct indices
         print *, "Passing SUBVAR to VAR", xpos, xposf, ypos, yposf
         VAR(ypos:yposf,xpos:xposf,:) = SUBVAR
         numlat=nlat
         numlon=nlon

         !Flip lat if needed:


         if(latord .ge. 0) then
         print *, "Lat is in increasing order", latord
         print *, "Data will be flipped in the lat direction"
         VAR = VAR(ubound(VAR,1):lbound(VAR,1):-1,:,:)
         endif
         END SUBROUTINE get_all_data 
      END MODULE snake_data_module
