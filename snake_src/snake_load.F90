!====================================================================
! Source code by Oliver Lopez
!
! Module name: snake_loadflm
! Purpose: Spherical harmonic analysis and synthesis using SPHEREPACK
! Called by:    snake_driver.F90
!====================================================================
      MODULE snake_load_module
         USE netcdf
         USE snake_common_module

         implicit none
         !netcdf ids
             integer ::  ncid, &
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
  


      END MODULE snake_load_module
