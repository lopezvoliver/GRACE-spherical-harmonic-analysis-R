!====================================================================
! Source code by Oliver Lopez
!
! Module name: snake_shout
! Purpose: Output results from the synthesis 
! Called by:    snake_driver.F90
!====================================================================
      MODULE snake_shout_module
         USE netcdf
         USE snake_common_module
         implicit none
         integer ncidc, londimid, latdimid, tcdimid, fvarid, fbvarid
         integer status
         integer lonvarid, latvarid, tcvarid
         CONTAINS

         SUBROUTINE SHOUT
         ! Output synthesis to netcdf
           print *, "SH synthesis on with output in ", FNAMES
           status = nf90_create(FNAMES, nf90_hdf5, ncidc)
            ! Define dimensions: lon, lat, t
         status = nf90_def_dim(ncidc, "lon", numlon, londimid)
         status = nf90_def_dim(ncidc, "lat", numlat, latdimid)
         status = nf90_def_dim(ncidc, "time", numt, tcdimid)
     ! Define variables: lon, lat, time, f, fback
         status = nf90_def_var(ncidc, "lon", nf90_double, &
                  (/londimid/), lonvarid)
         status = nf90_def_var(ncidc, "lat", nf90_double, &
                  (/latdimid/), latvarid)
         status = nf90_def_var(ncidc, "time", nf90_int, (/tcdimid/), &
                  tcvarid)
         status = nf90_def_var(ncidc, "f", nf90_double, &
               (/latdimid, londimid, tcdimid/), fvarid)
         status = nf90_def_var(ncidc, "fback", nf90_double, &
               (/latdimid, londimid, tcdimid/), fbvarid)
      ! Define fill values
         status = nf90_put_att(ncidc, fvarid, "_FillValue", &
               -9999.0d0)
         status = nf90_put_att(ncidc, fbvarid, "_FillValue", &
               -9999.0d0)
         status = nf90_enddef(ncidc)
      ! Put lon, lat and time values
         status = nf90_put_var(ncidc, lonvarid,lon)
         status = nf90_put_var(ncidc, latvarid,lat)
         status = nf90_put_var(ncidc, tcvarid, time)
        do k=1,numt
            !     Output: f and dreal(fback) 
!      print *, "Writing synthesis output"
         status = nf90_put_var(ncidc, fvarid, VAR(:,:,k), & 
               (/1,1,k/),(/numlat, numlon, 1/)) 
         status = nf90_put_var(ncidc, fbvarid, VARBACK(:,:,k), & 
               (/1,1,k/),(/numlat, numlon, 1/))  

         end do

        status = nf90_close(ncidc)

      END SUBROUTINE SHOUT


      END MODULE snake_shout_module

