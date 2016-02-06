!====================================================================
! Source code by Oliver Lopez, 2015
!
! Program name: snake_driver.F90
! Purpose: Perform spherical harmonic analysis/synthesis on data
!          read from a netcdf file.
!          Output analysis coefficients and synthesis to separate
!          netcdf files
! Algorithm source: (1) Wang et al. (2006) program dragon.F90
!                   Modified here for netcdf retrieval
!                   and netcdf output
!                   As well as modifications suggested in the same
!                   paper for faster computation of coefficients 
!                   on a 3D variable
!                   (2) SPHEREPACK library
!====================================================================
      PROGRAM snake_main
         USE snake_data_module
         USE dragon_module, only: DRGSHA, DRGSHS
         USE snake_shout_module, only: SHOUT
         USE snake_flmout_module, only: FLMOUT
      
         implicit none
         integer :: arg1, arg2
         timepatch = 0
        print *, 'Hello world.'
         
         call getarg(1, FNAMEI)
         call getarg(2, FNAMES)

         FNAMEI = TRIM(FNAMEI)
         FNAMES = TRIM(FNAMES)
         
         !Check arguments passed:
         arg1 = LEN_TRIM(FNAMEI)
         if (arg1 .eq. 0) then
            print *, "Program usage:"
            print *, "./visGRACE FNAMEI FNAMEO"
            print *, "FNAMEI: Name of input netcdf file"
            print *, "FNAMEO: Name of output netcdf file"
            stop
         endif

         ! PROGRAM: visGRACE  
         ! Description: Convert GRACE-type coefficients before using
         ! dragon synthesis
         CALL get_grace_data(FNAMEI)
         CALL DRGSHS
         CALL SHOUT
        
      END PROGRAM snake_main
