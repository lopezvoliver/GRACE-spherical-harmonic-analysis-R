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
         integer :: arg1, arg2, arg3, arg4
         character (len =100) :: LMAXARG
         timepatch = 0
        print *, 'Hello world.'
         
         call getarg(1, VARNAME)
         call getarg(2, FNAMEI)
         call getarg(3, FNAMEO)
         call getarg(4, LMAXARG) 

         VARNAME = TRIM(VARNAME)
         FNAMEI = TRIM(FNAMEI)
         FNAMEO = TRIM(FNAMEO)
         LMAXARG = TRIM(LMAXARG)
         
         !Check arguments passed:
         arg1 = LEN_TRIM(VARNAME)
         arg2 = LEN_TRIM(FNAMEI)
         arg3 = LEN_TRIM(FNAMEO)
         if (arg1 .eq. 0) then
            print *, "Program usage:"
            print *, "./snake VARNAME FNAMEI FNAMEO [FNAMES]"
            print *, "VARNAME: Name of variable to look for in &
               netcdf file."
            print *, "FNAMEI: Name of input netcdf file (data or&
                     coefficient)"
            print *, "FNAMEO: Name of output coefficient netcdf file"
            print *, "LMAX: lmax, only for DRAGON analysis "
            stop
         endif

         ! PROGRAM: snakeSHA (only analysis)
         CALL get_geo_data(FNAMEI)
         read (LMAXARG, *) lmax
         CALL DRGSHA
         CALL FLMOUT

      END PROGRAM snake_main
