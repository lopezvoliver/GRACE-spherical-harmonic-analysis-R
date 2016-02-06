!====================================================================
! Source code by Oliver Lopez
! Module name: snake_common
! Purpose: Declares common data needed by the different modules
!	   for the main program snake_driver.F90
!          These variables are GLOBAL and thus are saved
!====================================================================
      MODULE snake_common_module
         implicit none
         ! Strings
         character (len = 100) , save :: FNAMEI, FNAMEO, FNAMES, VARNAME
         ! Counters and single variables
         integer :: i, j, k, lmax, pos_split
         ! Switches
         integer :: flip_var, &
                    check_lon, &
                    timepatch
         ! Dimensions
         integer , save ::  numlon, &
                     numlat, &
                     numt
         real :: latord, miss_value
         ! 1D vectors (allocatable)
         real, allocatable , save :: lon(:), &
                              lat(:), &
                              time(:)

!                              v(:), sti, stj
         real, allocatable, save :: VARBACK(:,:,:), &  !synthesis
                                    VAR(:,:,:)         !main variable
         double precision, dimension(:,:,:), allocatable :: &
            aout, bout !used by snake_SPK_module and shout
                       !Currently for spherepack only

         !Output for dragon coefficients:
         double complex, allocatable :: flm(:,:,:)

         real gausr     
      END MODULE snake_common_module
