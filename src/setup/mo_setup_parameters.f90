!>
!! @par Revision History
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!
MODULE mo_setup_parameters
!-------------------------------------------------------------------------

  USE mo_global_variables
  USE mo_namelist,            ONLY: open_nml,  close_nml, open_nml_output, close_nml_output, &
    & position_nml, positioned


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setup_run_parameters, setup_sleve



  CONTAINS

 !-------------------------------------------------------------------------
 !>
 !!               Initialization of variables that contain general information
 !!               about the model run. The configuration is read from
 !!               namelist 'run_ctl'.
 !!
 SUBROUTINE setup_run_parameters
 
  
  NAMELIST /run_ctl/  nlev,itopo, igradp_method,    &
      &                rayleigh_coeff,damp_height, ivctype,  &
      &                vwind_offctr, exner_expol,            &
      &                calendar,               &
      &                ini_year, ini_month,  ini_day,        &
      &                ini_hour, ini_minute, ini_second,     &
      &                end_year, end_month,  end_day,        &
      &                end_hour, end_minute, end_second,     &
      &                run_day,                              &
      &                run_hour, run_minute, run_second,     &
      &                dtime,                        &
      &                i_cell_type,                          &
      &                ltimer, msg_level, iforcing, itime_scheme
  
  NAMELIST /grid_ctl/ patch_file_name, vertical_coordinates_file_name

INTEGER  :: i_status
   CHARACTER(*), PARAMETER :: method_name='setup_run_parameters'
!-----------------------------------------------------------------------
! set up the default values

  iforcing= inoforcing  ! index of parameterized forcing package
  itime_scheme=SSPRK54
  igradp_method  = 1
  ! Coefficient for temporal extrapolation for Exner pressure
  exner_expol    = 0.5_wp
  ! Rayleigh damping properties (Klemp,Dudhia,Hassiotis:MWR136,pp.3987-4004)
  rayleigh_coeff = 0.2_wp
  damp_height    = 10000.0_wp
  ! Type of vertical coordinate (1: Gal-Chen, 2: SLEVE)
  ivctype        = 1
  ! offcentering for w in vertical solver
  vwind_offctr   = 0.03_wp

   ! dimensions for new, initialized experiments
   nlev           = 31  ! number of full levels
   itopo          = 0

!    lcorio         = .TRUE.
   ltimer         = .TRUE.

   ! initial date and time
   calendar       = proleptic_gregorian
   ini_year       = 2008
   ini_month      = 9
   ini_day        = 1
   ini_hour       = 0
   ini_minute     = 0
   ini_second     = 0.0_wp
   !
   ! end date and time
   end_year       = 2008
   end_month      = 9
   end_day        = 1
   end_hour       = 1
   end_minute     = 40
   end_second     = 0.0_wp
   
   !
   ! length of integration = (number of timesteps)*(length of timestep)
   ! - If nsteps is set to a non-zero positive value, then the end date is computed
   !   from the initial date and time, the time step dtime, and nsteps.
   ! - Else if run_day, run_hour, run_minute or run_second is set to a non-zero,
   !   positive value, then the initial date and time and the run_... variables are
   !   used to compute the end date and time and, using dtime, nsteps.
   !   Else nsteps is computed from the initial and end date and time and dtime.
   !
   ! initialize run_... variables with zero
   run_day        = 0
   run_hour       = 0
   run_minute     = 0
   run_second     = 0.0_wp
   !
   !
   i_cell_type=3
   ! inner loop length to optimize vectorization
   locean = .false.
  
  ! Settings for icell_type=6
   upstr_beta     = 1.0_wp
   
   l_zdiffu_t     = .FALSE. ! not used by default
   thslp_zdiffu   = 0.025_wp ! slope threshold 0.025
   thhgtd_zdiffu  = 200._wp ! threshold for height difference between adjacent grid points 200 m

   !-------------------------------------------
   ! Open the namelist file.
   CALL open_nml('NAMELIST_ICON_TESTBED')
   
   ! read run_ctl namelist
   CALL position_nml ('run_ctl', status=i_status)
   SELECT CASE (i_status)
   CASE (positioned)
     READ (nnml, run_ctl)
   CASE default
     CALL finish( method_name,'run_ctl not positioned')
   END SELECT

   ! read grid_ctl namelist
   CALL position_nml ('grid_ctl', status=i_status)
   SELECT CASE (i_status)
   CASE (positioned)
     READ (nnml, grid_ctl)
   CASE default
     CALL finish( method_name,'grid_ctl not positioned')
   END SELECT
   
  ! Close the namelist files.
  CALL close_nml

   
  SELECT CASE (itopo)
  CASE (0,2)
    ! ok
  CASE default
     CALL finish( method_name,'wrong topography specifier, itopo must be in {0,1,2}]')
  END SELECT


   ! ----------------------------------
   ! check the validity of the settings
   ! ----------------------------------

   IF (nlev < 1)  CALL finish(TRIM(method_name),'"nlev" must be positive')
   nlevp1 = nlev+1
   nvclev = nlevp1

   ! time step
!    IF (dtime  <= 0._wp) CALL finish(method_name,'"dtime" must be positive')
! 
!    ! initial date and time
   ini_datetime%calendar = calendar
   ini_datetime%year     = ini_year
   ini_datetime%month    = ini_month
   ini_datetime%day      = ini_day
   ini_datetime%hour     = ini_hour
   ini_datetime%minute   = ini_minute
   ini_datetime%second   = ini_second
   CALL date_to_time      (ini_datetime) ! fill date time structure
   IF (msg_level > 8) THEN
     CALL message(' ',' ')
     CALL message(method_name,'Initial date and time')
     CALL message(method_name,'---------------------')
     CALL print_datetime_all(ini_datetime)  ! print all date and time components
   ENDIF
! 
!    ! end and length of integration
!    IF (nsteps/=0) THEN
!      IF (nsteps < 0    ) CALL finish(method_name,'"nsteps" must not be negative')
!      length_sec   = REAL(nsteps,wp)*dtime
!      end_datetime = ini_datetime
!      CALL add_time(length_sec,0,0,0,end_datetime)
!      !
!    ELSE IF (run_day/=0 .OR. run_hour/=0 .OR. run_minute/=0 .OR. run_second/=0.0_wp) THEN
!      IF (run_day    < 0    ) CALL finish(method_name,'"run_day" must not be negative')
!      IF (run_hour   < 0    ) CALL finish(method_name,'"run_hour" must not be negative')
!      IF (run_minute < 0    ) CALL finish(method_name,'"run_minute" must not be negative')
!      IF (run_second < 0._wp) CALL finish(method_name,'"run_second" must not be negative')
!      !
!      end_datetime = ini_datetime
!      CALL add_time(run_second,run_minute,run_hour,run_day,end_datetime)
!      !
!      ini_datetime_calsec    = (REAL(ini_datetime%calday,wp)+ini_datetime%caltime) &
!        &                      *REAL(ini_datetime%daylen,wp)
!      end_datetime_calsec    = (REAL(end_datetime%calday,wp)+end_datetime%caltime) &
!        &                      *REAL(end_datetime%daylen,wp)
!      nsteps=INT((end_datetime_calsec-ini_datetime_calsec)/dtime)
!      !
!    ELSE
!      ! compute nsteps from ini_datetime, end_datetime and dtime
!      end_datetime%calendar = calendar
!      end_datetime%year     = end_year
!      end_datetime%month    = end_month
!      end_datetime%day      = end_day
!      end_datetime%hour     = end_hour
!      end_datetime%minute   = end_minute
!      end_datetime%second   = end_second
!      CALL date_to_time      (end_datetime) ! fill date time structure
!      !
!      ini_datetime_calsec    = (REAL(ini_datetime%calday,wp)+ini_datetime%caltime) &
!        &                      *REAL(ini_datetime%daylen,wp)
!      end_datetime_calsec    = (REAL(end_datetime%calday,wp)+end_datetime%caltime) &
!        &                      *REAL(end_datetime%daylen,wp)
!      IF (end_datetime_calsec < ini_datetime_calsec) &
!        & CALL finish(method_name,'The end date and time must not be before the initial date and time')
!      !
!      nsteps=INT((end_datetime_calsec-ini_datetime_calsec)/dtime)
!      !
!    END IF
!    !
!    CALL message(' ',' ')
!    CALL message(method_name,'End date and time')
!    CALL message(method_name,'-----------------')
!    CALL print_datetime_all(end_datetime)  ! print all date and time components
!    !
!    CALL message(' ',' ')
!    CALL message(method_name,'Length of run')
!    CALL message(method_name,'-------------')
!    WRITE (message_text,'(a,f7.2)') 'dtime [s] :',dtime
!    CALL message(method_name,message_text)
!    WRITE (message_text,'(a,i7)')   'nsteps    :',nsteps
!    CALL message(method_name,message_text)
!    CALL message(' ',' ')

  SELECT CASE (i_cell_type)
  CASE (3,6)
     ! ok
  CASE default
     CALL finish( TRIM(method_name),'wrong cell type specifier, i_cell_type must be 3 or 6')
  END SELECT

   IF (upstr_beta > 1.0_wp .OR. upstr_beta < 0.0_wp) THEN
     CALL finish(method_name, 'upstr_beta out of range 0..1')
   ENDIF

 END SUBROUTINE setup_run_parameters
 !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
 !>
 !!  Initialization of the SLEVE coordinate namelist
 !!
 !! @par Revision History
 SUBROUTINE setup_sleve

  NAMELIST /sleve_ctl/ min_lay_thckn, top_height, decay_scale_1, &
                       decay_scale_2, decay_exp, flat_height

! !local variable
  INTEGER :: i_status

!-----------------------------------------------------------------------

  ! set up the default values

  ! a) Parameters determining the distribution of model layers
  !    (if not read in from a table)
  min_lay_thckn   = 50._wp      ! Layer thickness of lowermost layer
  top_height      = 23500._wp   ! Height of model top

  ! b) Parameters setting up the decay function of the topographic signal
  decay_scale_1   = 4000._wp    ! Decay scale of large-scale topography component
  decay_scale_2   = 2500._wp    ! Decay scale of small-scale topography component
  decay_exp       = 1.2_wp      ! Exponent for decay function
  flat_height     = 16000._wp   ! Height above which the coordinate surfaces are flat

  ! read namelist

  CALL position_nml ('sleve_ctl', status=i_status)
  SELECT CASE (i_status)
  CASE (positioned)
     READ (nnml, sleve_ctl)
  END SELECT

  ! write the contents of the namelist to an ASCII file

  IF(p_pe == p_io) WRITE(nnml_output,nml=sleve_ctl)

 END SUBROUTINE setup_sleve
 !-------------------------------------------------------------------------

END MODULE mo_setup_parameters
