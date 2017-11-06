!>
!!                 Contains global variables used by the shallow water model.
!!
!!                 More comments on many of the variable defined here
!!                 to be found in the <i>user_introduction</i>.
!!
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
MODULE mo_global_variables
!-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants
  USE mo_physical_constants, ONLY: grav
  USE mo_io_units,           ONLY: nnml, nnml_output, filename_max
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_datetime,           ONLY: t_datetime, proleptic_gregorian,          &
    &                              date_to_time, add_time, print_datetime_all

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC

  CHARACTER(len=*), PARAMETER :: modelname    = 'icon_testbed'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

  LOGICAL, PARAMETER :: ltestcase=.true.

  LOGICAL :: l_zdiffu_t     ! .true.: apply truly horizontal temperature diffusion over steep slopes
  REAL(wp):: thslp_zdiffu   ! threshold slope above which temperature diffusion is applied
  REAL(wp):: thhgtd_zdiffu  ! threshold height difference between adjacent model grid points
                            ! above which temperature diffusion is applied
  
!-------------------------------------------------------------
! !variables setting up the basic configuration of a model run
!-------------------------------------------------------------
  ! parameterized forcing (right hand side) of dynamics and transport
  ! parameterized forcing (right hand side) of dynamics and transport
  ! auxiliary variables
  INTEGER, PARAMETER :: inoforcing  =  0
  INTEGER, PARAMETER :: iheldsuarez =  1
  INTEGER, PARAMETER :: iecham      =  2
  INTEGER, PARAMETER :: inwp        =  3
  INTEGER, PARAMETER :: impiom      = -1
  INTEGER  :: iforcing= inoforcing  ! index of parameterized forcing package
                         ! - positive: for atmosphere
                         ! - negative for ocean
                         ! for non-zero value, read the namelist of
                         ! the specified forcing package and compute
                         ! the dynamical tendencies
                         !  0: no forcing
                         !  1: Held Suarez forcing
                         !  2: forcing by ECHAM parameterizations
                         !  3: forcing by NWP parameterizations
                         ! -1: forcing by MPIOM parameterizations
  INTEGER :: igradp_method = 1 ! Method for computing the horizontal presure gradient
  REAL(wp):: exner_expol    ! Temporal extrapolation of Exner for computation of
                            ! horizontal pressure gradient
  REAL(wp):: rayleigh_coeff ! Rayleigh damping coefficient in w-equation
  REAL(wp):: damp_height    ! height at which damping starts

  LOGICAL :: ltransport=.false.
  
  INTEGER :: ivctype        ! Type of vertical coordinate (Gal-Chen / SLEVE)
!--------------------------------------------------------------------
! flags for the SLEVE coordinate
!--------------------------------------------------------------------

  ! a) Parameters specifying the distrubution of the coordinate surfaces
  !     (the initializations are a workaround for a NEC compiler bug)
  REAL(wp):: min_lay_thckn = 1._wp  ! Layer thickness of lowermost level
  REAL(wp):: top_height    = 1._wp  ! Height of model top

  ! b) Parameters for SLEVE definition
  REAL(wp):: decay_scale_1 = 1._wp  ! Decay scale for large-scale topography component
  REAL(wp):: decay_scale_2 = 1._wp  ! Decay scale for small-scale topography component
  REAL(wp):: decay_exp     = 1._wp  ! Exponent for decay function
  REAL(wp):: flat_height   = 1._wp  ! Height above which the coordinate surfaces are exactly flat
                            ! additional feature not available in the standard SLEVE definition
  REAL(wp):: vwind_offctr   ! Off-centering in vertical wind solver

! time stepping scheme ---------------------------------------------------------

  LOGICAL  :: ltwotime   ! if .TRUE., two time level discretizations are used,
                         ! if .FALSE., three time level discretizations are used

  INTEGER  :: itime_scheme = leapfrog_si ! parameter used to select the time stepping scheme
                            ! = 1, explicit 2 time level scheme
                            ! = 2, semi implicit 2 time level scheme
                            ! = 3, explicit leapfrog
                            ! = 4, leapfrog with semi implicit correction
                            ! = 5, 4-stage Runge-Kutta method
                            ! = 6, SSPRK(5,4) (Runge-Kutta) method
  
  ! time information
  ! ----------------
  !
  ! calendar type
  INTEGER          :: calendar
  !
  ! initial date and time
  ! - namelist variables
  INTEGER          :: ini_year, ini_month, ini_day
  INTEGER          :: ini_hour, ini_minute
  REAL(wp)         :: ini_second
  ! - data and time structure
  TYPE(t_datetime) :: ini_datetime
  !
  ! end date and time
  ! - namelist variables
  INTEGER          :: end_year, end_month, end_day
  INTEGER          :: end_hour, end_minute
  REAL(wp)         :: end_second
  ! - data and time structure
  TYPE(t_datetime) :: end_datetime
  !
  ! run length
  ! - in day,hr,min,sec
  INTEGER          :: run_day
  INTEGER          :: run_hour, run_minute
  REAL(wp)         :: run_second
  ! - in time steps
  REAL(wp)         :: dtime      ! [s] length of a time step
  REAL(wp)         :: dtrk(3)    ! [s] Runge Kutta 3 time steps [s]

  ! computing setup
  ! ---------------

  INTEGER  :: nproma     ! vector length

  ! number of levels (atmosphere or ocean)
  ! ----------------
  INTEGER  :: nlev       ! number of full levels = number of layers
                         ! 0=default: determined by initial/restart file
  INTEGER  :: nlevp1     ! number of half levels = nlev+1

  INTEGER :: nvclev      ! no. of levels at which the coeffs A, B are given

  ! cell geometry
  ! -------------
  INTEGER  :: i_cell_type! =3 if cells are triangles
                         ! =6 if cells are hexagons/pentagons
                         ! =4 if cells are quadrilaterals
  !mag! Here we should consider to use separate grid files for tri and
  !mag! pen/hex cell grids, so that the objects of the grid files
  !mag! have always the same meaning, hence can be use in the same
  !mag! manner. Then i_cell_type could be read directly from the
  !mag! grid file (or later the initial or restart files) and there
  !mag! would be no need to differentiate different interpretations
  !mag! of cell centers,edges and vertices. This would simplify the
  !mag! code and make it at the same time more flexible.

  LOGICAL :: locean
  ! time integration scheme
  ! -----------------------
!!$  INTEGER  :: itimescheme! time integration scheme
!!$                         ! 1: ...
!!$                         ! 2: ...
  !mag! eventually we should add here a control mechanism for selecting
  !mag! the time integration scheme. If different schemes are used at
  !mag! different places, this should assure that only schemes consistent
  !mag! with each other are used.
  !mag! Currently the schemes are selected separately for the dynamics,
  !mag! in &dynamics_ctl, and for transport, in &transport_ctl

  ! processes
  ! ---------
  ! resolved dynamics
  INTEGER  :: itopo      ! flag for topography handling
                         ! 0: corresponds to analytical topography,
                         ! 1: corresponds to ETOPO2_topography, to be extended
                         ! 2: currently corresponds to netcdf files provided by
                         !    Herrmann Ascensio, to be extended

  LOGICAL  :: ltimer     ! if .TRUE.,  the timer is switched on

  ! dump/restore
  !-------------

  INTEGER  ::        &
    & idiv_method = 1     ! 1: Hydrostatic atmospheric model: 
                       !    Gauss integral with original normal velocity components
                       ! 1: Non-Hydrostatic atmospheric model: 
                       !    Gauss integral with averged normal velocity components
                       ! Thus, in a linear equilateral grid, methods 1 and 2 for
                       ! the non-hydrostatic model are the same.
                       ! 2: divergence averaging with bilinear averaging

  REAL(wp) :: divavg_cntrwgt  ! weight of central cell for divergence averaging


  INTEGER  :: n_zlev        ! number of ocean levels
  REAL(wp) :: dzlev_m(100)  ! namelist input of layer thickness
  REAL(wp) :: upstr_beta    ! =1 for 3rd order upstream, =0 for 4th order centered
                            ! theta advection

  INTEGER  :: msg_level=1  ! Determines how much printout is generated during runtime
  


  ! physics parameters
  INTEGER  ::     igbm      = 0
  INTEGER  ::     iwtr      = 1
  INTEGER  ::     iice      = 2
  INTEGER  ::     ilnd      = 3
  INTEGER  ::     nsfc_type = 3
  INTEGER, PARAMETER :: iqv    = 1    !> water vapour
  INTEGER, PARAMETER :: iqc    = 2    !! cloud water
  INTEGER, PARAMETER :: iqi    = 3    !! ice
  INTEGER, PARAMETER :: iqr    = 4    !! rain water
  INTEGER, PARAMETER :: iqs    = 5    !! snow
  INTEGER, PARAMETER :: icc    = 4    !! cloud cover
  ! - other species
  INTEGER, PARAMETER :: io3  = 6     !< O3
  INTEGER, PARAMETER :: iqcond = iqs  !! index of last hydrometeor to ease
                                      !! summation over all of them
  INTEGER, PARAMETER :: ico2 = 7     !< CO2
!   INTEGER  :: ntracer=iqcond
  INTEGER  :: ntracer=io3

  
  CHARACTER(filename_max) :: patch_file_name
  CHARACTER(filename_max) :: vertical_coordinates_file_name



END MODULE mo_global_variables
