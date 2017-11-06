!>
!! Driver for the nonhydrostatic solver.
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
MODULE mo_nh_driver
!-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_namelist,            ONLY: open_nml, close_nml, &
    & open_nml_output, close_nml_output, &
    & position_nml, positioned
  USE mo_io_units,            ONLY: nnml, filename_max, find_next_free_unit
  USE mo_nonhydro_state,      ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics, &
                                    construct_nh_state, destruct_nh_state,   &
                                    nnow, nnew, nnow_rcf, nnew_rcf, nsav1, nsav2
  USE mo_nh_global_variables
  USE mo_ext_data,            ONLY: ext_data
  USE mo_model_domain,        ONLY: t_patch, n_dom, n_dom_start,  &
    & global_patch_array, global_int_state_array
  USE mo_interpolation,       ONLY: t_int_state, cells2edges_scalar
                                   
  USE mo_interpolation,       ONLY: t_int_state
  USE mo_datetime,            ONLY: t_datetime, print_datetime, add_time
  USE mo_timer,               ONLY: timer_solve_nh, timer_nh_driver, timer_start, timer_stop
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH,iphysproc
  USE mo_physical_constants,  ONLY: cvd, cpd, cvd_o_rd, rd_o_cpd, p0ref, grav, rd, re, &
    & omega
  USE mo_math_constants,      ONLY: pi, pi_2
  USE mo_math_utilities,      ONLY: gc2cc, t_cartesian_coordinates, &
                                   t_geographical_coordinates, &
                                   arc_length
  USE mo_math_operators,      ONLY: nabla2_scalar, nabla2_vec, div_avg
!   USE mo_vector_operations,   ONLY: covariant_velocities,       &
!                                     contravariant_vorticities,  &
!                                     orthogonal_vorticities,     &
!                                     helicity_tendencies,        &
!                                     kinetic_energy
  USE mo_solve_nonhydro,      ONLY: solve_nh
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_mpi,                 ONLY: p_nprocs
!   USE mo_sync,                ONLY: global_sum_array, sync_patch_array_mult, SYNC_C, &
!                                     push_glob_comm, pop_glob_comm, global_max
!   USE mo_communication,       ONLY: start_delayed_exchange, do_delayed_exchange
  USE mo_vertical_grid,       ONLY: set_nh_metrics, init_sleve_coord, init_hybrid_coord
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e

  IMPLICIT NONE

  PRIVATE

! !DEFINED PARAMETERS for jablonowski williamson:  
  REAL(wp), PARAMETER :: eta0  = 0.252_wp ! 
  REAL(wp), PARAMETER :: etat  = 0.2_wp   ! tropopause
  REAL(wp), PARAMETER :: ps0   = 1.e5_wp  ! surface pressure (Pa)
  REAL(wp), PARAMETER :: u0    = 35._wp   ! maximum zonal wind (m/s)
  REAL(wp), PARAMETER :: temp0 = 288._wp  ! horizontal-mean temperature 
                                          ! at surface (K)
  REAL(wp), PARAMETER :: gamma = 0.005_wp ! temperature elapse rate (K/m)
  REAL(wp), PARAMETER :: dtemp = 4.8e5_wp ! empirical temperature difference (K)

  REAL(wp), PARAMETER :: lonC  = pi/9._wp ! longitude of the perturb. centre 
  REAL(wp), PARAMETER :: latC  = 2._wp*lonC !latitude of the perturb. centre
  
! !DEFINED PARAMETERS for mountain induced Rossby wave train:
  REAL(wp), PARAMETER :: pres_sp  = 93000.0_wp  !pressure surface at the south pole
  REAL(wp), PARAMETER :: temp_mrw = 288._wp     !temperature of isothermal atmosphere

  REAL(wp) :: mount_height           ! (m)
  REAL(wp) :: layer_thickness        ! (m)
  REAL(wp) :: nh_brunt_vais          ! (1/s)
  REAL(wp) :: nh_u0                  ! (m/s)
  REAL(wp) :: nh_t0                  ! (K)
  REAL(wp) :: jw_up                  ! amplitude of the u-perturbation (m/s), jabw  
  REAL(wp) :: u0_mrw                 ! (m/s) wind speed for mrw case 
  REAL(wp) :: mount_height_mrw       ! (m) maximum mount height in mrw and mwbr
  REAL(wp) :: mount_half_width       ! (m) half width of mountain in mrw, mwbr and bell
  REAL(wp) :: mount_lonctr_mrw_deg   ! (deg) lon of mountain center in mrw and mwbr
  REAL(wp) :: mount_latctr_mrw_deg   ! (deg) lat of mountain center in mrw and mwbr
  REAL(wp) :: p_int_mwbr_const       ! pressure at interface in mwbr_const test case
  REAL(wp) :: temp_i_mwbr_const      ! temp in isothermal lower layer in mwbr_const
  REAL(wp) :: bruntvais_u_mwbr_const ! brunt vaisala freq in upper layer in mwbr_const
  REAL(wp) :: u0_mwbr_const          ! (m/s) wind speed for mwbr_const case
  REAL(wp) :: rotate_axis_deg        ! (deg) rotation angle
  REAL(wp) :: torus_domain_length    ! (m) length of domain the slice (torus) grid
  LOGICAL  :: lhs_nh_vn_ptb          ! if true, random noise is added to vn in HS_nh test case
  LOGICAL  :: lhs_fric_heat          ! if true, frictional heating is switched on in HS_nh
  REAL(wp) :: hs_nh_vn_ptb_scale     ! amplitude of the random noise
  REAL(wp) :: rh_at_1000hpa          ! relative humidity at 1000 hPa [%]
  REAL(wp) :: qv_max                 ! limit of maximum specific humidity in the tropics [kg/kg]
  INTEGER :: n_flat_level

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  ! additional flow control variables that need to be dimensioned with the
  ! number of model domains
  LOGICAL, ALLOCATABLE :: lstep_adv(:)   ! determines whether tracer continuity equations
                                         ! should be integrated (.true.) or not (.false.)

  LOGICAL, ALLOCATABLE :: lcall_phy(:,:) ! contains information which physics package
                                         ! must be called at the current timestep
                                         ! and on the current domain.

  REAL(wp), ALLOCATABLE :: t_elapsed_phy(:,:)  ! time (in s) since the last call of
                                               ! the corresponding physics package

  LOGICAL, ALLOCATABLE :: linit_slowphy(:) ! determines whether slow physics has already been initialized

  LOGICAL, ALLOCATABLE :: linit_dyn(:)  ! determines whether dynamics has already been initialized

  INTEGER :: ntl_module

  TYPE(t_nh_state),   TARGET :: patches_nh_state(n_dom)
  
  PUBLIC ::  nh_solver_driver, construct_nh_solver, destruct_nh_solver

  INTEGER :: nh_iterations
  
  CONTAINS

  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE construct_nh_solver ()
    

    INTEGER :: ntl

    !-----------------------------------------------------------------------
    ! for the split explict scheme, ntl is always 2
    ntl_module = 2
    ntl = ntl_module

    CALL setup_nh_parameters()
    
    CALL init_hybrid_coord()
    
    CALL construct_nh_state(global_patch_array, patches_nh_state, ntl)

    IF (ltestcase) THEN
      CALL init_nh_testtopo(global_patch_array)    ! set analytic topography
    ENDIF

    CALL set_nh_metrics(global_patch_array, patches_nh_state, global_int_state_array)

    IF (ltestcase) THEN
      CALL init_nh_testcase(global_patch_array, patches_nh_state, global_int_state_array, ntl)
    ENDIF


  END SUBROUTINE construct_nh_solver
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE destruct_nh_solver ()
    
    CALL destruct_nh_state(patches_nh_state)

  END SUBROUTINE destruct_nh_solver
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE nh_solver_driver(datetime)
    TYPE(t_datetime), INTENT(inout)  :: datetime

    INTEGER :: i
    
    IF (ltimer) CALL timer_start(timer_nh_driver)
    DO i=1,nh_iterations
!       write(0,*) "perform_nh_stepping ", i
      CALL perform_nh_stepping(global_patch_array, &
                             global_int_state_array, &
                             patches_nh_state, &
                             datetime)
    ENDDO
    IF (ltimer) CALL timer_stop(timer_nh_driver)
    
  END SUBROUTINE nh_solver_driver
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !>
  !! Organizes nonhydrostatic time stepping.
  !!
  !! Organizes nonhydrostatic time stepping
  !! Currently we assume to have only one grid level.
  SUBROUTINE perform_nh_stepping (p_patch, p_int_state,  p_nh_state, &
                                  datetime)
!
  TYPE(t_patch), TARGET, INTENT(IN)            :: p_patch(n_dom_start:n_dom)
  TYPE(t_int_state), TARGET, INTENT(IN)        :: p_int_state(n_dom_start:n_dom)

  TYPE(t_nh_state), TARGET, INTENT(INOUT):: p_nh_state(n_dom)
  TYPE(t_datetime), INTENT(INOUT)      :: datetime

  REAL(wp)                             :: sim_time(n_dom)
  INTEGER                              :: jfile, jstep, jb, nlen, jg
  INTEGER                              :: ist
  REAL(wp)                             :: vnmax, wmax
  REAL(wp)                             :: vn_aux(p_patch(1)%nblks_int_e)
  REAL(wp)                             :: w_aux(p_patch(1)%nblks_int_c)
  REAL(wp), DIMENSION(:,:,:), POINTER  :: p_vn, p_w
  INTEGER :: nblks_c, blk, nn

  CHARACTER(len=MAX_CHAR_LENGTH) :: method_name = 'perform_nh_stepping'
  
!-----------------------------------------------------------------------

  sim_time(:) = 0._wp

  jfile = 1

  ! allocate flow control variables for transport and slow physics calls
  ALLOCATE(lstep_adv(n_dom),lcall_phy(n_dom,iphysproc),linit_slowphy(n_dom), &
    &      linit_dyn(n_dom),t_elapsed_phy(n_dom,iphysproc), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( 'mo_nh_stepping: perform_nh_stepping',           &
    &      'allocation for flow control variables failed' )
  ENDIF
  ! initialize
  t_elapsed_phy(:,:) = 0._wp
  linit_slowphy(:)   = .TRUE.
  linit_dyn(:)       = .TRUE.

!   TIME_LOOP: DO jstep = 1, nsteps

    CALL add_time(nh_dtime,0,0,0,datetime)

    ! Store first old exner pressure
    ! (to prepare some kind of divergence damping, or to account for
    ! physically based 'implicit weights' in forward backward time stepping)
    nn=nnow(1)
    IF (nn <= 0) CALL finish(TRIM(method_name),"nnow(1) <= 0")
!     write(0,*) "nnow(1)=",nnow(1)
    IF (jstep == 1) THEN
      DO jg = 1, n_dom
        nblks_c = p_patch(jg)%nblks_c
!$OMP PARALLEL
!$OMP DO PRIVATE(blk)
        DO blk=1,nblks_c
          p_nh_state(jg)%diag%exner_old(:,:,blk)=&
           & p_nh_state(jg)%prog(nn)%exner(:,:,blk)
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ENDDO
    ENDIF

!     CALL message(method_name, "Call dynamics_integration...")
    ! dynamics stepping
    ! GZ: will include physics tendencies in the future, so once we should rename this routine
    CALL dynamics_integration(p_nh_state, p_patch, p_int_state,  &
                              1, jstep, nh_dtime, sim_time, 1)

!     CALL message(method_name, "Call dynamics_integration done")
    ! output of results
    ! note: nnew has been replaced by nnow here because the update

!   ENDDO TIME_LOOP

  !
  ! deallocate flow control variables
  !
  DEALLOCATE( lstep_adv, lcall_phy, linit_slowphy, linit_dyn, t_elapsed_phy, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( 'mo_nh_stepping: perform_nh_stepping',          &
      &    'deallocation for lstep_adv, lcall_phy,' //            &
      &    't_elapsed_phy failed' )
  ENDIF


  END SUBROUTINE perform_nh_stepping

!-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !>
  !! dynamics_integration
  !!
  !! Performs dynamics time stepping:  Rotational modes (helicity bracket) and
  !! divergent modes (Poisson bracket) are splitted using Strang splitting.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-08-25)
  !! Adaptation for grid refinement by Guenther Zaengl, DWD (2010-02-09)
  !! Modification by Daniel Reinert, DWD (2010-04-15)
  !!  - Implementation of tracer transport
  !! Modification by Daniel Reinert, DWD (2010-07-23)
  !!  - optional reduced calling frequency for transport and physics
  !!
  RECURSIVE SUBROUTINE dynamics_integration (p_nh_state, p_patch, p_int_state,  &
  &        jg, nstep_global, dt_loc, sim_time, num_steps)

    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch(n_dom_start:n_dom)    !< patch
    TYPE(t_int_state),TARGET,INTENT(in)  :: p_int_state(n_dom_start:n_dom)!< interpolation state
    TYPE(t_nh_state), TARGET, INTENT(inout) :: p_nh_state(n_dom) !< nonhydrostatic state

    INTEGER, INTENT(IN)     :: jg           !< current grid level
    INTEGER, INTENT(IN)     :: nstep_global !< counter of global time step
    INTEGER, INTENT(IN)     :: num_steps    !< number of time steps to be executed
    REAL(wp), INTENT(IN)    :: dt_loc       !< time step applicable to local grid level
    REAL(wp), INTENT(INOUT) :: sim_time(n_dom) !< elapsed simulation time on each
                                               !< grid level

    ! Local variables

    ! Time levels
    INTEGER :: n_now_grf, n_now, n_new, n_save, n_temp
    INTEGER :: n_now_rcf, n_new_rcf, n_upt_rcf  ! accounts for reduced calling frequencies (rcf)

    INTEGER :: jstep, jgp, jgc, jn
    INTEGER :: nsteps_nest ! number of time steps executed in nested domain

    REAL(wp):: dt_sub, rdt_loc
    REAL(wp)::  rdt_rcf   ! time step for advection and fast physics

    REAL(wp), DIMENSION(:,:,:), POINTER  :: p_vn   => NULL()

    LOGICAL, PARAMETER :: l_straka=.FALSE.
    LOGICAL :: l_predictor
    LOGICAL :: l_bdy_nudge
    LOGICAL :: lclean_mflx   ! for reduced calling freqency: determines whether
                             ! mass-fluxes and trajectory-velocities are reset to zero
                             ! i.e. for starting new integration sweep

    ! Switch to determine manner of OpenMP parallelization in interpol_scal_grf
    LOGICAL :: lpar_fields=.FALSE.

    ! Switch to determine if nested domains are called at a given time step
    LOGICAL :: l_call_nests

  CHARACTER(len=MAX_CHAR_LENGTH) :: method_name = 'dynamics_integration'
!$  INTEGER :: num_threads_omp, omp_get_max_threads

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    ! settings for calling frequency for slow physics
    !--------------------------------------------------------------------------

    IF (jg == 1 .AND. linit_dyn(jg)) THEN

      nsav1(1:n_dom) = nnow(1:n_dom)
    ENDIF

    !--------------------------------------------------------------------------

!$  num_threads_omp = omp_get_max_threads()

    ! Determine parent domain ID
    IF ( jg > 1) THEN
      jgp = p_patch(jg)%parent_id
    ELSE IF (n_dom_start == 0) THEN
      jgp = 0
    ELSE
      jgp = 1
    ENDIF

    ! This executes one time step for the global domain and two steps for nested domains
    DO jstep = 1, num_steps


      ! update several switches which decide upon
      ! - calling transport running with rcf (lstep_adv)
      ! - re-initializing temporary transport fields for rcf (lclean_mflx)
      ! - switching order of operators in case of Marchuk-splitting



      ! Set local variables for time levels
      n_now  = nnow(jg)
      n_new  = nnew(jg)

      ! Set local variable for rcf-time levels
      n_now_rcf = nnow_rcf(jg)
      n_new_rcf = nnew_rcf(jg)
      ! the next time level is essential for physics packages, which are not
      ! synchronized with transport (i.e. which are not called for each advection
      ! step or for each ith advection step). Those unsynchronized physics-routines
      ! need to read from and update to timelevel n_upt_rcf and NOT n_new_rcf !!
      IF (lstep_adv(jg) ) THEN
        n_upt_rcf = nnew_rcf(jg)
      ELSE
        n_upt_rcf = nnow_rcf(jg)
      ENDIF


      l_bdy_nudge = .FALSE.
!       CALL message(method_name, "Call solve_nh...")
      IF (ltimer) CALL timer_start(timer_solve_nh)
      CALL solve_nh(p_nh_state(jg), p_patch(jg), p_int_state(jg),   &
                    n_now, n_new, linit_dyn(jg), l_bdy_nudge, dt_loc)
      IF (ltimer) CALL timer_stop(timer_solve_nh)
!       CALL message(method_name, "Call solve_nh done")
      
      ! counter for simulation time in seconds
      sim_time(jg) = sim_time(jg) + dt_loc

      ! Finally, switch between time levels now and new for next time step
      n_temp   = nnow(jg)
      nnow(jg) = nnew(jg)
      nsav1(jg) = nnow(jg)
      nnew(jg) = n_temp

      ! Special treatment for processes (i.e. advection) which can be treated with
      ! reduced calling frequency. Switch between time levels now and new immediately
      ! AFTER the last transport timestep.
      IF (lstep_adv(jg)) THEN
        n_temp       = nnow_rcf(jg)
        nnow_rcf(jg) = nnew_rcf(jg)
        nnew_rcf(jg) = n_temp
      ENDIF

    ENDDO

  END SUBROUTINE dynamics_integration


  !-------------------------------------------------------------------------
  !>
  !! Initialize topography for nonhydrostatic artificial testcases.
  !! 
  !! Initialize topography for nonhydrostatic artificial testcases
  !! 
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-03-19)
  !! Modification by Daniel Reinert, DWD (2010-07-15)
  !! - moved initialization of topography into new subroutine 
  !!   init_nh_testtopo, which is called after the domain-decomposition.
  !!   (because of possible conflicts with the external-data type)
  !! 
  SUBROUTINE init_nh_testtopo (p_patch)
!
! !INPUT VARIABLES:
  TYPE(t_patch),TARGET,  INTENT(INOUT) :: p_patch(n_dom)

  
  INTEGER        :: jg, jc, jv, jb, nlen
  REAL(wp)       :: z_lon, z_lat, z_dist
  REAL(wp)       :: zr, zexp
  TYPE(t_geographical_coordinates) :: z_x2_geo
  TYPE(t_cartesian_coordinates)    :: z_x1_cart, z_x2_cart
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &  routine = '(mo_nh_testcases) init_nh_testtopo:'
 
  REAL(wp)      :: zsiny, zcosy,tmp1,tmp2,tmp3
  REAL(wp)      :: z_lon_ctr, z_lat_ctr, z_fac1, z_fac2

!-----------------------------------------------------------------------

    ! default values

  ! Initialize topography to zero if idealized topo is used
  DO jg = 1, n_dom
    ext_data(jg)%atm%topography_c(1:nproma,1:p_patch(jg)%nblks_int_c) = 0.0_wp
    ext_data(jg)%atm%topography_v(1:nproma,1:p_patch(jg)%nblks_int_v) = 0.0_wp
  ENDDO

!  CASE ('bell')

    ! At present the mountain is at position lat=0,lon=0 (given in meters)
    z_x2_geo%lon = 0.0_wp
    z_x2_geo%lat = 0.0_wp

    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_int_c
        IF (jb /=  p_patch(jg)%nblks_int_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF
        DO jc = 1, nlen
            z_x1_cart    = gc2cc(p_patch(jg)%cells%center(jc,jb))
            z_x2_cart    = gc2cc(z_x2_geo)
            z_dist       = arc_length(z_x1_cart,z_x2_cart)
          ext_data(jg)%atm%topography_c(jc,jb) = mount_height/ &
                 (1.0_wp+ (z_dist/mount_half_width)**2)**1.5_wp
        ENDDO
      ENDDO 
    ENDDO 
    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_int_v
        IF (jb /=  p_patch(jg)%nblks_int_v) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_v
        ENDIF
        DO jv = 1, nlen
            z_x1_cart    = gc2cc(p_patch(jg)%verts%vertex(jv,jb))
            z_x2_cart    = gc2cc(z_x2_geo)
            z_dist       = arc_length(z_x1_cart,z_x2_cart)
          ext_data(jg)%atm%topography_v(jv,jb) = mount_height/ &
                 (1.0_wp+ (z_dist/mount_half_width)**2)**1.5_wp
        ENDDO
      ENDDO 
    ENDDO



  END SUBROUTINE init_nh_testtopo

  
  !-------------
  !>
  !! SUBROUTINE nh_prog_add_random
  !! Add random perturbation to the normal wind.
  !!
  !! @par Revision History
  !!  Based in  hydro_state_prog_add_random (by Hui)
  !!  Adapted to the Non-Hyd. core by Pilar Ripodas, DWD, 2010-09
  SUBROUTINE nh_prog_add_random(p_patch, & ! in
                                         p_nh_prog,  & ! inout
                                         pscale, nproma) ! input

    TYPE(t_patch)            :: p_patch
    TYPE(t_nh_prog) :: p_nh_prog

    REAL(wp), INTENT(IN) :: pscale ! magnitude of the perturbation

    ! LOCAL VARIABLES

    INTEGER :: jk
    INTEGER :: nproma, npromz, nblks

    INTEGER :: DateTimeArray(8)    ! Holds the date and time

    INTEGER :: seed_size, js, seed_trigger, ist

    INTEGER, ALLOCATABLE :: seed_array(:)

    REAL(wp) :: zrand(nproma,p_patch%nblks_e)

    CHARACTER(len=MAX_CHAR_LENGTH) :: string
    !-----
    CALL message('','=========== generating random number =============')

    !-----------------------------------------------------------
    ! 1. prepare memory for the seed
    !-----------------------------------------------------------

    CALL RANDOM_SEED(SIZE=seed_size)
    WRITE(string,*) 'The size of the intrinsic seed is', seed_size
    CALL message('',TRIM(string))

    ALLOCATE( seed_array(seed_size), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish('random number:','allocation of seed_array failed')
    ENDIF

    !----------------------------------------------------------
    ! 2. get seed
    !----------------------------------------------------------
    ! First inquire the current date and time.
    ! The output of the DATE_AND_TIME routine contains 8 elements:
    !   1: year (e.g. 2008)
    !   2: month
    !   3: day
    !   4: time difference with UTC in minutes
    !   5: hour
    !   6: minute
    !   7: seconds
    !   8: milliseconds

    CALL DATE_AND_TIME( VALUES=DateTimeArray )

    seed_trigger =   DateTimeArray(2)*100000000 + DateTimeArray(3)*1000000 &
                 & + DateTimeArray(5)*10000     + DateTimeArray(6)*100     &
                 & + DateTimeArray(7)

    WRITE(string,*) 'the seed trigger is', seed_trigger
    CALL message('',TRIM(string))

    DO js=1,seed_size
       seed_array(js)=ABS(seed_trigger)+(js-1)
    ENDDO

    CALL message('','Seed generated')

    !-----------------------------------------------------------
    ! 3. generate random numbers and perturb the normal wind
    !-----------------------------------------------------------

    CALL RANDOM_SEED( PUT=seed_array )
    CALL RANDOM_NUMBER( zrand )

    zrand = (zrand - 0.5_wp)*pscale

    nblks = p_patch%nblks_e
    npromz = p_patch%npromz_e

    ! add the same random noise to all vertical layers
    DO jk=1,nlev
       p_nh_prog%vn(:,jk,1:nblks-1)    = p_nh_prog%vn(:,jk,1:nblks-1)    + zrand(:,1:nblks-1)
       p_nh_prog%vn(1:npromz,jk,nblks) = p_nh_prog%vn(1:npromz,jk,nblks) + zrand(1:npromz,nblks)
    ENDDO

    !-----------------------------------------------------------
    ! 4. clean up
    !-----------------------------------------------------------

    DEALLOCATE( seed_array, STAT=ist)
    IF(ist/=SUCCESS) CALL finish('random number:','deallocation of seed_array failed')
    CALL message('','=========================================')

  END SUBROUTINE nh_prog_add_random

  !-------------------------------------------------------------------------
  !>
  !! Defines nonhydrostatic artificial initial conditions.
  !! 
  !! Initializes meteorological fields
  !! 
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-04-14)
  !! 
  SUBROUTINE init_nh_testcase (p_patch, p_nh_state, p_int, ntl)
!
! !INPUT VARIABLES:
  TYPE(t_patch),TARGET,  INTENT(INOUT) :: p_patch(n_dom)
  TYPE(t_int_state), INTENT(IN) :: p_int(n_dom)
  INTEGER :: ntl
  TYPE(t_nh_state), TARGET, INTENT(INOUT):: p_nh_state(n_dom)

  INTEGER        :: jg, je, jc, jb, jk, jt, jjt, jn, ji, niter, &
                    nlen, nblks_e, npromz_e,  nblks_c, npromz_c
  INTEGER        :: i_startidx, i_endidx, i_startblk, icount

  REAL(wp), DIMENSION(nproma) ::  &
             z_lat,z_siny,z_cosy, z_fac1, z_fac2, zeta_old, zcoszetav, &
             zsinzetav, z_tavg, z_favg, z_geopot, z_temp, z_fun, z_fund, &
             zeta, zu, zv, z_lon, z_exp, za, zb, zc, z_temp_kp1, z_fac3
  REAL(wp), ALLOCATABLE :: zeta_v(:,:,:)
  REAL(wp), ALLOCATABLE :: zeta_v_e(:,:,:)
  REAL(wp), ALLOCATABLE :: z_wsfc_e(:,:,:), z_wsfc_c(:,:,:)
  REAL(wp), ALLOCATABLE :: z_int_c(:,:)    ! z at interface in mwbr_const
  TYPE(t_nh_state), POINTER       :: p_nhdom
  REAL(wp)              :: zlat, zlon, z_u, bruntvaissq, kappa, zhelp1, &
                           zcoslat, zhelp2, z_pres, z_sfc, z_nlev
  REAL(wp)              :: zhelp1_i, zhelp1_u, zhelp2_i, bruntvaissq_i, &
                           bruntvaissq_u, zhelp3, zhelp4, rkappa
  REAL(wp)              :: theta_v_int !potential temp at the interface in mwbr_const

  REAL(wp) :: z_help

  REAL(wp) :: zsqv
  
  REAL(wp) :: zrhf,z_1_o_rh
  
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine =  &
                                   '(mo_nh_testcases) init_nh_testcase:' 
! Tracer related variables
  CHARACTER(LEN=1) :: ctracer
!-----------------------------------------------------------------------


  IF (msg_level > 10 ) THEN
      CALL message(TRIM(routine),'Jablonowski test')
    IF ( iforcing == inwp ) THEN
      CALL message(TRIM(routine),' iforcing == inwp')
    ELSE
      CALL message(TRIM(routine),'Attention: iforcing /= inwp')
    ENDIF
  ENDIF
  
  DO jg = 1, n_dom

    ALLOCATE (zeta_v(nproma,nlev,p_patch(jg)%nblks_c), &
              zeta_v_e(nproma,nlev,p_patch(jg)%nblks_e) )
    
    zeta_v    = 0._wp
    zeta_v_e  = 0._wp
    nblks_c   = p_patch(jg)%nblks_int_c
    npromz_c  = p_patch(jg)%npromz_int_c
    nblks_e   = p_patch(jg)%nblks_int_e
    npromz_e  = p_patch(jg)%npromz_int_e
    p_nhdom  => p_nh_state(jg)
    p_nhdom%diag%pres_sfc(:,:) = 100000._wp     !set surface pressure to 1000. hPa


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jn,jt,z_lat,z_siny,z_cosy,z_fac1,z_fac2,z_exp,zeta_old,&
!$OMP            zcoszetav,zsinzetav,z_tavg,z_favg,z_geopot,z_temp,z_fun,z_fund,zeta )
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
!       write(0,*) "jb=",jb
      DO jk = nlev, 1, -1
        DO jc = 1, nlen
          z_lat(jc) = p_patch(jg)%cells%center(jc,jb)%lat
          z_siny(jc) = SIN(z_lat(jc))
          z_cosy(jc) = COS(z_lat(jc))
          z_fac1(jc) = 1.0_wp/6.3_wp-2.0_wp*(z_siny(jc)**6)*(z_cosy(jc)**2+1.0_wp/3.0_wp)
          z_fac2(jc) = (8.0_wp/5.0_wp*(z_cosy(jc)**3)*(z_siny(jc)**2+2.0_wp/3.0_wp)&
                       -0.25_wp*pi)*re*omega
          z_exp(jc)  = rd*gamma/grav
          zeta_old(jc) = 1.0e-7_wp
        ENDDO
        ! Newton iteration to determine zeta
!         CALL message(TRIM(routine),'do 1 ends')
        DO jn = 1, 100
          DO jc = 1, nlen
            zeta_v(jc,jk,jb) = (zeta_old(jc) - eta0)*pi_2
            zcoszetav(jc)= COS(zeta_v(jc,jk,jb))
            zsinzetav(jc)= SIN(zeta_v(jc,jk,jb))
            z_tavg(jc)   = temp0*(zeta_old(jc)**z_exp(jc))
            z_favg(jc)   = temp0*grav/gamma*(1.0_wp-zeta_old(jc)**z_exp(jc))
            IF (zeta_old(jc) < etat ) THEN
               z_tavg(jc) = z_tavg(jc)+dtemp*((etat-zeta_old(jc))**5)
               z_favg(jc) = z_favg(jc)-rd*dtemp*(                           &
                        (log(zeta_old(jc)/etat)+137.0_wp/60.0_wp)*(etat**5) &
                        -5.0_wp*(etat**4)*zeta_old(jc)                      &
                        +5.0_wp*(etat**3)*(zeta_old(jc)**2)                 &
                        -10.0_wp/3.0_wp*(etat**2)*(zeta_old(jc)**3)         &
                        +1.25_wp*etat*(zeta_old(jc)**4)-0.2_wp*(zeta_old(jc)**5))
            ENDIF
            z_geopot(jc) = z_favg(jc)+u0*(zcoszetav(jc)**1.5_wp)*&
                          (z_fac1(jc)*u0*(zcoszetav(jc)**1.5_wp)+z_fac2(jc)) 
            z_temp(jc)   = z_tavg(jc)+0.75_wp*zeta_old(jc)*pi*u0/rd*zsinzetav(jc)*&
                       SQRT(zcoszetav(jc))*(2.0_wp*u0*z_fac1(jc)*(zcoszetav(jc)**1.5_wp) &
                       + z_fac2(jc))
            z_fun(jc)    = z_geopot (jc)- p_nhdom%metrics%geopot(jc,jk,jb)
            z_fund(jc)   = -rd/zeta_old(jc)*z_temp(jc)
            zeta(jc) = zeta_old(jc) - z_fun(jc)/z_fund(jc)
            zeta_old(jc) = zeta(jc)
          ENDDO ! jc
        ENDDO !jn
!         CALL message(TRIM(routine),'do 2 ends')
        ! Final update for zeta_v
        DO jc = 1, nlen
          zeta_v(jc,jk,jb) = (zeta_old(jc) - eta0)*pi_2
        ENDDO
!         CALL message(TRIM(routine),'do 3 ends')
        ! Save results for all time levels
        ! Use analytic expressions at all model level
        DO jt = 1, ntl
          DO jc = 1, nlen
            p_nhdom%prog(jt)%exner(jc,jk,jb) = zeta_old(jc)**(rd/cpd)
            p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) = &
            &        p_nhdom%prog(jt)%exner(jc,jk,jb)**cvd_o_rd*p0ref/rd
            p_nhdom%prog(jt)%theta_v(jc,jk,jb) = z_temp(jc) & 
            &        /p_nhdom%prog(jt)%exner(jc,jk,jb)
            p_nhdom%prog(jt)%rho(jc,jk,jb) = &
            &        p_nhdom%prog(jt)%rhotheta_v(jc,jk,jb) &
            &        /p_nhdom%prog(jt)%theta_v(jc,jk,jb)
            IF (jt == 1 ) THEN
              p_nhdom%diag%pres(jc,jk,jb) = p0ref*p_nhdom%prog(jt)%exner(jc,jk,jb)**(cpd/rd) 
              p_nhdom%diag%temp(jc,jk,jb) = z_temp(jc)  
            ENDIF
          ENDDO !jc
        ENDDO !jt
!         CALL message(TRIM(routine),'do 4 ends')
      ENDDO !jk
!       CALL message(TRIM(routine),'jk do ends')
    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

!     IF (p_test_run) zeta_v_e = 0._wp
    CALL cells2edges_scalar(zeta_v,p_patch(jg),p_int(jg)%c_lin_e,zeta_v_e)
!     CALL sync_patch_array(SYNC_E,p_patch(jg),zeta_v_e)

    i_startblk = p_patch(jg)%edges%start_blk(2,1)
    
!     CALL message(TRIM(routine),'start omp 2')
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,jt,z_lat,z_lon,zu,z_fac1,z_fac2,zv)
    DO jb = i_startblk, nblks_e


      CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, i_startidx, i_endidx, 2)

!       write(0,*) "jb=",jb,"i_startidx, i_endidx=",i_startidx, i_endidx
      DO jt = 1, ntl
!         write(0,*) "jb=",jb, "jt=", jt
        DO jk = 1, nlev
!           write(0,*) "jb=",jb, "jt=", jt, "jk=", jk
!           write(0,*) "i_startidx, i_endidx=",i_startidx, i_endidx
          DO je = i_startidx, i_endidx
!             write(0,*) "jb=",jb, "jt=", jt, "jk=", jk, "je=",je
            z_lat(je) = p_patch(jg)%edges%center(je,jb)%lat
            z_lon(je) = p_patch(jg)%edges%center(je,jb)%lon
!             write(0,*) "  z_lat(je), z_lon(je)=", z_lat(je), z_lon(je)
            zu(je)    = u0*(COS(zeta_v_e(je,jk,jb))**1.5_wp)*(SIN(2.0_wp*z_lat(je))**2)
!             write(0,*) "  zu(je)=", zu(je)
             z_fac1(je)= SIN(latC)*SIN(z_lat(je))+COS(latC)*COS(z_lat(je))*COS(z_lon(je)-lonC)
             z_fac2(je)  = 10._wp*ACOS(z_fac1(je))
!              write(0,*) "  z_fac1(je),z_fac2(je)=", z_fac1(je),z_fac2(je)
             zu(je) = zu(je) + jw_up* EXP(-(z_fac2(je)*z_fac2(je)))
!             write(0,*) "  zu(je)=", zu(je)
            zv(je) = 0._wp
            p_nhdom%prog(jt)%vn(je,jk,jb) = &
                       zu(je) * p_patch(jg)%edges%primal_normal(je,jb)%v1   &
                   & + zv(je) * p_patch(jg)%edges%primal_normal(je,jb)%v2
!             write(0,*) "  p_nhdom%prog(jt)%vn(je,jk,jb)=", p_nhdom%prog(jt)%vn(je,jk,jb)
          ENDDO
!           write(0,*) "end do je"
        ENDDO
!         write(0,*) "end do jk"
      ENDDO !jt
!       write(0,*) "end do jt"
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
!     CALL message(TRIM(routine),'end omp 2')


    DEALLOCATE (zeta_v,zeta_v_e)

  ENDDO !jg

  IF (msg_level > 10) &
    CALL message(TRIM(routine),'End setup Jablonowski test')
  

 END SUBROUTINE init_nh_testcase

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE init_nh_step (p_nh_state)
    TYPE(t_nh_state), TARGET, INTENT(INOUT):: p_nh_state

    INTEGER :: jt
    
    DO jt = 1, ntl_module
      p_nh_state%prog(jt)%vn(:,:,:)         = 0.0_wp
      p_nh_state%prog(jt)%w(:,:,:)          = 0.0_wp
      p_nh_state%prog(jt)%rho(:,:,:)        = 0.0_wp
      p_nh_state%prog(jt)%exner(:,:,:)      = 0.0_wp
      p_nh_state%prog(jt)%rhotheta_v(:,:,:) = 0.0_wp
      p_nh_state%prog(jt)%theta_v(:,:,:)    = 0.0_wp
    ENDDO
     
    p_nh_state%diag%u(:,:,:)             = 0.0_wp
    p_nh_state%diag%v(:,:,:)             = 0.0_wp

    p_nh_state%diag%vt(:,:,:)            = 0.0_wp
    p_nh_state%diag%e_kin(:,:,:)         = 0.0_wp
    p_nh_state%diag%omega_z(:,:,:)       = 0.0_wp
    p_nh_state%diag%omega_z_c(:,:,:)     = 0.0_wp
    p_nh_state%diag%ddt_vn(:,:,:)        = 0.0_wp
    p_nh_state%diag%ddt_w(:,:,:)         = 0.0_wp
    p_nh_state%diag%ddt_exner(:,:,:)     = 0.0_wp
    p_nh_state%diag%exner_old(:,:,:)     = 0.0_wp
    p_nh_state%diag%w_con(:,:,:)         = 0.0_wp
    p_nh_state%diag%temp(:,:,:)          = 0.0_wp
    p_nh_state%diag%temp_ifc(:,:,:)      = 0.0_wp
    p_nh_state%diag%pres(:,:,:)          = 0.0_wp
    p_nh_state%diag%pres_ifc(:,:,:)      = 0.0_wp
    p_nh_state%diag%dpres_mc(:,:,:)      = 0.0_wp
    p_nh_state%diag%pres_sfc(:,:)        = 0.0_wp
    p_nh_state%diag%div(:,:,:)           = 0.0_wp
    p_nh_state%diag%mass_fl_e(:,:,:)     = 0.0_wp
    p_nh_state%diag%rho_ic(:,:,:)        = 0.0_wp
    p_nh_state%diag%w_concorr_c(:,:,:)   = 0.0_wp
    p_nh_state%diag%e_kinh(:,:,:)        = 0.0_wp
    p_nh_state%diag%vn_half(:,:,:)       = 0.0_wp
    p_nh_state%diag%theta_v_h(:,:,:)     = 0.0_wp
    p_nh_state%diag%ddt_vn_adv(:,:,:,:)  = 0.0_wp
    p_nh_state%diag%ddt_w_adv(:,:,:,:)   = 0.0_wp
    p_nh_state%diag%grf_tend_vn(:,:,:)   = 0.0_wp
    p_nh_state%diag%grf_tend_w(:,:,:)    = 0.0_wp
    p_nh_state%diag%grf_tend_rho(:,:,:)  = 0.0_wp
    p_nh_state%diag%grf_tend_thv(:,:,:)  = 0.0_wp
    p_nh_state%diag%ddt_vn_phy    (:,:,:)  = 0.0_wp
    p_nh_state%diag%ddt_exner_phy (:,:,:)  = 0.0_wp

    
  END SUBROUTINE init_nh_step 

 !-------------------------------------------------------------------------
 !>
 !!               Initialization of variables that contain general information
 !!               about the nh run. The configuration is read from
 !!               namelist 'nh_ctl'.
 !!
  !! Defines nonhydrostatic artificial initial conditions.
  !! 
  !! Reads namelist
  !! 
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-03-19)
  !! Modification by Daniel Reinert, DWD (2010-07-15)
  !! - moved initialization of topography into new subroutine 
  !!   init_nh_testtopo, which is called after the domain-decomposition.
  !!   (because of possible conflicts with the external-data type)
  !! 
 SUBROUTINE setup_nh_parameters


  CHARACTER(*), PARAMETER :: method_name='setup_nh_parameters'

  CHARACTER(len=max_char_length), PARAMETER :: &
            routine = 'setup_nh_parameters'

! !local variable
  INTEGER :: i_status

  NAMELIST /nh_ctl/ nh_dtime, nh_iterations
  
  NAMELIST/nh_testcase_ctl/  mount_height, torus_domain_length, &
                            nh_brunt_vais, nh_u0, nh_t0, layer_thickness,    &
                            jw_up, u0_mrw, mount_height_mrw,   &
                            mount_half_width, mount_lonctr_mrw_deg,          &
                            mount_latctr_mrw_deg, p_int_mwbr_const,          &
                            temp_i_mwbr_const, bruntvais_u_mwbr_const,       &
                            u0_mwbr_const,  hs_nh_vn_ptb_scale,               &
                            lhs_nh_vn_ptb,  n_flat_level,              &
                            rh_at_1000hpa, qv_max
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  ! set up the default values

  ! reduced calling frequency for transport
!   iadv_rcf = 1  ! no reduced calling frequency

  ! Type of vertical coordinate (1: Gal-Chen, 2: SLEVE)
  ivctype        = 1
  damp_height    = 10000.0_wp

  ! Settings for icell_type=3
  rayleigh_coeff = 0.2_wp
  vwind_offctr   = 0.03_wp
!   iadv_rhotheta  = 2
  igradp_method  = 1
!   l_open_ubc     = .FALSE.
!   l_nest_rcf     = .TRUE.
!   l_masscorr_nest= .FALSE.
  exner_expol    = 0.5_wp

  ! truly horizontal temperature diffusion (experimental!)
  l_zdiffu_t     = .FALSE. ! not used by default
  thslp_zdiffu   = 0.025_wp ! slope threshold 0.025
  thhgtd_zdiffu  = 200._wp ! threshold for height difference between adjacent grid points 200 m

  ! Settings for icell_type=6
!   ltheta_up_hori =.FALSE.
!   gmres_rtol_nh  = 1.0e-6_wp
  upstr_beta     = 1.0_wp
!   l_impl_vert_adv=.TRUE.
!   k2_updamp_coeff= 2.0e6_wp
   nh_iterations = 1
  !-----------------------------------------------------------------------
  ! read namelist
   CALL open_nml('NAMELIST_ICON_TESTBED')

   ! read nh_ctl namelist
   CALL position_nml ('nh_ctl', status=i_status)
   SELECT CASE (i_status)
   CASE (positioned)
     READ (nnml, nh_ctl)
   CASE default
     CALL finish( method_name,'run_ctl not positioned')
   END SELECT


  !-----------------------------------------------------------------------
    ! default topo values
    mount_height           = 100.0_wp
    layer_thickness        = -999.0_wp
    n_flat_level           = 2
    nh_u0                  = 0.0_wp
    nh_brunt_vais          = 0.01_wp
    nh_t0                  = 300.0_wp
    jw_up                  = 1.0_wp
    u0_mrw                 = 20.0_wp
    mount_height_mrw       = 2000.0_wp
    mount_half_width       = 1500000._wp
    mount_lonctr_mrw_deg   = 90.0_wp
    mount_latctr_mrw_deg   = 30.0_wp
    u0_mwbr_const          = 20.0_wp
    p_int_mwbr_const       = 70000._wp
    temp_i_mwbr_const      = 288._wp
    bruntvais_u_mwbr_const = 0.025_wp
    rotate_axis_deg        = 0.0_wp
    torus_domain_length    = 100000.0_wp
    lhs_nh_vn_ptb          = .TRUE.
    lhs_fric_heat          = .FALSE.
    hs_nh_vn_ptb_scale     = 1._wp  ! magnitude of the random noise
    rh_at_1000hpa          = 0.7_wp
    qv_max                 = 20.e-3_wp ! 20 g/kg


    CALL position_nml ('nh_testcase_ctl', status=i_status)
    SELECT CASE (i_status)
    CASE (POSITIONED)
      READ (nnml, nh_testcase_ctl)
    END SELECT

    n_flat_level=MAX(2,n_flat_level)
  !-----------------------------------------------------------------------
! Close the namelist files.
  CALL close_nml

  !-----------------------------------------------------------------------
  ! check the validity of the settings
  IF (upstr_beta > 1.0_wp .OR. upstr_beta < 0.0_wp) THEN
    CALL finish(TRIM(routine), 'upstr_beta out of range 0..1')
  ENDIF

  IF (msg_level > 10) THEN
    WRITE(message_text,'(a,i6, a)') " nh_solver will iterate ",  nh_iterations, " times."
    CALL message('non-hydro_atmos',message_text)
  ENDIF

 END SUBROUTINE setup_nh_parameters
 !-------------------------------------------------------------------------
  
  

END MODULE mo_nh_driver


