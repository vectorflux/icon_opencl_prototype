!>
!! This is the main driver for the ICON testbed.
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
PROGRAM icon_testbed

  USE mo_kind,                ONLY: wp

  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_namelist,            ONLY: open_nml, close_nml, &
    & open_nml_output, close_nml_output, &
    & position_nml, positioned
  USE mo_io_units,           ONLY: nnml, filename_max

  USE mo_mpi,                 ONLY: p_start, p_stop, p_pe, p_io, p_nprocs
  USE mo_timer,               ONLY: init_timer, cleanup_timer, print_timer, &
    &  timer_start, timer_stop, timer_total
  USE mo_io_units,            ONLY: filename_max
  USE mo_datetime,            ONLY: t_datetime
  USE mo_setup_parameters,    ONLY: setup_run_parameters

  USE mo_parallel_ctl,        ONLY: setup_parallel_run,   & ! process parallel run ctl. params.
     &                              p_test_pe,            & !    internal parameter
     &                              p_comm_work,          &
     &                              p_test_run,           &
     &                              p_io_pe0                ! Number of first I/O PE

  USE mo_global_variables,    ONLY: ltimer, ini_datetime, dtime, &
    & nproma, msg_level, modelname

  USE mo_model_domain_import, ONLY: import_domain, destruct_domain        !
  USE mo_profile

#ifdef __ICON_TESTBED_BASE__
  USE mo_nh_driver,        ONLY: nh_solver_driver, construct_nh_solver, destruct_nh_solver
#endif
 
#ifdef __ICON_TESTBED_RADIATION__
  USE mo_radiation_driver, ONLY: radiation_driver, construct_radiaton, destruct_radiaton
#endif
  IMPLICIT NONE

  TYPE(t_datetime)  :: datetime
  CHARACTER(*), PARAMETER :: method_name='icon_testbed_driver'

  ! the testbed parameters
  INTEGER ::  main_iterations=1 ! how many times we will run the tests from the main driver
  LOGICAL ::  run_radiation=.false.
  LOGICAL ::  run_nh_solver=.true.
  LOGICAL ::  run_math_operators=.false.
  INTEGER  :: nproma_set(60) ! set of nproma values
  
  NAMELIST /testbed_ctl/ nproma_set, main_iterations, &
    & run_radiation, run_nh_solver, run_math_operators
  
  INTEGER :: iter, i_status, nproma_i

!-------------------------------------------------
! declaration of OpenMP Runtime Library Routines:
!$  INTEGER omp_get_max_threads
!$  INTEGER omp_get_num_threads
!$  INTEGER omp_get_num_procs
!$  INTEGER omp_get_thread_num
!$  LOGICAL omp_get_dynamic

!$  INTEGER :: max_threads_omp, num_procs_omp
!$  LOGICAL :: ldynamic_omp

  !print out some information about OpenMP parallelization
!$  max_threads_omp  = omp_get_max_threads()
!$  num_procs_omp    = omp_get_num_procs()
!$  ldynamic_omp     = omp_get_dynamic()
!$  WRITE(message_text,'(A,I3,A,I3)')                &
!$    & "OpenMP:  MAX_THREADS = ", max_threads_omp,  &
!$    & ",  NUM_PROCS = ", num_procs_omp
!$  CALL message(TRIM(method_name),message_text)
!$  WRITE(message_text,'(A,L3)')  &
!$    & "OpenMP:  DYNAMIC = ", ldynamic_omp
!$  CALL message(TRIM(method_name),message_text)
!-------------------------------------------------

  !-------------------------------------------------------------------
  ! Initialize MPI, this should always be the first call
  CALL p_start(modelname)

  !-------------------------------------------------------------------
  ! inner loop length to optimize vectorization
  nproma_set(:) = 0
  nproma_set(1) = 1
  !-------------------------------------------------------------------
  ! Open and read the namelist file.
  CALL open_nml('NAMELIST_ICON_TESTBED')
  ! read testbed_ctl namelist
  CALL position_nml ('testbed_ctl', status=i_status)
  SELECT CASE (i_status)
  CASE (positioned)
     READ (nnml, testbed_ctl)
  CASE default
     CALL finish( method_name,'testbed_ctl not positioned')
  END SELECT
  ! Close the namelist files.
  CALL close_nml
  !-------------------------------------------------------------------
  
  !-------------------------------------------------------------------
  IF (msg_level > 0) CALL message(method_name, 'setup...')
  CALL setup_run_parameters()
  CALL init_timer()

  nproma_i=1
  DO WHILE (nproma_set(nproma_i) > 0) 
    nproma=nproma_set(nproma_i)
    nproma_i = nproma_i + 1
    
    IF (msg_level > 5) THEN
      write(0,*) "=============================="
      write(0,*) "Setup icon_testbed with nproma=", nproma
    ENDIF
    CALL cleanup_timer()
  ! Import domain
    CALL import_domain()
    IF (msg_level > 10) CALL message(method_name, 'import_domain done.')

    !-------------------------------------------------------------------
    ! Initialize date and time
    datetime = ini_datetime

    !------------------------------------------------------------------
    ! Setup tests
#ifdef __ICON_TESTBED_BASE__
    IF (run_nh_solver) CALL construct_nh_solver()
#endif
#ifdef __ICON_TESTBED_RADIATION__
    IF (run_radiation) CALL construct_radiaton()
#endif
    !------------------------------------------------------------------
    ! run tests
    IF (msg_level > 0) THEN
      write(0,*) "=============================="
      write(0,*) "Running icon_testbed with nproma=", nproma
    ENDIF

    WRITE(message_text,'(a,i4.4)') 'nproma.', nproma
    CALL trace_init(TRIM(message_text),nproma)

    CALL trace_start("total", 1)

	IF (ltimer) CALL timer_start(timer_total)
	CALL start_cpp_timer()

    DO iter=1,main_iterations    
#ifdef __ICON_TESTBED_BASE__
      IF (run_nh_solver) CALL nh_solver_driver(datetime)
#endif
#ifdef __ICON_TESTBED_RADIATION__
      IF (run_radiation) CALL radiation_driver(datetime)
#endif
    ENDDO

	CALL stop_cpp_timer()
    IF (ltimer) CALL timer_stop(timer_total)

    CALL trace_stop("total", 1)
    CALL trace_finalize(nproma)

    CALL print_timer()
	
	CALL print_profiling_info()
    !------------------------------------------------------------------
    !  cleaning up tests
    IF (msg_level > 5) CALL message(method_name,'start to clean up...')
#ifdef __ICON_TESTBED_BASE__
    IF (run_nh_solver) CALL destruct_nh_solver()
#endif
#ifdef __ICON_TESTBED_RADIATION__
    IF (run_radiation) CALL destruct_radiaton()
#endif
    CALL destruct_domain()
    IF (msg_level > 5) CALL message(method_name,'clean-up finished')
    !------------------------------------------------------------------
  
  ENDDO ! nproma=nproma_range(1),nproma_range(2),nproma_range(3)

  !------------------------------------------------------------------
  ! Shut down MPI
  CALL p_stop

END PROGRAM icon_testbed

