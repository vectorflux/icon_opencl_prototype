!>
!! mo_solve_nonhydro
!!
!! This module contains the nonhydrostatic dynamical core for the triangular version
!! Its routines were previously contained in mo_divergent_modes and mo_vector_operations
!! but have been extracted for better memory efficiency
!!
!! @author Guenther Zaengl, DWD
!!
!! @par Revision History
!! Initial release by Guenther Zaengl (2010-10-13) based on earlier work
!! by Almut Gassmann, MPI-M
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
MODULE mo_solve_nonhydro

  USE mo_kind,              ONLY: wp
  USE mo_nh_global_variables,ONLY: igradp_method, nproma, ltimer, &
     itime_scheme, idiv_method, nlev, nlevp1
  USE mo_model_domain,      ONLY: t_patch
!   USE mo_model_domain_import,ONLY: l_limited_area
  USE mo_interpolation,     ONLY: t_int_state, cells2edges_scalar, edges2cells_scalar,  &
                                  rbf_vec_interpol_edge, cell_avg
  USE mo_nonhydro_state,    ONLY: t_nh_state, t_nh_metrics, t_nh_diag, t_nh_prog
!                                   t_buffer_memory
  USE mo_physical_constants,ONLY: cpd, rd, cvd, cvd_o_rd, grav, rd_o_cpd, p0ref
  USE mo_math_operators,    ONLY: div, rot_vertex, div_avg
  USE mo_vertical_grid,     ONLY: nflat, nrdmax, nflat_gradp
  USE mo_loopindices,       ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_impl_constants,    ONLY: min_rlcell_int, min_rledge_int, min_rlvert_int, min_rlcell
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
!   USE mo_advection_hflux,   ONLY: upwind_hflux_miura, upwind_hflux_miura3
!   USE mo_sync,              ONLY: SYNC_E, SYNC_C, sync_patch_array, sync_patch_array_mult, &
!                                   sync_patch_array_gm
!   USE mo_parallel_ctl,      ONLY: p_test_run, itype_comm
  USE mo_mpi,               ONLY: p_nprocs
!   USE mo_timer,             ONLY: timer_solve_nh, timer_start, timer_stop

  USE mo_intp_data_strc,    ONLY: rbf_vec_dim_e
  USE mo_global_variables,  ONLY: i_cell_type
  USE mo_impl_constants,      ONLY: min_rledge

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  REAL(wp), PARAMETER :: rd_o_cvd = 1._wp / cvd_o_rd
  REAL(wp), PARAMETER :: cpd_o_rd = 1._wp / rd_o_cpd
  REAL(wp), PARAMETER :: rd_o_p0ref = rd / p0ref
  REAL(wp), PARAMETER :: grav_o_cpd = grav / cpd

  PUBLIC :: solve_nh

  CONTAINS


  !----------------------------------------------------------------------------
  !>
  !! velocity_tendencies
  !!
  !! Discretization of nonhydrostatic momentum equation similar to hydrostatic core
  !! In particular, the Lamb transformation is applied only to the horizontal
  !! equation of motion, whereas the vertical wind equation is discretized
  !! in advective form
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl (2010-02-03)
  !!
  SUBROUTINE velocity_tendencies (p_prog,p_patch,p_int,p_metrics,p_diag,ntnd,istep,l_init)

    ! Passed variables
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_int_state), TARGET, INTENT(IN):: p_int
    TYPE(t_nh_prog), INTENT(IN)          :: p_prog
    TYPE(t_nh_metrics), INTENT(IN)       :: p_metrics
    TYPE(t_nh_diag), INTENT(inout)       :: p_diag

    INTEGER, INTENT(IN)             :: ntnd ! time level of ddt_adv fields used to store tendencies
    INTEGER, INTENT(IN)             :: istep ! 1: predictor step, 2: corrector step
    LOGICAL,       INTENT(IN)       :: l_init

    ! Local variables
    INTEGER :: jb, jk, jc, je
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom, nameIndex
    INTEGER :: rl_start, rl_end
    REAL(wp):: z_concorr_e(nproma,nlevp1,p_patch%nblks_e)
    REAL(wp):: z_w_con_c_full(nproma,nlev,p_patch%nblks_c)
    REAL(wp):: z_vt_ie(nproma,nlevp1)
    REAL(wp):: z_kin_hor_e(nproma,nlev,p_patch%nblks_e)
    REAL(wp):: z_ddxn_ekin_e(nproma,nlev,p_patch%nblks_e)
    REAL(wp):: z_vnw(nproma,nlevp1,p_patch%nblks_e)
    REAL(wp):: z_hadv_w(nproma,nlevp1,p_patch%nblks_c)

    INTEGER,  DIMENSION(:,:,:),   POINTER :: icidx, icblk, ieidx, ieblk, ividx, ivblk
!     INTEGER  :: nlev, nlevp1          !< number of full and half levels
    ! Preparation for vertical nesting; switch will become a patch attribute later on
    LOGICAL :: l_vert_nested, l_onGPU
	
	! additional resources needed for C++/OpenCL code
	INTEGER :: i_rcstartlev, nblks_e, ist, i_rlstart, i_nnow
    INTEGER,  DIMENSION(:), ALLOCATABLE :: blockLocalStart, blockLocalEnd

!===================================================================	
#ifdef TEST
	REAL(wp):: output_rbf(nproma,nlev  ,p_patch%nblks_e)
	REAL(wp):: output_rot(nproma,nlev  ,p_patch%nblks_v)
	REAL(wp):: output_e2c(nproma,nlev  ,p_patch%nblks_c)
	REAL(wp):: output_e2c_p1(nproma,nlevp1,p_patch%nblks_c)
	REAL(wp):: output_avg(nproma,nlevp1,p_patch%nblks_c,3)
	REAL(wp):: output_p1a(nproma,nlevp1,p_patch%nblks_e)
	REAL(wp):: output_p1b(nproma,nlev  ,p_patch%nblks_e)
	REAL(wp):: output_p1c(nproma,nlevp1)
	REAL(wp):: output_p1d(nproma,nlevp1,p_patch%nblks_e)
	REAL(wp):: output_p2a(nproma,nlevp1,p_patch%nblks_e)
	REAL(wp):: output_p2b(nproma,nlev  ,p_patch%nblks_e)
	REAL(wp):: output_p3a(nproma,nlevp1,p_patch%nblks_c)
	REAL(wp):: output_p3b(nproma,nlevp1,p_patch%nblks_c)
	REAL(wp):: output_p3c(nproma,nlev  ,p_patch%nblks_c)
	REAL(wp):: output_p4( nproma,nlev  ,p_patch%nblks_e,3)
	REAL(wp):: output_p5( nproma,nlevp1,p_patch%nblks_c,3)
#endif
!===================================================================

    !--------------------------------------------------------------------------

!     IF ((lvert_nest) .AND. (p_patch%nshift > 0)) THEN  
!       l_vert_nested = .TRUE.
!     ELSE
      l_vert_nested = .FALSE.
!     ENDIF

    ! number of vertical levels
!     nlev   = nlev
!     nlevp1 = nlevp1

    ! Set pointers to neighbor cells/edges/vertices
    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    i_nchdom   = MAX(1,p_patch%n_childdom)
	
	IF (istep==1) THEN
	  i_nnow=1
	ELSE
	  i_nnow=0
	ENDIF
	

    ! Tangential wind component using RBF reconstruction
    ! Note: vt is also diagnosed in divergent_modes. Thus, computation is needed in predictor step only
    IF (istep == 1) THEN
      i_rcstartlev = 2
	  i_startblk = p_patch%edges%start_blk(i_rcstartlev,1)
	  nblks_e = p_patch%nblks_int_e

#ifndef TEST	  
	  IF (l_init) THEN
#endif
        ist=0
        ALLOCATE(blockLocalStart(nblks_e-i_startblk+1), stat=ist); blockLocalStart=0
        ALLOCATE(blockLocalEnd(nblks_e-i_startblk+1), stat=ist); blockLocalEnd=0
	
	    DO jb = i_startblk, nblks_e
          CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
                             blockLocalStart(jb), blockLocalEnd(jb), i_rcstartlev)
          blockLocalStart(jb) = blockLocalStart(jb) - 1
          blockLocalEnd(jb)   = blockLocalEnd(jb)
        ENDDO
#ifndef TEST
	  ENDIF
#endif	  
!===================================================================
#ifdef TEST
      output_rbf = p_diag%vt
#endif
!===================================================================	  

	  CALL ocl_rbf_vec_interpol_edge(p_prog%vn, 18+i_nnow, .FALSE., &
									 istep, &
	                                 p_int%rbf_vec_coeff_e, &
									 p_int%rbf_vec_idx_e, &
									 p_int%rbf_vec_blk_e, &
									 i_startblk, &
									 nproma, &
									 nlev, &
									 rbf_vec_dim_e, &
									 nblks_e, &
									 p_diag%vt, 1, .FALSE., &
									 blockLocalStart, &
									 blockLocalEnd, &
									 0, nlev)
									 
#ifndef TEST
	  IF (l_init) THEN
#endif
	    DEALLOCATE(blockLocalStart,stat=ist)
	    DEALLOCATE(blockLocalEnd,  stat=ist)
#ifndef TEST
	  ENDIF
#endif
!===================================================================
#ifdef TEST
      CALL rbf_vec_interpol_edge(p_prog%vn, p_patch, p_int, output_rbf)
	  
      WRITE(*,*) 'checking rbf...'
      CALL checkresults(output_rbf, p_diag%vt, &
                        0,SIZE(p_diag%vt,1),SIZE(p_diag%vt,1), &
                        0,SIZE(p_diag%vt,2),SIZE(p_diag%vt,2), &
                        0,SIZE(p_diag%vt,3),SIZE(p_diag%vt,3),0)
#endif
!===================================================================
    ENDIF

    ! Compute vertical vorticity component at vertices
    IF (i_cell_type==3) THEN
	  i_rlstart = 2
	  i_startblk = p_patch%verts%start_blk(i_rlstart,1)
	  i_endblk   = p_patch%verts%end_blk(min_rlvert_int-1,MAX(1,p_patch%n_childdom))
	  
#ifndef TEST      
	  IF (l_init) THEN
#endif
	    ist=0
        ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
        ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0

	    DO jb = i_startblk, i_endblk
          CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, &
                             blockLocalStart(jb), blockLocalEnd(jb), i_rlstart)
          blockLocalStart(jb) = blockLocalStart(jb) - 1
          blockLocalEnd(jb)   = blockLocalEnd(jb)
        ENDDO
#ifndef TEST      
	  ENDIF
#endif

!===================================================================
#ifdef TEST
	  output_rot = p_diag%omega_z
#endif
!===================================================================
	  IF (istep==1) THEN
	    CALL ocl_rot_vertex_atmos_tri(p_prog%vn, 19, .TRUE., &
		                              istep, &
									  p_int%geofac_rot, &
									  p_patch%verts%edge_idx, &
									  p_patch%verts%edge_blk, &
									  p_patch%verts%start_blk(i_rlstart,1), &
									  nproma, &
									  0, nlev, &
									  nlev, &
									  p_patch%nblks_int_v, &
									  p_patch%nblks_int_e, &
									  p_diag%omega_z, 2, .FALSE., &
									  blockLocalStart, &
									  blockLocalEnd)
	  ELSE
	    CALL ocl_rot_vertex_atmos_tri(p_prog%vn, 18, .FALSE., &
		                              istep, &
									  p_int%geofac_rot, &
									  p_patch%verts%edge_idx, &
									  p_patch%verts%edge_blk, &
									  p_patch%verts%start_blk(i_rlstart,1), &
									  nproma, &
									  0, nlev, &
									  nlev, &
									  p_patch%nblks_int_v, &
									  p_patch%nblks_int_e, &
									  p_diag%omega_z, 2, .FALSE., &
									  blockLocalStart, &
									  blockLocalEnd)
	  ENDIF

#ifndef TEST
	  IF (l_init) THEN
#endif
        DEALLOCATE(blockLocalStart,stat=ist)
	    DEALLOCATE(blockLocalEnd,  stat=ist)
#ifndef TEST
	  ENDIF
#endif
!===================================================================
#ifdef TEST
    CALL rot_vertex (p_prog%vn, p_patch, p_int, output_rot, opt_rlend=min_rlvert_int-1)
	  
      WRITE(*,*) 'checking rot...'
      CALL checkresults(output_rot, p_diag%omega_z, &
                        0,SIZE(p_diag%omega_z,1),SIZE(p_diag%omega_z,1), &
                        0,SIZE(p_diag%omega_z,2),SIZE(p_diag%omega_z,2), &
                        0,SIZE(p_diag%omega_z,3),SIZE(p_diag%omega_z,3),1)
#endif
!===================================================================
    ELSE
      WRITE(*,*) 'this should not happen here'
      CALL rot_vertex (p_prog%vn, p_patch, p_int, p_diag%omega_z, opt_rlend=min_rlvert_int-1)
    ENDIF

    rl_start = 3
    rl_end = min_rledge_int - 2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

#ifndef TEST	  
	IF (l_init) THEN
#endif
    ist=0
    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
    ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0

	DO jb = i_startblk, i_endblk
	  CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
      blockLocalStart(jb) = blockLocalStart(jb) - 1
      blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif

	
!===================================================================
#ifdef TEST
	output_p1a = p_diag%vn_half
	output_p1b = z_kin_hor_e
	output_p1d = z_concorr_e
#endif
!===================================================================
	  
    CALL ocl_part1(i_startblk, i_endblk, &
	               p_patch%nblks_int_e, p_patch%nblks_int_c, &
	               nlev, nlevp1, &
				   blockLocalStart, blockLocalEnd, &
				   nproma, &
				   istep, &
				   nflat, &
				   l_vert_nested, &
				   p_metrics%wgtfac_e, &
				   p_metrics%wgtfacq1_e, &
				   p_metrics%wgtfacq_e, &
				   p_prog%vn, 18+i_nnow, &
				   p_diag%vt, &
				   p_metrics%ddxn_z_half, &
				   p_metrics%ddxt_z_half, &
				   p_diag%vn_half, &
				   z_kin_hor_e, &
				   z_concorr_e)

#ifndef TEST	  
	  IF (l_init) THEN
#endif	
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,  stat=ist)
#ifndef TEST      
	  ENDIF
#endif

!===================================================================
#ifdef TEST
	DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (istep == 1) THEN

        ! Interpolate vn to interface levels and compute horizontal part of kinetic energy on edges
        DO jk = 2, nlev
          DO je = i_startidx, i_endidx
            output_p1a(je,jk,jb) =                                  &
              p_metrics%wgtfac_e(je,jk,jb)*p_prog%vn(je,jk,jb) +        &
             (1._wp - p_metrics%wgtfac_e(je,jk,jb))*p_prog%vn(je,jk-1,jb)
            output_p1b(je,jk,jb) = 0.5_wp*(p_prog%vn(je,jk,jb)*p_prog%vn(je,jk,jb) + &
              p_diag%vt(je,jk,jb)*p_diag%vt(je,jk,jb) )
          ENDDO
        ENDDO

        IF (.NOT. l_vert_nested) THEN
          ! Top and bottom levels
          DO je = i_startidx, i_endidx
            output_p1a(je,1,jb) =                                &
              p_metrics%wgtfacq1_e(je,1,jb)*p_prog%vn(je,1,jb) +   &
              p_metrics%wgtfacq1_e(je,2,jb)*p_prog%vn(je,2,jb) + &
              p_metrics%wgtfacq1_e(je,3,jb)*p_prog%vn(je,3,jb)
            output_p1b(je,1,jb) = 0.5_wp*(p_prog%vn(je,1,jb)*p_prog%vn(je,1,jb) + &
              p_diag%vt(je,1,jb)*p_diag%vt(je,1,jb) )
            output_p1a(je,nlevp1,jb) =                           &
              p_metrics%wgtfacq_e(je,1,jb)*p_prog%vn(je,nlev,jb) +   &
              p_metrics%wgtfacq_e(je,2,jb)*p_prog%vn(je,nlev-1,jb) + &
              p_metrics%wgtfacq_e(je,3,jb)*p_prog%vn(je,nlev-2,jb)
          ENDDO
        ELSE
          ! vn_half(jk=1) is interpolated from parent domain in this case
          DO je = i_startidx, i_endidx
            output_p1b(je,1,jb) = 0.5_wp*(p_prog%vn(je,1,jb)*p_prog%vn(je,1,jb) + &
              p_diag%vt(je,1,jb)*p_diag%vt(je,1,jb) )
            output_p1a(je,nlevp1,jb) =                           &
              p_metrics%wgtfacq_e(je,1,jb)*p_prog%vn(je,nlev,jb) +   &
              p_metrics%wgtfacq_e(je,2,jb)*p_prog%vn(je,nlev-1,jb) + &
              p_metrics%wgtfacq_e(je,3,jb)*p_prog%vn(je,nlev-2,jb)
          ENDDO

        ENDIF

      ELSE ! corrector step

        ! Compute only horizontal kinetic energy
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            output_p1b(je,jk,jb) = 0.5_wp*                &
             (p_prog%vn(je,jk,jb)*p_prog%vn(je,jk,jb) +    &
              p_diag%vt(je,jk,jb)*p_diag%vt(je,jk,jb) )
          ENDDO
        ENDDO

      ENDIF ! istep == 1

! Why is this step executed for all the blocks?
! a run with the last block only is sufficient as the data is overwritten nblocks times
      IF (istep == 1 .AND. l_init) THEN
        ! Interpolate vt to interface levels
        DO jk = nflat+1, nlev
          DO je = i_startidx, i_endidx
            output_p1c(je,jk) =                                     &
              p_metrics%wgtfac_e(je,jk,jb)*p_diag%vt(je,jk,jb) + &
             (1._wp - p_metrics%wgtfac_e(je,jk,jb))*p_diag%vt(je,jk-1,jb)
          ENDDO
        ENDDO

        ! Bottom level
        DO je = i_startidx, i_endidx
          output_p1c(je,nlevp1) =                                     &
            p_metrics%wgtfacq_e(je,1,jb)*p_diag%vt(je,nlev,jb) +   &
            p_metrics%wgtfacq_e(je,2,jb)*p_diag%vt(je,nlev-1,jb) + &
            p_metrics%wgtfacq_e(je,3,jb)*p_diag%vt(je,nlev-2,jb)
        ENDDO

        ! Compute contravariant correction for vertical velocity at interface levels
        ! (will be interpolated to cell centers below)
        DO jk = nflat+1, nlevp1
          DO je = i_startidx, i_endidx
            output_p1d(je,jk,jb) = &
              output_p1a(je,jk,jb)*p_metrics%ddxn_z_half(je,jk,jb) + &
              output_p1c(je,jk)   *p_metrics%ddxt_z_half(je,jk,jb)
          ENDDO
        ENDDO
      ENDIF

    ENDDO
	
      WRITE(*,*) 'checking p1...'
      CALL checkresults(output_p1a, p_diag%vn_half, &
                        0,SIZE(p_diag%vn_half,1),SIZE(p_diag%vn_half,1), &
                        0,SIZE(p_diag%vn_half,2),SIZE(p_diag%vn_half,2), &
                        0,SIZE(p_diag%vn_half,3),SIZE(p_diag%vn_half,3),2)
      CALL checkresults(output_p1b, z_kin_hor_e, &
                        0,SIZE(z_kin_hor_e,1),SIZE(z_kin_hor_e,1), &
                        0,SIZE(z_kin_hor_e,2),SIZE(z_kin_hor_e,2), &
                        0,SIZE(z_kin_hor_e,3),SIZE(z_kin_hor_e,3),3)
	  CALL checkresults(output_p1d, z_concorr_e, &
                        0,SIZE(z_concorr_e,1),SIZE(z_concorr_e,1), &
                        0,SIZE(z_concorr_e,2),SIZE(z_concorr_e,2), &
                        0,SIZE(z_concorr_e,3),SIZE(z_concorr_e,3),4)
#endif
!===================================================================
   

    IF (i_cell_type==3) THEN
	  i_startblk = p_patch%cells%start_blk(2,1)
	  i_endblk = p_patch%cells%end_blk(min_rlcell_int-1,MAX(1,p_patch%n_childdom))
  	  ist=0
      ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
      ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0

	  DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           blockLocalStart(jb), blockLocalEnd(jb), 2, min_rlcell_int-1)
        blockLocalStart(jb) = blockLocalStart(jb) - 1
        blockLocalEnd(jb)   = blockLocalEnd(jb)
      ENDDO

!===================================================================
#ifdef TEST
	  output_e2c = p_diag%e_kinh
#endif
!===================================================================

      CALL ocl_edges2cells_scalar_tri(z_kin_hor_e, 3, .TRUE., &
									  istep, &
	                                  p_int%e_bln_c_s, &
	                                  p_patch%cells%edge_idx, &
									  p_patch%cells%edge_blk, &
									  blockLocalStart, &
									  blockLocalEnd, &
									  nproma, &
									  1, nlev, &
									  nlev, &
									  i_startblk, &
									  i_endblk, &
									  p_patch%nblks_int_c, &
									  p_patch%nblks_int_e, &
									  p_diag%e_kinh, 4, .FALSE.)
	
!===================================================================
#ifdef TEST
	  CALL edges2cells_scalar(z_kin_hor_e,p_patch,p_int%e_bln_c_s,output_e2c,&
                              opt_rlstart=2, opt_rlend=min_rlcell_int-1)	  
      
	  WRITE(*,*) 'checking e2c...'
      CALL checkresults(output_e2c, p_diag%e_kinh, &
                        0,SIZE(p_diag%e_kinh,1),SIZE(p_diag%e_kinh,1), &
                        0,SIZE(p_diag%e_kinh,2),SIZE(p_diag%e_kinh,2), &
                        0,SIZE(p_diag%e_kinh,3),SIZE(p_diag%e_kinh,3),5)
#endif
!===================================================================

	  IF (istep == 1 .AND. l_init) THEN
!===================================================================
#ifdef TEST
  	    output_e2c_p1 = p_diag%w_concorr_c
#endif
!===================================================================
	    CALL ocl_edges2cells_scalar_tri(z_concorr_e, 5, .TRUE., &
									    istep, &
	                                    p_int%e_bln_c_s, &
	                                    p_patch%cells%edge_idx, &
						  			    p_patch%cells%edge_blk, &
						  			    blockLocalStart, &
									    blockLocalEnd, &
									    nproma, &
									    nflat+1, nlevp1, &
									    nlevp1, &
									    i_startblk, &
									    i_endblk, &
									    p_patch%nblks_int_c, &
									    p_patch%nblks_int_e, &
									    p_diag%w_concorr_c, 6, .FALSE.)
	
!===================================================================
#ifdef TEST
	  CALL edges2cells_scalar(z_concorr_e,p_patch,p_int%e_bln_c_s,output_e2c_p1,&
                                nflat+1, nlevp1, opt_rlstart=2, opt_rlend=min_rlcell_int-1)	  
      
	  WRITE(*,*) 'checking e2c...'
      CALL checkresults(output_e2c_p1, p_diag%w_concorr_c, &
                        0,SIZE(p_diag%w_concorr_c,1),SIZE(p_diag%w_concorr_c,1), &
                        0,SIZE(p_diag%w_concorr_c,2),SIZE(p_diag%w_concorr_c,2), &
                        0,SIZE(p_diag%w_concorr_c,3),SIZE(p_diag%w_concorr_c,3),5)
#endif
!===================================================================
      ENDIF

	  DEALLOCATE(blockLocalStart,stat=ist)
	  DEALLOCATE(blockLocalEnd,  stat=ist)

	ELSE
	  WRITE(*,*) 'If you see this, move it to GPU!'
      CALL edges2cells_scalar(z_kin_hor_e,p_patch,p_int%e_bln_c_s,p_diag%e_kinh,&
                              opt_rlstart=2, opt_rlend=min_rlcell_int-1)
	  IF (istep == 1 .AND. l_init) THEN
        CALL edges2cells_scalar(z_concorr_e,p_patch,p_int%e_bln_c_s,p_diag%w_concorr_c,&
                                nflat+1, nlevp1, opt_rlstart=2, opt_rlend=min_rlcell_int-1)
      ENDIF
    ENDIF

    rl_start = 4
    rl_end = min_rledge_int - 2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

#ifndef TEST	  
	  IF (l_init) THEN
#endif
  	ist=0
    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
	ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0
						
	DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
      blockLocalStart(jb) = blockLocalStart(jb) - 1
      blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif

!===================================================================
#ifdef TEST
	output_p2a = z_vnw
	output_p2b = z_ddxn_ekin_e
#endif
!===================================================================

	CALL ocl_part2(i_startblk, &
	               i_endblk, &
				   p_patch%nblks_int_e, &
				   p_patch%nblks_int_c, &
				   1, &
				   nlev, &
				   nlev, &
				   nlevp1, &
				   nproma, &
				   blockLocalStart, &
				   blockLocalEnd, &
				   istep, &
				   icidx, &
				   icblk, &
				   p_diag%vn_half, &
				   p_int%c_lin_e, &
				   p_prog%w, 20+i_nnow, &
				   p_patch%edges%inv_dual_edge_length, &
				   p_diag%e_kinh, &
				   z_vnw, &
				   z_ddxn_ekin_e)

#ifndef TEST	  
	  IF (l_init) THEN
#endif
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,  stat=ist)
#ifndef TEST      
	  ENDIF
#endif

!===================================================================
#ifdef TEST
!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
!CDIR UNROLL=3
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          ! Multiply vn_half with w interpolated to edges for divergence computation
          output_p2a(je,jk,jb) = p_diag%vn_half(je,jk,jb)*                            &
           ( p_int%c_lin_e(je,1,jb) * p_prog%w(icidx(je,jb,1),jk,icblk(je,jb,1)) &
           + p_int%c_lin_e(je,2,jb) * p_prog%w(icidx(je,jb,2),jk,icblk(je,jb,2)) )

          ! Compute horizontal gradient of horizontal kinetic energy
          output_p2b(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb) *  &
           (p_diag%e_kinh(icidx(je,jb,2),jk,icblk(je,jb,2)) -                    &
            p_diag%e_kinh(icidx(je,jb,1),jk,icblk(je,jb,1)) )
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
      
	  WRITE(*,*) 'checking p2...'
      CALL checkresults(output_p2a, z_vnw, &
                        0,SIZE(z_vnw,1),SIZE(z_vnw,1), &
                        0,SIZE(z_vnw,2),SIZE(z_vnw,2), &
                        0,SIZE(z_vnw,3),SIZE(z_vnw,3),6)
      CALL checkresults(output_p2b, z_ddxn_ekin_e, &
                        0,SIZE(z_ddxn_ekin_e,1),SIZE(z_ddxn_ekin_e,1), &
                        0,SIZE(z_ddxn_ekin_e,2),SIZE(z_ddxn_ekin_e,2), &
                        0,SIZE(z_ddxn_ekin_e,3),SIZE(z_ddxn_ekin_e,3),7)
#endif
!===================================================================

    rl_start = 3
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

#ifndef TEST	  
	  IF (l_init) THEN
#endif
  	ist=0
    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
	ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0
						
	DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
      blockLocalStart(jb) = blockLocalStart(jb) - 1
      blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif

!===================================================================
#ifdef TEST
	output_p3a = z_hadv_w
	output_p3b = p_diag%w_con
	output_p3c = z_w_con_c_full
#endif
!===================================================================
	
    CALL ocl_part3(i_startblk, &
	               i_endblk, &
				   p_patch%nblks_int_e, &
				   p_patch%nblks_int_c, &
				   nflat, &
				   nlev, &
				   nlevp1, &
				   nproma, &
				   blockLocalStart, &
				   blockLocalEnd, &
				   istep, &
				   ieidx, &
				   ieblk, &
				   p_diag%vn_half, &
				   p_int%geofac_div, &
				   p_prog%w, 20+i_nnow, &
				   z_vnw, &
				   z_hadv_w, &
				   p_diag%w_con, &
				   p_diag%w_concorr_c, &
				   z_w_con_c_full)
				   
#ifndef TEST	  
	  IF (l_init) THEN
#endif
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,  stat=ist)
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Compute horizontal advection of w: -(div(vn*w)-w*div(vn))
      ! (combined into one step for efficiency improvement)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
!CDIR UNROLL=5
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif
          output_p3a(jc,jk,jb) = p_prog%w(jc,jk,jb)* ( &
            p_diag%vn_half(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
            p_diag%vn_half(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
            p_diag%vn_half(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%geofac_div(jc,3,jb))- &
           (z_vnw(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))         *p_int%geofac_div(jc,1,jb) + &
            z_vnw(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))         *p_int%geofac_div(jc,2,jb) + &
            z_vnw(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))         *p_int%geofac_div(jc,3,jb))

          output_p3b(jc,jk,jb) = p_prog%w(jc,jk,jb)
        ENDDO
      ENDDO

!CDIR UNROLL=5
      ! Contravariant vertical velocity on w points and interpolation to full levels
      DO jk = nlev, nflat+1, -1
        DO jc = i_startidx, i_endidx
          output_p3b(jc,jk,jb) = p_diag%w_con(jc,jk,jb) - p_diag%w_concorr_c(jc,jk,jb)
          output_p3c(jc,jk,jb) = 0.5_wp*(p_diag%w_con(jc,jk,jb)+p_diag%w_con(jc,jk+1,jb))
        ENDDO
      ENDDO

      DO jk = 1, nflat
        DO jc = i_startidx, i_endidx
          output_p3c(jc,jk,jb) = 0.5_wp*(p_diag%w_con(jc,jk,jb)+p_diag%w_con(jc,jk+1,jb))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL
      
	  WRITE(*,*) 'checking p3...'
      CALL checkresults(output_p3a, z_hadv_w, &
                        0,SIZE(z_hadv_w,1),SIZE(z_hadv_w,1), &
                        0,SIZE(z_hadv_w,2),SIZE(z_hadv_w,2), &
                        0,SIZE(z_hadv_w,3),SIZE(z_hadv_w,3),8)
      CALL checkresults(output_p3b, p_diag%w_con, &
                        0,SIZE(p_diag%w_con,1),SIZE(p_diag%w_con,1), &
                        0,SIZE(p_diag%w_con,2),SIZE(p_diag%w_con,2), &
                        0,SIZE(p_diag%w_con,3),SIZE(p_diag%w_con,3),9)
      CALL checkresults(output_p3c, z_w_con_c_full, &
                        0,SIZE(z_w_con_c_full,1),SIZE(z_w_con_c_full,1), &
                        0,SIZE(z_w_con_c_full,2),SIZE(z_w_con_c_full,2), &
                        0,SIZE(z_w_con_c_full,3),SIZE(z_w_con_c_full,3),10)
#endif
!===================================================================
				   
    ! Apply cell averaging to the components of horizontal w advection:
    rl_start = 4;
	rl_end = min_rlcell_int;
	i_startblk = p_patch%cells%start_blk(rl_start,1)
	i_endblk = p_patch%cells%end_blk(rl_end,MAX(1,p_patch%n_childdom))

#ifndef TEST	  
	  IF (l_init) THEN
#endif    
	ist=0
    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
    ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
      blockLocalStart(jb) = blockLocalStart(jb) - 1
      blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
	  output_avg = p_diag%ddt_w_adv
#endif
!===================================================================
	
	nameIndex = 25+ntnd

    CALL ocl_cell_avg(z_hadv_w, 7, .TRUE., &
					  istep, &
				      p_int%c_bln_avg, &
					  p_patch%cells%neighbor_idx, &
					  p_patch%cells%neighbor_blk, &
					  blockLocalStart, &
					  blockLocalEnd, &
					  nproma, &
					  1, &
					  nlev, &
					  nlevp1, & !cell_avg wrongly (?) "declares" it as being nlev
					  i_startblk, &
					  i_endblk, &
					  p_patch%nblks_int_c, &
					  p_diag%ddt_w_adv(:,:,:,ntnd), nameIndex, .FALSE.)

#ifndef TEST	  
	  IF (l_init) THEN
#endif					  
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,  stat=ist)
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
    CALL cell_avg(z_hadv_w, p_patch, p_int%c_bln_avg, output_avg(:,:,:,ntnd),&
                  opt_rlstart=4, opt_rlend=min_rlcell_int)  
      
	nameIndex = 35+ntnd
	  WRITE(*,*) 'checking avg...'
      CALL checkresults4d(output_avg, p_diag%ddt_w_adv, &
                        0,SIZE(p_diag%ddt_w_adv,1),SIZE(p_diag%ddt_w_adv,1), &
                        0,SIZE(p_diag%ddt_w_adv,2),SIZE(p_diag%ddt_w_adv,2), &
                        0,SIZE(p_diag%ddt_w_adv,3),SIZE(p_diag%ddt_w_adv,3), &
                        0,SIZE(p_diag%ddt_w_adv,4),SIZE(p_diag%ddt_w_adv,4),nameIndex)
#endif
!===================================================================

    rl_start = grf_bdywidth_e+1
    rl_end = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

#ifndef TEST	  
	  IF (l_init) THEN
#endif
	ist=0
    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
    ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0

    DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
      blockLocalStart(jb) = blockLocalStart(jb) - 1
      blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif

!===================================================================
#ifdef TEST
	  output_p4 = p_diag%ddt_vn_adv
#endif
!===================================================================	
	
    CALL ocl_part4(i_startblk, &
				   i_endblk, &
				   p_patch%nblks_int_e, &
				   p_patch%nblks_int_c, &
				   p_patch%nblks_int_v, &
				   nlev, &
				   nlevp1, &
				   nproma, &
				   ntnd, &
				   blockLocalStart, &
				   blockLocalEnd, &
				   istep, &
				   ividx, &
				   ivblk, &
				   icidx, &
				   icblk, &
				   z_ddxn_ekin_e, &
				   p_diag%vt, &
				   p_patch%edges%f_e, &
				   p_diag%omega_z, &
				   p_int%c_lin_e, &
				   z_w_con_c_full, &
				   p_diag%vn_half, &
				   p_metrics%ddqz_z_full_e, &
				   p_diag%ddt_vn_adv)
				   
#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,  stat=ist)
#ifndef TEST      
	  ENDIF
#endif	
!===================================================================
#ifdef TEST
!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Sum up terms of horizontal wind advection:
      ! grad(Ekin_h) + vt*(f+relvort_e) + wcon_e*dv/dz
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
! unrolling with 6 would be even better, but only if nlev is a multiple of 6
!CDIR UNROLL=2
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          output_p4(je,jk,jb,ntnd) = - ( z_ddxn_ekin_e(je,jk,jb) +           &
            p_diag%vt(je,jk,jb) * ( p_patch%edges%f_e(je,jb) + 0.5_wp*               &
           (p_diag%omega_z(ividx(je,jb,1),jk,ivblk(je,jb,1))   +                      &
            p_diag%omega_z(ividx(je,jb,2),jk,ivblk(je,jb,2))) ) +                     &
           (p_int%c_lin_e(je,1,jb)*z_w_con_c_full(icidx(je,jb,1),jk,icblk(je,jb,1)) + &
            p_int%c_lin_e(je,2,jb)*z_w_con_c_full(icidx(je,jb,2),jk,icblk(je,jb,2)))* &
           (p_diag%vn_half(je,jk,jb) - p_diag%vn_half(je,jk+1,jb))/   &
            p_metrics%ddqz_z_full_e(je,jk,jb) )
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
      
	  WRITE(*,*) 'checking p4...'
      CALL checkresults4d(output_p4, p_diag%ddt_vn_adv, &
                        0,SIZE(p_diag%ddt_vn_adv,1),SIZE(p_diag%ddt_vn_adv,1), &
                        0,SIZE(p_diag%ddt_vn_adv,2),SIZE(p_diag%ddt_vn_adv,2), &
                        0,SIZE(p_diag%ddt_vn_adv,3),SIZE(p_diag%ddt_vn_adv,3), &
                        0,SIZE(p_diag%ddt_vn_adv,4),SIZE(p_diag%ddt_vn_adv,4),12)
#endif
!===================================================================

    rl_start = grf_bdywidth_c+1
    rl_end = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

#ifndef TEST	  
	  IF (l_init) THEN
#endif
	ist=0
    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
    ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
      blockLocalStart(jb) = blockLocalStart(jb) - 1
      blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
	  output_p5 = p_diag%ddt_w_adv
#endif
!===================================================================

    CALL ocl_part5(i_startblk, &
	               i_endblk, &
				   p_patch%nblks_int_c, &
				   nlev, &
				   nlevp1, &
				   nproma, &
				   ntnd, &
				   blockLocalStart, &
				   blockLocalEnd, &
				   istep, &
				   p_diag%w_con, &
				   p_prog%w, 20+i_nnow, &
				   p_metrics%inv_ddqz_z_half2, &
				   p_diag%ddt_w_adv)
#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,  stat=ist)
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Sum up terms of vertical wind advection
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          output_p5(jc,jk,jb,ntnd) = output_p5(jc,jk,jb,ntnd) &
            - p_diag%w_con(jc,jk,jb)*(p_prog%w(jc,jk-1,jb)-p_prog%w(jc,jk+1,jb)) &
            * p_metrics%inv_ddqz_z_half2(jc,jk,jb)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
      
	  WRITE(*,*) 'checking p5...'
      CALL checkresults4d(output_p5, p_diag%ddt_w_adv, &
                        0,SIZE(p_diag%ddt_w_adv,1),SIZE(p_diag%ddt_w_adv,1), &
                        0,SIZE(p_diag%ddt_w_adv,2),SIZE(p_diag%ddt_w_adv,2), &
                        0,SIZE(p_diag%ddt_w_adv,3),SIZE(p_diag%ddt_w_adv,3), &
                        0,SIZE(p_diag%ddt_w_adv,4),SIZE(p_diag%ddt_w_adv,4),13)
#endif
!===================================================================

  END SUBROUTINE velocity_tendencies


  !>
  !! solve_nh
  !!
  !! Version of divergent_modes used for the triangular grid
  !!
  !! @par Revision History
  !! Based on the initial release of divergent_modes by Almut Gassmann (2009-05-12)
  !! Modified by Guenther Zaengl starting on 2010-02-03
  !!
  SUBROUTINE solve_nh (p_nh, p_patch, p_int, nnow, nnew, l_init, l_bdy_nudge, dtime)

    TYPE(t_nh_state),  TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_int_state), TARGET, INTENT(IN)    :: p_int
    TYPE(t_patch),     TARGET, INTENT(IN)    :: p_patch
!     TYPE(t_buffer_memory),     INTENT(INOUT) :: bufr

    ! Initialization switch that has to be .TRUE. at the initial time step only (not for restart)
    LOGICAL,                   INTENT(INOUT) :: l_init
    ! Switch to determine if boundary nudging is executed
    LOGICAL,                   INTENT(IN)    :: l_bdy_nudge
    ! Time levels
    INTEGER,                   INTENT(IN)    :: nnow, nnew
    ! Time step
    REAL(wp),                  INTENT(IN)    :: dtime

    ! Local variables
    INTEGER  :: jb, jk, jc, je
!     INTEGER  :: nlev, nlevp1              !< number of full levels
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER  :: rl_start, rl_end, istep, ntl1, ntl2, nvar
    INTEGER  :: ic, ie

    REAL(wp) :: z_theta_v_fl_e  (nproma,nlev  ,p_patch%nblks_e), &
                z_theta_v_e     (nproma,nlev  ,p_patch%nblks_e), &
                z_rho_e         (nproma,nlev  ,p_patch%nblks_e), &
                z_gradh_exner   (nproma,nlev  ,p_patch%nblks_e), &
                z_concorr_e     (nproma,nlevp1,p_patch%nblks_e)

    REAL(wp), TARGET :: z_vn_avg (nproma,nlev,p_patch%nblks_e)

    REAL(wp) :: z_mass_fl_div   (nproma,nlev  ,p_patch%nblks_c), &
                z_theta_v_fl_div(nproma,nlev  ,p_patch%nblks_c), &
                z_th_ddz_exner_c(nproma,nlevp1,p_patch%nblks_c), &
                z_dexner_dz_c (2,nproma,nlev  ,p_patch%nblks_c), &
                z_exner_ex_pr   (nproma,nlev  ,p_patch%nblks_c), &
                z_exner_pr      (nproma,nlev  ,p_patch%nblks_c)

    REAL(wp) :: z_w_expl        (nproma,nlevp1),          &
                z_contr_w_fl_l  (nproma,nlevp1),          &
                z_alpha         (nproma,nlevp1),          &
                z_beta          (nproma,nlev  ),          &
                z_gamma         (nproma,nlev  ),          &
                z_a             (nproma,nlevp1),          &
                z_b             (nproma,nlevp1),          &
                z_c             (nproma,nlev  ),          &
                z_g             (nproma,nlev  ),          &
                z_q             (nproma,nlev  ),          &
                z_rho_expl      (nproma,nlev  ),          &
                z_exner_expl    (nproma,nlev  ),          &
                z_exner_ic      (nproma,nlevp1),          &
                z_theta_v_pr_mc (nproma,nlev  ),          &
                z_theta_v_pr_ic (nproma,nlevp1),          &
                z_vt_ie         (nproma,nlevp1),          &
                z_thermal_exp   (nproma,p_patch%nblks_c),         &
                z_hydro_corr    (nproma,p_patch%nblks_e)


    REAL(wp):: fac_ex2pres, z_aux(nproma), z_theta1, z_theta2
    INTEGER :: nproma_gradp, nblks_gradp, npromz_gradp, nlen_gradp
    LOGICAL :: lcompute, lcleanup

    ! Preparation for vertical nesting; switch will become a patch attribute later on
    LOGICAL :: l_vert_nested

    ! Pointers to cell indices
    INTEGER,  DIMENSION(:,:,:),   POINTER :: icidx, icblk
    ! Pointers to vertical neighbor indices for pressure gradient computation
    INTEGER,  DIMENSION(:,:,:,:),   POINTER :: ikidx
    ! Pointers to quad edge indices
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iqidx, iqblk
    ! Pointer to velocity field used for mass flux computation
    REAL(wp), DIMENSION(:,:,:),   POINTER :: ptr_vn
    ! Pointers needed for igradp_method = 3
!     INTEGER,  DIMENSION(:),   POINTER :: iplev, ipeidx, ipeblk

	INTEGER :: ist, i_nnow
    INTEGER,  DIMENSION(:), ALLOCATABLE :: blockLocalStart, blockLocalEnd
	
#ifdef TEST
	REAL(wp):: output_c2e(nproma,nlev  ,p_patch%nblks_e)
	REAL(wp):: output_rbf(nproma,nlev  ,p_patch%nblks_e)
	REAL(wp):: output_e2c(nproma,nlevp1,p_patch%nblks_c)
	REAL(wp):: output_div(nproma,nlev  ,p_patch%nblks_c)
	
	REAL(wp):: output_exner_a(nproma,nlev  ,p_patch%nblks_c)
	REAL(wp):: output_exner_b(nproma,nlev  ,p_patch%nblks_c)
	REAL(wp):: output_exner_c(nproma,nlev  ,p_patch%nblks_c)
	REAL(wp):: output_exner_e(2,nproma,nlev,p_patch%nblks_c)
	
	REAL(wp):: output_theta_a(nproma,nlev  ,p_patch%nblks_c)
	REAL(wp):: output_theta_b(nproma,nlevp1,p_patch%nblks_c)
	REAL(wp):: output_theta_c(nproma,nlevp1,p_patch%nblks_c)
	REAL(wp):: output_theta_d(nproma,nlevp1,p_patch%nblks_c)
	REAL(wp):: output_theta_e(nproma,nlevp1,p_patch%nblks_c)
	REAL(wp):: output_theta_f(2,nproma,nlev,p_patch%nblks_c)
	
	REAL(wp):: output_grad_a(nproma,nlev  ,p_patch%nblks_e)
	REAL(wp):: output_grad_b(nproma,p_patch%nblks_e)
	
	REAL(wp):: output_vn(nproma,nlev,p_patch%nblks_e)
	
	REAL(wp):: output_vn_avg(nproma,nlev,p_patch%nblks_e)
	
	REAL(wp):: output_vn_half_a(nproma,nlevp1,p_patch%nblks_e)
	REAL(wp):: output_vn_half_b(nproma,nlevp1,p_patch%nblks_e)
	
	REAL(wp):: output_fluxes_a(nproma,nlev,p_patch%nblks_e)
	REAL(wp):: output_fluxes_b(nproma,nlev,p_patch%nblks_e)
	
	REAL(wp):: output_end_a(nproma,nlevp1,p_patch%nblks_c)
	REAL(wp):: output_end_b(nproma,nlev  ,p_patch%nblks_c)
	REAL(wp):: output_end_c(nproma,nlev  ,p_patch%nblks_c)
	REAL(wp):: output_end_d(nproma,nlev  ,p_patch%nblks_c)
	REAL(wp):: output_end_e(nproma,nlev  ,p_patch%nblks_c)
#endif
    !-------------------------------------------------------------------

!     IF (ltimer) CALL timer_start(timer_solve_nh)

!     IF ((lvert_nest) .AND. (p_patch%nshift > 0)) THEN  
!       l_vert_nested = .TRUE.
!     ELSE
      l_vert_nested = .FALSE.
!     ENDIF

    ! number of vertical levels
!     nlev   = nlev
!     nlevp1 = nlevp1

    ! Set pointers to neighbor cells
    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk

    ! Set pointer to vertical neighbor indices for pressure gradient
    ikidx => p_nh%metrics%vertidx_gradp

    ! Set pointers to quad edges
    iqidx => p_patch%edges%quad_idx
    iqblk => p_patch%edges%quad_blk

    ! Factor needed to convert Exner pressure into "ordinary" pressure
    fac_ex2pres = p0ref/grav/rd_o_cpd

    ! Set pointer to velocity field that is used for mass flux computation
    IF (idiv_method == 1) THEN
      ptr_vn => z_vn_avg
    ELSE
      ptr_vn => p_nh%prog(nnew)%vn
    ENDIF

!     IF (p_test_run) THEN
!       z_rho_e     = 0._wp
!       z_theta_v_e = 0._wp
!     ENDIF

    ! Set time levels of ddt_adv fields for call to velocity_tendencies
    IF (itime_scheme == 6) THEN ! Velocity advection 2nd order in time
      ntl1 = nnow
      ntl2 = nnew
    ELSE                        ! Velocity advection 1st order in time
      ntl1 = 1
      ntl2 = 1
    ENDIF
	
	!WRITE(*,*) 'iterating...'

    i_nchdom   = MAX(1,p_patch%n_childdom)

    DO istep = 1, 2

      IF (istep == 1) THEN ! predictor step
        IF (.NOT.(itime_scheme == 3 .AND. .NOT. l_init)) THEN
          CALL velocity_tendencies(p_nh%prog(nnow),p_patch,p_int,p_nh%metrics,&
                                   p_nh%diag,ntl1,istep,l_init)
        ENDIF
        nvar = nnow
		i_nnow = 1
      ELSE                 ! corrector step
        CALL velocity_tendencies(p_nh%prog(nnew),p_patch,p_int,p_nh%metrics,&
                                 p_nh%diag,ntl2,istep,l_init)
        nvar = nnew
		i_nnow = 0
      ENDIF

      ! Compute rho and theta at edges
      IF (istep == 1) THEN
	    i_startblk = p_patch%edges%start_blk(2,1)
		i_endblk = p_patch%edges%end_blk(min_rledge,MAX(1,p_patch%n_childdom))

#ifndef TEST	  
	  IF (l_init) THEN
#endif		
		ist=0
        ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
        ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0

        DO jb = i_startblk, i_endblk
          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             blockLocalStart(jb), blockLocalEnd(jb), 2, min_rledge)
          blockLocalStart(jb) = blockLocalStart(jb) - 1
          blockLocalEnd(jb)   = blockLocalEnd(jb)
        ENDDO
#ifndef TEST      
	  ENDIF
#endif

!===================================================================
#ifdef TEST
	  output_c2e = z_rho_e
#endif
!===================================================================

        ! density at edges
	    CALL ocl_cells2edges_scalar(p_nh%prog(nnow)%rho, 9, .FALSE., &
									istep, &
		                            p_int%c_lin_e, &
									p_patch%edges%cell_idx, &
									p_patch%edges%cell_blk, &
									blockLocalStart, blockLocalEnd, nproma, &
									1, nlev, nlev, &
									i_startblk, i_endblk, p_patch%nblks_int_c, p_patch%nblks_int_e, &
									z_rho_e, 10, .TRUE.)						
	
!===================================================================
#ifdef TEST
	  CALL cells2edges_scalar(p_nh%prog(nnow)%rho,p_patch,p_int%c_lin_e,output_c2e)
      
	  WRITE(*,*) 'checking c2e...'
      CALL checkresults(output_c2e, z_rho_e, &
                        0,SIZE(z_rho_e,1),SIZE(z_rho_e,1), &
                        0,SIZE(z_rho_e,2),SIZE(z_rho_e,2), &
                        0,SIZE(z_rho_e,3),SIZE(z_rho_e,3),14)
						
	  output_c2e = z_theta_v_e
#endif
!===================================================================
									
        ! virtual potential temperature at edges
	    CALL ocl_cells2edges_scalar(p_nh%prog(nnow)%theta_v, 11, .FALSE., &
									istep, &
		                            p_int%c_lin_e, &
									p_patch%edges%cell_idx, &
									p_patch%edges%cell_blk, &
									blockLocalStart, blockLocalEnd, nproma, &
									1, nlev, nlev, &
									i_startblk, i_endblk, p_patch%nblks_int_c, p_patch%nblks_int_e, &
									z_theta_v_e, 12, .TRUE.)
!===================================================================
#ifdef TEST
      CALL cells2edges_scalar(p_nh%prog(nnow)%theta_v,p_patch,p_int%c_lin_e,output_c2e)
      
	  WRITE(*,*) 'checking c2e...'
      CALL checkresults(output_c2e, z_theta_v_e, &
                        0,SIZE(z_theta_v_e,1),SIZE(z_theta_v_e,1), &
                        0,SIZE(z_theta_v_e,2),SIZE(z_theta_v_e,2), &
                        0,SIZE(z_theta_v_e,3),SIZE(z_theta_v_e,3),14)
#endif
!===================================================================

#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
        DEALLOCATE(blockLocalStart,stat=ist)
        DEALLOCATE(blockLocalEnd,stat=ist)
#ifndef TEST      
	  ENDIF
#endif

!       ENDIF

      ENDIF ! istep = 1

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk,jb,i_startidx,i_endidx)

    rl_start = 3
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    ! Computations at cell points; to be executed in predictor step only
    IF (istep == 1) THEN
#ifndef TEST	  
	  IF (l_init) THEN
#endif
	  ist=0
      ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
      ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0
 	
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
        blockLocalStart(jb) = blockLocalStart(jb) - 1
        blockLocalEnd(jb)   = blockLocalEnd(jb)
	  ENDDO
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
	output_exner_a = z_exner_ex_pr
	output_exner_b = z_exner_pr
	output_exner_c = p_nh%diag%exner_old
	output_exner_e = z_dexner_dz_c
#endif
!===================================================================

      CALL ocl_exner(i_startblk, i_endblk, p_patch%nblks_c, &
	                 nflat, nlev, nlevp1, &
				     blockLocalStart, blockLocalEnd, nproma, &
					 istep, &
				     p_nh%metrics%exner_ref_mc, &
				     p_nh%metrics%exner_exfac, &
				     p_nh%prog(nnow)%exner, &
				     p_nh%diag%exner_old, &
				     z_exner_ex_pr, &
				     z_exner_pr, &
				     p_nh%metrics%wgtfacq1_c, &
				     p_nh%metrics%wgtfacq_c, &
				     p_nh%metrics%wgtfac_c, &
				     p_nh%metrics%inv_ddqz_z_full, &
				     z_dexner_dz_c)

#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
      DEALLOCATE(blockLocalStart,stat=ist)
      DEALLOCATE(blockLocalEnd,stat=ist)
#ifndef TEST      
	  ENDIF
#endif
!===================================================================
#ifdef TEST

!$OMP DO PRIVATE(jk,jc,z_exner_ic)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

!CDIR UNROLL=2
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            ! extrapolated perturbation Exner pressure (used for horizontal gradients only)
            output_exner_a(jc,jk,jb) = - p_nh%metrics%exner_ref_mc(jc,jk,jb) +              &
              (1._wp + p_nh%metrics%exner_exfac(jc,jk,jb))*p_nh%prog(nnow)%exner(jc,jk,jb) &
                     - p_nh%metrics%exner_exfac(jc,jk,jb) *output_exner_c(jc,jk,jb)

            ! non-extrapolated perturbation Exner pressure
            output_exner_b(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb) - &
              p_nh%metrics%exner_ref_mc(jc,jk,jb)

            ! Now save current time level in exner_old
            output_exner_c(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb)
          ENDDO
        ENDDO

        ! Perturbation Exner pressure on bottom and top half level
        IF (nflat == 1) THEN
          DO jc = i_startidx, i_endidx
          z_exner_ic(jc,1) =                                          &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*output_exner_a(jc,1,jb) + &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*output_exner_a(jc,2,jb) + &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*output_exner_a(jc,3,jb)
          ENDDO
        ENDIF
        DO jc = i_startidx, i_endidx
          z_exner_ic(jc,nlevp1) =                                         &
            p_nh%metrics%wgtfacq_c(jc,1,jb)*output_exner_a(jc,nlev  ,jb) + &
            p_nh%metrics%wgtfacq_c(jc,2,jb)*output_exner_a(jc,nlev-1,jb) + &
            p_nh%metrics%wgtfacq_c(jc,3,jb)*output_exner_a(jc,nlev-2,jb)
        ENDDO

!CDIR UNROLL=3
        DO jk = nlev, MAX(2,nflat), -1
          DO jc = i_startidx, i_endidx
            ! Exner pressure on remaining half levels for metric correction term
            z_exner_ic(jc,jk) =                                              &
              p_nh%metrics%wgtfac_c(jc,jk,jb)*output_exner_a(jc,jk,jb) +      &
              (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*output_exner_a(jc,jk-1,jb)

            ! First vertical derivative of perturbation Exner pressure
            output_exner_e(1,jc,jk,jb) =                    &
             (z_exner_ic(jc,jk) - z_exner_ic(jc,jk+1)) *   &
              p_nh%metrics%inv_ddqz_z_full(jc,jk,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
      
	  WRITE(*,*) 'checking exner...'
      CALL checkresults(output_exner_a, z_exner_ex_pr, &
                        0,SIZE(z_exner_ex_pr,1),SIZE(z_exner_ex_pr,1), &
                        0,SIZE(z_exner_ex_pr,2),SIZE(z_exner_ex_pr,2), &
                        0,SIZE(z_exner_ex_pr,3),SIZE(z_exner_ex_pr,3),15)
      CALL checkresults(output_exner_b, z_exner_pr, &
                        0,SIZE(z_exner_pr,1),SIZE(z_exner_pr,1), &
                        0,SIZE(z_exner_pr,2),SIZE(z_exner_pr,2), &
                        0,SIZE(z_exner_pr,3),SIZE(z_exner_pr,3),16)
      CALL checkresults(output_exner_c, p_nh%diag%exner_old, &
                        0,SIZE(p_nh%diag%exner_old,1),SIZE(p_nh%diag%exner_old,1), &
                        0,SIZE(p_nh%diag%exner_old,2),SIZE(p_nh%diag%exner_old,2), &
                        0,SIZE(p_nh%diag%exner_old,3),SIZE(p_nh%diag%exner_old,3), &
						17)
      CALL checkresults4d(output_exner_e, z_dexner_dz_c, &
                        0,SIZE(z_dexner_dz_c,1),SIZE(z_dexner_dz_c,1), &
                        0,SIZE(z_dexner_dz_c,2),SIZE(z_dexner_dz_c,2), &
                        0,SIZE(z_dexner_dz_c,3),SIZE(z_dexner_dz_c,3), &
                        0,SIZE(z_dexner_dz_c,4),SIZE(z_dexner_dz_c,4),18)
#endif
!===================================================================
	
	ENDIF ! istep = 1

    IF (igradp_method == 1) THEN
      rl_end   = min_rlcell_int
      i_endblk = p_patch%cells%end_blk(rl_end,i_nchdom)
    ENDIF

#ifndef TEST	  
	  IF (l_init) THEN
#endif	
	ist=0
    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
    ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0
	 
	DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
      blockLocalStart(jb) = blockLocalStart(jb) - 1
      blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
	output_theta_b = p_nh%diag%rho_ic
	output_theta_d = p_nh%diag%theta_v_h
	output_theta_e = z_th_ddz_exner_c
	output_theta_f = z_dexner_dz_c
#endif
!===================================================================
	
    CALL ocl_theta(i_startblk, &
	               i_endblk, &
				   p_patch%nblks_c, &
				   nflat_gradp(p_patch%id), &
				   nlev, &
				   nlevp1, &
				   blockLocalStart, &
				   blockLocalEnd, &
				   nproma, &
				   istep, &
				   p_nh%prog(nvar)%theta_v, 22+i_nnow, &
				   p_nh%prog(nnow)%theta_v, &
				   p_nh%metrics%theta_ref_mc, &
				   p_nh%metrics%rho_refcorr_ic, &
				   p_nh%metrics%wgtfac_c, &
				   p_nh%prog(nvar)%rho, 24+i_nnow, &
				   p_nh%prog(nnow)%rho, &
				   p_nh%metrics%vwind_expl_wgt, &
				   z_exner_pr, &
				   p_nh%metrics%ddqz_z_half, &
				   p_nh%metrics%d_exner_dz_ref_ic, &
				   p_nh%diag%rho_ic, &
				   z_th_ddz_exner_c, &
				   p_nh%metrics%wgtfacq_c, &
				   p_nh%metrics%theta_ref_ic, &
				   p_nh%diag%theta_v_h, &
				   p_nh%metrics%d_exner_dz_ref_mc, &
				   p_nh%metrics%d2_exner_dz2_ref_mc, &
				   z_dexner_dz_c)

#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,stat=ist)
#ifndef TEST      
	  ENDIF
#endif

!===================================================================
#ifdef TEST
!$OMP DO PRIVATE(jk,jc,z_theta_v_pr_mc,z_theta_v_pr_ic)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      z_theta_v_pr_mc(i_startidx:i_endidx,1) =  0.5_wp *     &
        (p_nh%prog(nnow)%theta_v(i_startidx:i_endidx,1,jb) + &
        p_nh%prog(nvar)%theta_v(i_startidx:i_endidx,1,jb)) - &
        p_nh%metrics%theta_ref_mc(i_startidx:i_endidx,1,jb)

!CDIR UNROLL=8
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          ! density at interface levels for vertical flux divergence computation
          output_theta_b(jc,jk,jb) = p_nh%metrics%rho_refcorr_ic(jc,jk,jb)+0.5_wp*( &
            p_nh%metrics%wgtfac_c(jc,jk,jb)*(p_nh%prog(nnow)%rho(jc,jk,jb) +          &
            p_nh%prog(nvar)%rho(jc,jk,jb))+(1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*   &
            (p_nh%prog(nnow)%rho(jc,jk-1,jb)+p_nh%prog(nvar)%rho(jc,jk-1,jb)) )

          ! perturbation virtual potential temperature at main levels
          z_theta_v_pr_mc(jc,jk) = 0.5_wp*(p_nh%prog(nnow)%theta_v(jc,jk,jb) +     &
            p_nh%prog(nvar)%theta_v(jc,jk,jb)) - p_nh%metrics%theta_ref_mc(jc,jk,jb)

          ! perturbation virtual potential temperature at interface levels
          z_theta_v_pr_ic(jc,jk) = &
            p_nh%metrics%wgtfac_c(jc,jk,jb)*z_theta_v_pr_mc(jc,jk) +       &
            (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_theta_v_pr_mc(jc,jk-1)

          ! virtual potential temperature at interface levels
          output_theta_d(jc,jk,jb) = p_nh%metrics%theta_ref_ic(jc,jk,jb) + &
            z_theta_v_pr_ic(jc,jk)

          ! vertical pressure gradient * theta_v
          output_theta_e(jc,jk,jb) = p_nh%metrics%vwind_expl_wgt(jc,jb)* &
            output_theta_d(jc,jk,jb) * (z_exner_pr(jc,jk-1,jb)-       &
            z_exner_pr(jc,jk,jb)) / p_nh%metrics%ddqz_z_half(jc,jk,jb) +   &
           z_theta_v_pr_ic(jc,jk)*p_nh%metrics%d_exner_dz_ref_ic(jc,jk,jb)

        ENDDO
      ENDDO

      IF (istep == 1) THEN

        ! Perturbation theta at top and surface levels
        DO jc = i_startidx, i_endidx
          z_theta_v_pr_ic(jc,1)      = 0._wp
          z_theta_v_pr_ic(jc,nlevp1) =                                   &
            p_nh%metrics%wgtfacq_c(jc,1,jb)*z_theta_v_pr_mc(jc,nlev) +   &
            p_nh%metrics%wgtfacq_c(jc,2,jb)*z_theta_v_pr_mc(jc,nlev-1) + &
            p_nh%metrics%wgtfacq_c(jc,3,jb)*z_theta_v_pr_mc(jc,nlev-2)
            output_theta_d(jc,nlevp1,jb) =                                  &
              p_nh%metrics%theta_ref_ic(jc,nlevp1,jb) + z_theta_v_pr_ic(jc,nlevp1)
        ENDDO

!CDIR UNROLL=3
        DO jk = nflat_gradp(p_patch%id), nlev
          DO jc = i_startidx, i_endidx
            ! Second vertical derivative of perturbation Exner pressure (hydrostatic approximation)
            output_theta_f(2,jc,jk,jb) = -0.5_wp / p_nh%prog(nnow)%theta_v(jc,jk,jb) * &
             ((z_theta_v_pr_ic(jc,jk) - z_theta_v_pr_ic(jc,jk+1)) *                   &
              p_nh%metrics%d_exner_dz_ref_mc(jc,jk,jb) + z_theta_v_pr_mc(jc,jk)*      &
              p_nh%metrics%d2_exner_dz2_ref_mc(jc,jk,jb))
          ENDDO
        ENDDO

      ENDIF ! istep == 1

    ENDDO
!$OMP END DO
      
	  WRITE(*,*) 'checking theta...'
      CALL checkresults(output_theta_b, p_nh%diag%rho_ic, &
                        0,SIZE(p_nh%diag%rho_ic,1),SIZE(p_nh%diag%rho_ic,1), &
                        0,SIZE(p_nh%diag%rho_ic,2),SIZE(p_nh%diag%rho_ic,2), &
                        0,SIZE(p_nh%diag%rho_ic,3),SIZE(p_nh%diag%rho_ic,3),19)
      CALL checkresults(output_theta_d, p_nh%diag%theta_v_h, &
                        0,SIZE(p_nh%diag%theta_v_h,1),SIZE(p_nh%diag%theta_v_h,1), &
                        0,SIZE(p_nh%diag%theta_v_h,2),SIZE(p_nh%diag%theta_v_h,2), &
                        0,SIZE(p_nh%diag%theta_v_h,3),SIZE(p_nh%diag%theta_v_h,3),20)
      CALL checkresults(output_theta_e, z_th_ddz_exner_c, &
                        0,SIZE(z_th_ddz_exner_c,1),SIZE(z_th_ddz_exner_c,1), &
                        0,SIZE(z_th_ddz_exner_c,2),SIZE(z_th_ddz_exner_c,2), &
                        0,SIZE(z_th_ddz_exner_c,3),SIZE(z_th_ddz_exner_c,3),21)
      CALL checkresults4d(output_theta_f, z_dexner_dz_c, &
                        0,SIZE(z_dexner_dz_c,1),SIZE(z_dexner_dz_c,1), &
                        0,SIZE(z_dexner_dz_c,2),SIZE(z_dexner_dz_c,2), &
                        0,SIZE(z_dexner_dz_c,3),SIZE(z_dexner_dz_c,3), &
                        0,SIZE(z_dexner_dz_c,4),SIZE(z_dexner_dz_c,4),18)
#endif
!===================================================================

    ! Computations at edge points

    rl_start = grf_bdywidth_e + 1   ! boundary update follows below
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

#ifndef TEST	  
	  IF (l_init) THEN
#endif
	ist=0
    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
    ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0
	 
	DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
      blockLocalStart(jb) = blockLocalStart(jb) - 1
      blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif

    IF (istep == 1) THEN

!===================================================================
#ifdef TEST
	output_grad_a = z_gradh_exner
	output_grad_b = z_hydro_corr
#endif
!===================================================================
	!where is z_hydro_corr used?!
	  CALL ocl_grad(i_startblk, i_endblk, &
					p_patch%nblks_c, p_patch%nblks_e, &
					nflat, nflat_gradp(p_patch%id), &
					nlev, nlevp1, &
					blockLocalStart, blockLocalEnd, nproma, &
					igradp_method, istep, &
					icidx, icblk, ikidx, &
					grav_o_cpd, &
					p_patch%edges%inv_dual_edge_length, &
					z_exner_ex_pr, &
					z_dexner_dz_c, &
					p_nh%metrics%ddxn_z_full, &
					p_int%c_lin_e, &
					p_nh%metrics%zdiff_gradp, &
					p_nh%prog(nnow)%theta_v, &
					p_nh%diag%theta_v_h, &
					p_nh%metrics%inv_ddqz_z_full, &
					z_gradh_exner, &
					z_hydro_corr)
!===================================================================
#ifdef TEST
!$OMP DO PRIVATE(jk,je,z_theta1,z_theta2)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nflat-1
#else
!CDIR UNROLL=_URD-_URD2
        DO jk = 1, nflat-1
          DO je = i_startidx, i_endidx
#endif
            ! horizontal gradient of Exner pressure where coordinate surfaces are flat
            output_grad_a(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)* &
             (z_exner_ex_pr(icidx(je,jb,2),jk,icblk(je,jb,2)) -                  &
              z_exner_ex_pr(icidx(je,jb,1),jk,icblk(je,jb,1)) )
          ENDDO
        ENDDO

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = nflat, nflat_gradp(p_patch%id)
#else
!CDIR UNROLL=6
        DO jk = nflat, nflat_gradp(p_patch%id)
          DO je = i_startidx, i_endidx
#endif
            ! horizontal gradient of Exner pressure, including metric correction
            output_grad_a(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)*         &
             (z_exner_ex_pr(icidx(je,jb,2),jk,icblk(je,jb,2)) -                          &
              z_exner_ex_pr(icidx(je,jb,1),jk,icblk(je,jb,1)) ) -                        &
              p_nh%metrics%ddxn_z_full(je,jk,jb) *                                       &
             (p_int%c_lin_e(je,1,jb)*z_dexner_dz_c(1,icidx(je,jb,1),jk,icblk(je,jb,1)) + &
              p_int%c_lin_e(je,2,jb)*z_dexner_dz_c(1,icidx(je,jb,2),jk,icblk(je,jb,2)))

          ENDDO
        ENDDO

        IF (igradp_method >= 2) THEN
        ! remark: loop exchange is not beneficial here because of 3D indirect addressing
!CDIR UNROLL=5
          DO jk = nflat_gradp(p_patch%id)+1, nlev
            DO je = i_startidx, i_endidx

              ! horizontal gradient of Exner pressure, Taylor-expansion-based reconstruction
              output_grad_a(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)*         &
               (z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2)) +           &
                p_nh%metrics%zdiff_gradp(2,je,jk,jb)*                                      &
               (z_dexner_dz_c(1,icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2)) +         &
                p_nh%metrics%zdiff_gradp(2,je,jk,jb)*                                      &
                z_dexner_dz_c(2,icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2))) -        &
               (z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)) +           &
                p_nh%metrics%zdiff_gradp(1,je,jk,jb)*                                      &
               (z_dexner_dz_c(1,icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)) +         &
                p_nh%metrics%zdiff_gradp(1,je,jk,jb)*                                      &
                z_dexner_dz_c(2,icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)))))

            ENDDO
          ENDDO
        ENDIF

        IF (igradp_method == 3) THEN
        ! compute hydrostatically approximated correction term that replaces downward extrapolation
          
          DO je = i_startidx, i_endidx

            z_theta1 = &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb),icblk(je,jb,1)) + &
              p_nh%metrics%zdiff_gradp(1,je,nlev,jb)*                                      &
             (p_nh%diag%theta_v_h(icidx(je,jb,1),ikidx(1,je,nlev,jb),  icblk(je,jb,1)) -   &
              p_nh%diag%theta_v_h(icidx(je,jb,1),ikidx(1,je,nlev,jb)+1,icblk(je,jb,1))) *  &
              p_nh%metrics%inv_ddqz_z_full(icidx(je,jb,1),ikidx(1,je,nlev,jb),icblk(je,jb,1))

            z_theta2 = &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb),icblk(je,jb,2)) + &
              p_nh%metrics%zdiff_gradp(2,je,nlev,jb)*                                      &
             (p_nh%diag%theta_v_h(icidx(je,jb,2),ikidx(2,je,nlev,jb),  icblk(je,jb,2)) -   &
              p_nh%diag%theta_v_h(icidx(je,jb,2),ikidx(2,je,nlev,jb)+1,icblk(je,jb,2))) *  &
              p_nh%metrics%inv_ddqz_z_full(icidx(je,jb,2),ikidx(2,je,nlev,jb),icblk(je,jb,2))

            output_grad_b(je,jb) = grav_o_cpd*p_patch%edges%inv_dual_edge_length(je,jb)*    &
              (z_theta2-z_theta1)*4._wp/(z_theta1+z_theta2)**2

            ENDDO
        ENDIF

      ENDDO
!$OMP END DO
      
	  WRITE(*,*) 'checking grad...'
      CALL checkresults(output_grad_a, z_gradh_exner, &
                        0,SIZE(z_gradh_exner,1),SIZE(z_gradh_exner,1), &
                        0,SIZE(z_gradh_exner,2),SIZE(z_gradh_exner,2), &
                        0,SIZE(z_gradh_exner,3),SIZE(z_gradh_exner,3),22)
      CALL checkresults(output_grad_b, z_hydro_corr, &
                        0,SIZE(z_hydro_corr,1),SIZE(z_hydro_corr,1), &
                        0,SIZE(z_hydro_corr,2),SIZE(z_hydro_corr,2), &
						0,1,1,23)
#endif
!===================================================================
    ENDIF ! istep = 1

!===================================================================
#ifdef TEST
	output_vn = p_nh%prog(nnew)%vn
#endif
!===================================================================

    ! Update horizontal velocity field
	CALL ocl_vn(i_startblk, i_endblk, p_patch%nblks_e, &
				nlev, &
				blockLocalStart, blockLocalEnd, nproma, &
				itime_scheme, istep, &
				ntl1, ntl2, &
				cpd, dtime, &
				p_nh%prog(nnow)%vn, &
				p_nh%diag%ddt_vn_adv, &
				p_nh%diag%ddt_vn_phy, &
				z_theta_v_e, &
				z_gradh_exner, &
				p_nh%prog(nnew)%vn)

!===================================================================
#ifdef TEST
!$OMP DO PRIVATE(jk,je,ic)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (itime_scheme == 6 .AND. istep == 2) THEN
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            output_vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb)+ dtime     &
            & *(0.5_wp*(p_nh%diag%ddt_vn_adv(je,jk,jb,ntl1)                        &
            & +p_nh%diag%ddt_vn_adv(je,jk,jb,ntl2))+p_nh%diag%ddt_vn_phy(je,jk,jb) &
            & -cpd*z_theta_v_e(je,jk,jb)*z_gradh_exner(je,jk,jb))
          ENDDO
        ENDDO
      ELSE
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            output_vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb)+ dtime     &
            & *(p_nh%diag%ddt_vn_adv(je,jk,jb,ntl1)+p_nh%diag%ddt_vn_phy(je,jk,jb) &
            & -cpd*z_theta_v_e(je,jk,jb)*z_gradh_exner(je,jk,jb))
          ENDDO
        ENDDO
      ENDIF
    ENDDO
!$OMP END DO
      
	  WRITE(*,*) 'checking vn...'
      CALL checkresults(output_vn, p_nh%prog(nnew)%vn, &
                        0,SIZE(p_nh%prog(nnew)%vn,1),SIZE(p_nh%prog(nnew)%vn,1), &
                        0,SIZE(p_nh%prog(nnew)%vn,2),SIZE(p_nh%prog(nnew)%vn,2), &
                        0,SIZE(p_nh%prog(nnew)%vn,3),SIZE(p_nh%prog(nnew)%vn,3),24)
#endif
!===================================================================

#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,stat=ist)
#ifndef TEST      
	  ENDIF
#endif

    rl_start = 2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%nblks_int_e

#ifndef TEST	  
	  IF (l_init) THEN
#endif
	ist=0
    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
    ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0
	 
	DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         blockLocalStart(jb), blockLocalEnd(jb), rl_start)
      blockLocalStart(jb) = blockLocalStart(jb) - 1
      blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
	  output_rbf = p_nh%diag%vt
#endif
!===================================================================
	
    ! Tangential wind component using RBF reconstruction
	CALL ocl_rbf_vec_interpol_edge(p_nh%prog(nnew)%vn, 13, .TRUE., &
								   istep, &
	                               p_int%rbf_vec_coeff_e, &
								   p_int%rbf_vec_idx_e, &
								   p_int%rbf_vec_blk_e, &
								   i_startblk, &
								   nproma, &
								   nlev, &
								   rbf_vec_dim_e, &
								   p_patch%nblks_int_e, &
								   p_nh%diag%vt, 1, .FALSE., &
								   blockLocalStart, &
								   blockLocalEnd, &
								   0, nlev)
								   
!===================================================================
#ifdef TEST
    CALL rbf_vec_interpol_edge(p_nh%prog(nnew)%vn, p_patch, p_int, output_rbf)
      
	  WRITE(*,*) 'checking rbf...'
      CALL checkresults(output_rbf, p_nh%diag%vt, &
                        0,SIZE(p_nh%diag%vt,1),SIZE(p_nh%diag%vt,1), &
                        0,SIZE(p_nh%diag%vt,2),SIZE(p_nh%diag%vt,2), &
                        0,SIZE(p_nh%diag%vt,3),SIZE(p_nh%diag%vt,3),0)
#endif
!===================================================================

#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,stat=ist)
#ifndef TEST      
	  ENDIF
#endif

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

    IF (idiv_method == 1) THEN
	  rl_start = 2
      rl_end   = min_rledge_int - 2

      i_startblk = p_patch%edges%start_blk(rl_start,1)
      i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

#ifndef TEST	  
	  IF (l_init) THEN
#endif	
	  ist=0
      ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
      ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0
	 
	  DO jb = i_startblk, i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
        blockLocalStart(jb) = blockLocalStart(jb) - 1
        blockLocalEnd(jb)   = blockLocalEnd(jb)
	  ENDDO
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
	  output_vn_avg = z_vn_avg
#endif
!===================================================================

	  CALL ocl_vn_avg(i_startblk, i_endblk, p_patch%nblks_e, &
	  				  nlev, &
					  blockLocalStart, blockLocalEnd, nproma, &
					  istep, &
					  iqidx, iqblk, &
					  p_nh%prog(nnew)%vn, &
					  p_int%e_flx_avg, &
					  z_vn_avg)

#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
	  DEALLOCATE(blockLocalStart,stat=ist)
	  DEALLOCATE(blockLocalEnd,stat=ist)
#ifndef TEST      
	  ENDIF
#endif	  
!===================================================================
#ifdef TEST
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
!CDIR UNROLL=3
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          output_vn_avg(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)*p_int%e_flx_avg(je,1,jb)        &
            + p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1))*p_int%e_flx_avg(je,2,jb) &
            + p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2))*p_int%e_flx_avg(je,3,jb) &
            + p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3))*p_int%e_flx_avg(je,4,jb) &
            + p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))*p_int%e_flx_avg(je,5,jb)
         ENDDO
      ENDDO
    ENDDO
!$OMP END DO
      
	  WRITE(*,*) 'checking vn_avg...'
      CALL checkresults(output_vn_avg, z_vn_avg, &
                        0,SIZE(z_vn_avg,1),SIZE(z_vn_avg,1), &
                        0,SIZE(z_vn_avg,2),SIZE(z_vn_avg,2), &
                        0,SIZE(z_vn_avg,3),SIZE(z_vn_avg,3),25)
#endif
!===================================================================
    ENDIF

    rl_start = 3
    rl_end = min_rledge_int - 2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

#ifndef TEST	  
	  IF (l_init) THEN
#endif	
	ist=0
    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
    ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0
	 
	DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
      blockLocalStart(jb) = blockLocalStart(jb) - 1
      blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif

!===================================================================
#ifdef TEST
	output_vn_half_a = p_nh%diag%vn_half
	output_vn_half_b = z_concorr_e
#endif
!===================================================================

    CALL ocl_vn_half(i_startblk, i_endblk, p_patch%nblks_e, &
	                 nflat, nlev, nlevp1, &
					 blockLocalStart, blockLocalEnd, nproma, &
					 istep, l_vert_nested, &
					 p_nh%metrics%wgtfac_e, &
					 p_nh%metrics%wgtfacq_e, &
					 p_nh%metrics%wgtfacq1_e, &
					 p_nh%prog(nnew)%vn, &
					 p_nh%diag%vt, &
					 p_nh%diag%vn_half, &
					 p_nh%metrics%ddxn_z_half, &
					 p_nh%metrics%ddxt_z_half, &
					 z_concorr_e)

#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,stat=ist)
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,z_vt_ie)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Interpolate vn and vt to interface levels
!CDIR UNROLL=6
      DO jk = nflat, nlev
        DO je = i_startidx, i_endidx
          output_vn_half_a(je,jk,jb) = &
            p_nh%metrics%wgtfac_e(je,jk,jb)*p_nh%prog(nnew)%vn(je,jk,jb) +        &
           (1._wp - p_nh%metrics%wgtfac_e(je,jk,jb))*p_nh%prog(nnew)%vn(je,jk-1,jb)
          z_vt_ie(je,jk) =                                                  &
            p_nh%metrics%wgtfac_e(je,jk,jb)*p_nh%diag%vt(je,jk,jb) +        &
           (1._wp - p_nh%metrics%wgtfac_e(je,jk,jb))*p_nh%diag%vt(je,jk-1,jb)
        ENDDO
      ENDDO

        ! Reduced computations required for flat levels
      IF (istep == 1 .AND. .NOT. l_vert_nested) THEN
        DO je = i_startidx, i_endidx
            output_vn_half_a(je,1,jb) =                                     &
              p_nh%metrics%wgtfacq1_e(je,1,jb)*p_nh%prog(nnew)%vn(je,1,jb) + &
              p_nh%metrics%wgtfacq1_e(je,2,jb)*p_nh%prog(nnew)%vn(je,2,jb) + &
              p_nh%metrics%wgtfacq1_e(je,3,jb)*p_nh%prog(nnew)%vn(je,3,jb)
        ENDDO
      ENDIF

      IF (istep == 1) THEN
        DO jk = 2, nflat-1
          DO je = i_startidx, i_endidx
            output_vn_half_a(je,jk,jb) = &
              p_nh%metrics%wgtfac_e(je,jk,jb)*p_nh%prog(nnew)%vn(je,jk,jb) + &
             (1._wp - p_nh%metrics%wgtfac_e(je,jk,jb))*p_nh%prog(nnew)%vn(je,jk-1,jb)
          ENDDO
        ENDDO
      ENDIF

      ! Bottom level
      DO je = i_startidx, i_endidx
        output_vn_half_a(je,nlevp1,jb) =                                    &
          p_nh%metrics%wgtfacq_e(je,1,jb)*p_nh%prog(nnew)%vn(je,nlev,jb)   + &
          p_nh%metrics%wgtfacq_e(je,2,jb)*p_nh%prog(nnew)%vn(je,nlev-1,jb) + &
          p_nh%metrics%wgtfacq_e(je,3,jb)*p_nh%prog(nnew)%vn(je,nlev-2,jb)
        z_vt_ie(je,nlevp1) =                                           &
          p_nh%metrics%wgtfacq_e(je,1,jb)*p_nh%diag%vt(je,nlev,jb) +   &
          p_nh%metrics%wgtfacq_e(je,2,jb)*p_nh%diag%vt(je,nlev-1,jb) + &
          p_nh%metrics%wgtfacq_e(je,3,jb)*p_nh%diag%vt(je,nlev-2,jb)
      ENDDO

      ! Compute contravariant correction for vertical velocity at interface levels
      DO jk = nflat+1, nlevp1
        DO je = i_startidx, i_endidx
          output_vn_half_b(je,jk,jb) =                                            &
            output_vn_half_a(je,jk,jb)*p_nh%metrics%ddxn_z_half(je,jk,jb) + &
            z_vt_ie(je,jk)*p_nh%metrics%ddxt_z_half(je,jk,jb)
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
      
	  WRITE(*,*) 'checking vn_half...'
      CALL checkresults(output_vn_half_a, p_nh%diag%vn_half, &
                        0,SIZE(p_nh%diag%vn_half,1),SIZE(p_nh%diag%vn_half,1), &
                        0,SIZE(p_nh%diag%vn_half,2),SIZE(p_nh%diag%vn_half,2), &
                        0,SIZE(p_nh%diag%vn_half,3),SIZE(p_nh%diag%vn_half,3),26)
      CALL checkresults(output_vn_half_b, z_concorr_e, &
                        0,SIZE(z_concorr_e,1),SIZE(z_concorr_e,1), &
                        0,SIZE(z_concorr_e,2),SIZE(z_concorr_e,2), &
                        0,SIZE(z_concorr_e,3),SIZE(z_concorr_e,3),27)
#endif
!===================================================================

    rl_start = 7
    IF (idiv_method == 1) THEN
      rl_end = min_rledge_int - 2
    ELSE
      rl_end = min_rledge_int - 3
    ENDIF

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

#ifndef TEST	  
	  IF (l_init) THEN
#endif
	ist=0
    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
    ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0
	 
	DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
      blockLocalStart(jb) = blockLocalStart(jb) - 1
      blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
	output_fluxes_a = p_nh%diag%mass_fl_e
	output_fluxes_b = z_theta_v_fl_e
#endif
!===================================================================

	CALL ocl_fluxes(i_startblk, i_endblk, p_patch%nblks_e, &
	                nlev, &
					blockLocalStart, blockLocalEnd, nproma, &
					istep, &
					z_rho_e, &
					ptr_vn, idiv_method, &
					p_nh%metrics%ddqz_z_full_e, &
					z_theta_v_e, &
					p_nh%diag%mass_fl_e, &
					z_theta_v_fl_e)

#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,stat=ist)
#ifndef TEST      
	  ENDIF
#endif

!===================================================================
#ifdef TEST
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Fluxes at edges
      DO jk = 1,nlev
        DO je = i_startidx, i_endidx

          output_fluxes_a(je,jk,jb) = z_rho_e(je,jk,jb)         &
            * ptr_vn(je,jk,jb) * p_nh%metrics%ddqz_z_full_e(je,jk,jb)

          output_fluxes_b(je,jk,jb)= output_fluxes_a(je,jk,jb)   &
            * z_theta_v_e(je,jk,jb)

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
      
	  WRITE(*,*) 'checking fluxes...'
      CALL checkresults(output_fluxes_a, p_nh%diag%mass_fl_e, &
                        0,SIZE(p_nh%diag%mass_fl_e,1),SIZE(p_nh%diag%mass_fl_e,1), &
                        0,SIZE(p_nh%diag%mass_fl_e,2),SIZE(p_nh%diag%mass_fl_e,2), &
                        0,SIZE(p_nh%diag%mass_fl_e,3),SIZE(p_nh%diag%mass_fl_e,3),28)
      CALL checkresults(output_fluxes_b, z_theta_v_fl_e, &
                        0,SIZE(z_theta_v_fl_e,1),SIZE(z_theta_v_fl_e,1), &
                        0,SIZE(z_theta_v_fl_e,2),SIZE(z_theta_v_fl_e,2), &
                        0,SIZE(z_theta_v_fl_e,3),SIZE(z_theta_v_fl_e,3),29)
#endif
!===================================================================

	IF (i_cell_type == 3) THEN
	  rl_start = 3
	  rl_end = min_rlcell_int-1
	  
	  i_startblk = p_patch%cells%start_blk(rl_start,1)
	  i_endblk   = p_patch%cells%end_blk(rl_end,MAX(1,p_patch%n_childdom))
	  
	  ist=0
	  ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
	  ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0
	  
	  DO jb = i_startblk, i_endblk
	    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
		                   blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
        blockLocalStart(jb) = blockLocalStart(jb) - 1
        blockLocalEnd(jb)   = blockLocalEnd(jb)
	  ENDDO
	
!===================================================================
#ifdef TEST
	  output_e2c = p_nh%diag%w_concorr_c
#endif
!===================================================================

	  CALL ocl_edges2cells_scalar_tri(z_concorr_e, 5, .TRUE., &
									  istep, &
									  p_int%e_bln_c_s, &
									  p_patch%cells%edge_idx, &
									  p_patch%cells%edge_blk, &
									  blockLocalStart, &
									  blockLocalEnd, &
									  nproma, &
									  nflat+1, &
									  nlevp1, &
									  nlevp1, &
									  p_patch%cells%start_blk(rl_start,1), &
									  p_patch%cells%end_blk(rl_end,MAX(1,p_patch%n_childdom)), &
									  p_patch%nblks_int_c, &
									  p_patch%nblks_int_e, &
									  p_nh%diag%w_concorr_c, 6, .FALSE.)
									  
									  
!===================================================================
#ifdef TEST
      CALL edges2cells_scalar(z_concorr_e,p_patch,p_int%e_bln_c_s,output_e2c,&
                              nflat+1, nlevp1, opt_rlstart=3, opt_rlend=min_rlcell_int-1)
      
	  WRITE(*,*) 'checking e2c...'
      CALL checkresults(output_e2c, p_nh%diag%w_concorr_c, &
                        0,SIZE(p_nh%diag%w_concorr_c,1),SIZE(p_nh%diag%w_concorr_c,1), &
                        0,SIZE(p_nh%diag%w_concorr_c,2),SIZE(p_nh%diag%w_concorr_c,2), &
                        0,SIZE(p_nh%diag%w_concorr_c,3),SIZE(p_nh%diag%w_concorr_c,3),5)
#endif
!===================================================================
				   
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,stat=ist)
	
	ELSE
	  WRITE(*,*) 'TO MOVE ON GPU!'
      CALL edges2cells_scalar(z_concorr_e,p_patch,p_int%e_bln_c_s,p_nh%diag%w_concorr_c,&
                              nflat+1, nlevp1, opt_rlstart=3, opt_rlend=min_rlcell_int-1)
	ENDIF


    IF (idiv_method == 1) THEN ! use simple divergence based on averaged velocity

	  IF (i_cell_type == 3) THEN
	    rl_start = 4
	    rl_end = min_rlcell_int
	  
	    i_startblk = p_patch%cells%start_blk(rl_start,1)
	    i_endblk   = p_patch%cells%end_blk(rl_end,MAX(1,p_patch%n_childdom))

#ifndef TEST	  
	  IF (l_init) THEN
#endif	  
	    ist=0
	    ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
	    ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0

	    DO jb = i_startblk, i_endblk
	      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
		                     blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
          blockLocalStart(jb) = blockLocalStart(jb) - 1
          blockLocalEnd(jb)   = blockLocalEnd(jb)
	    ENDDO
#ifndef TEST      
	  ENDIF
#endif

!===================================================================
#ifdef TEST
	  output_div = z_mass_fl_div
#endif
!===================================================================
	
	    CALL ocl_div_tri(p_nh%diag%mass_fl_e, 14, .TRUE., &
					     istep, &
		                 p_int%geofac_div, &
						 p_patch%cells%edge_idx, &
						 p_patch%cells%edge_blk, &
						 i_startblk, &
						 i_endblk, &
						 nproma, &
						 nlev, &
						 p_patch%nblks_int_c, &
						 p_patch%nblks_int_e, &
						 z_mass_fl_div, 15, .TRUE., &
						 blockLocalStart, &
						 blockLocalEnd)
						 
!===================================================================
#ifdef TEST
	  CALL div(p_nh%diag%mass_fl_e, p_patch, p_int, output_div, &
                 opt_rlstart=4, opt_rlend=min_rlcell_int)
      
	  WRITE(*,*) 'checking div...'
      CALL checkresults(output_div, z_mass_fl_div, &
                        0,SIZE(z_mass_fl_div,1),SIZE(z_mass_fl_div,1), &
                        0,SIZE(z_mass_fl_div,2),SIZE(z_mass_fl_div,2), &
                        0,SIZE(z_mass_fl_div,3),SIZE(z_mass_fl_div,3),30)
						
	  output_div = z_theta_v_fl_div
#endif
!===================================================================
						 

		CALL ocl_div_tri(z_theta_v_fl_e, 16, .TRUE., &
					     istep, &
		                 p_int%geofac_div, &
						 p_patch%cells%edge_idx, &
						 p_patch%cells%edge_blk, &
						 i_startblk, &
						 i_endblk, &
						 nproma, &
						 nlev, &
						 p_patch%nblks_int_c, &
						 p_patch%nblks_int_e, &
						 z_theta_v_fl_div, 17, .TRUE., &
						 blockLocalStart, &
						 blockLocalEnd)
						 
!===================================================================
#ifdef TEST
	  CALL div(z_theta_v_fl_e, p_patch, p_int, output_div, &
                 opt_rlstart=4, opt_rlend=min_rlcell_int)
      
	  WRITE(*,*) 'checking div...'
      CALL checkresults(output_div, z_theta_v_fl_div, &
                        0,SIZE(z_theta_v_fl_div,1),SIZE(z_theta_v_fl_div,1), &
                        0,SIZE(z_theta_v_fl_div,2),SIZE(z_theta_v_fl_div,2), &
                        0,SIZE(z_theta_v_fl_div,3),SIZE(z_theta_v_fl_div,3),30)
#endif
!===================================================================
#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
	    DEALLOCATE(blockLocalStart,stat=ist)
	    DEALLOCATE(blockLocalEnd,stat=ist)
#ifndef TEST      
	  ENDIF
#endif	

	  ELSE
	    WRITE(*,*) 'TO MOVE ON GPU!'
		CALL div(p_nh%diag%mass_fl_e, p_patch, p_int, z_mass_fl_div, opt_in2=z_theta_v_fl_e, &
                 opt_out2=z_theta_v_fl_div, opt_rlstart=4, opt_rlend=min_rlcell_int)
	  ENDIF


    ELSE ! use averaged divergence
	  WRITE(*,*) 'TO MOVE ON GPU!'

      ! horizontal divergences of rho and rhotheta are processed in one step for efficiency
      CALL div_avg(p_nh%diag%mass_fl_e, p_patch, p_int, p_int%c_bln_avg, z_mass_fl_div, &
                   opt_in2=z_theta_v_fl_e, opt_out2=z_theta_v_fl_div, opt_rlstart=4,    &
                   opt_rlend=min_rlcell_int)
    ENDIF
	
    ! Vertical solution:
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

#ifndef TEST	  
	  IF (l_init) THEN
#endif	
	ist=0
	ALLOCATE(blockLocalStart(i_endblk-i_startblk+1), stat=ist); blockLocalStart=0
	ALLOCATE(blockLocalEnd(i_endblk-i_startblk+1), stat=ist); blockLocalEnd=0
	
	DO jb = i_startblk, i_endblk
	  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
						 blockLocalStart(jb), blockLocalEnd(jb), rl_start, rl_end)
	  blockLocalStart(jb) = blockLocalStart(jb) - 1
	  blockLocalEnd(jb)   = blockLocalEnd(jb)
	ENDDO
#ifndef TEST      
	  ENDIF
#endif
	
!===================================================================
#ifdef TEST
	output_end_a = p_nh%prog(nnew)%w
	output_end_b = p_nh%prog(nnew)%rho
	output_end_c = p_nh%prog(nnew)%exner
	output_end_d = p_nh%prog(nnew)%rhotheta_v
	output_end_e = p_nh%prog(nnew)%theta_v
#endif
!===================================================================

	CALL ocl_endphase(i_startblk, i_endblk, p_patch%nblks_c, &
	                  nlev, nlevp1, nrdmax, &
					  blockLocalStart, blockLocalEnd, nproma, &
					  ntl1, ntl2, istep, &
					  itime_scheme, &
					  dtime, cpd, cvd, rd, cvd_o_rd, &
					  p_nh%prog(nnow)%w, &
					  p_nh%prog(nnew)%w, &
					  p_nh%diag%ddt_w_adv, &
					  z_th_ddz_exner_c, &
					  p_nh%diag%rho_ic, &
					  p_nh%diag%w_concorr_c, &
					  p_nh%metrics%vwind_expl_wgt, &
					  p_nh%metrics%vwind_impl_wgt, &
					  p_nh%prog(nnow)%exner, &
					  p_nh%prog(nnew)%exner, &
					  p_nh%prog(nnow)%rhotheta_v, &
					  p_nh%prog(nnew)%rhotheta_v, &
					  p_nh%metrics%ddqz_z_full, &
					  p_nh%diag%theta_v_h, &
					  p_nh%metrics%ddqz_z_half, &
					  p_nh%prog(nnow)%rho, &
					  p_nh%prog(nnew)%rho, &
					  p_nh%metrics%inv_ddqz_z_full, &
					  z_mass_fl_div, &
					  z_exner_pr, &
					  z_theta_v_fl_div, &
					  p_nh%diag%ddt_exner, &
					  p_nh%diag%ddt_exner_phy, &
					  p_nh%metrics%rayleigh_w, &
					  p_nh%metrics%exner_ref_mc, &
					  p_nh%prog(nnew)%theta_v)

#ifndef TEST	  
	  IF (l_init) THEN
#endif				   
	DEALLOCATE(blockLocalStart,stat=ist)
	DEALLOCATE(blockLocalEnd,stat=ist)
#ifndef TEST      
	  ENDIF
#endif	
!===================================================================

#ifdef TEST
      !$OMP DO PRIVATE(jk,jc,z_w_expl,z_contr_w_fl_l,z_rho_expl,z_exner_expl,z_a,z_b,z_c,&
!$OMP            z_g,z_q,z_alpha,z_beta,z_gamma,z_aux,ic)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (istep == 2 .AND. itime_scheme == 6) THEN
!CDIR UNROLL=5
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx

            ! explicit part for w
            z_w_expl(jc,jk) = p_nh%prog(nnow)%w(jc,jk,jb) + dtime *                          &
              (0.5_wp*(p_nh%diag%ddt_w_adv(jc,jk,jb,ntl1)+p_nh%diag%ddt_w_adv(jc,jk,jb,ntl2))&
              -cpd*z_th_ddz_exner_c(jc,jk,jb))

            ! contravariant vertical velocity times density for explicit part
            z_contr_w_fl_l(jc,jk) = p_nh%diag%rho_ic(jc,jk,jb)*(-p_nh%diag%w_concorr_c(jc,jk,jb) &
              + p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,jk,jb) )
          ENDDO
        ENDDO
      ELSE
!CDIR UNROLL=5
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx

            ! explicit part for w
            z_w_expl(jc,jk) = p_nh%prog(nnow)%w(jc,jk,jb) + dtime *            &
              (p_nh%diag%ddt_w_adv(jc,jk,jb,ntl1)-cpd*z_th_ddz_exner_c(jc,jk,jb))

            ! contravariant vertical velocity times density for explicit part
            z_contr_w_fl_l(jc,jk) = p_nh%diag%rho_ic(jc,jk,jb)*(-p_nh%diag%w_concorr_c(jc,jk,jb) &
              + p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,jk,jb) )
          ENDDO
        ENDDO
      ENDIF

      ! Solver coefficients
!CDIR UNROLL=3
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          z_beta(jc,jk)=dtime*rd*p_nh%prog(nnow)%exner(jc,jk,jb)         &
          &                   /(cvd*p_nh%prog(nnow)%rhotheta_v(jc,jk,jb) &
          &                   *p_nh%metrics%ddqz_z_full(jc,jk,jb))

          z_alpha(jc,jk)= p_nh%metrics%vwind_impl_wgt(jc,jb)*         &
            &  p_nh%diag%theta_v_h(jc,jk,jb)*p_nh%diag%rho_ic(jc,jk,jb)
        ENDDO
      ENDDO

      z_alpha(:,nlevp1) = 0.0_wp
!CDIR UNROLL=2
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          z_gamma(jc,jk) =  dtime*cpd*p_nh%metrics%vwind_impl_wgt(jc,jb)* &
            p_nh%diag%theta_v_h(jc,jk,jb)/p_nh%metrics%ddqz_z_half(jc,jk,jb)

          z_a(jc,jk) = -z_gamma(jc,jk)*z_beta(jc,jk-1)*z_alpha(jc,jk-1)
          z_c(jc,jk) = -z_gamma(jc,jk)*z_beta(jc,jk  )*z_alpha(jc,jk+1)
          z_b(jc,jk) = 1.0_wp+z_gamma(jc,jk)*z_alpha(jc,jk) &
            *(z_beta(jc,jk-1)+z_beta(jc,jk))
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        z_q(jc,2) = -z_c(jc,2)/z_b(jc,2)
      ENDDO

!CDIR UNROLL=4
      DO jk = 3, nlev
        DO jc = i_startidx, i_endidx
          z_g(jc,jk) = 1.0_wp/(z_b(jc,jk)+z_a(jc,jk)*z_q(jc,jk-1))
          z_q(jc,jk) = - z_c(jc,jk)*z_g(jc,jk)
        ENDDO
      ENDDO

        output_end_a(:,1,jb) = 0._wp
        z_contr_w_fl_l(:,1)       = 0._wp

      ! lower boundary condition for w, consistent with contravariant correction
      DO jc = i_startidx, i_endidx
        output_end_a(jc,nlevp1,jb) = p_nh%diag%w_concorr_c(jc,nlevp1,jb)
        z_contr_w_fl_l(jc,nlevp1)       = 0.0_wp
      ENDDO


      ! other full level stuff
      ! Top level first
      DO jc = i_startidx, i_endidx
        z_rho_expl(jc,1)=        p_nh%prog(nnow)%rho(jc,1,jb) &
        &        -dtime*p_nh%metrics%inv_ddqz_z_full(jc,1,jb) &
        &                            *(z_mass_fl_div(jc,1,jb) &
        &                            +z_contr_w_fl_l(jc,1   ) &
        &                            -z_contr_w_fl_l(jc,2   ))

        z_exner_expl(jc,1)=             z_exner_pr(jc,1,jb)                    &
        &      -z_beta (jc,1)*(z_theta_v_fl_div(jc,1,jb)                       &
        & +p_nh%diag%theta_v_h(jc,1,jb)*z_contr_w_fl_l(jc,1)                   &
        & -p_nh%diag%theta_v_h(jc,2,jb)*z_contr_w_fl_l(jc,2))                  &
        & +dtime*(p_nh%diag%ddt_exner(jc,1,jb)+p_nh%diag%ddt_exner_phy(jc,1,jb))
      ENDDO

      ! Other levels
      DO jk = 2, nlev
!CDIR ON_ADB(output_end_a)
        DO jc = i_startidx, i_endidx
          z_rho_expl(jc,jk)=       p_nh%prog(nnow)%rho(jc,jk  ,jb) &
          &        -dtime*p_nh%metrics%inv_ddqz_z_full(jc,jk  ,jb) &
          &                            *(z_mass_fl_div(jc,jk  ,jb) &
          &                            +z_contr_w_fl_l(jc,jk     ) &
          &                             -z_contr_w_fl_l(jc,jk+1   ))

          z_exner_expl(jc,jk)=          z_exner_pr(jc,jk,jb) - z_beta(jc,jk)         &
          &                             *(z_theta_v_fl_div(jc,jk,jb)                 &
          &   +p_nh%diag%theta_v_h(jc,jk  ,jb)*z_contr_w_fl_l(jc,jk  )               &
          &   -p_nh%diag%theta_v_h(jc,jk+1,jb)*z_contr_w_fl_l(jc,jk+1))              &
          &   +dtime*(p_nh%diag%ddt_exner(jc,jk,jb)+p_nh%diag%ddt_exner_phy(jc,jk,jb))

          output_end_a(jc,jk,jb) = z_w_expl(jc,jk) - z_gamma(jc,jk)  &
          &      *(z_exner_expl(jc,jk-1)-z_exner_expl(jc,jk))
        ENDDO
      ENDDO

      ! Solve tridiagonal matrix for w
!CDIR ON_ADB(output_end_a)
      DO jc = i_startidx, i_endidx
        output_end_a(jc,2,jb)= output_end_a(jc,2,jb)/z_b(jc,2)
      ENDDO

      DO jk = 3, nlev
!CDIR ON_ADB(output_end_a)
        DO jc = i_startidx, i_endidx
          output_end_a(jc,jk,jb) = (output_end_a(jc,jk,jb)  &
            -z_a(jc,jk)*output_end_a(jc,jk-1,jb))*z_g(jc,jk)
        ENDDO
      ENDDO

      DO jk = nlev-1, 2, -1
!CDIR ON_ADB(output_end_a)
        DO jc = i_startidx, i_endidx
          output_end_a(jc,jk,jb) = output_end_a(jc,jk,jb)&
          &             +output_end_a(jc,jk+1,jb)*z_q(jc,jk)
        ENDDO
      ENDDO

      ! Rayleigh damping mechanism (Klemp,Dudhia,Hassiotis:MWR136,pp.3987-4004)
      DO jk = 2, nrdmax
!CDIR ON_ADB(output_end_a)
        DO jc = i_startidx, i_endidx
          output_end_a(jc,jk,jb) =                                &
            (output_end_a(jc,jk,jb)-output_end_a(jc,1,jb)) / &
            (1.0_wp+dtime*p_nh%metrics%rayleigh_w(jc,jk,jb))
        ENDDO
      ENDDO

      ! Results
      DO jk = 1, nlev
!CDIR ON_ADB(output_end_a)
        DO jc = i_startidx, i_endidx

          ! density
          output_end_b(jc,jk,jb) = z_rho_expl(jc,jk)              &
            - p_nh%metrics%vwind_impl_wgt(jc,jb)*dtime                   &
            * p_nh%metrics%inv_ddqz_z_full(jc,jk,jb)                     &
            *(p_nh%diag%rho_ic(jc,jk  ,jb)*output_end_a(jc,jk  ,jb) &
            - p_nh%diag%rho_ic(jc,jk+1,jb)*output_end_a(jc,jk+1,jb))

          ! exner
          output_end_c(jc,jk,jb) = z_exner_expl(jc,jk) &
            + p_nh%metrics%exner_ref_mc(jc,jk,jb)-z_beta(jc,jk) &
            *(z_alpha(jc,jk  )*output_end_a(jc,jk  ,jb)    &
            - z_alpha(jc,jk+1)*output_end_a(jc,jk+1,jb))

          ! rho*theta
          output_end_d(jc,jk,jb) = p_nh%prog(nnow)%rhotheta_v(jc,jk,jb)   &
            *( (output_end_c(jc,jk,jb)/p_nh%prog(nnow)%exner(jc,jk,jb)-1.0_wp) &
            *   cvd_o_rd+1.0_wp)

          ! theta
          output_end_e(jc,jk,jb) = &
            output_end_d(jc,jk,jb)/output_end_b(jc,jk,jb)

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO

	  WRITE(*,*) 'checking endphase...'
      CALL checkresults(output_end_a, p_nh%prog(nnew)%w, &
                        0,SIZE(p_nh%prog(nnew)%w,1),SIZE(p_nh%prog(nnew)%w,1), &
                        0,SIZE(p_nh%prog(nnew)%w,2),SIZE(p_nh%prog(nnew)%w,2), &
                        0,SIZE(p_nh%prog(nnew)%w,3),SIZE(p_nh%prog(nnew)%w,3),31)
      CALL checkresults(output_end_b, p_nh%prog(nnew)%rho, &
                        0,SIZE(p_nh%prog(nnew)%rho,1),SIZE(p_nh%prog(nnew)%rho,1), &
                        0,SIZE(p_nh%prog(nnew)%rho,2),SIZE(p_nh%prog(nnew)%rho,2), &
                        0,SIZE(p_nh%prog(nnew)%rho,3),SIZE(p_nh%prog(nnew)%rho,3),32)
      CALL checkresults(output_end_c, p_nh%prog(nnew)%exner, &
                        0,SIZE(p_nh%prog(nnew)%exner,1),SIZE(p_nh%prog(nnew)%exner,1), &
                        0,SIZE(p_nh%prog(nnew)%exner,2),SIZE(p_nh%prog(nnew)%exner,2), &
                        0,SIZE(p_nh%prog(nnew)%exner,3),SIZE(p_nh%prog(nnew)%exner,3),33)
      CALL checkresults(output_end_d, p_nh%prog(nnew)%rhotheta_v, &
                        0,SIZE(p_nh%prog(nnew)%rhotheta_v,1),SIZE(p_nh%prog(nnew)%rhotheta_v,1), &
                        0,SIZE(p_nh%prog(nnew)%rhotheta_v,2),SIZE(p_nh%prog(nnew)%rhotheta_v,2), &
                        0,SIZE(p_nh%prog(nnew)%rhotheta_v,3),SIZE(p_nh%prog(nnew)%rhotheta_v,3), &
						34)
      CALL checkresults(output_end_e, p_nh%prog(nnew)%theta_v, &
                        0,SIZE(p_nh%prog(nnew)%theta_v,1),SIZE(p_nh%prog(nnew)%theta_v,1), &
                        0,SIZE(p_nh%prog(nnew)%theta_v,2),SIZE(p_nh%prog(nnew)%theta_v,2), &
                        0,SIZE(p_nh%prog(nnew)%theta_v,3),SIZE(p_nh%prog(nnew)%theta_v,3),35)
#endif
!===================================================================
        l_init = .FALSE. ! should be .TRUE. only at initial predictor step

      ENDDO ! istep-loop
	
	
!     IF (ltimer) CALL timer_stop(timer_solve_nh)

    ! The remaining computations are needed for MPI-parallelized applications only
    IF (p_nprocs == 1) RETURN

	WRITE(*,*) 'TO PUT ON GPU!'
	
    rl_start = min_rlcell_int - 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

!$OMP DO PRIVATE(jk,jc)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
!           IF (p_nh%metrics%mask_prog_halo_c(jc,jb)) THEN

            p_nh%prog(nnew)%rhotheta_v(jc,jk,jb) = p_nh%prog(nnow)%rhotheta_v(jc,jk,jb)   &
              *( (p_nh%prog(nnew)%exner(jc,jk,jb)/p_nh%prog(nnow)%exner(jc,jk,jb)-1.0_wp) &
              *   cvd_o_rd+1.0_wp)

            p_nh%prog(nnew)%theta_v(jc,jk,jb) = &
              p_nh%prog(nnew)%rhotheta_v(jc,jk,jb)/p_nh%prog(nnew)%rho(jc,jk,jb)
!           ENDIF
        ENDDO
      ENDDO
!$OMP END DO
    ENDDO
!$OMP END PARALLEL


  END SUBROUTINE solve_nh

END MODULE mo_solve_nonhydro

