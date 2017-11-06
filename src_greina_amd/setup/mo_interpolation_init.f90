!>
!! Contains the implementation of interpolation and reconstruction
!! routines used by the shallow water model, including the RBF
!! reconstruction routines.
!!
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_interpolation_init
  !-------------------------------------------------------------------------

  USE mo_kind,               ONLY: wp
  USE mo_math_constants,      ONLY: pi2, pi_2
  USE mo_physical_constants,  ONLY: re
  USE mo_impl_constants,     ONLY: min_rlcell, min_rledge, min_rlvert, max_dom
  USE mo_model_domain,       ONLY: t_patch
  USE mo_math_utilities,     ONLY: gc2cc, cc2gc, gnomonic_proj,               &
                                & gvec2cvec, cvec2gvec,                      &
                                & t_cartesian_coordinates,                   &
                                & rotate_latlon, arc_length,                 &
                                & t_geographical_coordinates
  USE mo_global_variables,  ONLY: i_cell_type, nproma, msg_level
  USE mo_impl_constants,    ONLY: SUCCESS
  USE mo_intp_data_strc
  USE mo_loopindices,       ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_exception,         ONLY: message, finish
  USE mo_intp_rbf_coeffs

  !DR for testing purposes
!   USE mo_math_utilities,      ONLY: t_geographical_coordinates

  IMPLICIT NONE
  

  PRIVATE
  
  PUBLIC :: construct_int_state, destruct_int_state
  

CONTAINS


!------------------------------------------------------
SUBROUTINE construct_int_state(ptr_patch, ptr_int )
  TYPE(t_patch), INTENT(inout) :: ptr_patch
  TYPE(t_int_state), INTENT(inout) :: ptr_int

! LSQ reconstruction at cell center
llsq_high_consv  = .TRUE.   ! conservative reconstruction
lsq_high_ord     = 3        ! cubic polynomial
! RBF vector reconstruction at cell centers
rbf_vec_dim_c   = 9         ! use 2nd order reconstruction
rbf_vec_kern_c  = 1         ! use Gaussian kernel
! RBF vector reconstruction at vertices
rbf_vec_dim_v   = 6         ! use 6-point reconstruction
rbf_vec_kern_v  = 1         ! use Gaussian kernel
! RBF vector reconstruction at edge midpoints
rbf_vec_dim_e   = 4         ! use 4-point reconstruction
rbf_vec_kern_e  = 3         ! use Inverse multiquadric kernel
! Stencil size for reconstruction of gradient at cell midpoints
rbf_c2grad_dim  = 10

! Initialize namelist fields for scaling factors (dimension 1:depth) with dummy values
! A meaningful initialization follows after reading the namelist
rbf_vec_scale_c(:) = -1.0_wp
rbf_vec_scale_v(:) = -1.0_wp
rbf_vec_scale_e(:) = -1.0_wp

! Initialize the namelist for the method for the vorticity flux term
i_cori_method = 3
sick_a = 0.75_wp
sick_o = 1.0_wp-sick_a
l_corner_vort=.TRUE.

! Coefficients for lateral boundary nudging
nudge_max_coeff   = 0.02_wp  ! Maximum nudging coefficient
nudge_efold_width = 2._wp    ! e-folding width in units of cell rows
nudge_zone_width  = 8        ! Width of nudging zone in units of cell rows


   CALL allocate_int_state(ptr_patch, ptr_int)
  ! initializion of coefficients for averaging of scalars and kinetic energy
  !
   CALL scalar_int_coeff(ptr_patch, ptr_int)

   CALL complete_patchinfo( ptr_patch, ptr_int)
   CALL init_geo_factors(ptr_patch, ptr_int)
   CALL init_cellavg_wgt(ptr_patch, ptr_int)
   CALL bln_int_coeff_e2c( ptr_patch, ptr_int)
    ! ... at cell centers
    CALL rbf_vec_index_cell (ptr_patch, ptr_int)
    !
    CALL rbf_vec_compute_coeff_cell (ptr_patch, ptr_int)
    !
    ! ... at triangle vertices
    CALL rbf_vec_index_vertex (ptr_patch, ptr_int)
    !
    CALL rbf_vec_compute_coeff_vertex (ptr_patch, ptr_int)
    !
    ! ... at edge midpoints
    CALL rbf_vec_index_edge (ptr_patch, ptr_int)
    !
    CALL rbf_vec_compute_coeff_edge (ptr_patch, ptr_int)
    !
    ! Compute coefficients needed for gradient reconstruction at cell midpoints
    CALL rbf_c2grad_index (ptr_patch, ptr_int)
    CALL rbf_compute_coeff_c2grad (ptr_patch, ptr_int)

    !
    ! initialization of quadrature points and weights
    !
    CALL tri_quadrature_pts (ptr_patch, ptr_int)
   

END SUBROUTINE construct_int_state
!------------------------------------------------------

!------------------------------------------------------
SUBROUTINE init_geo_factors( ptr_patch, ptr_int )
  !
  !  patch on which computation is performed
  !
  TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

  ! Interpolation state
  TYPE(t_int_state), INTENT(inout):: ptr_int
  !

  INTEGER :: jc, jb, je, jv, je1, jn,jm, nincr
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom, nblks_e

  INTEGER :: ile, ibe, ilc1, ibc1, ilc2, ibc2, ifac, ic, ilnc, ibnc, &
            ilv1, ilv2, ibv1, ibv2, jr1, ilr, ibr
  TYPE(t_cartesian_coordinates)::z_pn_k,z_pn_j,z_pt_k,z_cart_no,z_cart_ea
  REAL(wp) :: z_proj, z_norm, z_lon, z_lat

  i_nchdom=1
  !-------------------------
  ! The calculation cannot be done for boundary edges
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk,nblks_e,ifac)
  i_startblk = ptr_patch%edges%start_blk(2,1)
  nblks_e  = ptr_patch%nblks_int_e
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ilc2,ibc1,ibc2,ilv1,ilv2,&
!$OMP            ibv1,ibv2)
  DO jb = i_startblk, nblks_e

    CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, 2)

    DO je = i_startidx, i_endidx

      ! inverse distance averaging (former subroutine cell2edge_lin_int_coeff)
      ! For the hexagonal grid, this is also the direct distance averaging
      ! because the edge is exactly half way between the cells
      ! (needed for proper bracket formalism)
      ptr_int%c_lin_e(je,1,jb) = ptr_patch%edges%edge_cell_length(je,jb,2)/&
                                           ptr_patch%edges%dual_edge_length(je,jb)
      ptr_int%c_lin_e(je,2,jb) = ptr_patch%edges%edge_cell_length(je,jb,1)/&
                                           ptr_patch%edges%dual_edge_length(je,jb)
    ENDDO

  ENDDO
!$OMP END DO
!-----------------------------------------------------------

! a) Geometrical factor for divergence
rl_start = 1
rl_end = min_rlcell

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      
!
! loop through all patch cells (and blocks)
!
!$OMP DO PRIVATE(jb,je,jc,i_startidx,i_endidx,ile,ibe)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je = 1, i_cell_type
    DO jc = i_startidx, i_endidx

      IF (je > ptr_patch%cells%num_edges(jc,jb)) CYCLE ! relevant for hexagons

      ile = ptr_patch%cells%edge_idx(jc,jb,je)
      ibe = ptr_patch%cells%edge_blk(jc,jb,je)

      ptr_int%geofac_div(jc,je,jb) = &
        ptr_patch%edges%primal_edge_length(ile,ibe) * &
        ptr_patch%cells%edge_orientation(jc,jb,je)  / &
        ptr_patch%cells%area(jc,jb)
      
    ENDDO !cell loop
  ENDDO

END DO !block loop
!$OMP END DO
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! b) Geometrical factor for rotation
rl_start = 2
rl_end = min_rlvert

! Vorticity should have the right sign
  ifac = 1
  ! values for the blocking
  i_startblk = ptr_patch%verts%start_blk(rl_start,1)
  i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)
  !
  ! loop through all patch cells (and blocks)
  !
!$OMP DO PRIVATE(jb,je,jv,i_startidx,i_endidx,ile,ibe)
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je = 1, 9-i_cell_type
      DO jv = i_startidx, i_endidx

        IF (je > ptr_patch%verts%num_edges(jv,jb)) CYCLE

        ile = ptr_patch%verts%edge_idx(jv,jb,je)
        ibe = ptr_patch%verts%edge_blk(jv,jb,je)

        ptr_int%geofac_rot(jv,je,jb) =                &
          ptr_patch%edges%dual_edge_length(ile,ibe) * &
          ptr_patch%verts%edge_orientation(jv,jb,je)/ &
          ptr_patch%verts%dual_area(jv,jb) * REAL(ifac,wp)

      ENDDO !vertex loop
    ENDDO

  END DO !block loop
!$OMP END DO


!-----------------------------------------------------------
! c) Geometrical factor for nabla2_scalar
rl_start = 2
rl_end = min_rlcell

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!
!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
!$OMP    ilc2,ibc2,ilnc,ibnc)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je = 1, i_cell_type
    DO jc = i_startidx, i_endidx

      ile = ptr_patch%cells%edge_idx(jc,jb,je)
      ibe = ptr_patch%cells%edge_blk(jc,jb,je)

      ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
      ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
      ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
      ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

      IF (jc == ilc1 .AND. jb == ibc1) THEN
        IF (i_cell_type == 3) THEN
          ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) - &
            ptr_int%geofac_div(jc,je,jb) /                            &
            ptr_patch%edges%dual_edge_length(ile,ibe)
        ELSE IF (i_cell_type == 6) THEN
          ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) - &
            ptr_int%geofac_div(jc,je,jb) /                            &
            ptr_patch%edges%dual_edge_length(ile,ibe)*                &
            ptr_patch%edges%system_orientation(ile,ibe)
        ENDIF
      ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
        IF (i_cell_type == 3) THEN
          ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) + &
            ptr_int%geofac_div(jc,je,jb) /                            &
            ptr_patch%edges%dual_edge_length(ile,ibe)
        ELSE IF (i_cell_type == 6) THEN
          ptr_int%geofac_n2s(jc,1,jb) = ptr_int%geofac_n2s(jc,1,jb) + &
            ptr_int%geofac_div(jc,je,jb) /                            &
            ptr_patch%edges%dual_edge_length(ile,ibe)*                &
            ptr_patch%edges%system_orientation(ile,ibe)
        ENDIF
      ENDIF
      DO ic = 1, i_cell_type
        ilnc = ptr_patch%cells%neighbor_idx(jc,jb,ic)
        ibnc = ptr_patch%cells%neighbor_blk(jc,jb,ic)
        IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
          IF (i_cell_type == 3) THEN
            ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) - &
              ptr_int%geofac_div(jc,je,jb) /                                  &
              ptr_patch%edges%dual_edge_length(ile,ibe)
          ELSE IF (i_cell_type == 6) THEN
            ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) - &
              ptr_int%geofac_div(jc,je,jb) /                                  &
              ptr_patch%edges%dual_edge_length(ile,ibe)*                      &
              ptr_patch%edges%system_orientation(ile,ibe)
          ENDIF
        ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
          IF (i_cell_type == 3) THEN
            ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) + &
              ptr_int%geofac_div(jc,je,jb) /                                  &
              ptr_patch%edges%dual_edge_length(ile,ibe)
          ELSE IF (i_cell_type == 6) THEN
            ptr_int%geofac_n2s(jc,ic+1,jb) = ptr_int%geofac_n2s(jc,ic+1,jb) + &
              ptr_int%geofac_div(jc,je,jb) /                                  &
              ptr_patch%edges%dual_edge_length(ile,ibe)*                      &
              ptr_patch%edges%system_orientation(ile,ibe)
          ENDIF
        ENDIF
      ENDDO

      ! To ensure that dummy edges have a factor of 0:
      IF (je > ptr_patch%cells%num_edges(jc,jb)) THEN
        ptr_int%geofac_n2s(jc,je+1,jb) = 0._wp
      ENDIF

    ENDDO !cell loop
  ENDDO

END DO !block loop
!$OMP END DO



!-----------------------------------------------------------
  rl_start = 2
  rl_end = min_rledge

  ! values for the blocking
  i_startblk = ptr_patch%edges%start_blk(rl_start,1)
  i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,je,je1,i_startidx,i_endidx,ile,ibe)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je1 = 1, 4
      DO je = i_startidx, i_endidx

      ile = ptr_patch%edges%quad_idx(je,jb,je1)
      ibe = ptr_patch%edges%quad_blk(je,jb,je1)

      ptr_int%geofac_qdiv(je,je1,jb) = &
        ptr_patch%edges%primal_edge_length(ile,ibe) * &
        ptr_patch%edges%quad_orientation(je,jb,je1)  / &
        ptr_patch%edges%quad_area(je,jb)

      ENDDO !edge loop
    ENDDO

  END DO !block loop
!$OMP END DO

! f) Geometrical factor for Green-Gauss gradient
rl_start = 2
rl_end = min_rlcell

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!
!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
!$OMP    ilc2,ibc2,ilnc,ibnc)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je = 1, i_cell_type
    DO jc = i_startidx, i_endidx

      ile = ptr_patch%cells%edge_idx(jc,jb,je)
      ibe = ptr_patch%cells%edge_blk(jc,jb,je)

      ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
      ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
      ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
      ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

      IF (jc == ilc1 .AND. jb == ibc1) THEN
        ptr_int%geofac_grg(jc,1,jb,1) = ptr_int%geofac_grg(jc,1,jb,1) +      &
          ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
          ptr_int%c_lin_e(ile,1,ibe)
        ptr_int%geofac_grg(jc,1,jb,2) = ptr_int%geofac_grg(jc,1,jb,2) +      &
          ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
          ptr_int%c_lin_e(ile,1,ibe)
      ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
        ptr_int%geofac_grg(jc,1,jb,1) = ptr_int%geofac_grg(jc,1,jb,1) +      &
          ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
          ptr_int%c_lin_e(ile,2,ibe)
        ptr_int%geofac_grg(jc,1,jb,2) = ptr_int%geofac_grg(jc,1,jb,2) +      &
          ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
          ptr_int%c_lin_e(ile,2,ibe)
      ENDIF
      DO ic = 1, i_cell_type
        ilnc = ptr_patch%cells%neighbor_idx(jc,jb,ic)
        ibnc = ptr_patch%cells%neighbor_blk(jc,jb,ic)
        IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
          ptr_int%geofac_grg(jc,ic+1,jb,1) = ptr_int%geofac_grg(jc,ic+1,jb,1)+ &
            ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
            ptr_int%c_lin_e(ile,1,ibe)
          ptr_int%geofac_grg(jc,ic+1,jb,2) = ptr_int%geofac_grg(jc,ic+1,jb,2)+ &
            ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
            ptr_int%c_lin_e(ile,1,ibe)
        ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
          ptr_int%geofac_grg(jc,ic+1,jb,1) = ptr_int%geofac_grg(jc,ic+1,jb,1)+ &
            ptr_int%primal_normal_ec(jc,jb,je,1)*ptr_int%geofac_div(jc,je,jb)* &
            ptr_int%c_lin_e(ile,2,ibe)
          ptr_int%geofac_grg(jc,ic+1,jb,2) = ptr_int%geofac_grg(jc,ic+1,jb,2)+ &
            ptr_int%primal_normal_ec(jc,jb,je,2)*ptr_int%geofac_div(jc,je,jb)* &
            ptr_int%c_lin_e(ile,2,ibe)
        ENDIF
      ENDDO

      ! To ensure that dummy edges have a factor of 0:
      IF (je > ptr_patch%cells%num_edges(jc,jb)) THEN
        ptr_int%geofac_grg(jc,je+1,jb,1:2) = 0._wp
      ENDIF

    ENDDO !cell loop
  ENDDO

END DO !block loop
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE init_geo_factors
!------------------------------------------------------

!-------------------------------------------------------------------------
!> Allocation of components of interpolation state.
SUBROUTINE allocate_int_state( ptr_patch, ptr_int)
!
  TYPE(t_patch), INTENT(IN) :: ptr_patch

  TYPE(t_int_state), INTENT(inout) :: ptr_int

  INTEGER :: nblks_c, nblks_e, nblks_v, nincr
  INTEGER :: ist, tot_status
  INTEGER :: idummy

!-----------------------------------------------------------------------

  !
  ! determine size of arrays, i.e.
  ! values for the blocking
  !
  nblks_c  = ptr_patch%nblks_c
  nblks_e  = ptr_patch%nblks_e
  nblks_v  = ptr_patch%nblks_v
  !
  ! allocate interpolation state
  !
  ! c_lin_e
  !
  ALLOCATE (ptr_int%c_lin_e(nproma,2,nblks_e), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for c_lin_e failed')
  ENDIF
  
    ALLOCATE (ptr_int%e_bln_c_s(nproma,3,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_bln_c_s failed')
    ENDIF
    !
    ! e_bln_c_u
    !
    ALLOCATE (ptr_int%e_bln_c_u(nproma,6,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_bln_c_u failed')
    ENDIF
    !
    ! e_bln_c_v
    !
    ALLOCATE (ptr_int%e_bln_c_v(nproma,6,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_bln_c_v failed')
    ENDIF
    !
    ! c_bln_avg
    !
    ALLOCATE (ptr_int%c_bln_avg(nproma,4,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for c_bln_avg failed')
    ENDIF
    !
    ! e_flx_avg
    !
    ALLOCATE (ptr_int%e_flx_avg(nproma,5,nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_flx_avg failed')
    ENDIF
  !
  ! v_1o2_e
  !
  ALLOCATE (ptr_int%v_1o2_e(nproma,2,nblks_e), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state', &
    &            'allocation for v_1o2_e failed')
  ENDIF
                                            
  ALLOCATE (ptr_int%geofac_div(nproma, i_cell_type, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_div failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_rot(nproma, 9-i_cell_type, nblks_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state', &
    &             'allocation for geofac_rot failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_n2s(nproma, i_cell_type+1, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_n2s failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_grg(nproma, i_cell_type+1, nblks_c, 2), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_grg failed')
  ENDIF

  ALLOCATE (ptr_int%cart_edge_coord(nproma, nblks_e, 3), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for cart_edge_coord failed')
  ENDIF

  ALLOCATE (ptr_int%cart_cell_coord(nproma, nblks_c, 3), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for cart_cell_coord failed')
  ENDIF

  ALLOCATE (ptr_int%primal_normal_ec(nproma, nblks_c, i_cell_type, 2), STAT=ist)
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for primal_normal_ec failed')
  ENDIF

  ALLOCATE (ptr_int%edge_cell_length(nproma, nblks_c, i_cell_type), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for edge_cell_length failed')
  ENDIF

  ALLOCATE (ptr_int%e_inn_c(nproma,i_cell_type,nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state', &
    &            'allocation for e_inn_c failed')
  ENDIF
  !
  !
  ! verts_aw_cells
  !
  ALLOCATE (ptr_int%verts_aw_cells(nproma,i_cell_type,nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                     &
      &             'allocation for verts_aw_cells failed')
  ENDIF
  !
  ! cells_aw_verts
  !
  ALLOCATE (ptr_int%cells_aw_verts(nproma,9-i_cell_type,nblks_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                     &
      &             'allocation for cells_aw_verts failed')
  ENDIF
    !
    ! rbf_vec_idx_c, rbf_vec_blk_c
    !
    ALLOCATE (ptr_int%rbf_vec_idx_c(rbf_vec_dim_c, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_idx_c failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_vec_blk_c(rbf_vec_dim_c, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_blk_c failed')
    ENDIF
    !
    ! rbf_vec_stencil_c
    !
    ALLOCATE (ptr_int%rbf_vec_stencil_c(nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_stencil_c failed')
    ENDIF
    !
    ! rbf_vec_coeff_c
    !
    ALLOCATE (ptr_int%rbf_vec_coeff_c(rbf_vec_dim_c, 2, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_coeff_c failed')
    ENDIF
    !
    ! rbf_c2grad_idx, rbf_c2grad_blk
    !
    ALLOCATE (ptr_int%rbf_c2grad_idx(rbf_c2grad_dim, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_c2grad_idx failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_c2grad_blk(rbf_c2grad_dim, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_c2grad_blk failed')
    ENDIF
    !
    ! rbf_c2grad_coeff
    !
    ALLOCATE (ptr_int%rbf_c2grad_coeff(rbf_c2grad_dim, 2, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_c2grad_coeff failed')
    ENDIF
    !
    ! rbf_vec_idx_v, rbf_vec_blk_v
    !
    ALLOCATE (ptr_int%rbf_vec_idx_v(rbf_vec_dim_v, nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_idx_v failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_vec_blk_v(rbf_vec_dim_v, nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for rbf_vec_blk_v failed')
    ENDIF
    !
    ! rbf_vec_stencil_v
    !
    ALLOCATE (ptr_int%rbf_vec_stencil_v(nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_stencil_v failed')
    ENDIF
    !
    ! rbf_vec_coeff_v
    !
    ALLOCATE (ptr_int%rbf_vec_coeff_v(rbf_vec_dim_v, 2, nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_coeff_v failed')
    ENDIF
    !
    ! rbf_vec_idx_e, rbf_vec_blk_e
    !
    ALLOCATE (ptr_int%rbf_vec_idx_e(rbf_vec_dim_e, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_idx_e failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_vec_blk_e(rbf_vec_dim_e, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_blk_e failed')
    ENDIF
    !
    ! rbf_vec_stencil_e
    !
    ALLOCATE (ptr_int%rbf_vec_stencil_e(nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_stencil_e failed')
    ENDIF
    !
    ! rbf_vec_coeff_e
    !
!     write(0,*) "ALLOCATE (ptr_int%rbf_vec_coeff_e:",rbf_vec_dim_e, nproma, nblks_e
!     STOP
    ALLOCATE (ptr_int%rbf_vec_coeff_e(rbf_vec_dim_e, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_coeff_e failed')
    ENDIF

     !
    ALLOCATE (ptr_int%lsq_lin%lsq_dim_stencil(0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_dim_stencil failed')
    ENDIF
    !
    ! lsq_idx_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_idx_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_idx_c failed')
    ENDIF
    !
    ! lsq_blk_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_blk_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_blk_c failed')
    ENDIF
    !
    ! lsq_weights_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_weights_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_weights_c failed')
    ENDIF
    !
    ! lsq_qtmat_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_qtmat_c(0, 0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_qtmat_c failed')
    ENDIF
    !
    ! lsq_rmat_rdiag_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_rmat_rdiag_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_rdiag_c failed')
    ENDIF
    !
    ! lsq_rmat_utri_c
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_rmat_utri_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_utri_c failed')
    ENDIF
    !
    ! lsq_moments
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_moments(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for lsq_moments failed')
    ENDIF
    !
    ! lsq_moments_hat
    !
    ALLOCATE (ptr_int%lsq_lin%lsq_moments_hat(0, 0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                       &
        &             'allocation for lsq_moments_hat failed')
    ENDIF

    ! *** higher order ***
    !
    !
    ! lsq_dim_stencil
    !
    ALLOCATE (ptr_int%lsq_high%lsq_dim_stencil(0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_dim_stencil failed')
    ENDIF
    !
    ! lsq_idx_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_idx_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_idx_c failed')
    ENDIF
    !
    ! lsq_blk_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_blk_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_blk_c failed')
    ENDIF
    !
    ! lsq_weights_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_weights_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_weights_c failed')
    ENDIF
    !
    ! lsq_qtmat_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_qtmat_c(0, 0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_qtmat_c failed')
    ENDIF
    !
    ! lsq_rmat_rdiag_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_rmat_rdiag_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_rdiag_c failed')
    ENDIF
    !
    ! lsq_rmat_utri_c
    !
    ALLOCATE (ptr_int%lsq_high%lsq_rmat_utri_c(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for lsq_rmat_utri_c failed')
    ENDIF
    !
    ! lsq_moments
    !
    ALLOCATE (ptr_int%lsq_high%lsq_moments(0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for lsq_moments failed')
    ENDIF
    !
    ! lsq_moments_hat
    !
    ALLOCATE (ptr_int%lsq_high%lsq_moments_hat(0, 0, 0, 0), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for lsq_moments_hat failed')
    ENDIF


    ALLOCATE (ptr_int%geofac_qdiv(nproma, 4, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for geofac_qdiv failed')
    ENDIF
    ALLOCATE (ptr_int%nudgecoeff_c(nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for nudgecoeff_c failed')
    ENDIF
    ALLOCATE (ptr_int%nudgecoeff_e(nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for nudgecoeff_e failed')
    ENDIF

    !
    ! Quadrature points and weights for integration over triangular element
    !
    ALLOCATE (ptr_int%gquad%qpts_tri_l(nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for qpts_tri_l failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%qpts_tri_q(nproma, nblks_c,3), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for qpts_tri_q failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%qpts_tri_c(nproma, nblks_c,4), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for qpts_tri_c failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%weights_tri_q(3), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for weights_tri_q failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%weights_tri_c(4), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for weights_tri_c failed')
    ENDIF

    
    ptr_int%e_bln_c_s     = 0._wp
    ptr_int%e_bln_c_u     = 0._wp
    ptr_int%e_bln_c_v     = 0._wp
    ptr_int%c_bln_avg     = 0._wp
    ptr_int%e_flx_avg     = 0._wp
  
  ptr_int%v_1o2_e       = 0.5_wp
  ptr_int%c_lin_e       = 0._wp
  ptr_int%e_inn_c       = 0._wp
  ptr_int%verts_aw_cells= 0._wp
  ptr_int%cells_aw_verts= 0._wp
    
    ptr_int%rbf_vec_idx_c     = 0
    ptr_int%rbf_vec_blk_c     = 0
    ptr_int%rbf_vec_stencil_c = 0
    ptr_int%rbf_vec_coeff_c   = 0._wp

    ptr_int%rbf_c2grad_idx    = 0
    ptr_int%rbf_c2grad_blk    = 0
    ptr_int%rbf_c2grad_coeff  = 0._wp

    ptr_int%rbf_vec_idx_v     = 0
    ptr_int%rbf_vec_blk_v     = 0
    ptr_int%rbf_vec_stencil_v = 0
    ptr_int%rbf_vec_coeff_v   = 0._wp

    ptr_int%rbf_vec_idx_e     = 0
    ptr_int%rbf_vec_blk_e     = 0
    ptr_int%rbf_vec_stencil_e = 0
    ptr_int%rbf_vec_coeff_e   = 0._wp
    ptr_int%geofac_qdiv = 0._wp
    ptr_int%nudgecoeff_c = 0._wp
    ptr_int%nudgecoeff_e = 0._wp

    ptr_int%gquad%qpts_tri_l(:,:)%lat    = 0._wp
    ptr_int%gquad%qpts_tri_l(:,:)%lon    = 0._wp
    ptr_int%gquad%qpts_tri_q(:,:,:)%lat  = 0._wp
    ptr_int%gquad%qpts_tri_q(:,:,:)%lon  = 0._wp
    ptr_int%gquad%qpts_tri_c(:,:,:)%lat  = 0._wp
    ptr_int%gquad%qpts_tri_c(:,:,:)%lon  = 0._wp
    ptr_int%gquad%weights_tri_q(:)       = 0._wp
    ptr_int%gquad%weights_tri_c(:)       = 0._wp
  
  ptr_int%geofac_div = 0._wp
  ptr_int%geofac_rot = 0._wp
  ptr_int%geofac_n2s = 0._wp
  ptr_int%geofac_grg = 0._wp
  ptr_int%cart_edge_coord = 0._wp
  ptr_int%cart_cell_coord = 0._wp
  ptr_int%primal_normal_ec = 0._wp
  ptr_int%edge_cell_length = 0._wp
 
  IF (msg_level > 10) &
    & CALL message ('mo_intp_state:allocate_int_state','memory allocation finished')

END SUBROUTINE allocate_int_state


!-------------------------------------------------------------------------
!> Allocation of components of interpolation state.
SUBROUTINE destruct_int_state( ptr_int)
!

  TYPE(t_int_state), INTENT(inout) :: ptr_int

  INTEGER :: ist, tot_status
  INTEGER :: idummy

  !
  DEALLOCATE (ptr_int%c_lin_e)
  
  DEALLOCATE (ptr_int%e_bln_c_s)
    !
  DEALLOCATE (ptr_int%e_bln_c_u)
    !
  DEALLOCATE (ptr_int%e_bln_c_v)
    !
  DEALLOCATE (ptr_int%c_bln_avg)
    !
  DEALLOCATE (ptr_int%e_flx_avg)
  !
  DEALLOCATE (ptr_int%v_1o2_e)
                                            
  DEALLOCATE (ptr_int%geofac_div)

  DEALLOCATE (ptr_int%geofac_rot)

  DEALLOCATE (ptr_int%geofac_n2s)

  DEALLOCATE (ptr_int%geofac_grg)

  DEALLOCATE (ptr_int%cart_edge_coord)

  DEALLOCATE (ptr_int%cart_cell_coord)

  DEALLOCATE (ptr_int%primal_normal_ec)

  DEALLOCATE (ptr_int%edge_cell_length)

  DEALLOCATE (ptr_int%e_inn_c)
  !
  DEALLOCATE (ptr_int%verts_aw_cells)
  !
  DEALLOCATE (ptr_int%cells_aw_verts)
    !
    DEALLOCATE (ptr_int%rbf_vec_idx_c )
    DEALLOCATE (ptr_int%rbf_vec_blk_c )
    !
    DEALLOCATE (ptr_int%rbf_vec_stencil_c)
    !
    DEALLOCATE (ptr_int%rbf_vec_coeff_c )
    !
    DEALLOCATE (ptr_int%rbf_c2grad_idx)
    DEALLOCATE (ptr_int%rbf_c2grad_blk )
    !
    DEALLOCATE (ptr_int%rbf_c2grad_coeff )
    !
    DEALLOCATE (ptr_int%rbf_vec_idx_v )
    DEALLOCATE (ptr_int%rbf_vec_blk_v )
    !
    DEALLOCATE (ptr_int%rbf_vec_stencil_v)
    !
    DEALLOCATE (ptr_int%rbf_vec_coeff_v )
    !
    DEALLOCATE (ptr_int%rbf_vec_idx_e)
    DEALLOCATE (ptr_int%rbf_vec_blk_e)
    !
    DEALLOCATE (ptr_int%rbf_vec_stencil_e)
    DEALLOCATE (ptr_int%rbf_vec_coeff_e)
    DEALLOCATE (ptr_int%lsq_lin%lsq_dim_stencil)
   !
    DEALLOCATE (ptr_int%lsq_lin%lsq_idx_c)
    DEALLOCATE (ptr_int%lsq_lin%lsq_blk_c)
    DEALLOCATE (ptr_int%lsq_lin%lsq_weights_c)
    DEALLOCATE (ptr_int%lsq_lin%lsq_qtmat_c)
    DEALLOCATE (ptr_int%lsq_lin%lsq_rmat_rdiag_c)
    DEALLOCATE (ptr_int%lsq_lin%lsq_rmat_utri_c)
    DEALLOCATE (ptr_int%lsq_lin%lsq_moments)
    DEALLOCATE (ptr_int%lsq_lin%lsq_moments_hat)
    DEALLOCATE (ptr_int%lsq_high%lsq_dim_stencil )
    DEALLOCATE (ptr_int%lsq_high%lsq_idx_c)
    DEALLOCATE (ptr_int%lsq_high%lsq_blk_c)
    DEALLOCATE (ptr_int%lsq_high%lsq_weights_c)
    DEALLOCATE (ptr_int%lsq_high%lsq_qtmat_c)
    DEALLOCATE (ptr_int%lsq_high%lsq_rmat_rdiag_c)
    DEALLOCATE (ptr_int%lsq_high%lsq_rmat_utri_c)
    DEALLOCATE (ptr_int%lsq_high%lsq_moments)
    DEALLOCATE (ptr_int%lsq_high%lsq_moments_hat )
    DEALLOCATE (ptr_int%geofac_qdiv)
    DEALLOCATE (ptr_int%nudgecoeff_c )
    DEALLOCATE (ptr_int%nudgecoeff_e )
    DEALLOCATE (ptr_int%gquad%qpts_tri_l )
    DEALLOCATE (ptr_int%gquad%qpts_tri_q )
    DEALLOCATE (ptr_int%gquad%qpts_tri_c )
    DEALLOCATE (ptr_int%gquad%weights_tri_q )
    DEALLOCATE (ptr_int%gquad%weights_tri_c )

END SUBROUTINE destruct_int_state

!-------------------------------------------------------------------------
!>
!! This routine initializes the coefficients used.
!!
!! This routine initializes the coefficients used
!! for interpolations needed for scalars. The original routines were aw_int_coeff
!! and cell2edge_lin_int_coeff
!!
!! @par Revision History
!!  Original version by Hui Wan, MPI-M (2007-08-02)
!!  Modified by Almut Gassmann, MPI-M (2009-01-05)
!!  - added interpolation weights for edges to verts
!!  Renamed by Almut Gassmann, MPI-M (2009-01-27)
!!  - joining aw_int_coeff and cell2edge_lin_int_coeff to this routine
!!  - use of different edge-midpoint to cell distances
!!  - implement area weightings
!!
SUBROUTINE scalar_int_coeff( ptr_patch, ptr_int_state )
!

TYPE(t_patch), INTENT(inout) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int_state

INTEGER :: nlen, nblks_c, npromz_c, nblks_e, npromz_e, nblks_v, npromz_v
INTEGER :: jc, je, jb, jv  ! integer over edges, blocks and levels
INTEGER :: ile, ibe, &
           ilc, ibc, ilc1, ilc2, ibc1, ibc2, idx_ve,&
           ilv, ibv, ilv1, ilv2, ibv1, ibv2, idx_ce
INTEGER :: i_startblk                ! start block
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index

REAL(wp) :: z_sum


!--------------------------------------------------------------------

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_int_c
  npromz_c = ptr_patch%npromz_int_c
  nblks_e  = ptr_patch%nblks_int_e
  npromz_e = ptr_patch%npromz_int_e
  nblks_v  = ptr_patch%nblks_int_v
  npromz_v = ptr_patch%npromz_int_v

  ! a) the control volume associated to each edge is defined as the
  ! quadrilateral whose edges are the primal edge and the associated dual edge
  !----------------------------------------------------------------------------
  ! loop over all blocks and edges

!$OMP PARALLEL PRIVATE(i_startblk)
!$OMP DO PRIVATE(jb,je,nlen)
  DO jb = 1, nblks_e

    IF (jb /= nblks_e) THEN
      nlen = nproma
    ELSE
      nlen = npromz_e
    ENDIF

    DO je = 1, nlen
      ptr_patch%edges%area_edge(je,jb) =  &
        &    ptr_patch%edges%primal_edge_length(je,jb)  &
        &  * ptr_patch%edges%dual_edge_length(je,jb)
    END DO

  END DO
!$OMP END DO

  ! b1) cell to edge averages
  !-------------------------
  ! The calculation cannot be done for boundary edges
  i_startblk = ptr_patch%edges%start_blk(2,1)
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ilc2,ibc1,ibc2,ilv1,ilv2,&
!$OMP            ibv1,ibv2)
  DO jb = i_startblk, nblks_e

    CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, 2)

    DO je = i_startidx, i_endidx

      ! inverse distance averaging (former subroutine cell2edge_lin_int_coeff)
      ! For the hexagonal grid, this is also the direct distance averaging
      ! because the edge is exactly half way between the cells
      ! (needed for proper bracket formalism)
      ptr_int_state%c_lin_e(je,1,jb) = ptr_patch%edges%edge_cell_length(je,jb,2)/&
                                           ptr_patch%edges%dual_edge_length(je,jb)
      ptr_int_state%c_lin_e(je,2,jb) = ptr_patch%edges%edge_cell_length(je,jb,1)/&
                                           ptr_patch%edges%dual_edge_length(je,jb)

      IF (i_cell_type == 6) THEN
        ilv1 = ptr_patch%edges%vertex_idx(je,jb,1)
        ilv2 = ptr_patch%edges%vertex_idx(je,jb,2)
        ibv1 = ptr_patch%edges%vertex_blk(je,jb,1)
        ibv2 = ptr_patch%edges%vertex_blk(je,jb,2)
        ptr_int_state%tria_aw_rhom(je,1,jb)=ptr_patch%verts%dual_area(ilv1,ibv1)/&
              (ptr_patch%verts%dual_area(ilv1,ibv1)+ptr_patch%verts%dual_area(ilv2,ibv2))
        ptr_int_state%tria_aw_rhom(je,2,jb)=ptr_patch%verts%dual_area(ilv2,ibv2)/&
              (ptr_patch%verts%dual_area(ilv1,ibv1)+ptr_patch%verts%dual_area(ilv2,ibv2))
      ENDIF

    ENDDO

  ENDDO
!$OMP END DO



  ! b2) vert to edge averages
  !-------------------------
  IF (i_cell_type == 6) THEN
    ! The calculation cannot be done for boundary edges
    i_startblk = ptr_patch%edges%start_blk(2,1)
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx)
    DO jb = i_startblk, nblks_e

      CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, i_startidx, i_endidx, 2)

      DO je = i_startidx, i_endidx

        ! distance averaging
        ptr_int_state%v_1o2_e(je,1,jb) = ptr_patch%edges%edge_vert_length(je,jb,1)/&
                                             ptr_patch%edges%primal_edge_length(je,jb)
        ptr_int_state%v_1o2_e(je,2,jb) = ptr_patch%edges%edge_vert_length(je,jb,2)/&
                                             ptr_patch%edges%primal_edge_length(je,jb)

      ENDDO

    ENDDO
!$OMP END DO
  ENDIF

  ! c) vert to cell averagings, edge to cell inner product
  !-------------------------------------------------------
  ! loop over all blocks and cells

!$OMP DO PRIVATE(jb,jc,je,jv,nlen,ile,ibe,idx_ce,ilv1,ilv2,ibv1,ibv2,ilv,ibv)
  DO jb = 1, nblks_c
    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF

    DO jc = 1, nlen

       ptr_int_state%verts_aw_cells(jc,:,jb) = 0.0_wp

       IF (i_cell_type == 6) z_sum = 0.0_wp

       DO je = 1, ptr_patch%cells%num_edges(jc,jb)

          ile = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe = ptr_patch%cells%edge_blk(jc,jb,je)
          IF ( ptr_patch%edges%cell_idx(ile,ibe,1) == jc .AND. &
               ptr_patch%edges%cell_blk(ile,ibe,1) == jb ) THEN
               idx_ce = 1
          ELSE
               idx_ce = 2
          ENDIF

          ptr_int_state%e_inn_c(jc,je,jb) = &
              ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
              ptr_patch%edges%primal_edge_length(ile,ibe)/&
              ptr_patch%cells%area(jc,jb)

          IF (i_cell_type == 6) THEN
            ptr_int_state%e_aw_c(jc,je,jb) = 0.5_wp*&
              ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
              ptr_patch%edges%primal_edge_length(ile,ibe)/&
              ptr_patch%cells%area(jc,jb)
            ptr_int_state%r_aw_c(jc,je,jb) = &
              ptr_patch%edges%quad_area(ile,ibe)
            z_sum = z_sum+ptr_patch%edges%quad_area(ile,ibe)
          ENDIF

          ilv1 = ptr_patch%edges%vertex_idx(ile,ibe,1)
          ibv1 = ptr_patch%edges%vertex_blk(ile,ibe,1)
          ilv2 = ptr_patch%edges%vertex_idx(ile,ibe,2)
          ibv2 = ptr_patch%edges%vertex_blk(ile,ibe,2)

          DO jv = 1, ptr_patch%cells%num_edges(jc,jb)
            ilv = ptr_patch%cells%vertex_idx(jc,jb,jv)
            ibv = ptr_patch%cells%vertex_blk(jc,jb,jv)

            IF (ilv == ilv1 .AND. ibv == ibv1) THEN
              ptr_int_state%verts_aw_cells(jc,jv,jb) =   &
                ptr_int_state%verts_aw_cells(jc,jv,jb) + &
                0.5_wp/ptr_patch%cells%area(jc,jb) *             &
                ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
                ptr_patch%edges%edge_vert_length(ile,ibe,1)
            ENDIF
            IF (ilv == ilv2 .AND. ibv == ibv2) THEN
              ptr_int_state%verts_aw_cells(jc,jv,jb)  =  &
                ptr_int_state%verts_aw_cells(jc,jv,jb) + &
                0.5_wp/ptr_patch%cells%area(jc,jb) *             &
                ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
                ptr_patch%edges%edge_vert_length(ile,ibe,2)
            ENDIF

          ENDDO

       ENDDO
       IF (i_cell_type == 6) THEN
         ptr_int_state%r_aw_c(jc,:,jb)=ptr_int_state%r_aw_c(jc,:,jb)/z_sum
       ENDIF

    ENDDO !loop over all cells

  ENDDO   !loop over all blocks
!$OMP END DO

  ! d) cells to verts averagings, edge to verts averagings
  !-------------------------------------------------------
  ! loop over all blocks and verts

  i_startblk = ptr_patch%verts%start_blk(2,1)
!$OMP DO PRIVATE(jb,jc,je,jv,i_startidx,i_endidx,ile,ibe,idx_ve,ilc,ibc, &
!$OMP            ilc1,ilc2,ibc1,ibc2 )
  DO jb = i_startblk, nblks_v

    CALL get_indices_v(ptr_patch, jb, i_startblk, nblks_v, &
                       i_startidx, i_endidx, 2)

    DO jv = i_startidx, i_endidx

       ptr_int_state%cells_aw_verts(jv,:,jb) = 0.0_wp


       DO je = 1, ptr_patch%verts%num_edges(jv,jb)

          ile = ptr_patch%verts%edge_idx(jv,jb,je)
          ibe = ptr_patch%verts%edge_blk(jv,jb,je)
          IF ( ptr_patch%edges%vertex_idx(ile,ibe,1) == jv .AND. &
               ptr_patch%edges%vertex_blk(ile,ibe,1) == jb ) THEN
               idx_ve = 1
          ELSE
               idx_ve = 2
          ENDIF

          IF (i_cell_type == 6 ) THEN
            ptr_int_state%e_aw_v(jv,je,jb) = 0.5_wp*&
            & ptr_patch%edges%edge_vert_length(ile,ibe,idx_ve) &
            &*ptr_patch%edges%dual_edge_length(ile,ibe) &
            &/ptr_patch%verts%dual_area(jv,jb)
            ptr_int_state%e_inn_v(jv,je,jb) = &
            & ptr_patch%edges%edge_vert_length(ile,ibe,idx_ve) &
            &*ptr_patch%edges%dual_edge_length(ile,ibe) &
            &/ptr_patch%verts%dual_area(jv,jb)
          ENDIF

          ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
          ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
          ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
          ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

          DO jc = 1, ptr_patch%verts%num_edges(jv,jb)
            ilc = ptr_patch%verts%cell_idx(jv,jb,jc)
            ibc = ptr_patch%verts%cell_blk(jv,jb,jc)

            IF (ilc == ilc1 .AND. ibc == ibc1) THEN
              ptr_int_state%cells_aw_verts(jv,jc,jb) =   &
                ptr_int_state%cells_aw_verts(jv,jc,jb) + &
                0.5_wp/ptr_patch%verts%dual_area(jv,jb) *             &
                ptr_patch%edges%edge_vert_length(ile,ibe,idx_ve)*&
                ptr_patch%edges%edge_cell_length(ile,ibe,1)
            ENDIF
            IF (ilc == ilc2 .AND. ibc == ibc2) THEN
              ptr_int_state%cells_aw_verts(jv,jc,jb)  =  &
                ptr_int_state%cells_aw_verts(jv,jc,jb) + &
                0.5_wp/ptr_patch%verts%dual_area(jv,jb) *             &
                ptr_patch%edges%edge_vert_length(ile,ibe,idx_ve)*&
                ptr_patch%edges%edge_cell_length(ile,ibe,2)
            ENDIF

          ENDDO

       ENDDO

    ENDDO !loop over all cells

  ENDDO   !loop over all blocks
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE scalar_int_coeff

!-------------------------------------------------------------------------
!
!
!>
!! Computes the weighting coefficients for bilinear edge-to-cell interpolation.
!!
!! Results are stored in ptr_int_state\\%e_bln_c_s
!!
!! @par Revision History
!!  developed by Guenther Zaengl, 2009-01-06
!!
SUBROUTINE bln_int_coeff_e2c( ptr_patch, ptr_int_state )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(inout) :: ptr_int_state
!

INTEGER :: jc, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3

REAL(wp) :: xtemp,ytemp,wgt(3),xloc,yloc,x(3),y(3), &
            pollat,pollon

!-----------------------------------------------------------------------

!$OMP PARALLEL PRIVATE(rl_start, rl_end,i_nchdom,i_startblk,i_endblk)
rl_start = 1
rl_end = min_rlcell

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,yloc,xloc,pollat,pollon,ile1,ibe1,&
!$OMP            ile2,ibe2,ile3,ibe3,xtemp,ytemp,wgt,x,y)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc = i_startidx, i_endidx

    yloc = ptr_patch%cells%center(jc,jb)%lat
    xloc = ptr_patch%cells%center(jc,jb)%lon

    ! Rotate local point into the equator for better accuracy of bilinear weights
    IF (yloc >= 0._wp) THEN
      pollat = yloc - pi2/4._wp
    ELSE
      pollat = yloc + pi2/4._wp
    ENDIF
    pollon = xloc

    CALL rotate_latlon( yloc, xloc, pollat, pollon )

    !  get the line and block indices of the cell edges

    ile1 = ptr_patch%cells%edge_idx(jc,jb,1)
    ibe1 = ptr_patch%cells%edge_blk(jc,jb,1)
    ile2 = ptr_patch%cells%edge_idx(jc,jb,2)
    ibe2 = ptr_patch%cells%edge_blk(jc,jb,2)
    ile3 = ptr_patch%cells%edge_idx(jc,jb,3)
    ibe3 = ptr_patch%cells%edge_blk(jc,jb,3)

    ! x and y are the zonal and meridional distances from the local
    ! cell point (ignoring the earth's radius, which drops out anyway)

    xtemp = ptr_patch%edges%center(ile1,ibe1)%lon
    ytemp = ptr_patch%edges%center(ile1,ibe1)%lat
    CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

    y(1)  = ytemp-yloc
    x(1)  = xtemp-xloc
    ! This is needed when the date line is crossed
    IF (x(1) >  3.5_wp) x(1) = x(1) - pi2
    IF (x(1) < -3.5_wp) x(1) = x(1) + pi2

    xtemp = ptr_patch%edges%center(ile2,ibe2)%lon
    ytemp = ptr_patch%edges%center(ile2,ibe2)%lat
    CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

    y(2)  = ytemp-yloc
    x(2)  = xtemp-xloc
    ! This is needed when the date line is crossed
    IF (x(2) >  3.5_wp) x(2) = x(2) - pi2
    IF (x(2) < -3.5_wp) x(2) = x(2) + pi2

    xtemp = ptr_patch%edges%center(ile3,ibe3)%lon
    ytemp = ptr_patch%edges%center(ile3,ibe3)%lat
    CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

    y(3)  = ytemp-yloc
    x(3)  = xtemp-xloc
    ! This is needed when the date line is crossed
    IF (x(3) >  3.5_wp) x(3) = x(3) - pi2
    IF (x(3) < -3.5_wp) x(3) = x(3) + pi2

    ! The weighting factors are based on the requirement that sum(w(i)*x(i)) = 0
    ! and sum(w(i)*y(i)) = 0, which ensures that linear horizontal gradients
    ! are not aliased into a checkerboard pattern between upward- and downward
    ! directed cells. The third condition is sum(w(i)) = 1. Analytical elimination yields...

    wgt(3) = 1.0_wp/( (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))/(x(2)-x(1)) ) * &
                ( -y(1) + x(1)*(y(2)-y(1))/(x(2)-x(1)) )
    wgt(2) = (-x(1) - wgt(3)*(x(3)-x(1)))/(x(2)-x(1))
    wgt(1) = 1.0_wp - wgt(2) - wgt(3)

    ! Store results in ptr_int_state%e_bln_c_s
    ptr_int_state%e_bln_c_s(jc,1,jb) = wgt(1)
    ptr_int_state%e_bln_c_s(jc,2,jb) = wgt(2)
    ptr_int_state%e_bln_c_s(jc,3,jb) = wgt(3)

  ENDDO !cell loop

END DO !block loop
!$OMP END DO

! Now compute vector interpolation weights: These take the normal and tangential
! wind components at the edges and reconstruct u and v at the cell midpoints

rl_start = 2

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ile1,ibe1,ile2,ibe2,ile3,ibe3)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc = i_startidx, i_endidx

    !  get the line and block indices of the cell edges

    ile1 = ptr_patch%cells%edge_idx(jc,jb,1)
    ibe1 = ptr_patch%cells%edge_blk(jc,jb,1)
    ile2 = ptr_patch%cells%edge_idx(jc,jb,2)
    ibe2 = ptr_patch%cells%edge_blk(jc,jb,2)
    ile3 = ptr_patch%cells%edge_idx(jc,jb,3)
    ibe3 = ptr_patch%cells%edge_blk(jc,jb,3)

    IF (ptr_patch%edges%cell_idx(ile1,ibe1,1) == jc .AND. &
        ptr_patch%edges%cell_blk(ile1,ibe1,1) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v1
      ptr_int_state%e_bln_c_u(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,1)%v1

      ptr_int_state%e_bln_c_v(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v2
      ptr_int_state%e_bln_c_v(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,1)%v2

    ELSE IF (ptr_patch%edges%cell_idx(ile1,ibe1,2) == jc .AND. &
        ptr_patch%edges%cell_blk(ile1,ibe1,2) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v1
      ptr_int_state%e_bln_c_u(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,2)%v1

      ptr_int_state%e_bln_c_v(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v2
      ptr_int_state%e_bln_c_v(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,2)%v2

    ENDIF

    IF (ptr_patch%edges%cell_idx(ile2,ibe2,1) == jc .AND. &
        ptr_patch%edges%cell_blk(ile2,ibe2,1) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,1)%v1
      ptr_int_state%e_bln_c_u(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,1)%v1

      ptr_int_state%e_bln_c_v(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,1)%v2
      ptr_int_state%e_bln_c_v(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,1)%v2

    ELSE IF (ptr_patch%edges%cell_idx(ile2,ibe2,2) == jc .AND. &
        ptr_patch%edges%cell_blk(ile2,ibe2,2) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,2)%v1
      ptr_int_state%e_bln_c_u(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,2)%v1

      ptr_int_state%e_bln_c_v(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,2)%v2
      ptr_int_state%e_bln_c_v(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,2)%v2

    ENDIF

    IF (ptr_patch%edges%cell_idx(ile3,ibe3,1) == jc .AND. &
        ptr_patch%edges%cell_blk(ile3,ibe3,1) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,1)%v1
      ptr_int_state%e_bln_c_u(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,1)%v1

      ptr_int_state%e_bln_c_v(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,1)%v2
      ptr_int_state%e_bln_c_v(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,1)%v2

    ELSE IF (ptr_patch%edges%cell_idx(ile3,ibe3,2) == jc .AND. &
        ptr_patch%edges%cell_blk(ile3,ibe3,2) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,2)%v1
      ptr_int_state%e_bln_c_u(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,2)%v1

      ptr_int_state%e_bln_c_v(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,2)%v2
      ptr_int_state%e_bln_c_v(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,2)%v2

    ENDIF
  ENDDO !cell loop

END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE bln_int_coeff_e2c


!-------------------------------------------------------------------------
!>
!! Computes the local orientation of the edge primal normal and dual normal.
!!
!! Computes the local orientation of the edge primal normal and dual normal
!! at the location of the cell centers and vertices.
!! Moreover, the Cartesian orientation vectors of the edge primal normals
!! are stored for use in the RBF initialization routines, and inverse
!! primal and dual edge lengths are computed
!!
!! @par Revision History
!!  developed by Guenther Zaengl, 2009-03-31
!!
SUBROUTINE complete_patchinfo( ptr_patch, ptr_int )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation state
TYPE(t_int_state),     INTENT(inout) :: ptr_int
!

INTEGER :: jb, je, jc
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: ilc1, ibc1, ilv1, ibv1, ilc2, ibc2, ilv2, ibv2, &
           ilv3, ibv3, ilv4, ibv4, ile1, ibe1

REAL(wp) :: z_nu, z_nv, z_lon, z_lat, z_nx1(3), z_nx2(3), z_norm

TYPE(t_cartesian_coordinates) :: cc_edge, cc_cell, cc_ev3, cc_ev4

!-----------------------------------------------------------------------

i_nchdom   = MAX(1,ptr_patch%n_childdom)

!$OMP PARALLEL  PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
rl_start = 1
rl_end = min_rlcell

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! First step: compute Cartesian coordinates of cell centers on full domain
! this is needed for least squares gradient reconstruction;
!
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,cc_cell)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc =  i_startidx, i_endidx

    ! compute Cartesian coordinates (needed for RBF initialization)

    cc_cell = gc2cc(ptr_patch%cells%center(jc,jb))
    ptr_int%cart_cell_coord(jc,jb,1:3)= cc_cell%x(1:3)

  ENDDO

END DO !block loop
!$OMP END DO


rl_start = 1
rl_end = min_rledge

! values for the blocking
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)
!
! First step: compute Cartesian coordinates and Cartesian vectors on full domain
! this is needed to vectorize RBF initialization; the existing field carrying
! the Cartesian orientation vectors (primal_cart_normal) did not work for that
! because it is a derived data type
! In addition, the fields for the inverse primal and dual edge lengths are
! initialized here.
!
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,cc_edge)
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je =  i_startidx, i_endidx

    ! compute Cartesian coordinates (needed for RBF initialization)

    cc_edge = gc2cc(ptr_patch%edges%center(je,jb))
    ptr_int%cart_edge_coord(je,jb,1:3)= cc_edge%x(1:3)

    ! finally, compute inverse primal edge length
    ! (dual follows below in the rl_start=2 section)

    ptr_patch%edges%inv_primal_edge_length(je,jb) = &
      1._wp/ptr_patch%edges%primal_edge_length(je,jb)

  ENDDO

END DO !block loop
!$OMP END DO

rl_start = 2
rl_end = min_rledge

! Second step: computed projected orientation vectors and related information
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)
!
! loop through all patch edges
!
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,ilc1,ibc1,ilv1,ibv1,ilc2,ibc2,ilv2, &
!$OMP            ibv2,ilv3,ibv3,ilv4,ibv4,z_nu,z_nv,z_lon,z_lat,z_nx1,z_nx2,   &
!$OMP            cc_ev3,cc_ev4,z_norm)
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO je =  i_startidx, i_endidx

    ! compute inverse dual edge length (undefined for refin_ctrl=1)

    ptr_patch%edges%inv_dual_edge_length(je,jb) = &
      1._wp/ptr_patch%edges%dual_edge_length(je,jb)

    ! compute edge-vertex indices (and blocks) 3 and 4, which
    ! are the outer vertices of cells 1 and 2, respectively,
    ! and the inverse length bewtween vertices 3 and 4

    ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
    ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
    ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
    ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

    ilv1 = ptr_patch%edges%vertex_idx(je,jb,1)
    ibv1 = ptr_patch%edges%vertex_blk(je,jb,1)
    ilv2 = ptr_patch%edges%vertex_idx(je,jb,2)
    ibv2 = ptr_patch%edges%vertex_blk(je,jb,2)

    IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,1) /= &
         ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
         ptr_patch%cells%vertex_blk(ilc1,ibc1,1) /= &
         ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
        (ptr_patch%cells%vertex_idx(ilc1,ibc1,1) /= &
         ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
         ptr_patch%cells%vertex_blk(ilc1,ibc1,1) /= &
         ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,1)
      ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,1)

    ELSE IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,2) /= &
              ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
              ptr_patch%cells%vertex_blk(ilc1,ibc1,2) /= &
              ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
             (ptr_patch%cells%vertex_idx(ilc1,ibc1,2) /= &
              ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
              ptr_patch%cells%vertex_blk(ilc1,ibc1,2) /= &
              ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,2)
      ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,2)

    ELSE IF ((ptr_patch%cells%vertex_idx(ilc1,ibc1,3) /= &
              ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
              ptr_patch%cells%vertex_blk(ilc1,ibc1,3) /= &
              ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
             (ptr_patch%cells%vertex_idx(ilc1,ibc1,3) /= &
              ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
              ptr_patch%cells%vertex_blk(ilc1,ibc1,3) /= &
              ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,3) = ptr_patch%cells%vertex_idx(ilc1,ibc1,3)
      ptr_patch%edges%vertex_blk(je,jb,3) = ptr_patch%cells%vertex_blk(ilc1,ibc1,3)

    ENDIF

    IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,1) /= &
         ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
         ptr_patch%cells%vertex_blk(ilc2,ibc2,1) /= &
         ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
        (ptr_patch%cells%vertex_idx(ilc2,ibc2,1) /= &
         ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
         ptr_patch%cells%vertex_blk(ilc2,ibc2,1) /= &
         ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,1)
      ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,1)

    ELSE IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,2) /= &
              ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
              ptr_patch%cells%vertex_blk(ilc2,ibc2,2) /= &
              ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
             (ptr_patch%cells%vertex_idx(ilc2,ibc2,2) /= &
              ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
              ptr_patch%cells%vertex_blk(ilc2,ibc2,2) /= &
              ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,2)
      ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,2)

    ELSE IF ((ptr_patch%cells%vertex_idx(ilc2,ibc2,3) /= &
              ptr_patch%edges%vertex_idx(je,jb,1) .OR.  &
              ptr_patch%cells%vertex_blk(ilc2,ibc2,3) /= &
              ptr_patch%edges%vertex_blk(je,jb,1)) .AND.  &
             (ptr_patch%cells%vertex_idx(ilc2,ibc2,3) /= &
              ptr_patch%edges%vertex_idx(je,jb,2) .OR.  &
              ptr_patch%cells%vertex_blk(ilc2,ibc2,3) /= &
              ptr_patch%edges%vertex_blk(je,jb,2)) )        THEN

      ptr_patch%edges%vertex_idx(je,jb,4) = ptr_patch%cells%vertex_idx(ilc2,ibc2,3)
      ptr_patch%edges%vertex_blk(je,jb,4) = ptr_patch%cells%vertex_blk(ilc2,ibc2,3)

    ENDIF

    ilv3 = ptr_patch%edges%vertex_idx(je,jb,3)
    ibv3 = ptr_patch%edges%vertex_blk(je,jb,3)
    ilv4 = ptr_patch%edges%vertex_idx(je,jb,4)
    ibv4 = ptr_patch%edges%vertex_blk(je,jb,4)

    cc_ev3 = gc2cc(ptr_patch%verts%vertex(ilv3,ibv3))
    cc_ev4 = gc2cc(ptr_patch%verts%vertex(ilv4,ibv4))

    ! inverse length bewtween vertices 3 and 4
    IF (i_cell_type == 3 ) THEN
      ptr_patch%edges%inv_vert_vert_length(je,jb) = 1._wp/(re*arc_length(cc_ev3,cc_ev4))
    ENDIF

    ! next step: compute projected orientation vectors for cells and vertices
    ! bordering to each edge (incl. vertices 3 and 4 intorduced above)

    ! transform orientation vectors at local edge center to Cartesian space
    z_lon = ptr_patch%edges%center(je,jb)%lon
    z_lat = ptr_patch%edges%center(je,jb)%lat

    ! transform primal normal to cartesian vector z_nx1
    z_nx1(:)=ptr_patch%edges%primal_cart_normal(je,jb)%x(:)

    ! transform dual normal to cartesian vector z_nx2
    z_nx2(:)=ptr_patch%edges%dual_cart_normal(je,jb)%x(:)

    ! get location of cell 1

    z_lon = ptr_patch%cells%center(ilc1,ibc1)%lon
    z_lat = ptr_patch%cells%center(ilc1,ibc1)%lat

    ! compute local primal and dual normals at cell 1

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_cell(je,jb,1)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_cell(je,jb,1)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_cell(je,jb,1)%v2 = z_nv/z_norm

    ! get location of cell 2

    z_lon = ptr_patch%cells%center(ilc2,ibc2)%lon
    z_lat = ptr_patch%cells%center(ilc2,ibc2)%lat

    ! compute local primal and dual normals at cell 2

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_cell(je,jb,2)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_cell(je,jb,2)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_cell(je,jb,2)%v2 = z_nv/z_norm

    ! get location of vertex 1

    z_lon = ptr_patch%verts%vertex(ilv1,ibv1)%lon
    z_lat = ptr_patch%verts%vertex(ilv1,ibv1)%lat

    ! compute local primal and dual normals at vertex 1

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_vert(je,jb,1)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_vert(je,jb,1)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_vert(je,jb,1)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_vert(je,jb,1)%v2 = z_nv/z_norm

    ! get location of vertex 2

    z_lon = ptr_patch%verts%vertex(ilv2,ibv2)%lon
    z_lat = ptr_patch%verts%vertex(ilv2,ibv2)%lat

    ! compute local primal and dual normals at vertex 2

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_vert(je,jb,2)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_vert(je,jb,2)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_vert(je,jb,2)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_vert(je,jb,2)%v2 = z_nv/z_norm

    ! get location of vertex 3

    z_lon = ptr_patch%verts%vertex(ilv3,ibv3)%lon
    z_lat = ptr_patch%verts%vertex(ilv3,ibv3)%lat

    ! compute local primal and dual normals at vertex 3

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_vert(je,jb,3)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_vert(je,jb,3)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_vert(je,jb,3)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_vert(je,jb,3)%v2 = z_nv/z_norm

    ! get location of vertex 4

    z_lon = ptr_patch%verts%vertex(ilv4,ibv4)%lon
    z_lat = ptr_patch%verts%vertex(ilv4,ibv4)%lat

    ! compute local primal and dual normals at vertex 2

    CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%primal_normal_vert(je,jb,4)%v1 = z_nu/z_norm
    ptr_patch%edges%primal_normal_vert(je,jb,4)%v2 = z_nv/z_norm

    CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
    z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

    ptr_patch%edges%dual_normal_vert(je,jb,4)%v1 = z_nu/z_norm
    ptr_patch%edges%dual_normal_vert(je,jb,4)%v2 = z_nv/z_norm

  ENDDO

END DO !block loop
!$OMP END DO


rl_start = 2
rl_end = min_rlcell

! Final step: store primal_normal_cell also with respect to cell points
! in order to reduce indirect addressing during runtime
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells
!
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,ile1,ibe1)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc =  i_startidx, i_endidx

    DO je = 1, ptr_patch%cells%num_edges(jc,jb)

      ile1 = ptr_patch%cells%edge_idx(jc,jb,je)
      ibe1 = ptr_patch%cells%edge_blk(jc,jb,je)


      IF ((ptr_patch%edges%cell_idx(ile1,ibe1,1) == jc) .AND. &
          (ptr_patch%edges%cell_blk(ile1,ibe1,1) == jb)) THEN

        ptr_int%primal_normal_ec(jc,jb,je,1) = &
          ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v1
        ptr_int%primal_normal_ec(jc,jb,je,2) = &
          ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v2
        ptr_int%edge_cell_length(jc,jb,je) = &
          ptr_patch%edges%edge_cell_length(ile1,ibe1,1)

      ELSE IF ((ptr_patch%edges%cell_idx(ile1,ibe1,2) == jc) .AND. &
               (ptr_patch%edges%cell_blk(ile1,ibe1,2) == jb)) THEN

        ptr_int%primal_normal_ec(jc,jb,je,1) = &
          ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v1
        ptr_int%primal_normal_ec(jc,jb,je,2) = &
          ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v2
        ptr_int%edge_cell_length(jc,jb,je) = &
          ptr_patch%edges%edge_cell_length(ile1,ibe1,2)

      ENDIF

    ENDDO

  ENDDO

END DO !block loop
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE complete_patchinfo

!-------------------------------------------------------------------------
!>
!! Computes the weighting coefficients for cell averaging with.
!!
!! Computes the weighting coefficients for cell averaging with
!! variable interpolation factors. Results are stored in ptr_patch%cells%avg_wgt
!!
!! @par Revision History
!!  developed by Guenther Zaengl, 2008-12-05
!! @par
!!  modification by Guenther Zaengl, 2009-09-02
!!  revised weights to achieve mass conservation
!!
SUBROUTINE init_cellavg_wgt( ptr_patch, ptr_int )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(inout):: ptr_int
!

INTEGER :: jc, je, jb, iter, niter
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: ilc1, ibc1, ilc2, ibc2, ilc3, ibc3, inb1, inb2, inb3, ie4, ie5
INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3, ile4, ibe4
INTEGER, DIMENSION(nproma) :: iie1, iie2, iie3, iie4
REAL(wp), DIMENSION (nproma,3) :: z_nx1, z_nx2, z_nx3, z_nx4, z_nx5
REAL(wp) :: checksum(nproma,ptr_patch%nblks_e)

REAL(wp) :: xtemp,ytemp,wgt(3),xloc,yloc,x(3),y(3), &
            pollat,pollon,relax_coeff,wgt_loc,maxwgt_loc,minwgt_loc

REAL(wp), DIMENSION(nproma,ptr_patch%nblks_c)  :: wgt_loc_sum, resid
INTEGER, DIMENSION(nproma,ptr_patch%nblks_c,3) :: inv_neighbor_id

#ifdef DEBUG_COEFF
REAL(wp) :: sum1
#endif
!-----------------------------------------------------------------------

! Number of iterations for computation of bilinear weights
niter = 1000

! Relaxation coefficient for adaptation of local weight (empirically determined)
relax_coeff = 0.46_wp

! Initial weighting factor of the local grid point
! wgt_loc = divavg_cntrwgt
wgt_loc = 0.5_wp

! Maximum/minimum  weighting factors of the local grid point
maxwgt_loc = wgt_loc + 0.003_wp
minwgt_loc = wgt_loc - 0.003_wp

! Initialization of the residuum  field
resid(:,:) = 0._wp

! values for the blocking
rl_start = 2
rl_end = min_rlcell

i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! Compute inverse neighbor ID's
! The inverse neigbor ID of a neighbor cell (ilc1,ibc1) is the neighbor ID
! the local cell (jc,jb) has from the point of view of the neighbor cell

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,ilc3,ibc3)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc = i_startidx, i_endidx

    ! line and block indices of the neighbouring cells

    ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
    ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
    ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
    ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
    ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
    ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)

    IF ( (ilc1>0) .AND. (ibc1>0) ) THEN
      IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,1) == jc) .AND. &
          (ptr_patch%cells%neighbor_blk(ilc1,ibc1,1) == jb)) THEN
        inv_neighbor_id(jc,jb,1) = 1
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,2) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc1,ibc1,2) == jb)) THEN
        inv_neighbor_id(jc,jb,1) = 2
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc1,ibc1,3) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc1,ibc1,3) == jb)) THEN
        inv_neighbor_id(jc,jb,1) = 3
      ENDIF
    ENDIF
    IF ( (ilc2>0) .AND. (ibc2>0) ) THEN
      IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,1) == jc) .AND. &
          (ptr_patch%cells%neighbor_blk(ilc2,ibc2,1) == jb)) THEN
        inv_neighbor_id(jc,jb,2)  = 1
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,2) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc2,ibc2,2) == jb)) THEN
        inv_neighbor_id(jc,jb,2)  = 2
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc2,ibc2,3) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc2,ibc2,3) == jb)) THEN
        inv_neighbor_id(jc,jb,2)  = 3
      ENDIF
    ENDIF
    IF ( (ilc3>0) .AND. (ibc3>0) ) THEN
      IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,1) == jc) .AND. &
          (ptr_patch%cells%neighbor_blk(ilc3,ibc3,1) == jb)) THEN
        inv_neighbor_id(jc,jb,3)  = 1
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,2) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc3,ibc3,2) == jb)) THEN
        inv_neighbor_id(jc,jb,3)  = 2
      ELSE IF ((ptr_patch%cells%neighbor_idx(ilc3,ibc3,3) == jc) .AND. &
               (ptr_patch%cells%neighbor_blk(ilc3,ibc3,3) == jb)) THEN
        inv_neighbor_id(jc,jb,3)  = 3
      ENDIF
    ENDIF

  ENDDO !cell loop

END DO !block loop
!$OMP END DO

! Compute coefficients for bilinear interpolation

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,yloc,xloc,pollat,pollon,&
!$OMP            ilc1,ibc1,ilc2,ibc2,ilc3,ibc3,xtemp,ytemp,wgt,x,y)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx

      yloc = ptr_patch%cells%center(jc,jb)%lat
      xloc = ptr_patch%cells%center(jc,jb)%lon

      ! Rotate local point into the equator for better accuracy of bilinear weights
      IF (yloc >= 0._wp) THEN
        pollat = yloc - pi2/4._wp
      ELSE
        pollat = yloc + pi2/4._wp
      ENDIF
      pollon = xloc

      CALL rotate_latlon( yloc, xloc, pollat, pollon )

      ! line and block indices of the neighbouring cells

      ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
      ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
      ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
      ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
      ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
      ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)

      ! x and y are the zonal and meridional distances from the local
      ! cell point (ignoring the earth's radius, which drops out anyway)

      xtemp = ptr_patch%cells%center(ilc1,ibc1)%lon
      ytemp = ptr_patch%cells%center(ilc1,ibc1)%lat
      CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

      y(1)  = ytemp-yloc
      x(1)  = xtemp-xloc
      ! This is needed when the date line is crossed
      IF (x(1) >  3.5_wp) x(1) = x(1) - pi2
      IF (x(1) < -3.5_wp) x(1) = x(1) + pi2

      xtemp = ptr_patch%cells%center(ilc2,ibc2)%lon
      ytemp = ptr_patch%cells%center(ilc2,ibc2)%lat
      CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

      y(2)  = ytemp-yloc
      x(2)  = xtemp-xloc
      ! This is needed when the date line is crossed
      IF (x(2) >  3.5_wp) x(2) = x(2) - pi2
      IF (x(2) < -3.5_wp) x(2) = x(2) + pi2

      xtemp = ptr_patch%cells%center(ilc3,ibc3)%lon
      ytemp = ptr_patch%cells%center(ilc3,ibc3)%lat
      CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

      y(3)  = ytemp-yloc
      x(3)  = xtemp-xloc
      ! This is needed when the date line is crossed
      IF (x(3) >  3.5_wp) x(3) = x(3) - pi2
      IF (x(3) < -3.5_wp) x(3) = x(3) + pi2

      ! The weighting factors are based on the requirement that sum(w(i)*x(i)) = 0
      ! and sum(w(i)*y(i)) = 0, which ensures that linear horizontal gradients
      ! are not aliased into a checkerboard pattern between upward- and downward
      ! directed cells. The third condition is sum(w(i)) = 1., and the weight
      ! of the local point is 0.5 (see above). Analytical elimination yields...

      wgt(3) = 1._wp/( (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))/(x(2)-x(1)) ) * &
                  (1._wp-wgt_loc)*( -y(1) + x(1)*(y(2)-y(1))/(x(2)-x(1)) )
      wgt(2) = (-(1._wp-wgt_loc)*x(1) - wgt(3)*(x(3)-x(1)))/(x(2)-x(1))
      wgt(1) = 1._wp - wgt_loc - wgt(2) - wgt(3)

      ! Store results in ptr_patch%cells%avg_wgt
      ptr_int%c_bln_avg(jc,1,jb) = wgt_loc
      ptr_int%c_bln_avg(jc,2,jb) = wgt(1)
      ptr_int%c_bln_avg(jc,3,jb) = wgt(2)
      ptr_int%c_bln_avg(jc,4,jb) = wgt(3)

    ENDDO !cell loop

  END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

! The coefficients for bilinear interpolation are now iteratively modified
! in order to obtain mass conservation.
! The criterion for conservation is that the three-point divergence
! calculated for any given grid point is used with a total factor of 1

DO iter = 1, niter

  ! Compute sum of weighting coefficients with which
  ! each local divergence value is used
  ! Note: the summation needs to be split into 4 loops in order to
  ! allow for vectorization and parallelization

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

  rl_start = 2
  rl_end = min_rlcell
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,ilc3,ibc3,inb1,inb2,inb3)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx

      ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
      ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
      ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
      ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
      ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
      ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)
      inb1 = inv_neighbor_id(jc,jb,1) + 1
      inb2 = inv_neighbor_id(jc,jb,2) + 1
      inb3 = inv_neighbor_id(jc,jb,3) + 1

      wgt_loc_sum(jc,jb) = &
        ptr_int%c_bln_avg(jc,1,jb)*ptr_patch%cells%area(jc,jb)          + &
        ptr_int%c_bln_avg(ilc1,inb1,ibc1)*ptr_patch%cells%area(ilc1,ibc1)  + &
        ptr_int%c_bln_avg(ilc2,inb2,ibc2)*ptr_patch%cells%area(ilc2,ibc2)  + &
        ptr_int%c_bln_avg(ilc3,inb3,ibc3)*ptr_patch%cells%area(ilc3,ibc3)

    ENDDO !cell loop

  END DO !block loop
!$OMP END DO

  rl_start = 3
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx

      ! For mass conservation, wgt_loc_sum/area should be 1 for each cell
      ! The deviation therefrom is termed residuum here.

      resid(jc,jb) = wgt_loc_sum(jc,jb)/ptr_patch%cells%area(jc,jb)-1._wp

    ENDDO !cell loop

  END DO !block loop
!$OMP END DO
IF (iter < niter) THEN ! Apply iterative correction to weighting coefficients
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,ilc3,ibc3)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx

      ! line and block indices of the neighbouring cells

      ilc1 = ptr_patch%cells%neighbor_idx(jc,jb,1)
      ibc1 = ptr_patch%cells%neighbor_blk(jc,jb,1)
      ilc2 = ptr_patch%cells%neighbor_idx(jc,jb,2)
      ibc2 = ptr_patch%cells%neighbor_blk(jc,jb,2)
      ilc3 = ptr_patch%cells%neighbor_idx(jc,jb,3)
      ibc3 = ptr_patch%cells%neighbor_blk(jc,jb,3)

      ! Modify weighting coefficients

      ptr_int%c_bln_avg(jc,1,jb) = ptr_int%c_bln_avg(jc,1,jb) - relax_coeff*resid(jc,jb)
      ptr_int%c_bln_avg(jc,2,jb) = ptr_int%c_bln_avg(jc,2,jb) - relax_coeff*resid(ilc1,ibc1)
      ptr_int%c_bln_avg(jc,3,jb) = ptr_int%c_bln_avg(jc,3,jb) - relax_coeff*resid(ilc2,ibc2)
      ptr_int%c_bln_avg(jc,4,jb) = ptr_int%c_bln_avg(jc,4,jb) - relax_coeff*resid(ilc3,ibc3)

      wgt_loc_sum(jc,jb) = SUM(ptr_int%c_bln_avg(jc,1:4,jb)) - 1._wp

      ptr_int%c_bln_avg(jc,1,jb) = ptr_int%c_bln_avg(jc,1,jb) - 0.25_wp*wgt_loc_sum(jc,jb)
      ptr_int%c_bln_avg(jc,2,jb) = ptr_int%c_bln_avg(jc,2,jb) - 0.25_wp*wgt_loc_sum(jc,jb)
      ptr_int%c_bln_avg(jc,3,jb) = ptr_int%c_bln_avg(jc,3,jb) - 0.25_wp*wgt_loc_sum(jc,jb)
      ptr_int%c_bln_avg(jc,4,jb) = ptr_int%c_bln_avg(jc,4,jb) - 0.25_wp*wgt_loc_sum(jc,jb)

      ! To be safe: Avoid runaway of central weight
      ptr_int%c_bln_avg(jc,1,jb) = MAX(ptr_int%c_bln_avg(jc,1,jb),minwgt_loc)
      ptr_int%c_bln_avg(jc,1,jb) = MIN(ptr_int%c_bln_avg(jc,1,jb),maxwgt_loc)

    ENDDO !cell loop

  END DO !block loop
!$OMP END DO
ELSE ! In the last iteration, enforce the mass conservation condition
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx

      ! Modify weighting coefficients

      ptr_int%c_bln_avg(jc,1,jb) = ptr_int%c_bln_avg(jc,1,jb) - resid(jc,jb)

    ENDDO !cell loop

  END DO !block loop
!$OMP END DO

! Compute coefficients needed to reconstruct averaged mass fluxes
! for approximately mass-consistent transport with divergence-averaging
! They can alternatively be used to average the velocity going into the divergence
! computation (without div averaging), yielding exact mass consistency but somewhat
! larger discretization errors for divergence

  rl_start = 4
  rl_end   = min_rledge
  i_startblk = ptr_patch%edges%start_blk(rl_start,1)
  i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,inb1,inb2,inb3,ie4,ie5)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je = i_startidx, i_endidx

      ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
      ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
      ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
      ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

      IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,1) .AND. &
          jb == ptr_patch%cells%edge_blk(ilc1,ibc1,1)) THEN

        inb1 = inv_neighbor_id(ilc1,ibc1,1)
        ie4  = MOD(inb1,  3)+1
        ie5  = MOD(inb1+1,3)+1

!         write(0,*) je,jb,ilc2,inb1,ibc2,ilc1,ibc1,ilc2,inb1,ibc2
!         write(0,*) ptr_int%c_bln_avg(ilc2,inb1+1,ibc2)
!         write(0,*) ptr_int%geofac_div(ilc1,2,ibc1)
!         write(0,*) ptr_int%geofac_div(ilc2,inb1,ibc2)

        ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb1+1,ibc2) &
          *ptr_int%geofac_div(ilc1,2,ibc1)/ptr_int%geofac_div(ilc2,inb1,ibc2)

        ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb1+1,ibc2) &
          *ptr_int%geofac_div(ilc1,3,ibc1)/ptr_int%geofac_div(ilc2,inb1,ibc2)

        ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,2,ibc1) &
          *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,1,ibc1)

        ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,2,ibc1) &
          *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,1,ibc1)


      ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,2) .AND. &
               jb == ptr_patch%cells%edge_blk(ilc1,ibc1,2)) THEN

        inb2 = inv_neighbor_id(ilc1,ibc1,2)
        ie4  = MOD(inb2  ,3)+1
        ie5  = MOD(inb2+1,3)+1

        ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb2+1,ibc2) &
          *ptr_int%geofac_div(ilc1,3,ibc1)/ptr_int%geofac_div(ilc2,inb2,ibc2)

        ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb2+1,ibc2) &
          *ptr_int%geofac_div(ilc1,1,ibc1)/ptr_int%geofac_div(ilc2,inb2,ibc2)

        ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,3,ibc1) &
          *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,2,ibc1)

        ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,3,ibc1) &
          *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,2,ibc1)


      ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,3) .AND. &
               jb == ptr_patch%cells%edge_blk(ilc1,ibc1,3)) THEN

        inb3 = inv_neighbor_id(ilc1,ibc1,3)
        ie4  = MOD(inb3  ,3)+1
        ie5  = MOD(inb3+1,3)+1

        ptr_int%e_flx_avg(je,2,jb) = ptr_int%c_bln_avg(ilc2,inb3+1,ibc2) &
          *ptr_int%geofac_div(ilc1,1,ibc1)/ptr_int%geofac_div(ilc2,inb3,ibc2)

        ptr_int%e_flx_avg(je,3,jb) = ptr_int%c_bln_avg(ilc2,inb3+1,ibc2) &
          *ptr_int%geofac_div(ilc1,2,ibc1)/ptr_int%geofac_div(ilc2,inb3,ibc2)

        ptr_int%e_flx_avg(je,4,jb) = ptr_int%c_bln_avg(ilc1,4,ibc1) &
          *ptr_int%geofac_div(ilc2,ie4,ibc2)/ptr_int%geofac_div(ilc1,3,ibc1)

        ptr_int%e_flx_avg(je,5,jb) = ptr_int%c_bln_avg(ilc1,4,ibc1) &
          *ptr_int%geofac_div(ilc2,ie5,ibc2)/ptr_int%geofac_div(ilc1,3,ibc1)

      ENDIF

    ENDDO !edge loop

  END DO !block loop
!$OMP END DO

  rl_start = 5
  i_startblk = ptr_patch%edges%start_blk(rl_start,1)

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,inb1,inb2,inb3,ie4,ie5, &
!$OMP            ile1,ibe1,ile2,ibe2,ile3,ibe3,ile4,ibe4,iie1,iie2,iie3,iie4)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je = i_startidx, i_endidx

      ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
      ibc1 = ptr_patch%edges%cell_blk(je,jb,1)
      ilc2 = ptr_patch%edges%cell_idx(je,jb,2)
      ibc2 = ptr_patch%edges%cell_blk(je,jb,2)

      ile1 = ptr_patch%edges%quad_idx(je,jb,1)
      ibe1 = ptr_patch%edges%quad_blk(je,jb,1)
      ile2 = ptr_patch%edges%quad_idx(je,jb,2)
      ibe2 = ptr_patch%edges%quad_blk(je,jb,2)
      ile3 = ptr_patch%edges%quad_idx(je,jb,3)
      ibe3 = ptr_patch%edges%quad_blk(je,jb,3)
      ile4 = ptr_patch%edges%quad_idx(je,jb,4)
      ibe4 = ptr_patch%edges%quad_blk(je,jb,4)

      IF (ptr_patch%edges%cell_idx(ile1,ibe1,1) == ilc1 .AND. &
          ptr_patch%edges%cell_blk(ile1,ibe1,1) == ibc1 ) THEN
        iie1(je) = 3
      ELSE IF (ptr_patch%edges%cell_idx(ile1,ibe1,2) == ilc1 .AND. &
               ptr_patch%edges%cell_blk(ile1,ibe1,2) == ibc1 ) THEN
        iie1(je) = 5
      ENDIF
      IF (ptr_patch%edges%cell_idx(ile2,ibe2,1) == ilc1 .AND. &
          ptr_patch%edges%cell_blk(ile2,ibe2,1) == ibc1 ) THEN
        iie2(je) = 2
      ELSE IF (ptr_patch%edges%cell_idx(ile2,ibe2,2) == ilc1 .AND. &
               ptr_patch%edges%cell_blk(ile2,ibe2,2) == ibc1 ) THEN
        iie2(je) = 4
      ENDIF
      IF (ptr_patch%edges%cell_idx(ile3,ibe3,1) == ilc2 .AND. &
          ptr_patch%edges%cell_blk(ile3,ibe3,1) == ibc2 ) THEN
        iie3(je) = 3
      ELSE IF (ptr_patch%edges%cell_idx(ile3,ibe3,2) == ilc2 .AND. &
               ptr_patch%edges%cell_blk(ile3,ibe3,2) == ibc2 ) THEN
        iie3(je) = 5
      ENDIF
      IF (ptr_patch%edges%cell_idx(ile4,ibe4,1) == ilc2 .AND. &
          ptr_patch%edges%cell_blk(ile4,ibe4,1) == ibc2 ) THEN
        iie4(je) = 2
      ELSE IF (ptr_patch%edges%cell_idx(ile4,ibe4,2) == ilc2 .AND. &
               ptr_patch%edges%cell_blk(ile4,ibe4,2) == ibc2 ) THEN
        iie4(je) = 4
      ENDIF

      IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,1) .AND. &
          jb == ptr_patch%cells%edge_blk(ilc1,ibc1,1)) THEN

        inb1 = inv_neighbor_id(ilc1,ibc1,1)
        ie4  = MOD(inb1  ,3)+1
        ie5  = MOD(inb1+1,3)+1

        ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
          ( ptr_int%geofac_div(ilc1,1,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
          + ptr_int%geofac_div(ilc2,inb1,ibc2)*ptr_int%c_bln_avg(ilc1,2,ibc1)         &
          - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,2,ibc1)     &
          - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,3,ibc1) )   &
          / ptr_int%geofac_div(ilc1,1,ibc1)   +                                       &
          ( ptr_int%geofac_div(ilc2,inb1,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
          + ptr_int%geofac_div(ilc1,1,ibc1)*ptr_int%c_bln_avg(ilc2,inb1+1,ibc2)       &
          - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
          - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
          / ptr_int%geofac_div(ilc2,inb1,ibc2) )


      ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,2) .AND. &
               jb == ptr_patch%cells%edge_blk(ilc1,ibc1,2)) THEN

        inb2 = inv_neighbor_id(ilc1,ibc1,2)
        ie4  = MOD(inb2  ,3)+1
        ie5  = MOD(inb2+1,3)+1

        ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
          ( ptr_int%geofac_div(ilc1,2,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
          + ptr_int%geofac_div(ilc2,inb2,ibc2)*ptr_int%c_bln_avg(ilc1,3,ibc1)         &
          - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,3,ibc1)     &
          - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,1,ibc1) )   &
          / ptr_int%geofac_div(ilc1,2,ibc1)   +                                       &
          ( ptr_int%geofac_div(ilc2,inb2,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
          + ptr_int%geofac_div(ilc1,2,ibc1)*ptr_int%c_bln_avg(ilc2,inb2+1,ibc2)       &
          - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
          - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
          / ptr_int%geofac_div(ilc2,inb2,ibc2) )


      ELSE IF (je == ptr_patch%cells%edge_idx(ilc1,ibc1,3) .AND. &
               jb == ptr_patch%cells%edge_blk(ilc1,ibc1,3)) THEN

        inb3 = inv_neighbor_id(ilc1,ibc1,3)
        ie4  = MOD(inb3  ,3)+1
        ie5  = MOD(inb3+1,3)+1

        ptr_int%e_flx_avg(je,1,jb) = 0.5_wp *(                                        &
          ( ptr_int%geofac_div(ilc1,3,ibc1)*ptr_int%c_bln_avg(ilc1,1,ibc1)            &
          + ptr_int%geofac_div(ilc2,inb3,ibc2)*ptr_int%c_bln_avg(ilc1,4,ibc1)         &
          - ptr_int%e_flx_avg(ile1,iie1(je),ibe1)*ptr_int%geofac_div(ilc1,1,ibc1)     &
          - ptr_int%e_flx_avg(ile2,iie2(je),ibe2)*ptr_int%geofac_div(ilc1,2,ibc1) )   &
          / ptr_int%geofac_div(ilc1,3,ibc1)   +                                       &
          ( ptr_int%geofac_div(ilc2,inb3,ibc2)*ptr_int%c_bln_avg(ilc2,1,ibc2)         &
          + ptr_int%geofac_div(ilc1,3,ibc1)*ptr_int%c_bln_avg(ilc2,inb3+1,ibc2)       &
          - ptr_int%e_flx_avg(ile3,iie3(je),ibe3)*ptr_int%geofac_div(ilc2,ie4,ibc2)   &
          - ptr_int%e_flx_avg(ile4,iie4(je),ibe4)*ptr_int%geofac_div(ilc2,ie5,ibc2) ) &
          / ptr_int%geofac_div(ilc2,inb3,ibc2) )

      ENDIF

    ENDDO !edge loop

  END DO !block loop
!$OMP END DO

  ! Finally, the weighting coefficients are scaled in order to
  ! yield the right result for a constant wind field

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ile1,ibe1,ile2,ibe2,ile3,ibe3,ile4,ibe4, &
!$OMP            z_nx1,z_nx2,z_nx3,z_nx4,z_nx5)
  DO jb = i_startblk, i_endblk

    CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    DO je = i_startidx, i_endidx

      ile1 = ptr_patch%edges%quad_idx(je,jb,1)
      ibe1 = ptr_patch%edges%quad_blk(je,jb,1)
      ile2 = ptr_patch%edges%quad_idx(je,jb,2)
      ibe2 = ptr_patch%edges%quad_blk(je,jb,2)
      ile3 = ptr_patch%edges%quad_idx(je,jb,3)
      ibe3 = ptr_patch%edges%quad_blk(je,jb,3)
      ile4 = ptr_patch%edges%quad_idx(je,jb,4)
      ibe4 = ptr_patch%edges%quad_blk(je,jb,4)

      z_nx1(je,1:3) = ptr_patch%edges%primal_cart_normal(je,jb)%x(1:3)
      z_nx2(je,1:3) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(1:3)
      z_nx3(je,1:3) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(1:3)
      z_nx4(je,1:3) = ptr_patch%edges%primal_cart_normal(ile3,ibe3)%x(1:3)
      z_nx5(je,1:3) = ptr_patch%edges%primal_cart_normal(ile4,ibe4)%x(1:3)

      ! The sum of the coefficients - multiplied by the projection factors -
      ! is enforced to be 1 so that a constant vector field is processed correctly

      checksum(je,jb) = ptr_int%e_flx_avg(je,1,jb)                            &
        + DOT_PRODUCT(z_nx1(je,1:3),z_nx2(je,1:3))*ptr_int%e_flx_avg(je,2,jb) &
        + DOT_PRODUCT(z_nx1(je,1:3),z_nx3(je,1:3))*ptr_int%e_flx_avg(je,3,jb) &
        + DOT_PRODUCT(z_nx1(je,1:3),z_nx4(je,1:3))*ptr_int%e_flx_avg(je,4,jb) &
        + DOT_PRODUCT(z_nx1(je,1:3),z_nx5(je,1:3))*ptr_int%e_flx_avg(je,5,jb)

      ptr_int%e_flx_avg(je,1,jb) = ptr_int%e_flx_avg(je,1,jb)/checksum(je,jb)
      ptr_int%e_flx_avg(je,2,jb) = ptr_int%e_flx_avg(je,2,jb)/checksum(je,jb)
      ptr_int%e_flx_avg(je,3,jb) = ptr_int%e_flx_avg(je,3,jb)/checksum(je,jb)
      ptr_int%e_flx_avg(je,4,jb) = ptr_int%e_flx_avg(je,4,jb)/checksum(je,jb)
      ptr_int%e_flx_avg(je,5,jb) = ptr_int%e_flx_avg(je,5,jb)/checksum(je,jb)

    ENDDO !edge loop

  END DO !block loop
!$OMP END DO

ENDIF ! end of last iteration
!$OMP END PARALLEL
ENDDO ! iteration loop


END SUBROUTINE init_cellavg_wgt

 !----------------------------------------------------------------------------
  !>
  !! Primal cell quadrature points and weights
  !!
  !! Computes quadrature points and weights for triangular grid cells.
  !! Quadrature points and weights are provided for accurately integrating 
  !! linear, quadratic and cubic functions. This is necessary for initializing 
  !! idealized testcases.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-11-16)
  !!
  !! @par Literature
  !! Numerical Methods in Engineering with Python, Jaan Kiusalaas (2005), 
  !! 233-247
  SUBROUTINE tri_quadrature_pts (ptr_patch, ptr_int)

    TYPE(t_patch), INTENT(IN) :: ptr_patch  !< patch

    TYPE(t_int_state), INTENT(INOUT) :: ptr_int  !< interpolation state

    REAL(wp) ::  alpha_l(3),        & !< area coordinates for quadrature up to 
      &  alpha_q(3,3), alpha_c(3,4)   !< fourth order
                                      !< (n_area_coords,n_pts))

    TYPE(t_cartesian_coordinates)    :: z_vert_cc(3) ! cell vertices in cartesian
                                                     ! coordinates
    TYPE(t_cartesian_coordinates)    :: z_quad_cc  ! triangle quadrature point in cartesian
                                                   ! coordinates
    TYPE(t_geographical_coordinates) :: z_quad_gg  ! triangle quadrature point in geographical
                                                   ! coordinates

    INTEGER :: ilv, ibv               !< line and block indices of cell vertices
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rcstartlev
    INTEGER :: nq, nv                 !< loop index for quadrature points and
                                      !< cell vertices
    INTEGER :: jc, jb                 !< loop index for cells

  !-------------------------------------------------------------------------

    IF (msg_level > 10) &
      CALL message('mo_interpolation:tri_quadrature_pts', '')

    ! set area coordinates
    !
    ! linear
    alpha_l(1) = 1._wp/3._wp
    alpha_l(2) = 1._wp/3._wp
    alpha_l(3) = 1._wp/3._wp

    !
    ! quadratic
    ! 
    alpha_q(1,1) = 0.5_wp 
    alpha_q(2,1) = 0._wp
    alpha_q(3,1) = 0.5_wp

    alpha_q(1,2) = 0.5_wp
    alpha_q(2,2) = 0.5_wp
    alpha_q(3,2) = 0._wp

    alpha_q(1,3) = 0._wp
    alpha_q(2,3) = 0.5_wp
    alpha_q(3,3) = 0.5_wp
 
    !
    ! cubic
    !
    alpha_c(1,1) = 1._wp/3._wp 
    alpha_c(2,1) = 1._wp/3._wp
    alpha_c(3,1) = 1._wp/3._wp

    alpha_c(1,2) = 1._wp/5._wp
    alpha_c(2,2) = 1._wp/5._wp
    alpha_c(3,2) = 3._wp/5._wp

    alpha_c(1,3) = 3._wp/5._wp
    alpha_c(2,3) = 1._wp/5._wp
    alpha_c(3,3) = 1._wp/5._wp

    alpha_c(1,4) = 1._wp/5._wp
    alpha_c(2,4) = 3._wp/5._wp
    alpha_c(3,4) = 1._wp/5._wp

    ! note that the linear weighting factor is 1 (not stored)
    ptr_int%gquad%weights_tri_q(1:3) = 1._wp/3._wp
    ptr_int%gquad%weights_tri_c(1)   = -27._wp/48._wp
    ptr_int%gquad%weights_tri_c(2:4) = 25._wp/48._wp


    i_rcstartlev = 2

    ! start and end block
    i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)
    i_endblk   = ptr_patch%nblks_int_c

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,nv,nq,i_startidx,i_endidx,ilv,ibv,z_vert_cc,z_quad_cc, &
!$OMP           z_quad_gg)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rcstartlev)

      DO jc = i_startidx, i_endidx

       ! loop over triangle vertices
!CDIR EXPAND=3
       DO nv=1,3
         ! get line and block indices of cell vertices
         ilv= ptr_patch%cells%vertex_idx(jc,jb,nv)
         ibv= ptr_patch%cells%vertex_blk(jc,jb,nv)

         ! Transform geographical coordinates to cartesian coordinates for vertices
          z_vert_cc(nv)=gc2cc(ptr_patch%verts%vertex(ilv,ibv))
       ENDDO

       !
       ! Linear
       !
       ! Compute quadrature point in cartesian coordinates (= triangle centroid)
       ! i.e. map area coordinates into cartesian coordinates
       z_quad_cc%x(1)= alpha_l(1)*z_vert_cc(1)%x(1)  &
         &           + alpha_l(2)*z_vert_cc(2)%x(1)  &
         &           + alpha_l(3)*z_vert_cc(3)%x(1)

       z_quad_cc%x(2)= alpha_l(1)*z_vert_cc(1)%x(2)  &
         &           + alpha_l(2)*z_vert_cc(2)%x(2)  &
         &           + alpha_l(3)*z_vert_cc(3)%x(2)

       z_quad_cc%x(3)= alpha_l(1)*z_vert_cc(1)%x(3)  &
         &           + alpha_l(2)*z_vert_cc(2)%x(3)  &
         &           + alpha_l(3)*z_vert_cc(3)%x(3)


       ! Transform back to geographical coordinates
       z_quad_gg = cc2gc(z_quad_cc)

       ! store
       ptr_int%gquad%qpts_tri_l(jc,jb)%lat = z_quad_gg%lat
       ptr_int%gquad%qpts_tri_l(jc,jb)%lon = z_quad_gg%lon


       !
       ! quadratic
       !
       ! Loop over quadrature points
!CDIR EXPAND=3
       DO nq=1,3
         ! map area coordinates into cartesian coordinates
         z_quad_cc%x(1)= alpha_q(1,nq)*z_vert_cc(1)%x(1)  &
           &           + alpha_q(2,nq)*z_vert_cc(2)%x(1)  &
           &           + alpha_q(3,nq)*z_vert_cc(3)%x(1)

         z_quad_cc%x(2)= alpha_q(1,nq)*z_vert_cc(1)%x(2)  &
           &           + alpha_q(2,nq)*z_vert_cc(2)%x(2)  &
           &           + alpha_q(3,nq)*z_vert_cc(3)%x(2)

         z_quad_cc%x(3)= alpha_q(1,nq)*z_vert_cc(1)%x(3)  &
           &           + alpha_q(2,nq)*z_vert_cc(2)%x(3)  &
           &           + alpha_q(3,nq)*z_vert_cc(3)%x(3)


         ! Transform back to geographical coordinates
         z_quad_gg = cc2gc(z_quad_cc)

         ! store
         ptr_int%gquad%qpts_tri_q(jc,jb,nq)%lat = z_quad_gg%lat
         ptr_int%gquad%qpts_tri_q(jc,jb,nq)%lon = z_quad_gg%lon
       ENDDO


       !
       ! cubic
       !
       ! Loop over quadrature points
!CDIR EXPAND=4
       DO nq=1,4
         ! map area coordinates into cartesian coordinates
         z_quad_cc%x(1)= alpha_c(1,nq)*z_vert_cc(1)%x(1)  &
           &           + alpha_c(2,nq)*z_vert_cc(2)%x(1)  &
           &           + alpha_c(3,nq)*z_vert_cc(3)%x(1)

         z_quad_cc%x(2)= alpha_c(1,nq)*z_vert_cc(1)%x(2)  &
           &           + alpha_c(2,nq)*z_vert_cc(2)%x(2)  &
           &           + alpha_c(3,nq)*z_vert_cc(3)%x(2)

         z_quad_cc%x(3)= alpha_c(1,nq)*z_vert_cc(1)%x(3)  &
           &           + alpha_c(2,nq)*z_vert_cc(2)%x(3)  &
           &           + alpha_c(3,nq)*z_vert_cc(3)%x(3)


         ! Transform back to geographical coordinates
         z_quad_gg = cc2gc(z_quad_cc)

         ! store
         ptr_int%gquad%qpts_tri_c(jc,jb,nq)%lat = z_quad_gg%lat
         ptr_int%gquad%qpts_tri_c(jc,jb,nq)%lon = z_quad_gg%lon
       ENDDO

      END DO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE tri_quadrature_pts

END MODULE mo_interpolation_init
