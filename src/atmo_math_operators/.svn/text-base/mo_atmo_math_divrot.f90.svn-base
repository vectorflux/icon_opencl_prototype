!>
!!   Contains the implementation of the div,rot,recon mathematical operators.
!!
!!   Contains the implementation of the mathematical operators
!!   employed by the shallow water prototype.
!!
!! @par Revision History
!!  Developed  by Luca Bonaventura and Will Sawyer (2002-4).
!!  Modified to ProTeX-style by  Luca Bonaventura and Thomas Heinze (2004).
!!  Adapted to new data structure by Thomas Heinze,
!!  Peter Korn and Luca Bonaventura (2005).
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Subroutine for divergence multiplied by area added by P.Korn (2006).
!!  Modification by Peter Korn, MPI-M, (2006-11-23):
!!  - replacements in TYPE patch: ic by l2g_c, ie by l2g_e, iv by l2g_v,
!!    iic by g2l_c, iie by g2l_e, iiv by g2l_v
!!  - replaced edge_index by edge_idx
!!  - replaced vertex_index by vertex_idx
!!  - replaced cell_index by cell_idx
!!  - replaced neighbor_index by neighbor_idx
!!  - replaced child_index by child_idx
!!  Modified by P Ripodas (2007-02):
!!  - include the system orientation factor in the vorticity term of nabla2_vec
!!  - solved errors in nabla4_vec and nabla4_scalar
!!  Modification by Peter Korn, MPI-M, (2006/2007):
!!  -operator overloading of curl operator and nabla2vec to handle atmosphere and ocean version
!!  -change of input/output arguments of subroutines: arrays of fixed size are
!!   changed to pointers, to avoid occurence of not-initialized numbers
!!  Modified by Almut Gassmann, MPI-M, (2007-04)
!!  - removed references to unused halo_verts
!!  - summing over all halos corresponding to different parallel patches
!!  Modified by Hui Wan, MPI-M, (2007-11)
!!  - added subroutine cell_avg
!!  Modified by Hui Wan, MPI-M, (2008-04-04)
!!  - control variable loce renamed locean
!!  Modification by Jochen Foerstner, DWD, (2008-05-05)
!!  - div and div_times_area are now generic subroutines
!!  - the divergence can now be computed either
!!    using the midpoint rule
!!    (div_midpoint, div_midpoint_times_area) or
!!    using the Simpson's rule
!!    (div_simpson, div_simpson_times_area)
!!  Modification by Jochen Foerstner, DWD, (2008-07-16)
!!  - introduction of several new operators (to be) used in combination with
!!    the tracer advection:
!!    grad_green_gauss_cell, grad_green_gauss_edge and div_quad_twoadjcells
!!    (for the new div operator there is again a version using the midpoint
!!    and a version using the Simpson's rule).
!!    The first operator is used to calculate a cell centered value of the
!!    gradient for the piecewise linear reconstruction. The second and third
!!    will be used in combination with the MPDATA scheme. Both deal with the
!!    quadrilateral control volumes formed by two adjacent triangles.
!!  Modification by Marco Restelli, MPI (2008-07-17)
!!  - included subroutine dtan.
!!  Modification by Jochen Foerstner, DWD (2008-09-12)
!!  - moved SUBROUTINE ravtom_normgrad2 from mo_interpolation to this module
!!    because of conflicting use statements.
!!  Modification by Jochen Foerstner, DWD (2008-09-16)
!!  - removed SUBROUTINE ravtom_normgrad2 (not used)
!!  Modification by Daniel Reinert, DWD (2009-07-20)
!!  - added subroutine grad_lsq_cell for gradient reconstruction via the
!!    least-squares method and grad_green_gauss_gc_cell for Green-Gauss
!!    gradient in geographical coordinates
!!  Modification by Daniel Reinert, DWD (2009-12-14)
!!  - renamed grad_lsq_cell -> recon_lsq_cell_l
!!  Modification by Leonidas Linardakis, MPI-M (2010-21-01)
!!  - splitted mo_math_operators into submodules
!!  Modification by Daniel Reinert, DWD (2010-04-12)
!!  - added subroutine recon_lsq_cell_q for third order accurate least-squares
!!    reconstruction of an arbitrary field. Based on cell centered values.
!!  Modification by Daniel Reinert, DWD (2010-10-14)
!!  - added subroutine recon_lsq_cell_c for fitting a cubic polynomial in a least
!!    squares sense. Based on cell centered values.
!!
!! @par To Do
!! Boundary exchange, nblks in presence of halos and dummy edge
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
MODULE mo_atmo_math_divrot
!-------------------------------------------------------------------------
!
USE mo_kind,               ONLY: wp
USE mo_impl_constants,     ONLY: SEA, min_rlcell, min_rledge, min_rlvert
USE mo_interpolation,      ONLY: t_int_state
USE mo_model_domain,       ONLY: t_patch
! USE mo_model_domain_import, ONLY: l_limited_area
USE mo_global_variables,   ONLY: nproma, i_cell_type, nlev !, ltimer
USE mo_exception,          ONLY: finish
! USE mo_timer,              ONLY: timer_start, timer_stop, timer_div
USE mo_loopindices,        ONLY: get_indices_c, get_indices_e, get_indices_v

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'


PUBLIC :: div, div_avg
PUBLIC :: div_quad_twoadjcells
PUBLIC :: rot_vertex
PUBLIC :: rot_vertex_atmos


INTERFACE rot_vertex

  MODULE PROCEDURE rot_vertex_atmos

END INTERFACE


CONTAINS



!-------------------------------------------------------------------------
!>
!! Computes discrete divergence of a vector field.
!!
!! Computes discrete divergence of a vector field
!! given by its components in the directions normal to triangle edges.
!! The midpoint rule is used for quadrature.
!! input:  lives on edges (velocity points)
!! output: lives on centers of triangles
!!
!! @par Revision History
!! Developed  by  Luca Bonaventura, MPI-M (2002-5).
!! Changes according to programming guide by Thomas Heinze, DWD (2006-08-18).
!! Modification by Thomas Heinze, DWD (2006-09-11):
!! - loop only over the inner cells of a patch, not any more over halo cells
!! Modifications by P. Korn, MPI-M(2007-2)
!! -Switch fom array arguments to pointers
!! Modification by Almut Gassmann, MPI-M (2007-04-20)
!! - abandon grid for the sake of patch
!! Modification by Guenther Zaengl, DWD (2009-03-17)
!! - vector optimization
!! Modification by Guenther Zaengl, DWD (2010-08-20) :
!! - Option for processing two fields at once for efficiency optimization
!!
SUBROUTINE div( vec_e, ptr_patch, ptr_int, div_vec_c, &
  &             opt_slev, opt_elev, opt_in2, opt_out2, opt_rlstart, opt_rlend )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
! edge based variable of which divergence
! is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

! optional second input field for more efficient processing in NH core
REAL(wp), OPTIONAL, INTENT(in) ::  &
  &  opt_in2(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! cell based variable in which divergence is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  div_vec_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

! optional second output field
REAL(wp), OPTIONAL, INTENT(inout) ::  &
  &  opt_out2(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: nlen, npromz_c, nblks_c

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
LOGICAL :: l2fields

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = nlev
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF

IF ( PRESENT(opt_in2) .AND. PRESENT(opt_out2)) THEN
  l2fields = .TRUE.
ELSE
  l2fields = .FALSE.
ENDIF

iidx => ptr_patch%cells%edge_idx
iblk => ptr_patch%cells%edge_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

! loop through all patch cells (and blocks)
!

!IF(ltimer) CALL timer_start(timer_div)

!$OMP PARALLEL


SELECT CASE (i_cell_type)

CASE (3) ! (i_cell_type == 3)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    ! original comment for divergence computation;
    ! everything that follows in this explanation has been combined into geofac_div

    ! compute the discrete divergence for cell jc by finite volume
    ! approximation (see Bonaventura and Ringler MWR 2005);
    ! multiplication of the normal vector component vec_e at the edges
    ! by the appropriate cell based edge_orientation is required to
    ! obtain the correct value for the application of Gauss theorem
    ! (which requires the scalar product of the vector field with the
    ! OUTWARD pointing unit vector with respect to cell jc; since the
    ! positive direction for the vector components is not necessarily
    ! the outward pointing one with respect to cell jc, a correction
    ! coefficient (equal to +-1) is necessary, given by
    ! ptr_patch%grid%cells%edge_orientation)

    IF (l2fields) THEN

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
!CDIR UNROLL=5
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif

          div_vec_c(jc,jk,jb) =  &
            vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

          opt_out2(jc,jk,jb) =  &
            opt_in2(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            opt_in2(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            opt_in2(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        END DO
      END DO

    ELSE

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
!CDIR UNROLL=6
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif

          div_vec_c(jc,jk,jb) =  &
            vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        END DO
      END DO

    ENDIF
  END DO
!$OMP END DO

CASE (6) ! (i_cell_type == 6)

  ! no grid refinement in hexagonal model
  nblks_c   = ptr_patch%nblks_int_c
  npromz_c  = ptr_patch%npromz_int_c

!$OMP DO PRIVATE(jb,nlen,jc,jk)
  DO jb = 1, nblks_c

    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF

#ifdef __LOOP_EXCHANGE
    DO jc = 1, nlen
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jc = 1, nlen
#endif

          div_vec_c(jc,jk,jb) =   &
          &   vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) &
          & + vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) &
          & + vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb) &
          & + vec_e(iidx(jc,jb,4),jk,iblk(jc,jb,4)) * ptr_int%geofac_div(jc,4,jb) &
          & + vec_e(iidx(jc,jb,5),jk,iblk(jc,jb,5)) * ptr_int%geofac_div(jc,5,jb) &
          & + vec_e(iidx(jc,jb,6),jk,iblk(jc,jb,6)) * ptr_int%geofac_div(jc,6,jb)

      END DO
    END DO
  END DO
!$OMP END DO

END SELECT
!$OMP END PARALLEL

!IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div


!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field.
!!
!! Computes discrete divergence of a vector field
!! given by its components in the directions normal to triangle edges,
!! followed by bilinear averaging to remove checkerboard noise
!! (Combines div_midpoint and cell_avg_varwgt to increase computing efficiency)
!!
!! @par Revision History
!! Developed by Guenther Zaengl, DWD (2009-03-30)
!! Modification by Guenther Zaengl, DWD (2010-04-20) :
!! - Option for processing two fields at once for efficiency optimization
!!
SUBROUTINE div_avg( vec_e, ptr_patch, ptr_int, avg_coeff, div_vec_c, &
  &                 opt_in2, opt_out2, opt_rlstart, opt_rlend )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int

!  averaging coefficients
REAL(wp), INTENT(in) :: avg_coeff(:,:,:) ! dim: (nproma,nlev,nblks_c)
!
! edge based variable of which divergence
! is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

! optional second input field for more efficient processing in NH core
REAL(wp), OPTIONAL, INTENT(in) ::  &
  &  opt_in2(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! cell based variable in which divergence is stored
!
REAL(wp), INTENT(inout) ::  &
  &  div_vec_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

! optional second output field
REAL(wp), OPTIONAL, INTENT(inout) ::  &
  &  opt_out2(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end, rl_start_l2
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

REAL(wp), DIMENSION (nproma,nlev,ptr_patch%nblks_c) :: aux_c, aux_c2

INTEGER,  DIMENSION(:,:,:),   POINTER :: inidx, inblk, ieidx, ieblk
LOGICAL :: l2fields

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF
IF ( PRESENT(opt_in2) .AND. PRESENT(opt_out2)) THEN
  l2fields = .TRUE.
ELSE
  l2fields = .FALSE.
ENDIF

rl_start_l2 = rl_start + 1

inidx => ptr_patch%cells%neighbor_idx
inblk => ptr_patch%cells%neighbor_blk
ieidx => ptr_patch%cells%edge_idx
ieblk => ptr_patch%cells%edge_blk


! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)

! First compute divergence
!
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), SCHEDULE(runtime)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  IF (l2fields) THEN

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = 1, nlev
#else
!CDIR UNROLL=5
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
#endif

        aux_c(jc,jk,jb) =  &
          vec_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
          vec_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
          vec_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        aux_c2(jc,jk,jb) =  &
          opt_in2(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
          opt_in2(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
          opt_in2(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

      END DO
    END DO
  ELSE

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = 1, nlev
#else
!CDIR UNROLL=6
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
#endif
        aux_c(jc,jk,jb) =  &
          vec_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
          vec_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
          vec_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

      END DO
    END DO
  ENDIF

END DO
!$OMP END DO

! IF (l_limited_area .OR. ptr_patch%id > 1) THEN
!   ! Fill div_vec_c along the lateral boundaries of nests
! 
!   i_startblk = ptr_patch%cells%start_blk(rl_start,1)
!   i_endblk   = ptr_patch%cells%end_blk(rl_start_l2,1)
! !
!   IF (l2fields) THEN
! !$OMP WORKSHARE
!      div_vec_c(:,:,i_startblk:i_endblk) =  aux_c (:,:,i_startblk:i_endblk)
!      opt_out2 (:,:,i_startblk:i_endblk) =  aux_c2(:,:,i_startblk:i_endblk)
! !$OMP END WORKSHARE
!   ELSE
! !$OMP WORKSHARE
!      div_vec_c(:,:,i_startblk:i_endblk) =  aux_c(:,:,i_startblk:i_endblk)
! !$OMP END WORKSHARE
!   ENDIF
! ENDIF

!
! Now do averaging with weights given by avg_coeff

! values for the blocking
i_startblk = ptr_patch%cells%start_blk(rl_start_l2,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk), SCHEDULE(runtime)
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start_l2, rl_end)

  IF (l2fields) THEN

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = 1, nlev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
#endif

        !  calculate the weighted average
        div_vec_c(jc,jk,jb) =  &
          &    aux_c(jc,jk,jb)                         * avg_coeff(jc,1,jb) &
          &  + aux_c(inidx(jc,jb,1),jk,inblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + aux_c(inidx(jc,jb,2),jk,inblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + aux_c(inidx(jc,jb,3),jk,inblk(jc,jb,3)) * avg_coeff(jc,4,jb)

        opt_out2(jc,jk,jb) =  &
          &    aux_c2(jc,jk,jb)                         * avg_coeff(jc,1,jb) &
          &  + aux_c2(inidx(jc,jb,1),jk,inblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + aux_c2(inidx(jc,jb,2),jk,inblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + aux_c2(inidx(jc,jb,3),jk,inblk(jc,jb,3)) * avg_coeff(jc,4,jb)

      END DO !cell loop
    END DO !vertical levels loop
  ELSE

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = 1, nlev
#else
#ifdef _URD
!CDIR UNROLL=_URD
#endif
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
#endif

        !  calculate the weighted average
        div_vec_c(jc,jk,jb) =  &
          &    aux_c(jc,jk,jb)                         * avg_coeff(jc,1,jb) &
          &  + aux_c(inidx(jc,jb,1),jk,inblk(jc,jb,1)) * avg_coeff(jc,2,jb) &
          &  + aux_c(inidx(jc,jb,2),jk,inblk(jc,jb,2)) * avg_coeff(jc,3,jb) &
          &  + aux_c(inidx(jc,jb,3),jk,inblk(jc,jb,3)) * avg_coeff(jc,4,jb)

      END DO !cell loop
    END DO !vertical levels loop
  ENDIF

END DO !block loop
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE div_avg

!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field for the quadrilateral.
!!
!! Computes discrete divergence of a vector field for the quadrilateral
!! control volume formed by two adjacent triangles
!! The vector field is given by its components in the directions normal
!! to the edges.
!! The midpoint rule is used for quadrature.
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
!! @par Revision History
!! Developed  by  Jochen Foerstner, DWD (2008-07-16).
!!
SUBROUTINE div_quad_twoadjcells( vec_e, ptr_patch, ptr_int, div_vec_e,  &
  &                              opt_slev, opt_elev,                    &
                                 opt_rlstart, opt_rlend )
!
!
! patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
! edge based variable of which divergence is computed
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! edge based variable in which divergence is stored
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  div_vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

!

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: je, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER,  DIMENSION(:,:,:),   POINTER :: iqidx, iqblk

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = nlev
END IF


IF ( PRESENT(opt_rlstart) ) THEN
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_math_operators:div_quad_twoadjcells_midpoint',  &
          &      'opt_rlstart must not be equal to 1')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF

! Set pointers to required index lists
iqidx => ptr_patch%edges%quad_idx
iqblk => ptr_patch%edges%quad_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%edges%start_blk(rl_start,1)
i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!
! loop through all patch edges (and blocks)
!
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx)
DO jb = i_startblk, i_endblk

  CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
  DO je = i_startidx, i_endidx
    DO jk = slev, elev
#else
!CDIR UNROLL=3
  DO jk = slev, elev
    DO je = i_startidx, i_endidx
#endif

      div_vec_e(je,jk,jb) =  &
        &    vec_e(iqidx(je,jb,1),jk,iqblk(je,jb,1))*ptr_int%geofac_qdiv(je,1,jb) &
        &  + vec_e(iqidx(je,jb,2),jk,iqblk(je,jb,2))*ptr_int%geofac_qdiv(je,2,jb) &
        &  + vec_e(iqidx(je,jb,3),jk,iqblk(je,jb,3))*ptr_int%geofac_qdiv(je,3,jb) &
        &  + vec_e(iqidx(je,jb,4),jk,iqblk(je,jb,4))*ptr_int%geofac_qdiv(je,4,jb)

    END DO
  END DO

END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE div_quad_twoadjcells

!-------------------------------------------------------------------------
!
!>
!! Computes discrete rotation.
!!
!! Computes discrete rotation at
!! (i) vertices of triangle cells (centers of dual grid cells) --or--
!! (ii) edges of hexagonal cells (centers of dual grid cells)
!! from a vector field given by its components in the directions normal
!! to triangle edges.
!! input:  lives on edges (velocity points)
!! output: lives on dual of cells (vertices for triangular grid, edges for
!!         hexagonal grid)
!!
!! @par Revision History
!! Developed and tested  by L.Bonaventura  (2002-4).
!! Adapted to new grid and patch structure by P. Korn
!! and L. Bonaventura (2005).
!! Modifications by P. Korn, MPI-M(2007-2)
!! - Switch fom array arguments to pointers
!! Modifications by Almut Gassmann, MPI-M (2007-04-20)
!! - abandon grid for the sake of patch
!! Modification by Almut Gassmann, MPI-M (2009-12-17)
!! - vorticity of hexagonal grid lives on rhombi
!!
SUBROUTINE rot_vertex_atmos( vec_e, ptr_patch, ptr_int, rot_vec, &
  &                          opt_slev, opt_elev, opt_rlstart, opt_rlend )
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  edge based variable of which rotation is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  vertex based variable in which rotation is stored
!
REAL(wp), INTENT(inout) ::  &
  &  rot_vec(:,:,:) ! dim: (nproma,nlev,nblks_v or nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jv, jk, jb

INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: nlen, npromz_v, nblks_v

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = nlev
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_math_operators:rot_vertex_atmos',  &
          &      'opt_rlstart must not be equal to 1')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF

!
!  loop through over all patch vertices (and blocks)
!
! The special treatment of 2D fields is essential for efficiency on the NEC

SELECT CASE (i_cell_type)

CASE (3)

  iidx => ptr_patch%verts%edge_idx
  iblk => ptr_patch%verts%edge_blk

  ! values for the blocking
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%verts%start_blk(rl_start,1)
  i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk), SCHEDULE(runtime)
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

    !
    ! compute the discrete rotation for vertex jv by
    ! finite volume approximation
    ! (see Bonaventura and Ringler MWR 2005);
    ! multiplication of the vector component vec_e by
    ! the appropriate dual cell based verts%edge_orientation
    ! is required to obtain the correct value for the
    ! application of Stokes theorem (which requires the scalar
    ! product of the vector field with the tangent unit vectors
    ! going around dual cell jv COUNTERCLOKWISE;
    ! since the positive direction for the vec_e components is
    ! not necessarily the one yelding counterclockwise rotation
    ! around dual cell jv, a correction coefficient (equal to +-1)
    ! is necessary, given by g%verts%edge_orientation
    !

#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
#endif
        !
        ! calculate rotation, i.e.
        ! add individual edge contributions to rotation
        !

        rot_vec(jv,jk,jb) =   &
          vec_e(iidx(jv,jb,1),jk,iblk(jv,jb,1)) * ptr_int%geofac_rot(jv,1,jb) + &
          vec_e(iidx(jv,jb,2),jk,iblk(jv,jb,2)) * ptr_int%geofac_rot(jv,2,jb) + &
          vec_e(iidx(jv,jb,3),jk,iblk(jv,jb,3)) * ptr_int%geofac_rot(jv,3,jb) + &
          vec_e(iidx(jv,jb,4),jk,iblk(jv,jb,4)) * ptr_int%geofac_rot(jv,4,jb) + &
          vec_e(iidx(jv,jb,5),jk,iblk(jv,jb,5)) * ptr_int%geofac_rot(jv,5,jb) + &
          vec_e(iidx(jv,jb,6),jk,iblk(jv,jb,6)) * ptr_int%geofac_rot(jv,6,jb)

      END DO

    END DO

  ENDDO
!$OMP END DO
!$OMP END PARALLEL

CASE (6) ! (i_cell_type == 6)

  iidx => ptr_patch%verts%edge_idx
  iblk => ptr_patch%verts%edge_blk

  ! values for the blocking
  ! no grid refinement in hexagonal model
  nblks_v   = ptr_patch%nblks_int_v
  npromz_v  = ptr_patch%npromz_int_v

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jv,jk)
  DO jb = 1, nblks_v
    IF (jb /= nblks_v) THEN
      nlen = nproma
    ELSE
      nlen = npromz_v
    ENDIF
    !
    ! Compute the discrete rotation for a triangle by
    ! application of Stokes theorem (which requires the scalar
    ! product of the vector field with the tangent unit vectors
    ! going around dual cell jv counterclockwise)
    !
#ifdef __LOOP_EXCHANGE
    DO jv = 1, nlen
      DO jk = slev, elev
#else
    DO jk = slev, elev
      DO jv = 1, nlen
#endif
        !
        ! add individual edge contributions to rotation
        !
        rot_vec(jv,jk,jb) =   &
          vec_e(iidx(jv,jb,1),jk,iblk(jv,jb,1)) * ptr_int%geofac_rot(jv,1,jb) + &
          vec_e(iidx(jv,jb,2),jk,iblk(jv,jb,2)) * ptr_int%geofac_rot(jv,2,jb) + &
          vec_e(iidx(jv,jb,3),jk,iblk(jv,jb,3)) * ptr_int%geofac_rot(jv,3,jb)

      END DO

    END DO

  ENDDO
!$OMP END DO
!$OMP END PARALLEL
END SELECT

END SUBROUTINE rot_vertex_atmos

END MODULE mo_atmo_math_divrot
