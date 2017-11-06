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
MODULE mo_oce_math_divrot
!-------------------------------------------------------------------------
!
USE mo_kind,               ONLY: wp
USE mo_impl_constants,     ONLY: SEA, min_rlcell, min_rledge, min_rlvert
USE mo_interpolation,      ONLY: t_int_state
USE mo_model_domain,       ONLY: t_patch
USE mo_global_variables,   ONLY: nproma, i_cell_type, nlev !, ltimer
USE mo_exception,          ONLY: finish
! USE mo_timer,              ONLY: timer_start, timer_stop, timer_div
USE mo_loopindices,        ONLY: get_indices_c, get_indices_e, get_indices_v

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'


PUBLIC :: div_oce

CONTAINS




!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field in presence of lateral boundaries as in ocean setting.
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
!!Boundary treatment by Peter Korn (2009)
!!
SUBROUTINE div_oce( vec_e, ptr_patch, ptr_int, div_vec_c, &
  &                 opt_slev, opt_elev, opt_rlstart, opt_rlend )
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
! Attention - vec_e is changed at boundaries - better use scratch array #slo2010/04
REAL(wp), INTENT(inout) ::  &
!REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

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

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

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
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF

iidx => ptr_patch%cells%edge_idx
iblk => ptr_patch%cells%edge_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

! loop through all patch cells (and blocks)
!

! IF(ltimer) CALL timer_start(timer_div)

!$OMP PARALLEL

! The special treatment of 2D fields is essential for efficiency on the NEC

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __SX__
!CDIR UNROLL=6
#endif
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
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
        !
        !Distinghuish: case of a land cell (put div to zero), and
        !cases where one of the triangle edges are boundary or land
        !(put corresponding velocity to zero).
        IF(ptr_patch%patch_oce%lsm_oce_c(jc,jk,jb)/=SEA)THEN
          div_vec_c(jc,jk,jb) = 0.0_wp
        ELSE
          IF(ptr_patch%patch_oce%lsm_oce_e(iidx(jc,jb,1),jk,iblk(jc,jb,1))/=SEA)THEN

            vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1))=0.0_wp

          ELSEIF(ptr_patch%patch_oce%lsm_oce_e(iidx(jc,jb,2),jk,iblk(jc,jb,2))/=SEA)THEN

            vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2))=0.0_wp

          ELSEIF(ptr_patch%patch_oce%lsm_oce_e(iidx(jc,jb,3),jk,iblk(jc,jb,3))/=SEA)THEN

            vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3))=0.0_wp

          ENDIF

          div_vec_c(jc,jk,jb) =  &
            vec_e(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            vec_e(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            vec_e(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)

        ENDIF
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

! IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div_oce



END MODULE mo_oce_math_divrot
