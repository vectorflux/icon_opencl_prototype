!>
!! Contains the implementation of interpolation and reconstruction.
!!
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
MODULE mo_intp_rbf
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                ONLY: wp
USE mo_impl_constants,      ONLY: min_rlvert
USE mo_model_domain,        ONLY: t_patch
USE mo_global_variables,    ONLY: nlev
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v

USE mo_intp_data_strc

IMPLICIT NONE

PRIVATE

PUBLIC :: rbf_vec_interpol_cell, rbf_interpol_c2grad,    &
        & rbf_vec_interpol_vertex, rbf_vec_interpol_edge

CONTAINS
!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
!
!
!>
!! Performs vector RBF reconstruction at cell center.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each cell center.
!!
!! @par Revision History
!! Developed  by L.Bonaventura  (2004).
!! Enhanced efficiency with  direct solve for coefficients of data,
!! introduced by Will Sawyer (2004).
!! Adapted to new data structure by L.Bonaventura and P.Korn (2006).
!! This a new version of the subroutine with the cartesian
!!   components of the reconstructed
!!   vector as output (Pilar Ripodas November 2006)
!! Modifications by Tobias Ruppert, DWD (2007-02-08):
!! - replaced rbf_vec_dim by rbf_mat_dim
!! Modifications by Almut Gassmann, MPI-M (2007-04-30)
!! - abandon grid for the sake of patch
!! Modification by Pilar Ripodas, DWD (2007-07):
!! - substruct the outgoing component of the reconstructed
!!   vector
!! Modification by Guenther Zaengl, DWD (2009-02-13)
!! - change to direct reconstruction of vector components
!!
SUBROUTINE rbf_vec_interpol_cell( p_vn_in, ptr_patch, ptr_int, p_u_out,  &
  &                               p_v_out, opt_slev, opt_elev, opt_rlstart )

! !INPUT PARAMETERS
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input normal components of (velocity) vectors at edge midpoints
REAL(wp), INTENT(IN) ::  &
  &  p_vn_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start value of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart

! reconstructed x-component (u) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_u_out(:,:,:) ! dim: (nproma,nlev,nblks_c)

! reconstructed y-component (v) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_v_out(:,:,:) ! dim: (nproma,nlev,nblks_c)

! !LOCAL VARIABLES
INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: nblks_c
INTEGER :: jc, jk, jb                ! integer over cells, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: i_rcstartlev    ! refinement control start level


INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

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
  i_rcstartlev = opt_rlstart
ELSE
  i_rcstartlev = 2
END IF


iidx => ptr_int%rbf_vec_idx_c
iblk => ptr_int%rbf_vec_blk_c

ptr_coeff => ptr_int%rbf_vec_coeff_c

! values for the blocking
nblks_c  = ptr_patch%nblks_int_c

! The start block depends on the width of the stencil
i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), SCHEDULE(runtime)

DO jb = i_startblk, nblks_c

  CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                     i_startidx, i_endidx, i_rcstartlev)

#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = slev, elev
#else
!CDIR UNROLL=2
  DO jk = slev, elev
    DO jc = i_startidx, i_endidx
#endif

      p_u_out(jc,jk,jb) =  &
        ptr_coeff(1,1,jc,jb)*p_vn_in(iidx(1,jc,jb),jk,iblk(1,jc,jb)) + &
        ptr_coeff(2,1,jc,jb)*p_vn_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,1,jc,jb)*p_vn_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,1,jc,jb)*p_vn_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,1,jc,jb)*p_vn_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,1,jc,jb)*p_vn_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,1,jc,jb)*p_vn_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,1,jc,jb)*p_vn_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,1,jc,jb)*p_vn_in(iidx(9,jc,jb),jk,iblk(9,jc,jb))
      p_v_out(jc,jk,jb) =  &
        ptr_coeff(1,2,jc,jb)*p_vn_in(iidx(1,jc,jb),jk,iblk(1,jc,jb)) + &
        ptr_coeff(2,2,jc,jb)*p_vn_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,2,jc,jb)*p_vn_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,2,jc,jb)*p_vn_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,2,jc,jb)*p_vn_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,2,jc,jb)*p_vn_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,2,jc,jb)*p_vn_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,2,jc,jb)*p_vn_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,2,jc,jb)*p_vn_in(iidx(9,jc,jb),jk,iblk(9,jc,jb))

    ENDDO
  ENDDO

ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE rbf_vec_interpol_cell

SUBROUTINE rbf_interpol_c2grad( p_cell_in, ptr_patch, ptr_int, grad_x,  &
  &                             grad_y, opt_slev, opt_elev, opt_rlstart )

! !INPUT PARAMETERS
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input cell-based variable for which gradient at cell center is computed
REAL(wp), INTENT(IN) ::  &
  &  p_cell_in(:,:,:) ! dim: (nproma,nlev,nblks_c)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start value of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart

! reconstructed zonal (x) component of gradient vector
REAL(wp),INTENT(INOUT) ::  &
  &  grad_x(:,:,:) ! dim: (nproma,nlev,nblks_c)

! reconstructed zonal (x) component of gradient vector
REAL(wp),INTENT(INOUT) ::  &
  &  grad_y(:,:,:) ! dim: (nproma,nlev,nblks_c)

! !LOCAL VARIABLES
INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: nblks_c
INTEGER :: jc, jk, jb                ! integer over cells, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: i_rcstartlev    ! refinement control start level


INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

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
  i_rcstartlev = opt_rlstart
ELSE
  i_rcstartlev = 2
END IF


iidx => ptr_int%rbf_c2grad_idx
iblk => ptr_int%rbf_c2grad_blk

ptr_coeff => ptr_int%rbf_c2grad_coeff

! values for the blocking
nblks_c  = ptr_patch%nblks_int_c

! The start block depends on the width of the stencil
i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), SCHEDULE(runtime)

DO jb = i_startblk, nblks_c

  CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                     i_startidx, i_endidx, i_rcstartlev)

#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = slev, elev
#else
!CDIR UNROLL=3
  DO jk = slev, elev
    DO jc = i_startidx, i_endidx
#endif

      grad_x(jc,jk,jb) =  &
        ptr_coeff(1,1,jc,jb)*p_cell_in(jc,jk,jb) + &
        ptr_coeff(2,1,jc,jb)*p_cell_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,1,jc,jb)*p_cell_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,1,jc,jb)*p_cell_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,1,jc,jb)*p_cell_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,1,jc,jb)*p_cell_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,1,jc,jb)*p_cell_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,1,jc,jb)*p_cell_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,1,jc,jb)*p_cell_in(iidx(9,jc,jb),jk,iblk(9,jc,jb)) + &
        ptr_coeff(10,1,jc,jb)*p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))
      grad_y(jc,jk,jb) =  &
        ptr_coeff(1,2,jc,jb)*p_cell_in(jc,jk,jb) + &
        ptr_coeff(2,2,jc,jb)*p_cell_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
        ptr_coeff(3,2,jc,jb)*p_cell_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
        ptr_coeff(4,2,jc,jb)*p_cell_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
        ptr_coeff(5,2,jc,jb)*p_cell_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
        ptr_coeff(6,2,jc,jb)*p_cell_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
        ptr_coeff(7,2,jc,jb)*p_cell_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
        ptr_coeff(8,2,jc,jb)*p_cell_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
        ptr_coeff(9,2,jc,jb)*p_cell_in(iidx(9,jc,jb),jk,iblk(9,jc,jb)) + &
        ptr_coeff(10,2,jc,jb)*p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))

    ENDDO
  ENDDO

ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE rbf_interpol_c2grad

!-------------------------------------------------------------------------
!
!
!>
!! Performs vector RBF reconstruction at triangle vertices.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each vertex.
!!
!! @par Revision History
!! Developed  by Jochen Foerstner, DWD (2008-05-02)
!! Modification by Guenther Zaengl, DWD (2009-02-13)
!! - change to direct reconstruction of vector components
!! Modification by Almut Gassmann, MPI-M (2009-11-06)
!! - include switch which distinguishes normal or tangential components as input
!!
SUBROUTINE rbf_vec_interpol_vertex( p_e_in, ptr_patch, ptr_int, &
                                    p_u_out, p_v_out,           &
                                    opt_slev, opt_elev, opt_rlstart, opt_rlend )
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input components of velocity or horizontal vorticity vectors at edge midpoints
REAL(wp), INTENT(IN) ::  &
  &  p_e_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend

! reconstructed x-component (u) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_u_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

! reconstructed y-component (v) of velocity vector
REAL(wp),INTENT(INOUT) ::  &
  &  p_v_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

! !LOCAL VARIABLES

INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: jv, jk, jb                ! integer over vertices, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_endblk        ! end block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: i_nchdom        ! number of child domains
INTEGER :: i_rcstartlev    ! refinement control start level
INTEGER :: i_rcendlev      ! refinement control end level

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

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
  i_rcstartlev = opt_rlstart
ELSE
  i_rcstartlev = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  i_rcendlev = opt_rlend
ELSE
  i_rcendlev = min_rlvert
END IF

iidx => ptr_int%rbf_vec_idx_v
iblk => ptr_int%rbf_vec_blk_v

ptr_coeff => ptr_int%rbf_vec_coeff_v

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(i_rcstartlev,1)
i_endblk   = ptr_patch%verts%end_blk(i_rcendlev,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jv), SCHEDULE(runtime)
DO jb = i_startblk, i_endblk

  CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, i_rcstartlev, i_rcendlev)

#ifdef __LOOP_EXCHANGE
  DO jv = i_startidx, i_endidx
    DO jk = slev, elev
#else
!CDIR UNROLL=6
  DO jk = slev, elev
    DO jv = i_startidx, i_endidx
#endif

      p_u_out(jv,jk,jb) =  &
        ptr_coeff(1,1,jv,jb)*p_e_in(iidx(1,jv,jb),jk,iblk(1,jv,jb)) + &
        ptr_coeff(2,1,jv,jb)*p_e_in(iidx(2,jv,jb),jk,iblk(2,jv,jb)) + &
        ptr_coeff(3,1,jv,jb)*p_e_in(iidx(3,jv,jb),jk,iblk(3,jv,jb)) + &
        ptr_coeff(4,1,jv,jb)*p_e_in(iidx(4,jv,jb),jk,iblk(4,jv,jb)) + &
        ptr_coeff(5,1,jv,jb)*p_e_in(iidx(5,jv,jb),jk,iblk(5,jv,jb)) + &
        ptr_coeff(6,1,jv,jb)*p_e_in(iidx(6,jv,jb),jk,iblk(6,jv,jb))
      p_v_out(jv,jk,jb) =  &
        ptr_coeff(1,2,jv,jb)*p_e_in(iidx(1,jv,jb),jk,iblk(1,jv,jb)) + &
        ptr_coeff(2,2,jv,jb)*p_e_in(iidx(2,jv,jb),jk,iblk(2,jv,jb)) + &
        ptr_coeff(3,2,jv,jb)*p_e_in(iidx(3,jv,jb),jk,iblk(3,jv,jb)) + &
        ptr_coeff(4,2,jv,jb)*p_e_in(iidx(4,jv,jb),jk,iblk(4,jv,jb)) + &
        ptr_coeff(5,2,jv,jb)*p_e_in(iidx(5,jv,jb),jk,iblk(5,jv,jb)) + &
        ptr_coeff(6,2,jv,jb)*p_e_in(iidx(6,jv,jb),jk,iblk(6,jv,jb))

      ENDDO
    ENDDO

ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE rbf_vec_interpol_vertex

!-------------------------------------------------------------------------
!
!
!>
!! Performs vector RBF reconstruction at edge midpionts.
!!
!! Theory described in Narcowich and Ward (Math Comp. 1994) and
!! Bonaventura and Baudisch (Mox Report n. 75).
!! It takes edge based variables as input and combines them
!! into three dimensional cartesian vectors at each edge.
!!
!! @par Revision History
!! Developed  by Jochen Foerstner, DWD (2008-07-15)
!! Modification by Guenther Zaengl, DWD (2009-02-13)
!! - change to direct reconstruction of tangential vector component
!!
SUBROUTINE rbf_vec_interpol_edge( p_vn_in, ptr_patch, ptr_int, p_vt_out,  &
  &                               opt_slev, opt_elev, opt_rlstart )
!
TYPE(t_patch), TARGET, INTENT(in) ::  &
  &  ptr_patch

! input normal components of velocity vectors at edge midpoints
REAL(wp), INTENT(IN) ::  &
  &  p_vn_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(IN)  ::  ptr_int

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

! start value of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL :: opt_rlstart

! reconstructed tangential velocity component
REAL(wp),INTENT(INOUT) ::  &
  &  p_vt_out(:,:,:) ! dim: (nproma,nlev,nblks_e)


INTEGER :: slev, elev                ! vertical start and end level
INTEGER :: nblks_e
INTEGER :: je, jk, jb                ! integer over edges, levels, and blocks,

INTEGER :: i_startblk      ! start block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: i_rcstartlev    ! refinement control start level

INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
REAL(wp), DIMENSION(:,:,:), POINTER :: ptr_coeff

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

iidx => ptr_int%rbf_vec_idx_e
iblk => ptr_int%rbf_vec_blk_e

ptr_coeff => ptr_int%rbf_vec_coeff_e
IF ( PRESENT(opt_rlstart) ) THEN
  i_rcstartlev = opt_rlstart
ELSE
  i_rcstartlev = 2
END IF

! values for the blocking
nblks_e  = ptr_patch%nblks_int_e

i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je)
  DO jb = i_startblk, nblks_e

    CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, i_rcstartlev)

#ifdef __LOOP_EXCHANGE
    DO je = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=3
    DO jk = slev, elev
      DO je = i_startidx, i_endidx
#endif

        p_vt_out(je,jk,jb) =  &
          ptr_coeff(1,je,jb)*p_vn_in(iidx(1,je,jb),jk,iblk(1,je,jb)) + &
          ptr_coeff(2,je,jb)*p_vn_in(iidx(2,je,jb),jk,iblk(2,je,jb)) + &
          ptr_coeff(3,je,jb)*p_vn_in(iidx(3,je,jb),jk,iblk(3,je,jb)) + &
          ptr_coeff(4,je,jb)*p_vn_in(iidx(4,je,jb),jk,iblk(4,je,jb))

      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE rbf_vec_interpol_edge



END MODULE mo_intp_rbf
