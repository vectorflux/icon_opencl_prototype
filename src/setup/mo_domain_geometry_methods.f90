!>
!! @par Revision History
!!!
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
MODULE mo_domain_geometry_methods
!-------------------------------------------------------------------------

USE mo_kind,               ONLY: wp
USE mo_impl_constants,      ONLY: SUCCESS, &
     &                           MAX_CHAR_LENGTH,  &
     &                           min_rlcell, max_rlcell, &
     &                           min_rledge, max_rledge, &
     &                           min_rlvert, max_rlvert
USE mo_exception
USE mo_model_domain
USE mo_physical_constants, ONLY: omega
USE mo_global_variables,   ONLY: nproma, i_cell_type
USE mo_math_utilities,     ONLY: gvec2cvec, t_cartesian_coordinates
USE mo_math_constants,     ONLY: rad2deg, pi_2
USE mo_loopindices,        ONLY: get_indices_e
USE mo_mpi,                ONLY: p_pe, p_io
USE mo_reshape_arrays

IMPLICIT NONE

PRIVATE
#if defined(NOMPI) && !defined(STANDALONE)
INCLUDE 'netcdf.inc'
#endif

CHARACTER(len=*), PARAMETER :: version = '$Id$'
!subroutines
PUBLIC :: calculate_cart_normal, &
        & init_quad_twoadjcells, lcorio, init_coriolis 
LOGICAL :: lcorio = .true.
!-------------------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------------
!>
!! This routine calculates the Cartesian components of the edge primal normals.
!!
!! Therefore the given zonal and meridional
!! components of the edge primal normals are used.
!! This information should be provided by the grid generator and
!! read in the routine read_patch.
!!
SUBROUTINE calculate_cart_normal( p_patch )
!
TYPE(t_patch), TARGET, INTENT(inout) :: p_patch  ! patch on a specific level

TYPE(t_cartesian_coordinates) :: z_vec
REAL(wp) :: z_lon, z_lat, z_u, z_v  ! location and components of normal
REAL(wp) :: z_norm                  ! norm of Cartesian normal

INTEGER :: nlen, nblks_e, npromz_e
INTEGER :: jb, je                   ! loop indices

!-----------------------------------------------------------------------

! values for the blocking
nblks_e  = p_patch%nblks_e
npromz_e = p_patch%npromz_e


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,nlen,z_lon,z_lat,z_u,z_v,z_norm,z_vec)
DO jb = 1, nblks_e

  IF (jb /= nblks_e) THEN
    nlen = nproma
  ELSE
    nlen = npromz_e
  ENDIF

  ! loop over edges
  DO je = 1, nlen

    ! location of edge midpoint
      z_lon = p_patch%edges%center(je,jb)%lon
      z_lat = p_patch%edges%center(je,jb)%lat

    ! zonal and meridional component of primal normal
    z_u = p_patch%edges%primal_normal(je,jb)%v1
    z_v = p_patch%edges%primal_normal(je,jb)%v2

    ! calculate Cartesian components of primal normal
    CALL gvec2cvec( z_u, z_v, z_lon, z_lat, z_vec%x(1), z_vec%x(2), z_vec%x(3) )

    ! compute unit normal to edge je
    z_norm = SQRT( DOT_PRODUCT(z_vec%x(1:3),z_vec%x(1:3)) )
    z_vec%x(1:3) = 1._wp / z_norm * z_vec%x(1:3)

    ! save the values in the according type structure of the patch
    p_patch%edges%primal_cart_normal(je,jb)%x(1) = z_vec%x(1)
    p_patch%edges%primal_cart_normal(je,jb)%x(2) = z_vec%x(2)
    p_patch%edges%primal_cart_normal(je,jb)%x(3) = z_vec%x(3)

    ! zonal and meridional component of dual normal
    z_u = p_patch%edges%dual_normal(je,jb)%v1
    z_v = p_patch%edges%dual_normal(je,jb)%v2

    ! calculate Cartesian components of primal normal
    CALL gvec2cvec( z_u, z_v, z_lon, z_lat, z_vec%x(1), z_vec%x(2), z_vec%x(3) )

    ! compute unit normal to edge je
    z_norm = SQRT( DOT_PRODUCT(z_vec%x(1:3),z_vec%x(1:3)) )
    z_vec%x(1:3) = 1._wp / z_norm * z_vec%x(1:3)

    ! save the values in the according type structure of the patch
    p_patch%edges%dual_cart_normal(je,jb)%x(1) = z_vec%x(1)
    p_patch%edges%dual_cart_normal(je,jb)%x(2) = z_vec%x(2)
    p_patch%edges%dual_cart_normal(je,jb)%x(3) = z_vec%x(3)

  END DO

END DO
!$OMP END DO
!$OMP END PARALLEL


END SUBROUTINE calculate_cart_normal

!-------------------------------------------------------------------------
!>
!! This routine initializes the data for the quadrilateral cells
!! formed by the two adjacent cells of an edge.
!! The four edge indices of this qudrilateral and its area are stored
!! in the derived type for the edges.
!! This information should be provided by the grid generator and
!! read in the routine read_patch.
!!
!! Also suitable for hexagons: Almut Gassmann (2009-10-01)
!!
SUBROUTINE init_quad_twoadjcells( p_patch )
!
TYPE(t_patch), TARGET, INTENT(inout) :: p_patch  ! patch on a specific level

INTEGER :: nblks_e, npromz_e
INTEGER :: je, jb                   ! loop index
INTEGER :: iie, iie1(nproma), iie2(nproma), ierror(nproma)
INTEGER :: ilc1, ibc1, ilc2, ibc2                   ! cell indices
INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3       ! edge indices
INTEGER :: ilv1, ibv1, ilv2, ibv2, ilev, ibev, jjev ! vertex indices
INTEGER :: i_startblk, i_startidx, i_endidx


!-----------------------------------------------------------------------

! values for the blocking
nblks_e  = p_patch%nblks_e
npromz_e = p_patch%npromz_e

! Quad cells cannot be computed along the lateral edges of a nested domain
i_startblk = p_patch%edges%start_blk(2,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,iie,ilc1,ibc1,ilc2,ibc2,ile1,ibe1,ile2,ibe2,&
!$OMP            ile3,ibe3,iie1,iie2,ierror,ilv1,ibv1,ilv2,ibv2,ilev,ibev,jjev)
DO jb = i_startblk, nblks_e

  CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
                     i_startidx, i_endidx, 2)

  IF (i_cell_type == 3) THEN ! vectorized code for triangular grid

    ierror = 0

    DO je = i_startidx, i_endidx

      !
      ! get global indices of the edges of the two neighboring cells
      !
      ilc1 = p_patch%edges%cell_idx(je,jb,1)
      ibc1 = p_patch%edges%cell_blk(je,jb,1)
      ilc2 = p_patch%edges%cell_idx(je,jb,2)
      ibc2 = p_patch%edges%cell_blk(je,jb,2)

      ! sum area of the two adjacent cells
      p_patch%edges%quad_area(je,jb) = &
        p_patch%cells%area(ilc1,ibc1) + p_patch%cells%area(ilc2,ibc2)

      ! Indices of edges adjacent to cell 1
      ile1 = p_patch%cells%edge_idx(ilc1,ibc1,1)
      ibe1 = p_patch%cells%edge_blk(ilc1,ibc1,1)
      ile2 = p_patch%cells%edge_idx(ilc1,ibc1,2)
      ibe2 = p_patch%cells%edge_blk(ilc1,ibc1,2)
      ile3 = p_patch%cells%edge_idx(ilc1,ibc1,3)
      ibe3 = p_patch%cells%edge_blk(ilc1,ibc1,3)

      IF (je == ile1 .AND. jb == ibe1) THEN
        iie1(je) = 2
        iie2(je) = 3
      ELSE IF (je == ile2 .AND. jb == ibe2) THEN
        iie1(je) = 3
        iie2(je) = 1
      ELSE IF (je == ile3 .AND. jb == ibe3) THEN
        iie1(je) = 1
        iie2(je) = 2
      ELSE
        ierror(je) = ierror(je) + 1
      ENDIF

      p_patch%edges%quad_idx(je,jb,1) = p_patch%cells%edge_idx(ilc1,ibc1,iie1(je))
      p_patch%edges%quad_blk(je,jb,1) = p_patch%cells%edge_blk(ilc1,ibc1,iie1(je))
      p_patch%edges%quad_orientation(je,jb,1) =  &
        &  p_patch%cells%edge_orientation(ilc1,ibc1,iie1(je))

      p_patch%edges%quad_idx(je,jb,2) = p_patch%cells%edge_idx(ilc1,ibc1,iie2(je))
      p_patch%edges%quad_blk(je,jb,2) = p_patch%cells%edge_blk(ilc1,ibc1,iie2(je))
      p_patch%edges%quad_orientation(je,jb,2) =  &
        &  p_patch%cells%edge_orientation(ilc1,ibc1,iie2(je))

      ! Indices of edges adjacent to cell 2
      ile1 = p_patch%cells%edge_idx(ilc2,ibc2,1)
      ibe1 = p_patch%cells%edge_blk(ilc2,ibc2,1)
      ile2 = p_patch%cells%edge_idx(ilc2,ibc2,2)
      ibe2 = p_patch%cells%edge_blk(ilc2,ibc2,2)
      ile3 = p_patch%cells%edge_idx(ilc2,ibc2,3)
      ibe3 = p_patch%cells%edge_blk(ilc2,ibc2,3)

      IF (je == ile1 .AND. jb == ibe1) THEN
        iie1(je) = 2
        iie2(je) = 3
      ELSE IF (je == ile2 .AND. jb == ibe2) THEN
        iie1(je) = 3
        iie2(je) = 1
      ELSE IF (je == ile3 .AND. jb == ibe3) THEN
        iie1(je) = 1
        iie2(je) = 2
      ELSE
        ierror(je) = ierror(je) + 1
      ENDIF

      p_patch%edges%quad_idx(je,jb,3) = p_patch%cells%edge_idx(ilc2,ibc2,iie1(je))
      p_patch%edges%quad_blk(je,jb,3) = p_patch%cells%edge_blk(ilc2,ibc2,iie1(je))
      p_patch%edges%quad_orientation(je,jb,3) =  &
        &  p_patch%cells%edge_orientation(ilc2,ibc2,iie1(je))

      p_patch%edges%quad_idx(je,jb,4) = p_patch%cells%edge_idx(ilc2,ibc2,iie2(je))
      p_patch%edges%quad_blk(je,jb,4) = p_patch%cells%edge_blk(ilc2,ibc2,iie2(je))
      p_patch%edges%quad_orientation(je,jb,4) =  &
        &  p_patch%cells%edge_orientation(ilc2,ibc2,iie2(je))

    ENDDO
    IF (MAXVAL(ierror) > 0) THEN
      WRITE(0,*) "Number of errors", SUM(ierror(1:nproma))
      CALL finish ('mo_model_domain_import:init_quad_twoadjcells',  &
        &          'edge-cell index relationships are apparently incorrect')
    ENDIF

  ELSE ! i_cell_type == 6

    DO je = i_startidx, i_endidx

      iie = 0
      !
      ! get global indices of the edges of the two neighboring verts
      !
      ilv1 = p_patch%edges%vertex_idx(je,jb,1)
      ibv1 = p_patch%edges%vertex_blk(je,jb,1)
      ilv2 = p_patch%edges%vertex_idx(je,jb,2)
      ibv2 = p_patch%edges%vertex_blk(je,jb,2)

      ! sum area of the two adjacent verts
      p_patch%edges%quad_area(je,jb) = &
        p_patch%verts%dual_area(ilv1,ibv1) + p_patch%verts%dual_area(ilv2,ibv2)

      DO jjev = 1, 3

        ilev = p_patch%verts%edge_idx(ilv1,ibv1,jjev)
        ibev = p_patch%verts%edge_blk(ilv1,ibv1,jjev)

        ! test if edge is not the current one, then store line and
        ! block indices and orientation of quad edge
        IF ( ilev /= je .OR. ibev /= jb ) THEN
          iie = iie+1
          p_patch%edges%quad_idx(je,jb,iie) = ilev
          p_patch%edges%quad_blk(je,jb,iie) = ibev
          p_patch%edges%quad_orientation(je,jb,iie) =  &
            &  p_patch%verts%edge_orientation(ilv1,ibv1,jjev)
        ENDIF

      ENDDO

      DO jjev = 1, 3

        ilev = p_patch%verts%edge_idx(ilv2,ibv2,jjev)
        ibev = p_patch%verts%edge_blk(ilv2,ibv2,jjev)

        ! test if edge is not the current one, then store line and
        ! block indices and orientation of quad edge
        IF ( ilev /= je .OR. ibev /= jb ) THEN
          iie = iie+1
          p_patch%edges%quad_idx(je,jb,iie) = ilev
          p_patch%edges%quad_blk(je,jb,iie) = ibev
          p_patch%edges%quad_orientation(je,jb,iie) =  &
            &  p_patch%verts%edge_orientation(ilv2,ibv2,jjev)
        ENDIF

      ENDDO

      IF (iie /= 4) THEN
        PRINT *, "ERROR ==>  iie = ", iie,  " /= 4 "
        CALL finish ('mo_model_domain_import:init_quad_twoadjcells',  &
          &             'wrong number of edge indices for quad')
      ENDIF

    END DO ! end edge loop

  ENDIF ! i_cell_type

END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE init_quad_twoadjcells

!-------------------------------------------------------------------------
!>
!! Initializes the Coriolis components of the grid with analytical values.
SUBROUTINE init_coriolis( p_patch )
!

TYPE(t_patch), TARGET, INTENT(inout) :: p_patch ! patch on specific level

INTEGER :: nlen,                         &
  &        nblks_c,  nblks_e,  nblks_v,  &
  &        npromz_c, npromz_e, npromz_v
INTEGER :: jc, je, jv, jb

REAL(wp) :: zlat

!-----------------------------------------------------------------------

! values for the blocking
nblks_c  = p_patch%nblks_c
npromz_c = p_patch%npromz_c
nblks_e  = p_patch%nblks_e
npromz_e = p_patch%npromz_e
nblks_v  = p_patch%nblks_v
npromz_v = p_patch%npromz_v


IF (lcorio ) THEN

  DO jb = 1, nblks_c
    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF
    DO jc = 1, nlen
      zlat = p_patch%cells%center(jc,jb)%lat
      p_patch%cells%f_c(jc,jb) = 2._wp*omega*SIN(zlat)
    END DO
  END DO

  DO jb = 1, nblks_e
    IF (jb /= nblks_e) THEN
      nlen = nproma
    ELSE
      nlen = npromz_e
    ENDIF
    DO je = 1, nlen
      zlat = p_patch%edges%center(je,jb)%lat
      p_patch%edges%f_e(je,jb) = 2._wp*omega*SIN(zlat)
    END DO
  END DO

  DO jb = 1, nblks_v
    IF (jb /= nblks_v) THEN
      nlen = nproma
    ELSE
      nlen = npromz_v
    ENDIF
    DO jv = 1, nlen
      zlat = p_patch%verts%vertex(jv,jb)%lat
      p_patch%verts%f_v(jv,jb) = 2._wp*omega*SIN(zlat)
    END DO
  END DO


ELSE

  p_patch%cells%f_c(:,:) = 0.0_wp
  p_patch%edges%f_e(:,:) = 0.0_wp
  p_patch%verts%f_v(:,:) = 0.0_wp

ENDIF

END SUBROUTINE init_coriolis


END MODULE mo_domain_geometry_methods
