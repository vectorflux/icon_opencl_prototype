!>
!! Contains the implementation of interpolation and reconstruction.
!!
!! Contains the implementation of interpolation and reconstruction
!! routines used by the shallow water model, including the RBF
!! reconstruction routines.
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
MODULE mo_intp_data_strc
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------

USE mo_kind,                ONLY: wp
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert, max_dom
!DR for testing purposes
USE mo_math_utilities,      ONLY: t_geographical_coordinates

IMPLICIT NONE


! NOTE: The variables will be use along the mo_interpolation sub-modules
!       They are declared to be public
PUBLIC

INTEGER  :: rbf_vec_kern_c,   & ! parameter determining the type
            rbf_vec_kern_v,   & ! of vector rbf kernel
            rbf_vec_kern_e

INTEGER  :: i_cori_method       ! Identifier for the method with wich the tangential
                                ! wind reconstruction in Coriolis force is computed,
                                ! if the Thuburn method is used. (To be
                                ! implemented for triangles, currently only for
                                ! hexagons)
                                ! i_cori_method = 1 : Almut's method for reconstruction
                                !                     but TRSK method for PV
                                ! i_cori_method = 2 : Thuburn/Ringler/Skamarock/Klemp
                                ! i_cori_method = 3 : Almut's method for reconstruction
                                !                     Almut's method also for PV
REAL(wp) :: sick_a, sick_o      ! if i_cori_method = 3: To avoid the SICK instability
                                ! (Symmetric Instability of Computational Kind or
                                !  Hollingsworth instability), an average of the kinetic
                                ! energy and thus of the mass flux must be defined.
                                ! sick_a is to be given in the namelist as the weight of
                                ! the fully averaged kinetic energy, sick_o=1-sick_a.
LOGICAL :: l_corner_vort        ! yields for i_cori_method=3
                                ! Decision wheter the hexagon vector reconstruction is
                                ! combined with either of the two vorticities :
                                ! .TRUE. : Three rhombi are combined to the corner
                                !          and afterwards averaged to the hexagon center
                                ! .FALSE.: 6 rhombi are directly averaged to the
                                !          hexagon center (original method).  
                                ! After the writing of the paper to be published in JCP 
                                ! it seems that l_corner_vort=.TRUE. should be the right way.
                                ! Possibly the 'sick' handling must be retuned in that case.

INTEGER  :: rbf_vec_dim_c,    & ! parameter determining the size
            rbf_vec_dim_v,    & ! of vector rbf stencil
            rbf_vec_dim_e,    & !
            rbf_c2grad_dim      ! ... and for cell-to-gradient reconstruction

! Parameter fields determining the scale factor used by the vector rbf
! interpolator.
! Note: these fields are defined on each grid level; to allow the namelist input
! going from 1 to depth (rather than from start_lev to end_lev), the namelist input
! fields defined here differ from those used in the model
REAL(wp) :: rbf_vec_scale_c(max_dom),  &
            rbf_vec_scale_v(max_dom),  &
            rbf_vec_scale_e(max_dom)

! Namelist variables setting up the lateral boundary nudging (applicable to limited-area
! runs and one-way nesting). The nudging coefficients start with nudge_max_coeff in
! the cell row bordering to the boundary interpolation zone, and decay exponentially
! with nudge_efold_width (in units of cell rows)
REAL(wp) :: nudge_max_coeff, nudge_efold_width
INTEGER  :: nudge_zone_width     ! total width of nudging zone in units of cell rows


LOGICAL  :: llsq_high_consv     ! flag to determine whether the high order least 
                                ! squares reconstruction should be conservative

INTEGER  :: lsq_high_ord        ! specific order for higher order lsq


TYPE t_lsq_set
  LOGICAL :: l_consv             ! flag to determine whether the least squares
                                 ! reconstruction should be conservative

  INTEGER :: dim_c,        &     ! parameter determining the size of the lsq stencil
    &        dim_unk             ! parameter determining the dimension of the solution
                                 ! vector (== number of unknowns) of the lsq system
END TYPE t_lsq_set



TYPE(t_lsq_set) :: lsq_lin_set, &! Parameter settings for linear and higher order  
  &                lsq_high_set  ! least squares

TYPE t_lsq
  ! fields related to weighted least squares polynomial reconstruction
  !---------------------------------------------------------------------
  INTEGER, ALLOCATABLE  :: lsq_dim_stencil(:,:)  ! stencil size as a function of jc and jb
                                                 ! necessary in order to account for pentagon
                                                 ! points.
  INTEGER, ALLOCATABLE  :: lsq_idx_c(:,:,:)      ! index array defining the stencil for
                                                 ! lsq reconstruction (nproma,nblks_c,lsq_dim_c)
  INTEGER, ALLOCATABLE  :: lsq_blk_c(:,:,:)      ! dito for the blocks
  REAL(wp), ALLOCATABLE :: lsq_weights_c(:,:,:)  ! weights for lsq reconstruction
                                                 ! (nproma,lsq_dim_c,nblks_c)
  REAL(wp), ALLOCATABLE :: lsq_qtmat_c(:,:,:,:)  ! transposed Q of QR-factorization for
                                                 ! lsq reconstruction
                                                 ! (nproma,lsq_dim_unk,lsq_dim_c,nblks_c)
  REAL(wp), ALLOCATABLE :: lsq_rmat_rdiag_c(:,:,:)! reciprocal diagonal elements of R-matrix
                                                  ! resulting from QR-decomposition
                                                  ! (nproma,lsq_dim_unk,nblks_c)
  REAL(wp), ALLOCATABLE :: lsq_rmat_utri_c(:,:,:)! upper triangular elements without diagonal
                                                 ! elements of R-matrix (starting from the bottom
                                                 ! right)
                                                 ! (nproma,(lsq_dim_unk^2-lsq_dim_unk)/2,nblks_c)
  REAL(wp), ALLOCATABLE :: lsq_moments(:,:,:)      ! Moments (x^ny^m)_{i} for control volume
                                                   ! (nproma,nblks_c,lsq_dim_unk)
  REAL(wp), ALLOCATABLE :: lsq_moments_hat(:,:,:,:)! Moments (\hat{x^ny^m}_{ij}) for control volume
                                                   ! (nproma,nblks_c,lsq_dim_c,lsq_dim_unk)
END TYPE t_lsq


TYPE t_gauss_quad
  !
  ! quadrature points for intergration over triangular element
  ! 
  TYPE(t_geographical_coordinates), ALLOCATABLE :: & !< linear (nproma, nblks_c)  
    &  qpts_tri_l(:,:)     
  TYPE(t_geographical_coordinates), ALLOCATABLE :: & !< quadratic (nproma, nblks_c,3)
    &  qpts_tri_q(:,:,:)     
  TYPE(t_geographical_coordinates), ALLOCATABLE :: & !< cubic (nproma, nblks_c,4)
    &  qpts_tri_c(:,:,:)     
  REAL(wp), ALLOCATABLE :: weights_tri_q(:)      !< quadratic weight (3)
  REAL(wp), ALLOCATABLE :: weights_tri_c(:)      !< cubic weights (4)
END TYPE t_gauss_quad


TYPE t_int_state

  ! a) weights which are inconsistent with the Hamiltonian viewpoint
  !-----------------------------------------------------------------

  REAL(wp), ALLOCATABLE :: c_lin_e(:,:,:)   ! coefficient for interpolation
                                            ! from adjacent cells onto edge
                                            ! (nproma,2,nblks_e)

  REAL(wp), ALLOCATABLE :: e_bln_c_s(:,:,:) ! coefficient for bilinear
                                            ! interpolation from edges to cells
                                            ! for scalar quantities

  REAL(wp), ALLOCATABLE :: e_bln_c_u(:,:,:) ! coefficient for bilinear interpolation
                                            ! from edges to cells for vector components
                                            ! (input: v_t, v_n, output: u)

  REAL(wp), ALLOCATABLE :: e_bln_c_v(:,:,:) ! coefficient for bilinear interpolation
                                            ! from edges to cells for vector components
                                            ! (input: v_t, v_n, output: v)

  REAL(wp), ALLOCATABLE :: c_bln_avg(:,:,:) ! coefficients for bilinear divergence
                                            ! averaging (nproma,4,nblks_c)

  REAL(wp), ALLOCATABLE :: e_flx_avg(:,:,:) ! coefficients for related velocity or mass flux
                                            ! averaging (nproma,5,nblks_e)

  REAL(wp), ALLOCATABLE :: v_1o2_e(:,:,:)   ! coefficient for interpolation
                                            ! from vertices onto edges by 1/2
                                            ! weighting (nproma,2,nblks_e),

  ! b) weights which are consistent with the Hamiltonian viewpoint
  !---------------------------------------------------------------
  ! The following weights are needed for the mass and theta brackets

  REAL(wp), ALLOCATABLE :: e_inn_c(:,:,:)   ! coefficient for inner product
                                            ! of 2 vector components
                                            ! from edges to cells

  REAL(wp), ALLOCATABLE :: e_inn_v(:,:,:)   ! coefficient for inner product
                                            ! of 2 vector components
                                            ! from edges to verts

  REAL(wp), ALLOCATABLE :: e_aw_c(:,:,:)    ! coefficient for scalar interp
                                            ! from edges to cells

  REAL(wp), ALLOCATABLE :: r_aw_c(:,:,:)    ! coefficient for scalar interp
                                            ! from rhombi to cells

  REAL(wp), ALLOCATABLE :: e_aw_v(:,:,:)    ! coefficient for scalar interp
                                            ! from edges to vertices

  REAL(wp), ALLOCATABLE :: e_1o3_v(:,:,:)   ! coefficient for hexagonal grid
                                            ! averaging from edges to vertices

  REAL(wp), ALLOCATABLE :: tria_aw_rhom(:,:,:)! coefficient for interpolation
                                            ! from triangles to rhombi by area
                                            ! weighting

  REAL(wp), ALLOCATABLE :: verts_aw_cells(:,:,:)! coefficient for interpolation
                                            ! from vertices to cells by
                                            ! area weighting

  REAL(wp), ALLOCATABLE :: cells_aw_verts(:,:,:)! coefficient for interpolation
                                            ! from cells to verts by
                                            ! area weighting

  ! c) RBF related fields
  !----------------------
  INTEGER, ALLOCATABLE  :: rbf_vec_idx_c(:,:,:)  ! index array defining the
                                            ! stencil of surrounding edges for
                                            ! vector rbf interpolation at each
                                            ! cell center
                                            ! (rbf_vec_dim_c,nproma,nblks_c)
  INTEGER, ALLOCATABLE  :: rbf_vec_blk_c(:,:,:)  ! ... dito for the blocks

  INTEGER, ALLOCATABLE  :: rbf_vec_stencil_c(:,:) ! array defining number of
                                            ! surrounding edges in the stencil
                                            ! for vector rbf interpolation at
                                            ! each cell center
                                            ! (nproma,nblks_c)
  REAL(wp), ALLOCATABLE :: rbf_vec_coeff_c(:,:,:,:) ! array containing the
                                            ! coefficients used for
                                            ! vector rbf interpolation
                                            ! at each cell center
                                            ! (rbf_vec_dim_c,2,nproma,nblks_c)

  INTEGER, ALLOCATABLE  :: rbf_c2grad_idx(:,:,:)  ! index array defining the
                                            ! stencil of surrounding cells for
                                            ! 2D gradient reconstruction at each
                                            ! cell center
                                            ! (rbf_c2grad_dim,nproma,nblks_c)
  INTEGER, ALLOCATABLE  :: rbf_c2grad_blk(:,:,:)  ! ... dito for the blocks

  REAL(wp), ALLOCATABLE :: rbf_c2grad_coeff(:,:,:,:) ! array containing the
                                            ! coefficients used for
                                            ! 2D gradient reconstruction
                                            ! at each cell center
                                            ! (rbf_c2grad_dim,2,nproma,nblks_c)

  INTEGER, ALLOCATABLE  :: rbf_vec_idx_v(:,:,:) ! index array defining the
                                            ! stencil of surrounding edges for
                                            ! vector rbf interpolation at each
                                            ! triangle vertex
                                            ! (rbf_vec_dim_v,nproma,nblks_v)
  INTEGER, ALLOCATABLE  :: rbf_vec_blk_v(:,:,:) ! ... dito for the blocks

  INTEGER, ALLOCATABLE  :: rbf_vec_stencil_v(:,:) ! array defining number of
                                            ! surrounding edges in the stencil
                                            ! for vector rbf interpolation at
                                            ! each triangle vertex
                                            ! (nproma,nblks_v)

  REAL(wp), ALLOCATABLE :: rbf_vec_coeff_v(:,:,:,:) ! array containing the
                                            ! coefficients used for vector rbf
                                            ! interpolation at each tringle
                                            ! vertex (input is normal component)
                                            ! (rbf_vec_dim_v,2,nproma,nblks_v)

  INTEGER, ALLOCATABLE  :: rbf_vec_idx_e(:,:,:) ! index array defining the
                                            ! stencil of surrounding edges for
                                            ! vector rbf interpolation at each
                                            ! triangle edge
                                            ! (rbf_vec_dim_e,nproma,nblks_e)
  INTEGER, ALLOCATABLE  :: rbf_vec_blk_e(:,:,:) ! ... dito for the blocks

  INTEGER, ALLOCATABLE  :: rbf_vec_stencil_e(:,:) ! array defining number of
                                            ! surrounding edges in the stencil
                                            ! for vector rbf interpolation at
                                            ! each triangle edge
                                            ! (nproma,nblks_e)

  REAL(wp), ALLOCATABLE :: rbf_vec_coeff_e(:,:,:) ! array containing the
                                            ! coefficients used for rbf inter-
                                            ! polation of the tangential velo-
                                            ! city component (from the
                                            ! surrounding normals) at each
                                            ! triangle edge
                                            ! (rbf_vec_dim_e,nproma,nblks_e)

  ! d) fields needed for reconstructions in hexagonal model
  !--------------------------------------------------------
  REAL(wp), ALLOCATABLE :: heli_coeff(:,:,:)  ! coefficients for lamb_rot computation
  INTEGER, ALLOCATABLE  :: heli_vn_idx(:,:,:) ! len indices for vn in lamb_rot
  INTEGER, ALLOCATABLE  :: heli_vn_blk(:,:,:) ! blk indices for vn in lamb_rot
  REAL(wp), ALLOCATABLE :: hex_north(:,:,:)   ! coeffs for north vector in hexagon center
  REAL(wp), ALLOCATABLE :: hex_east(:,:,:)    ! coeffs for east vector in hexagon center
  REAL(wp), ALLOCATABLE :: tria_north(:,:,:)  ! coeffs for north vector in triangle center
  REAL(wp), ALLOCATABLE :: tria_east(:,:,:)   ! coeffs for east vector in triangle center
  REAL(wp), ALLOCATABLE :: quad_north(:,:,:)  ! coeffs for north vector in rhombus center
  REAL(wp), ALLOCATABLE :: quad_east(:,:,:)   ! coeffs for east vector in rhomus center
  REAL(wp), ALLOCATABLE :: cno_en(:,:,:)      ! coeffs for north vector in edge center
  REAL(wp), ALLOCATABLE :: cea_en(:,:,:)      ! coeffs for east vector in edge center

  ! e) precomputed geometrical factors for mathematical operators (for efficiency)
  !------------------------------------------------------------------------------
  REAL(wp), ALLOCATABLE :: geofac_div(:,:,:)    ! factor for divergence (nproma,i_cell_type,nblks_c)
  REAL(wp), ALLOCATABLE :: geofac_qdiv(:,:,:)   ! factor for quad-cell divergence (nproma,4,nblks_e)
  REAL(wp), ALLOCATABLE :: geofac_rot(:,:,:)    ! factor for divergence (nproma,9-i_cell_type,nblks_v)
  REAL(wp), ALLOCATABLE :: geofac_n2s(:,:,:)    ! factor for nabla2-scalar (nproma,i_cell_type+1,nblks_c)
  REAL(wp), ALLOCATABLE :: geofac_grg(:,:,:,:)  ! factor for Green-Gauss gradient (nproma,4,nblks_c,2)

  ! f) precomputed Cartesian orientation and location vectors of edge midpoints
  !    and location of cell centers(for efficiency)
  !------------------------------------------------------------------------------
  REAL(wp), ALLOCATABLE :: cart_edge_coord(:,:,:)    ! Cartesian edge coordinates (nproma,nblks_e)
  REAL(wp), ALLOCATABLE :: cart_cell_coord(:,:,:)    ! Cartesian cell coordinates (nproma,nblks_c)

  ! g) patch elements restored from edges to cells to reduce frequency of indirect addressing
  !------------------------------------------------------------------------------
  REAL(wp), ALLOCATABLE :: primal_normal_ec(:,:,:,:) ! p_patch%edges%primal_normal_cell stored on
                                                     ! the cell data type (nproma,nblks_c,3,2)
  REAL(wp), ALLOCATABLE :: edge_cell_length(:,:,:)   ! p_patch%edges%edge_cell_length stored on
                                                     ! the cell data type (nproma,nblks_c,3)

  ! h) fields related to calculation of backward trajectories on local plane
  !    tangential to the edge midpoint
  !------------------------------------------------------------------------------
  REAL(wp), ALLOCATABLE :: pos_on_tplane_e(:,:,:,:)  ! positions of various points on local plane
                                                     ! tangential to the edge midpoint.
                                                     ! currently: (nproma,nblks_e,8,2)
  REAL(wp), ALLOCATABLE :: tplane_e_dotprod(:,:,:,:) ! Dot product between unit vectors at inner
                                                     ! edge midpoint and edge midpoints of
                                                     ! corresponding quadrilateral cell.
                                                     ! (nproma,nblks_e,4,4)
                                                     

  ! i) fields related to weighted least squares polynomial reconstruction
  !------------------------------------------------------------------------------
  TYPE(t_lsq) :: lsq_lin,  &  ! coefficients for linear lsq-reconstruction
    &            lsq_high     ! coefficients for higher order lsq-reconstruction


  ! j) fields related to third order advection and Smagorinski diffusion (on hexagons):
  ! directional laplacian ( as du/dx > xx, dv/dy > yy, du/dx=dv/dy > xy )
  !----------------------------------------------------------------------
  INTEGER, ALLOCATABLE :: dir_gradh_i1(:,:,:) &
  & ; !< index array for edges of neighbor hexagon cell 1 of considered edge (6,nproma,nblks_e)
  INTEGER, ALLOCATABLE :: dir_gradh_b1(:,:,:) &
  & ; !< block array for edges of neighbor hexagon cell 1 of considered edge (6,nproma,nblks_e)
  REAL(wp), ALLOCATABLE :: dir_gradhxx_c1(:,:,:) &
  & ; !< coeff array for edges of neighbor hexagon cell 1 of considered edge (6,nproma,nblks_e)
  INTEGER, ALLOCATABLE :: dir_gradh_i2(:,:,:) &
  & ; !< index array for edges of neighbor hexagon cell 2 of considered edge (6,nproma,nblks_e)
  INTEGER, ALLOCATABLE :: dir_gradh_b2(:,:,:) &
  & ; !< block array for edges of neighbor hexagon cell 2 of considered edge (6,nproma,nblks_e)
  REAL(wp), ALLOCATABLE :: dir_gradhxx_c2(:,:,:) &
  & ; !< coeff array for edges of neighbor hexagon cell 2 of considered edge (6,nproma,nblks_e)

  INTEGER, ALLOCATABLE :: dir_gradt_i1(:,:,:) &
  & ; !< index array for edges of neighbor triangle cell 1 of considered edge (12,nproma,nblks_e)
  INTEGER, ALLOCATABLE :: dir_gradt_b1(:,:,:) &
  & ; !< block array for edges of neighbor triangle cell 1 of considered edge (12,nproma,nblks_e)
  REAL(wp), ALLOCATABLE :: dir_gradtxy_v1(:,:,:) &
  & ; !< coeff array for edges of neighbor triangle cell 1 of considered edge (12,nproma,nblks_e)
  INTEGER, ALLOCATABLE :: dir_gradt_i2(:,:,:) &
  & ; !< index array for edges of neighbor triangle cell 2 of considered edge (12,nproma,nblks_e)
  INTEGER, ALLOCATABLE :: dir_gradt_b2(:,:,:) &
  & ; !< block array for edges of neighbor triangle cell 2 of considered edge (12,nproma,nblks_e)
  REAL(wp), ALLOCATABLE :: dir_gradtxy_v2(:,:,:) &
  & ; !< coeff array for edges of neighbor triangle cell 2 of considered edge (12,nproma,nblks_e)

  ! k) Nudging coefficients used for 1-way nesting and limited-area mode (defined here
  !    rather than in grf_state because the limited-area mode may be used without nesting)
  !------------------------------------------------------------------------------
  REAL(wp), ALLOCATABLE :: nudgecoeff_c(:,:)  !< Nudging coefficient for cells
  REAL(wp), ALLOCATABLE :: nudgecoeff_e(:,:)  !< Nudging coefficient for cells


  ! l) Quadrature points and weights for integration over triangular element
  !--------------------------------------------------------------------------
  TYPE(t_gauss_quad) ::gquad

END TYPE t_int_state




END MODULE mo_intp_data_strc
