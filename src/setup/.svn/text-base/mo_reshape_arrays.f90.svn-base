!>
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
MODULE mo_reshape_arrays
!-------------------------------------------------------------------------

USE mo_kind,               ONLY: wp
USE mo_impl_constants,      ONLY: SUCCESS, &
     &                           MAX_CHAR_LENGTH,  &
     &                           min_rlcell, max_rlcell, &
     &                           min_rledge, max_rledge, &
     &                           min_rlvert, max_rlvert
USE mo_exception,          ONLY: finish
USE mo_model_domain
USE mo_physical_constants, ONLY: omega
USE mo_global_variables,   ONLY: nproma, i_cell_type
USE mo_loopindices,        ONLY: get_indices_e
USE mo_mpi,                ONLY: p_pe, p_io

IMPLICIT NONE

PUBLIC

CONTAINS



!-------------------------------------------------------------------------
!>
!!               reshape an integer array for the blocking.
SUBROUTINE reshape_int( p_int_array_in, nblks, npromz,  &
  &                     p_reshaped_int_array_out )

! input array
INTEGER, INTENT(in) :: p_int_array_in(:)
! number of blocks
INTEGER, INTENT(in) :: nblks
! chunk length
INTEGER, INTENT(in) :: npromz

! output array
INTEGER, INTENT(inout) :: p_reshaped_int_array_out(:,:)

INTEGER :: nlen
INTEGER :: jl, jb
INTEGER :: il

!-----------------------------------------------------------------------

DO jb = 1, nblks

  IF (jb /= nblks) THEN
    nlen = nproma
  ELSE
    nlen = npromz
    DO jl = npromz+1, nproma
      p_reshaped_int_array_out(jl,nblks) = 0
    END DO
  END IF

  DO jl = 1, nlen
    il = jl + (jb-1)*nproma
    p_reshaped_int_array_out(jl,jb) = p_int_array_in(il)
  END DO

END DO

END SUBROUTINE reshape_int

!-------------------------------------------------------------------------
!>
!!               reshape a real array for the blocking.
SUBROUTINE reshape_real( p_real_array_in, nblks, npromz,  &
  &                      p_reshaped_real_array_out )



! input array
REAL(wp), INTENT(in):: p_real_array_in(:)
! number of blocks
INTEGER, INTENT(in) :: nblks
! chunk length
INTEGER, INTENT(in) :: npromz

! output array
REAL(wp), INTENT(inout) :: p_reshaped_real_array_out(:,:)

INTEGER :: nlen
INTEGER :: jl, jb
INTEGER :: il

!-----------------------------------------------------------------------

DO jb = 1, nblks

  IF (jb /= nblks) THEN
    nlen = nproma
  ELSE
    nlen = npromz
    DO jl = npromz+1, nproma
      p_reshaped_real_array_out(jl,nblks) = 0._wp
    END DO
  END IF

  DO jl = 1, nlen
    il = jl + (jb-1)*nproma
    p_reshaped_real_array_out(jl,jb) = p_real_array_in(il)
  END DO

END DO

END SUBROUTINE reshape_real

!-------------------------------------------------------------------------
!>
!!               reshape an integer index array for the blocking.
SUBROUTINE reshape_idx( p_idx_array_in, nblks, npromz,  &
  &                     p_reshaped_idx_array_out,       &
  &                     p_reshaped_blk_array_out )



! input index array
INTEGER, INTENT(in) :: p_idx_array_in(:)
! number of blocks
INTEGER, INTENT(in) :: nblks
! chunk length
INTEGER, INTENT(in) :: npromz

! output line index array
INTEGER, INTENT(inout) :: p_reshaped_idx_array_out(:,:)
! output block index array
INTEGER, INTENT(inout) :: p_reshaped_blk_array_out(:,:)

INTEGER :: nlen
INTEGER :: jl, jb
INTEGER :: il, idx_in, idx, blk

!-----------------------------------------------------------------------

DO jb = 1, nblks

  IF (jb /= nblks) THEN
    nlen = nproma
  ELSE
    nlen = npromz
    DO jl = npromz+1, nproma
      p_reshaped_idx_array_out(jl,nblks) = 0
      p_reshaped_blk_array_out(jl,nblks) = 0
    END DO
  END IF

  DO jl = 1, nlen
    il  = jl + ( jb - 1 )*nproma
    idx_in = p_idx_array_in(il)
    blk = ( ABS(idx_in) - 1 ) / nproma + 1
    idx = SIGN( ABS(idx_in) - ( blk - 1 )*nproma, idx_in )
    p_reshaped_idx_array_out(jl,jb) = idx
    p_reshaped_blk_array_out(jl,jb) = blk
  END DO

END DO

END SUBROUTINE reshape_idx
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!>
!!               reshape an index list array for the blocking.
SUBROUTINE reshape_indlist( indlist_in,                     &
  &                         start_idx_out, start_blk_out,   &
  &                         end_idx_out,   end_blk_out )

! input index array
INTEGER, INTENT(in) :: indlist_in(:,:)

! output line index array
INTEGER, INTENT(inout) :: start_idx_out(:,:), end_idx_out(:,:)
! output block index array
INTEGER, INTENT(inout) :: start_blk_out(:,:), end_blk_out(:,:)

!-----------------------------------------------------------------------

start_blk_out(:,1) = (indlist_in(:,1) - 1) / nproma + 1
start_idx_out(:,1) =  indlist_in(:,1) - (start_blk_out(:,1) - 1) * nproma

end_blk_out(:,1) = (indlist_in(:,2) - 1) / nproma + 1
end_idx_out(:,1) =  indlist_in(:,2) - (end_blk_out(:,1) - 1) * nproma


END SUBROUTINE reshape_indlist


!-------------------------------------------------------------------------
!>
!!               reshape the start_idx / end_idx fields for blocking.
!!
SUBROUTINE reshape_idx_list( indlist_in, idx_out, blk_out )
! input index array
INTEGER, INTENT(in) :: indlist_in(:,:)
! output line index array
INTEGER, INTENT(inout) :: idx_out(:,:)
! output block index array
INTEGER, INTENT(inout) :: blk_out(:,:)

! write(0,*) 'reshape_idx_list, global_max_childdom:', global_max_childdom
!-----------------------------------------------------------------------
blk_out(:,1:global_max_childdom) = (indlist_in(:,1:global_max_childdom) - 1) / nproma + 1
idx_out(:,1:global_max_childdom) =  indlist_in(:,1:global_max_childdom) - &
                            (blk_out(:,1:global_max_childdom) - 1) * nproma

END SUBROUTINE reshape_idx_list
!-------------------------------------------------------------------------


END MODULE mo_reshape_arrays
