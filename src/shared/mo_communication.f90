!>
!!               This module provides the communication routines.
!!
!!               This module provides the communication routines
!! for parallel runs
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
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
MODULE mo_communication
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!

USE mo_kind,               ONLY: wp
USE mo_exception,          ONLY: finish
USE mo_global_variables,   ONLY: nproma
USE mo_mpi,                ONLY: p_pe, p_nprocs, p_send, p_recv, p_irecv, p_wait
USE mo_parallel_ctl,       ONLY: p_pe_work, p_test_pe, p_n_work, p_comm_work

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

!modules interface-------------------------------------------
!subroutines
PUBLIC :: blk_no, idx_no, idx_1d
PUBLIC :: setup_comm_pattern, exchange_data, exchange_data_reverse, exchange_data_hydro_c, &
          exchange_data_mult, exchange_data_grf, exchange_data_reverse_mult, &
          start_delayed_exchange, do_delayed_exchange
!
!variables

TYPE t_comm_pattern

   ! Number of points we receive in communication,
   ! this is the same as recv_limits

   INTEGER :: n_recv ! Number of points we receive from other PEs
   INTEGER :: n_pnts ! Number of points we output into local array;
                     ! this may be bigger than n_recv due to
                     ! duplicate entries

   INTEGER, ALLOCATABLE :: recv_limits(:)

   INTEGER, ALLOCATABLE :: recv_src(:)
   INTEGER, ALLOCATABLE :: recv_dst_blk(:)
   INTEGER, ALLOCATABLE :: recv_dst_idx(:)

   INTEGER :: n_send

   INTEGER, ALLOCATABLE :: send_limits(:)

   INTEGER, ALLOCATABLE :: send_src_blk(:)
   INTEGER, ALLOCATABLE :: send_src_idx(:)

END TYPE

PUBLIC t_comm_pattern

TYPE t_buffer

   INTEGER :: nelems
   REAL(wp), ALLOCATABLE :: buf(:)

END TYPE

TYPE(t_buffer), ALLOCATABLE :: send_bufs(:)


INTEGER, PARAMETER :: max_delayed_requests = 10000

INTEGER :: n_delayed_requests = 0

TYPE t_request

   TYPE(t_comm_pattern), POINTER :: p_pat
   INTEGER :: ndim2
   LOGICAL :: reverse
   REAL(wp), POINTER :: recv2(:,:)
   REAL(wp), POINTER :: recv3(:,:,:)
   REAL(wp), POINTER :: add2(:,:)
   REAL(wp), POINTER :: add3(:,:,:)

END TYPE

TYPE(t_request) :: delayed_request(max_delayed_requests)

LOGICAL :: use_exchange_delayed = .FALSE.


INTERFACE exchange_data
   MODULE PROCEDURE exchange_data_3
   MODULE PROCEDURE exchange_data_2
END INTERFACE

INTERFACE exchange_data_reverse
   MODULE PROCEDURE exchange_data_reverse_3
   MODULE PROCEDURE exchange_data_reverse_2
END INTERFACE

!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
! The following functions are for conversion of 1D to 2D indices and vice versa
!
! Treatment of 0 (important for empty patches) and negative numbers:
!
! Converting 1D => 2D:
!
! 0 always is mapped to blk_no = 1, idx_no = 0
! negative numbers: Convert usings ABS(j) and negate idx_no
!
! Thus: blk_no >= 1 always!
!       idx_no > 0  for j > 0
!       idx_no = 0  for j = 0
!       idx_no < 0  for j < 0
!
! This mimics mostly the behaviour of reshape_idx in mo_model_domimp_patches
! with a difference for nproma=1 and j=0 (where reshape_idx returns blk_no=0, idx_no=1)
!
! The consisten treatment of 0 in the above way is very important for empty patches
! where start_index=1, end_index=0
!
! Converting 2D => 1D:
! Trying to invert the above and catching cases with blk_no < 1
!
!-------------------------------------------------------------------------
ELEMENTAL INTEGER FUNCTION blk_no(j)
  INTEGER, INTENT(in) :: j
  blk_no = MAX((ABS(j)-1)/nproma + 1, 1) ! i.e. also 1 for j=0, nproma=1
END FUNCTION blk_no
!-------------------------------------------------------------------------
ELEMENTAL INTEGER FUNCTION idx_no(j)
  INTEGER, INTENT(in) :: j
  IF(j==0) THEN
    idx_no = 0
  ELSE
    idx_no = SIGN(MOD(ABS(j)-1,nproma)+1, j)
  ENDIF
END FUNCTION idx_no
!-------------------------------------------------------------------------
ELEMENTAL INTEGER FUNCTION idx_1d(jl,jb)
  INTEGER, INTENT(in) :: jl, jb
  IF(jb<=0) THEN
    idx_1d = 0 ! This covers the special case nproma==1,jb=0,jl=1
               ! All other cases are invalid and get also a 0 returned
  ELSE
    idx_1d = SIGN((jb-1)*nproma + ABS(jl), jl)
  ENDIF
END FUNCTION idx_1d
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!
!

!>
!! Sets up a communication pattern for exchanging data.
!!
!! n_points       Total number of points in the RECEIVER array,
!!                not every point is necessarily set during exchange
!!                (see owner!)
!!
!! owner          Owner PE number of every point in the RECEIVER array,
!!                if owner(.) == -1, this point will not be set during exchange.
!!                If owner(.) == p_pe, this point will be exchanged,
!!                this is necessary if sender and receiver arrays are
!!                different (e.g. feedback, gather, scatter)
!!
!! global_index   Global index of of every point in the RECEIVER array
!!                There may be more tha 1 point with the same global index,
!!                in this case the point is exchanged only once and
!!                locally distributed.
!!
!! local_index    Local index in the SENDER array.
!!                This array must have the local index at the global
!!                position for points which are local and a negative
!!                value at every position which is not owned by the local PE.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
SUBROUTINE setup_comm_pattern(n_points, owner, global_index, local_index, p_pat)

!

   INTEGER, INTENT(IN) :: n_points         ! Total number of points
   INTEGER, INTENT(IN) :: owner(:)         ! Owner of every point
   INTEGER, INTENT(IN) :: global_index(:)  ! Global index of every point
   INTEGER, INTENT(IN) :: local_index(:)   ! Array mapping global indices to local ones
                                           ! valid on the remote side!

   TYPE(t_comm_pattern), INTENT(OUT) :: p_pat



   INTEGER, ALLOCATABLE :: icnt(:), flag(:), global_recv_index(:), send_src(:)

   INTEGER :: i, n, np, nr, num_recv, irs, ire, num_send, iss, ise, max_glb

!-----------------------------------------------------------------------

   if(p_nprocs == 1 .OR. p_pe == p_test_pe) &
      CALL finish('setup_comm_pattern','must not be called on single PE/test PE')

   ALLOCATE(icnt(0:p_n_work-1))
   max_glb = MAX(MAXVAL(ABS(global_index(1:n_points)),mask=(owner(1:n_points)>=0)),1)
   ALLOCATE(flag(max_glb))

   ! Count the number of points we want to receive from every PE
   ! and the total number of points to output

   icnt(:) = 0
   flag(:) = 0

   p_pat%n_pnts = 0

   DO i = 1, n_points
      IF(owner(i)>=0) THEN
         p_pat%n_pnts = p_pat%n_pnts + 1 ! Count total number of points we output
         IF(flag(ABS(global_index(i)))==0) THEN
            icnt(owner(i)) = icnt(owner(i))+1 ! Number to get from owner(i)
            flag(ABS(global_index(i))) = 1 ! Flag that this global point is already on the list
         ENDIF
      ENDIF
   ENDDO

   ! Allocate and set up the recv_limits array

   ALLOCATE(p_pat%recv_limits(0:p_n_work))

   p_pat%recv_limits(0) = 0
   DO np = 0, p_n_work-1
      p_pat%recv_limits(np+1) = p_pat%recv_limits(np) + icnt(np)
   ENDDO

   ! The last entry in recv_limits is the total number of points we receive

   p_pat%n_recv = p_pat%recv_limits(p_n_work)

   ! Allocate and set up the recv_src array

   ALLOCATE(p_pat%recv_src(p_pat%n_pnts))
   ALLOCATE(p_pat%recv_dst_blk(p_pat%n_pnts))
   ALLOCATE(p_pat%recv_dst_idx(p_pat%n_pnts))
   ALLOCATE(global_recv_index(p_pat%n_recv))

   DO np = 0, p_n_work-1
      icnt(np) = p_pat%recv_limits(np)
   ENDDO

   flag(:) = 0
   n = 0 ! Counts total number of local points

   DO i = 1, n_points
      IF(owner(i)>=0) THEN
         n = n+1
         IF(flag(ABS(global_index(i)))==0) THEN
            icnt(owner(i)) = icnt(owner(i)) + 1    ! Current index in recv array
            global_recv_index(icnt(owner(i))) = ABS(global_index(i))
                                                   ! Global index of points in receive array
            p_pat%recv_src(n) = icnt(owner(i))     ! From where in the receive array we get
                                                   ! the local point
            p_pat%recv_dst_blk(n) = blk_no(i)      ! Where to put the local point
            p_pat%recv_dst_idx(n) = idx_no(i)      ! Where to put the local point
            flag(ABS(global_index(i))) = icnt(owner(i)) ! Store from where to get duplicates
         ELSE
            p_pat%recv_src(n) = flag(ABS(global_index(i)))
            p_pat%recv_dst_blk(n) = blk_no(i)
            p_pat%recv_dst_idx(n) = idx_no(i)
         ENDIF
      ENDIF
   ENDDO


   ! Exchange the number of points we want to receive with the respective senders

   DO np = 0, p_n_work-1 ! loop over PEs where to send the data

      num_recv = p_pat%recv_limits(np+1) - p_pat%recv_limits(np)
      irs = p_pat%recv_limits(np)+1 ! Start index in global_recv_index
      ire = p_pat%recv_limits(np+1) ! End   index in global_recv_index

      IF(np/=p_pe_work) THEN
         ! Just send the global index of the points we need from PE np
         CALL p_send(num_recv, np, 1, comm=p_comm_work)
         if(num_recv>0) &
           CALL p_send(global_recv_index(irs:ire), np, 2, comm=p_comm_work)
      ELSE

         ! In this turn, we receive the info about which points are needed from us

         ! First get the number of points only
         DO nr = 0, p_n_work-1
            IF(nr /= p_pe_work) THEN
               CALL p_recv(icnt(nr), nr, 1,  comm=p_comm_work)
            ELSE
               icnt(nr) = num_recv
            ENDIF
         ENDDO

         ! Allocate and set up the send_limits array

         ALLOCATE(p_pat%send_limits(0:p_n_work))

         p_pat%send_limits(0) = 0
         DO nr = 0, p_n_work-1
            p_pat%send_limits(nr+1) = p_pat%send_limits(nr) + icnt(nr)
         ENDDO

         ! The last entry in send_limits is the total number of points we receive

         p_pat%n_send = p_pat%send_limits(p_n_work)

         ! Allocate and set up the send_src array

         ALLOCATE(send_src(p_pat%n_send))

         DO nr = 0, p_n_work-1
            num_send = p_pat%send_limits(nr+1) - p_pat%send_limits(nr)
            iss = p_pat%send_limits(nr)+1 ! Start index in send_src
            ise = p_pat%send_limits(nr+1) ! End   index in send_src
            IF(nr /= p_pe_work) THEN
               IF(num_send>0) &
                  CALL p_recv(send_src(iss:ise), nr, 2, comm=p_comm_work)
            ELSE
               IF(num_send>0) &
                  send_src(iss:ise) = global_recv_index(irs:ire)
            ENDIF
         ENDDO

      ENDIF

   ENDDO

   ALLOCATE(p_pat%send_src_blk(p_pat%n_send))
   ALLOCATE(p_pat%send_src_idx(p_pat%n_send))

   ! The indices in p_pat%send_src are global, convert ot local

   DO i = 1, p_pat%n_send
      np = send_src(i)
      if(np<1 .or. np>SIZE(local_index)) CALL finish('setup_comm_pattern','Got illegal index')
      np = local_index(np)
      if(np<0) CALL finish('setup_comm_pattern','Got illegal index')
      p_pat%send_src_blk(i) = blk_no(np)
      p_pat%send_src_idx(i) = idx_no(np)
   ENDDO

   DEALLOCATE(icnt, flag, global_recv_index, send_src)

END SUBROUTINE setup_comm_pattern

!-------------------------------------------------------------------------
!
!
!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Modified by Guenther Zaengl for vectorization
!!
SUBROUTINE exchange_data_3(p_pat, recv, send, add, send_lbound3)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   REAL(wp), INTENT(INOUT), TARGET        :: recv(:,:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
   INTEGER, OPTIONAL :: send_lbound3


   REAL(wp) :: send_buf(SIZE(recv,2),p_pat%n_send), &
               recv_buf(SIZE(recv,2),p_pat%n_recv)

   REAL(wp), POINTER :: send_ptr(:,:,:)

   INTEGER :: i, k, np, irs, ire, iss, ise, ndim2, lbound3

!-----------------------------------------------------------------------

   IF(p_nprocs == 1 .OR. p_pe == p_test_pe) &
     CALL finish('exchange_data','must not be called on single PE/test PE')

   IF(SIZE(recv,1) /= nproma) THEN
     CALL finish('exchange_data','Illegal first dimension of data array')
   ENDIF

   ndim2 = SIZE(recv,2)

   IF(PRESENT(send) .AND. PRESENT(send_lbound3)) THEN
     lbound3 = send_lbound3
   ELSE
     lbound3 = 1
   ENDIF

   ! Set up send buffer

   IF(PRESENT(send)) THEN
     send_ptr => send
   ELSE
     send_ptr => recv
   ENDIF


   IF (ndim2 == 1) THEN
     DO i = 1, p_pat%n_send
       send_buf(1,i) = send_ptr(p_pat%send_src_idx(i),1,p_pat%send_src_blk(i)-lbound3+1)
     ENDDO
   ELSE
#ifdef __SX__
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_send
         send_buf(k,i) = send_ptr(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)-lbound3+1)
       ENDDO
     ENDDO
#else
     DO i = 1, p_pat%n_send
       send_buf(1:ndim2,i) = send_ptr(p_pat%send_src_idx(i),1:ndim2,   &
         &                            p_pat%send_src_blk(i)-lbound3+1)
     ENDDO
#endif
   ENDIF

   IF(use_exchange_delayed) THEN
     CALL add_delayed_request(p_pat, ndim2, .FALSE.)
     delayed_request(n_delayed_requests)%recv3 => recv
     IF(PRESENT(add)) delayed_request(n_delayed_requests)%add3 => add

     CALL buffer_data(p_pat%send_limits, send_buf)
     RETURN
   ENDIF


   ! Set up irecv's for receive buffers

   DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

      ! If the compiler makes a copy of recv_buf, we are lost !!!!!!!!!!!!!!!!!!

      irs = p_pat%recv_limits(np)+1
      ire = p_pat%recv_limits(np+1)
      IF(ire >= irs) &
         CALL p_irecv(recv_buf(1,irs), np, 1, p_count=(ire-irs+1)*ndim2, comm=p_comm_work)

   ENDDO

   ! Send our data

   DO np = 0, p_n_work-1 ! loop over PEs where to send the data

      iss = p_pat%send_limits(np)+1
      ise = p_pat%send_limits(np+1)

      IF(ise >= iss) &
         CALL p_send(send_buf(1,iss), np, 1, p_count=(ise-iss+1)*ndim2, comm=p_comm_work)

   ENDDO

   ! Wait for all outstanding requests to finish

   CALL p_wait

   ! Fill in receive buffer

   IF(PRESENT(add)) THEN
     IF (ndim2 == 1) THEN
       k = 1
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
           recv_buf(k,p_pat%recv_src(i)) +                       &
           add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
       ENDDO
     ELSE
#ifdef __SX__
!CDIR UNROLL=6
       DO k = 1, ndim2
         DO i = 1, p_pat%n_pnts
           recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
             recv_buf(k,p_pat%recv_src(i)) +                       &
             add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
         ENDDO
       ENDDO
#else
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),:,p_pat%recv_dst_blk(i)) =   &
           recv_buf(:,p_pat%recv_src(i)) +                       &
           add(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i))
       ENDDO
#endif
     ENDIF
   ELSE
     IF (ndim2 == 1) THEN
       k = 1
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
           recv_buf(k,p_pat%recv_src(i))
       ENDDO
     ELSE
#ifdef __SX__
!CDIR UNROLL=6
       DO k = 1, ndim2
         DO i = 1, p_pat%n_pnts
           recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
             recv_buf(k,p_pat%recv_src(i))
         ENDDO
       ENDDO
#else
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),:,p_pat%recv_dst_blk(i)) = &
           recv_buf(:,p_pat%recv_src(i))
       ENDDO
#endif
     ENDIF
   ENDIF

END SUBROUTINE exchange_data_3

!--------------------------------------------------------------------------------------------------
!>
!! Does a reverse data exchange according to a communication pattern (in p_pat),
!! i.e. if exchange_data transfers data from source point x to destination point y
!! this routine transfers from data from y to x.
!!
!! Please note:
!! The communication pattern used should have exactly one destination point
!! for every source point, otherwise it is undefined which of the different
!! destination data points is transferred to the source.
!!
!! !REVISION HISTORY:
!! Initial version by Rainer Johanni, Dec 2009
!!
SUBROUTINE exchange_data_reverse_3(p_pat, recv, send)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   REAL(wp), INTENT(INOUT), TARGET        :: recv(:,:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)

! !LOCAL VARIABLES:

   ! Note reversed second dimensions
   REAL(wp) :: send_buf(SIZE(recv,2),p_pat%n_recv), &
               recv_buf(SIZE(recv,2),p_pat%n_send)

   REAL(wp), POINTER :: send_ptr(:,:,:)

   INTEGER :: i, k, jb, jl, np, irs, ire, iss, ise, ndim2

!EOP
!-----------------------------------------------------------------------
!BOC

   IF(p_nprocs == 1 .OR. p_pe == p_test_pe) &
      CALL finish('exchange_data_reverse','must not be called on single PE/test PE')

   IF(SIZE(recv,1) /= nproma) THEN
      CALL finish('exchange_data_reverse','Illegal first dimension of data array')
   ENDIF

   ndim2 = SIZE(recv,2)

   ! Set up send buffer

   IF(PRESENT(send)) THEN
     send_ptr => send
   ELSE
     send_ptr => recv
   ENDIF

   IF (ndim2 == 1) THEN
     DO i = 1, p_pat%n_pnts
       send_buf(1,p_pat%recv_src(i)) = send_ptr(p_pat%recv_dst_idx(i),1,p_pat%recv_dst_blk(i))
     ENDDO
   ELSE
#ifndef __SX__
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_pnts
         send_buf(k,p_pat%recv_src(i)) = send_ptr(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
       ENDDO
     ENDDO
#else
     DO i = 1, p_pat%n_pnts
       send_buf(1:ndim2,p_pat%recv_src(i)) = &
         & send_ptr(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i))
     ENDDO
#endif
   ENDIF

   IF(use_exchange_delayed) THEN
     CALL add_delayed_request(p_pat, ndim2, .TRUE.)
     delayed_request(n_delayed_requests)%recv3 => recv

     CALL buffer_data(p_pat%recv_limits, send_buf)
     RETURN
   ENDIF


   ! Set up irecv's for receive buffers

   DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

      ! If the compiler makes a copy of recv_buf, we are lost !!!!!!!!!!!!!!!!!!

      irs = p_pat%send_limits(np)+1
      ire = p_pat%send_limits(np+1)
      IF(ire >= irs) &
         CALL p_irecv(recv_buf(1,irs), np, 1, p_count=(ire-irs+1)*ndim2, comm=p_comm_work)

   ENDDO

   ! Send our data

   DO np = 0, p_n_work-1 ! loop over PEs where to send the data

      iss = p_pat%recv_limits(np)+1
      ise = p_pat%recv_limits(np+1)

      IF(ise >= iss) &
         CALL p_send(send_buf(1,iss), np, 1, p_count=(ise-iss+1)*ndim2, comm=p_comm_work)

   ENDDO

   ! Wait for all outstanding requests to finish

   CALL p_wait

   ! Fill in receive buffer

#ifndef __SX__
   DO i = 1, p_pat%n_send
      jb = p_pat%send_src_blk(i)
      jl = p_pat%send_src_idx(i)
      recv(jl,:,jb) = recv_buf(:,i)
   ENDDO
#else
   IF (ndim2 == 1) THEN
     k = 1
     DO i = 1, p_pat%n_send
       recv(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)) = recv_buf(k,i)
     ENDDO
   ELSE
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_send
         recv(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)) = recv_buf(k,i)
       ENDDO
     ENDDO
   ENDIF
#endif

END SUBROUTINE exchange_data_reverse_3

!-------------------------------------------------------------------------
!
!

!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process surface pressure and temperature
!! in one step
!!
SUBROUTINE exchange_data_hydro_c(p_pat, recv2d, recv3d, send2d, send3d, add2d, add3d)

!

   TYPE(t_comm_pattern), INTENT(IN) :: p_pat
   REAL(wp), INTENT(IN), OPTIONAL :: send2d(:,:), send3d(:,:,:)
   REAL(wp), INTENT(IN), OPTIONAL :: add2d(:,:),  add3d(:,:,:)

   REAL(wp), INTENT(INOUT) :: recv2d(:,:), recv3d(:,:,:)



   REAL(wp) :: send_buf(SIZE(recv3d,2)+1,p_pat%n_send), &
               recv_buf(SIZE(recv3d,2)+1,p_pat%n_recv)

   INTEGER :: i, k, np, irs, ire, iss, ise, ndim2, np1

!-----------------------------------------------------------------------

   ndim2 = SIZE(recv3d,2)
   np1   = ndim2 + 1

!   IF ( PRESENT(send3d) .NEQV. PRESENT(send2d) ) &
!     CALL finish('exchange_data_hydro_c','inconsistent argument list')
!  IF ( PRESENT(add3d) .NEQV. PRESENT(add2d) ) &
!     CALL finish('exchange_data_hydro_c','inconsistent argument list')

   IF(use_exchange_delayed) THEN
     IF(PRESENT(send3d)) THEN
       IF(PRESENT(add3d)) THEN
         CALL exchange_data_3(p_pat, recv3d, send3d, add3d)
         CALL exchange_data_2(p_pat, recv2d, send2d, add2d)
       ELSE
         CALL exchange_data_3(p_pat, recv3d, send3d)
         CALL exchange_data_2(p_pat, recv2d, send2d)
       ENDIF
     ELSE
       IF(PRESENT(add3d)) THEN
         CALL exchange_data_3(p_pat, recv3d, add=add3d)
         CALL exchange_data_2(p_pat, recv2d, add=add2d)
       ELSE
         CALL exchange_data_3(p_pat, recv3d)
         CALL exchange_data_2(p_pat, recv2d)
       ENDIF
     ENDIF
     RETURN
   ENDIF

   ! Set up send buffer

   IF ( PRESENT(send3d) ) THEN

!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_send
         send_buf(k,i) = send3d(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i))
       ENDDO
     ENDDO

     DO i = 1, p_pat%n_send
       send_buf(np1,i) = send2d(p_pat%send_src_idx(i),p_pat%send_src_blk(i))
     ENDDO

   ELSE
     ! Send and receive arrays are identical (for boundary exchange)

!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_send
         send_buf(k,i) = recv3d(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i))
       ENDDO
     ENDDO

     DO i = 1, p_pat%n_send
       send_buf(np1,i) = recv2d(p_pat%send_src_idx(i),p_pat%send_src_blk(i))
     ENDDO

   ENDIF


   ! Set up irecv's for receive buffers

   DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

      ! If the compiler makes a copy of recv_buf, we are lost !!!!!!!!!!!!!!!!!!

      irs = p_pat%recv_limits(np)+1
      ire = p_pat%recv_limits(np+1)
      IF(ire >= irs) &
         CALL p_irecv(recv_buf(1,irs), np, 1, p_count=(ire-irs+1)*np1, comm=p_comm_work)

   ENDDO

   ! Send our data

   DO np = 0, p_n_work-1 ! loop over PEs where to send the data

      iss = p_pat%send_limits(np)+1
      ise = p_pat%send_limits(np+1)

      IF(ise >= iss) &
         CALL p_send(send_buf(1,iss), np, 1, p_count=(ise-iss+1)*np1, comm=p_comm_work)

   ENDDO

   ! Wait for all outstanding requests to finish

   CALL p_wait

   ! Fill in receive buffer

   IF(PRESENT(add3d)) THEN
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_pnts
         recv3d(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
           recv_buf(k,p_pat%recv_src(i)) +                         &
           add3d(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
       ENDDO
     ENDDO

     DO i = 1, p_pat%n_pnts
       recv2d(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i)) =   &
         recv_buf(np1,p_pat%recv_src(i)) +                     &
         add2d(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i))
     ENDDO

   ELSE
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_pnts
         recv3d(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
           recv_buf(k,p_pat%recv_src(i))
       ENDDO
     ENDDO

     DO i = 1, p_pat%n_pnts
       recv2d(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i)) = &
         recv_buf(np1,p_pat%recv_src(i))
     ENDDO

   ENDIF

END SUBROUTINE exchange_data_hydro_c


!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process 4D fields or up to five 3D fields
!! in one step
!!
SUBROUTINE exchange_data_mult(p_pat, nfields, ndim2tot, recv1, send1, add1, recv2, send2, &
                              add2, recv3, send3, add3, recv4, send4, add4, recv5, send5, &
                              add5, recv6, send6, add6, recv4d, send4d, add4d, lpar)

   TYPE(t_comm_pattern), INTENT(IN) :: p_pat

   REAL(wp), INTENT(INOUT), TARGET, OPTIONAL ::  &
     recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4(:,:,:), recv5(:,:,:), recv6(:,:,:), &
     recv4d(:,:,:,:)
   REAL(wp), INTENT(IN   ), TARGET, OPTIONAL ::  &
     send1(:,:,:), send2(:,:,:), send3(:,:,:), send4(:,:,:), send5(:,:,:), send6(:,:,:), &
     send4d(:,:,:,:)
   REAL(wp), INTENT(IN   ), TARGET, OPTIONAL ::  &
     add1(:,:,:), add2(:,:,:), add3(:,:,:), add4(:,:,:), add5(:,:,:), add6(:,:,:),       &
     add4d(:,:,:,:)

   INTEGER, INTENT(IN)           :: nfields, ndim2tot
   LOGICAL, OPTIONAL, INTENT(IN) :: lpar

   TYPE t_fieldptr
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr

   TYPE(t_fieldptr) :: recv(nfields), send(nfields), add(nfields)
   INTEGER        :: ndim2(nfields), noffset(nfields)

   REAL(wp) :: send_buf(ndim2tot,p_pat%n_send),recv_buf(ndim2tot,p_pat%n_recv)

   INTEGER :: i, k, ik, jb, jl, n, np, irs, ire, iss, ise
   LOGICAL :: lsend, ladd, l_par

!-----------------------------------------------------------------------

   lsend     = .FALSE.
   ladd      = .FALSE.

   ! Set pointers to input fields
   IF (PRESENT(recv4d)) THEN
     DO i = 1, nfields
       recv(i)%fld => recv4d(:,:,:,i)
     ENDDO
     IF (PRESENT(send4d)) THEN ! all 4D fields must have the same dimensions
       DO i = 1, nfields
         send(i)%fld => send4d(:,:,:,i)
       ENDDO
       lsend = .TRUE.
     ENDIF
     IF (PRESENT(add4d)) THEN
       DO i = 1, nfields
         add(i)%fld => add4d(:,:,:,i)
       ENDDO
       ladd = .TRUE.
     ENDIF
   ELSE
     IF (PRESENT(recv1)) THEN
       recv(1)%fld => recv1
       IF (PRESENT(send1)) THEN
         send(1)%fld => send1
         lsend = .TRUE.
       ENDIF
       IF (PRESENT(add1)) THEN
         add(1)%fld => add1
         ladd = .TRUE.
       ENDIF
       IF (PRESENT(recv2)) THEN
         recv(2)%fld => recv2
         IF (lsend) send(2)%fld => send2
         IF (ladd)  add(2)%fld  => add2
         IF (PRESENT(recv3)) THEN
           recv(3)%fld => recv3
           IF (lsend) send(3)%fld => send3
           IF (ladd)  add(3)%fld  => add3
           IF (PRESENT(recv4)) THEN
             recv(4)%fld => recv4
             IF (lsend) send(4)%fld => send4
             IF (ladd)  add(4)%fld  => add4
             IF (PRESENT(recv5)) THEN
               recv(5)%fld => recv5
               IF (lsend) send(5)%fld => send5
               IF (ladd)  add(5)%fld  => add5
               IF (PRESENT(recv6)) THEN
                 recv(6)%fld => recv6
                 IF (lsend) send(6)%fld => send6
                 IF (ladd)  add(6)%fld  => add6
               ENDIF
             ENDIF
           ENDIF
         ENDIF
       ENDIF
     ENDIF
   ENDIF

   IF(use_exchange_delayed) THEN
     DO n = 1, nfields
       IF(lsend) THEN
         IF(ladd) THEN
           CALL exchange_data_3(p_pat, recv(n)%fld(:,:,:), send(n)%fld(:,:,:), add(n)%fld(:,:,:))
         ELSE
           CALL exchange_data_3(p_pat, recv(n)%fld(:,:,:), send(n)%fld(:,:,:))
         ENDIF
       ELSE
         IF(ladd) THEN
           CALL exchange_data_3(p_pat, recv(n)%fld(:,:,:), add=add(n)%fld(:,:,:))
         ELSE
           CALL exchange_data_3(p_pat, recv(n)%fld(:,:,:))
         ENDIF
       ENDIF
     ENDDO
     RETURN
   ENDIF

   IF (PRESENT(lpar)) THEN
     l_par = lpar
   ELSE
     l_par = .FALSE.
   ENDIF

   noffset(1) = 0
   ndim2(1)   = SIZE(recv(1)%fld,2)
   DO n = 2, nfields
     noffset(n) = noffset(n-1)+ndim2(n-1)
     ndim2(n)   = SIZE(recv(n)%fld,2)
   ENDDO

   ! Set up send buffer
#ifdef __SX__
   IF ( lsend .AND. l_par ) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(n,k,i)
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_send
           send_buf(k+noffset(n),i) = &
             send(n)%fld(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i))
         ENDDO
       ENDDO
     ENDDO
!$OMP END DO
!$OMP END PARALLEL
   ELSE IF ( lsend ) THEN
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_send
           send_buf(k+noffset(n),i) = &
             send(n)%fld(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ELSE
       ! Send and receive arrays are identical (for boundary exchange)
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_send
           send_buf(k+noffset(n),i) = &
             recv(n)%fld(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ENDIF
#else
   DO i = 1, p_pat%n_send
     jb = p_pat%send_src_blk(i)
     jl = p_pat%send_src_idx(i)
     IF ( lsend ) THEN
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),i) = send(n)%fld(jl,k,jb)
         ENDDO
       ENDDO
     ELSE
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),i) = recv(n)%fld(jl,k,jb)
         ENDDO
       ENDDO
     ENDIF
   ENDDO
#endif


   ! Set up irecv's for receive buffers

   DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

      ! If the compiler makes a copy of recv_buf, we are lost !!!!!!!!!!!!!!!!!!

      irs = p_pat%recv_limits(np)+1
      ire = p_pat%recv_limits(np+1)
      IF(ire >= irs) &
         CALL p_irecv(recv_buf(1,irs), np, 1, p_count=(ire-irs+1)*ndim2tot, comm=p_comm_work)

   ENDDO

   ! Send our data

   DO np = 0, p_n_work-1 ! loop over PEs where to send the data

      iss = p_pat%send_limits(np)+1
      ise = p_pat%send_limits(np+1)

      IF(ise >= iss) &
         CALL p_send(send_buf(1,iss), np, 1, p_count=(ise-iss+1)*ndim2tot, comm=p_comm_work)

   ENDDO

   ! Wait for all outstanding requests to finish

   CALL p_wait

   ! Fill in receive buffer

#ifdef __SX__
   IF (ladd .AND. l_par) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(n,k,i)
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           recv(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
             recv_buf(k+noffset(n),p_pat%recv_src(i)) + &
             add(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
         ENDDO
       ENDDO
     ENDDO
!$OMP END DO
!$OMP END PARALLEL
   ELSE IF (ladd) THEN
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           recv(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
             recv_buf(k+noffset(n),p_pat%recv_src(i)) + &
             add(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ELSE
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           recv(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
             recv_buf(k+noffset(n),p_pat%recv_src(i))
         ENDDO
       ENDDO
     ENDDO
   ENDIF
#else
   DO i = 1, p_pat%n_pnts
     jb = p_pat%recv_dst_blk(i)
     jl = p_pat%recv_dst_idx(i)
     ik  = p_pat%recv_src(i)
     IF (ladd) THEN
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           recv(n)%fld(jl,k,jb) = recv_buf(k+noffset(n),ik)+add(n)%fld(jl,k,jb)
         ENDDO
       ENDDO
     ELSE
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           recv(n)%fld(jl,k,jb) = recv_buf(k+noffset(n),ik)
         ENDDO
       ENDDO
     ENDIF
   ENDDO
#endif

END SUBROUTINE exchange_data_mult

!>
!! Does reverse data exchange according to a communication pattern (in p_pat).
!!
!! Please note:
!! The communication pattern used should have exactly one destination point
!! for every source point, otherwise it is undefined which of the different
!! destination data points is transferred to the source.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process 4D fields or up to five 3D fields
!! in one step
!!
SUBROUTINE exchange_data_reverse_mult(p_pat, nfields, ndim2tot, recv1, send1, recv2, send2,    &
                                      recv3, send3, recv4, send4, recv5, send5, recv4d, send4d )

   TYPE(t_comm_pattern), INTENT(IN) :: p_pat

   REAL(wp), INTENT(INOUT), TARGET, OPTIONAL ::  &
     recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4(:,:,:), recv5(:,:,:), recv4d(:,:,:,:)
   REAL(wp), INTENT(IN   ), TARGET, OPTIONAL ::  &
     send1(:,:,:), send2(:,:,:), send3(:,:,:), send4(:,:,:), send5(:,:,:), send4d(:,:,:,:)

   INTEGER, INTENT(IN)           :: nfields, ndim2tot

   TYPE t_fieldptr
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr

   TYPE(t_fieldptr) :: recv(nfields), send(nfields)
   INTEGER        :: ndim2(nfields), noffset(nfields)

   ! Note reversed second dimensions for reverse communication
   REAL(wp) :: send_buf(ndim2tot,p_pat%n_recv),recv_buf(ndim2tot,p_pat%n_send)

   INTEGER :: i, k, jb, jl, n, np, irs, ire, iss, ise
   LOGICAL :: lsend

!-----------------------------------------------------------------------

   lsend     = .FALSE.

   ! Set pointers to input fields
   IF (PRESENT(recv4d)) THEN
     DO i = 1, nfields
       recv(i)%fld => recv4d(:,:,:,i)
     ENDDO
     IF (PRESENT(send4d)) THEN ! all 4D fields must have the same dimensions
       DO i = 1, nfields
         send(i)%fld => send4d(:,:,:,i)
       ENDDO
       lsend = .TRUE.
     ENDIF
   ELSE
     IF (PRESENT(recv1)) THEN
       recv(1)%fld => recv1
       IF (PRESENT(send1)) THEN
         send(1)%fld => send1
         lsend = .TRUE.
       ENDIF
       IF (PRESENT(recv2)) THEN
         recv(2)%fld => recv2
         IF (lsend) send(2)%fld => send2
         IF (PRESENT(recv3)) THEN
           recv(3)%fld => recv3
           IF (lsend) send(3)%fld => send3
           IF (PRESENT(recv4)) THEN
             recv(4)%fld => recv4
             IF (lsend) send(4)%fld => send4
             IF (PRESENT(recv5)) THEN
               recv(5)%fld => recv5
               IF (lsend) send(5)%fld => send5
             ENDIF
           ENDIF
         ENDIF
       ENDIF
     ENDIF
   ENDIF

   IF(use_exchange_delayed) THEN
     DO n = 1, nfields
       IF(lsend) THEN
         CALL exchange_data_reverse_3(p_pat, recv(n)%fld(:,:,:), send(n)%fld(:,:,:))
       ELSE
         CALL exchange_data_reverse_3(p_pat, recv(n)%fld(:,:,:))
       ENDIF
     ENDDO
     RETURN
   ENDIF

   noffset(1) = 0
   ndim2(1)   = SIZE(recv(1)%fld,2)
   DO n = 2, nfields
     noffset(n) = noffset(n-1)+ndim2(n-1)
     ndim2(n)   = SIZE(recv(n)%fld,2)
   ENDDO

   ! Set up send buffer
#ifdef __SX__
   IF ( lsend ) THEN
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           send_buf(k+noffset(n),p_pat%recv_src(i)) = &
             send(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ELSE
       ! Send and receive arrays are identical (for boundary exchange)
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           send_buf(k+noffset(n),p_pat%recv_src(i)) = &
             recv(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ENDIF
#else
   DO i = 1, p_pat%n_pnts
     jb = p_pat%recv_dst_blk(i)
     jl = p_pat%recv_dst_idx(i)
     IF ( lsend ) THEN
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),p_pat%recv_src(i)) = send(n)%fld(jl,k,jb)
         ENDDO
       ENDDO
     ELSE
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),p_pat%recv_src(i)) = recv(n)%fld(jl,k,jb)
         ENDDO
       ENDDO
     ENDIF
   ENDDO
#endif


   ! Set up irecv's for receive buffers

   DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

      ! If the compiler makes a copy of recv_buf, we are lost !!!!!!!!!!!!!!!!!!

      irs = p_pat%send_limits(np)+1
      ire = p_pat%send_limits(np+1)
      IF(ire >= irs) &
         CALL p_irecv(recv_buf(1,irs), np, 1, p_count=(ire-irs+1)*ndim2tot, comm=p_comm_work)

   ENDDO

   ! Send our data

   DO np = 0, p_n_work-1 ! loop over PEs where to send the data

      iss = p_pat%recv_limits(np)+1
      ise = p_pat%recv_limits(np+1)

      IF(ise >= iss) &
         CALL p_send(send_buf(1,iss), np, 1, p_count=(ise-iss+1)*ndim2tot, comm=p_comm_work)

   ENDDO

   ! Wait for all outstanding requests to finish

   CALL p_wait

   ! Fill in receive buffer

#ifdef __SX__
   DO n = 1, nfields
!CDIR UNROLL=6
     DO k = 1, ndim2(n)
       DO i = 1, p_pat%n_send
         recv(n)%fld(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)) =   &
           recv_buf(k+noffset(n),i)
       ENDDO
     ENDDO
   ENDDO
#else
   DO i = 1, p_pat%n_send
     jb = p_pat%send_src_blk(i)
     jl = p_pat%send_src_idx(i)
     DO n = 1, nfields
       DO k = 1, ndim2(n)
         recv(n)%fld(jl,k,jb) = recv_buf(k+noffset(n),i)
       ENDDO
     ENDDO
   ENDDO
#endif

END SUBROUTINE exchange_data_reverse_mult


!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process 4D fields or up to three 3D fields
!! for an array-sized communication pattern (as needed for boundary interpolation) in one step
!!
SUBROUTINE exchange_data_grf(p_pat, nfields, ndim2tot, nsendtot, nrecvtot, recv1, send1, &
                             recv2, send2, recv3, send3, recv4d, send4d, send_lbound3)

   TYPE(t_comm_pattern), INTENT(IN) :: p_pat(:)

   REAL(wp), INTENT(INOUT), TARGET, OPTIONAL ::  &
     recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4d(:,:,:,:)
   ! Note: the send fields have one additional dimension in this case because
   ! the fourth index corresponds to the dimension of p_pat
   REAL(wp), INTENT(IN   ), TARGET, OPTIONAL ::  &
     send1(:,:,:,:), send2(:,:,:,:), send3(:,:,:,:), send4d(:,:,:,:,:)

   INTEGER, INTENT(IN)           :: nfields  ! total number of input fields
   INTEGER, INTENT(IN)           :: ndim2tot ! sum of vertical levels of input fields
   INTEGER, INTENT(IN)           :: nsendtot ! total number of send points
   INTEGER, INTENT(IN)           :: nrecvtot ! total number of receive points
   INTEGER, OPTIONAL, INTENT(IN) :: send_lbound3

   TYPE t_fieldptr
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr

   TYPE(t_fieldptr) :: recv(nfields), send(nfields*SIZE(p_pat))
   INTEGER        :: ndim2(nfields), noffset(nfields),            &
                     ioffset_s(SIZE(p_pat)), ioffset_r(SIZE(p_pat))

   REAL(wp) :: send_buf(ndim2tot,nsendtot),recv_buf(ndim2tot,nrecvtot), &
               auxs_buf(ndim2tot,nsendtot),auxr_buf(ndim2tot,nrecvtot)

   INTEGER :: i, k, ik, jb, jl, n, np, irs, ire, iss, ise, &
              npats, isum, ioffset, isum1

!-----------------------------------------------------------------------

    npats = SIZE(p_pat)  ! Number of communication patterns provided on input


   IF(use_exchange_delayed) THEN
     IF(.NOT.PRESENT(send_lbound3)) CALL finish('exchange_data_grf','Missing send_lbound3')
     IF (PRESENT(recv4d)) THEN
       DO n = 1, UBOUND(recv4d, 4)
         DO np = 1, npats
           CALL exchange_data_3(p_pat(np), recv4d(:,:,:,n), send4d(:,:,:,np,n), &
             &                  send_lbound3=send_lbound3)
         ENDDO
       ENDDO
     ENDIF
     IF (PRESENT(recv1)) THEN
       DO np = 1, npats
         CALL exchange_data_3(p_pat(np), recv1(:,:,:), send1(:,:,:,np),send_lbound3=send_lbound3)
       ENDDO
     ENDIF
     IF (PRESENT(recv2)) THEN
       DO np = 1, npats
         CALL exchange_data_3(p_pat(np), recv2(:,:,:), send2(:,:,:,np),send_lbound3=send_lbound3)
       ENDDO
     ENDIF
     IF (PRESENT(recv3)) THEN
       DO np = 1, npats
         CALL exchange_data_3(p_pat(np), recv3(:,:,:), send3(:,:,:,np),send_lbound3=send_lbound3)
       ENDDO
     ENDIF

     RETURN
   ENDIF

   ! Set pointers to input fields
   IF (PRESENT(recv4d)) THEN
     DO n = 1, nfields
       recv(n)%fld => recv4d(:,:,:,n)
       DO np = 1, npats
         send(np+(n-1)*npats)%fld => send4d(:,:,:,np,n)
       ENDDO
     ENDDO
   ELSE
     IF (PRESENT(recv1)) THEN
       recv(1)%fld => recv1
       DO np = 1, npats
         send(np)%fld => send1(:,:,:,np)
       ENDDO
     ENDIF
     IF (PRESENT(recv2)) THEN
       recv(2)%fld => recv2
       DO np = 1, npats
         send(np+npats)%fld => send2(:,:,:,np)
       ENDDO
     ENDIF
     IF (PRESENT(recv3)) THEN
       recv(3)%fld => recv3
       DO np = 1, npats
         send(np+2*npats)%fld => send3(:,:,:,np)
       ENDDO
     ENDIF
   ENDIF

   noffset(1) = 0
   ndim2(1)   = SIZE(recv(1)%fld,2)
   DO n = 2, nfields
     noffset(n) = noffset(n-1)+ndim2(n-1)
     ndim2(n)   = SIZE(recv(n)%fld,2)
   ENDDO

   ioffset_r(1) = 0
   ioffset_s(1) = 0
   DO np = 2, npats
     ioffset_r(np) = ioffset_r(np-1) + p_pat(np-1)%n_recv
     ioffset_s(np) = ioffset_s(np-1) + p_pat(np-1)%n_send
   ENDDO

   ! Set up send buffer
#ifdef __SX__
   DO n = 1, nfields
     DO np = 1, npats
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
!CDIR NODEP
         DO i = 1, p_pat(np)%n_send
           send_buf(k+noffset(n),i+ioffset_s(np)) =                                            &
             & send(np+(n-1)*npats)%fld(p_pat(np)%send_src_idx(i),k,p_pat(np)%send_src_blk(i)- &
             & send_lbound3+1)
         ENDDO
       ENDDO
     ENDDO
   ENDDO
#else
   DO np = 1, npats
     DO i = 1, p_pat(np)%n_send
       jb = p_pat(np)%send_src_blk(i)
       jl = p_pat(np)%send_src_idx(i)
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),i+ioffset_s(np)) = &
             send(np+(n-1)*npats)%fld(jl,k,jb-send_lbound3+1)
         ENDDO
       ENDDO
     ENDDO
   ENDDO
#endif

   ! Set up irecv's for receive buffers
   ioffset = 0
   DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

     ! Sum up receive points over all communication patterns to be processed
     isum = ioffset
     DO n = 1, npats
       isum = isum + p_pat(n)%recv_limits(np+1) - p_pat(n)%recv_limits(np)
     ENDDO

     IF(isum > ioffset) &
       CALL p_irecv(auxr_buf(1,ioffset+1), np, 1, p_count=(isum-ioffset)*ndim2tot, &
                    comm=p_comm_work)

     ioffset = isum

   ENDDO

   ! Send our data
   ioffset = 0
   DO np = 0, p_n_work-1 ! loop over PEs where to send the data

     ! Copy send points for all communication patterns into one common send buffer
     isum = ioffset
     DO n = 1, npats
       iss = p_pat(n)%send_limits(np)+1 + ioffset_s(n)
       ise = p_pat(n)%send_limits(np+1) + ioffset_s(n)
       isum1 = ise - iss + 1
       IF (isum1 > 0) THEN
!CDIR COLLAPSE
         auxs_buf(:,isum+1:isum+isum1) = send_buf(:,iss:ise)
         isum = isum+isum1
       ENDIF
     ENDDO

     IF(isum > ioffset) &
       CALL p_send(auxs_buf(1,ioffset+1), np, 1, p_count=(isum-ioffset)*ndim2tot, &
                    comm=p_comm_work)

     ioffset = isum

   ENDDO

   ! Wait for all outstanding requests to finish
   CALL p_wait

   ! Copy exchanged data back to receive buffer
   ioffset = 0
   DO np = 0, p_n_work-1

     isum = ioffset
     DO n = 1, npats
       irs = p_pat(n)%recv_limits(np)+1 + ioffset_r(n)
       ire = p_pat(n)%recv_limits(np+1) + ioffset_r(n)
       isum1 = ire - irs + 1
       IF (isum1 > 0) THEN
!CDIR COLLAPSE
         recv_buf(:,irs:ire) = auxr_buf(:,isum+1:isum+isum1)
         isum = isum + isum1
       ENDIF
     ENDDO

     ioffset = isum

   ENDDO

   ! Fill in receive buffer

#ifdef __SX__
   DO n = 1, nfields
     DO np = 1, npats
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat(np)%n_pnts
           recv(n)%fld(p_pat(np)%recv_dst_idx(i),k,p_pat(np)%recv_dst_blk(i)) =   &
             recv_buf(k+noffset(n),p_pat(np)%recv_src(i)+ioffset_r(np))
         ENDDO
       ENDDO
     ENDDO
   ENDDO
#else
   DO np = 1, npats
     DO i = 1, p_pat(np)%n_pnts
       jb = p_pat(np)%recv_dst_blk(i)
       jl = p_pat(np)%recv_dst_idx(i)
       ik  = p_pat(np)%recv_src(i)+ioffset_r(np)
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           recv(n)%fld(jl,k,jb) = recv_buf(k+noffset(n),ik)
         ENDDO
       ENDDO
     ENDDO
   ENDDO
#endif


END SUBROUTINE exchange_data_grf


!-------------------------------------------------------------------------
!
!

!>
!! Interface for 2D arrays for exchange_data.
!!
!! Just reshapes the arrays and calls exchange_data_3.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
SUBROUTINE exchange_data_2(p_pat, recv, send, add, send_lbound2)

!

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   REAL(wp), INTENT(INOUT), TARGET        :: recv(:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
   INTEGER, OPTIONAL :: send_lbound2

   REAL(wp) :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))
   REAL(wp) :: send_buf(1,p_pat%n_send)
   INTEGER :: i, lbound2

!-----------------------------------------------------------------------

   IF(use_exchange_delayed) THEN

     CALL add_delayed_request(p_pat, 1, .FALSE.)

     delayed_request(n_delayed_requests)%recv2 => recv
     IF(PRESENT(add)) delayed_request(n_delayed_requests)%add2 => add

     IF(PRESENT(send) .AND. PRESENT(send_lbound2)) THEN
       lbound2 = send_lbound2
     ELSE
       lbound2 = 1
     ENDIF

     IF ( PRESENT(send)) THEN
       DO i = 1, p_pat%n_send
         send_buf(1,i) = send(p_pat%send_src_idx(i),p_pat%send_src_blk(i)-lbound2+1)
       ENDDO
     ELSE
       DO i = 1, p_pat%n_send
         send_buf(1,i) = recv(p_pat%send_src_idx(i),p_pat%send_src_blk(i))
       ENDDO
     ENDIF

     CALL buffer_data(p_pat%send_limits, send_buf)
     RETURN

   ENDIF

   tmp_recv(:,1,:) = recv(:,:)

   IF(PRESENT(send) .AND. PRESENT(send_lbound2)) THEN
      IF(PRESENT(add)) THEN
         CALL exchange_data_3(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
                              add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)), &
                              send_lbound3=send_lbound2 )
      ELSE
         CALL exchange_data_3(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
                              send_lbound3=send_lbound2)
      ENDIF
   ELSE IF(PRESENT(send)) THEN
      IF(PRESENT(add)) THEN
         CALL exchange_data_3(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
                              add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data_3(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
      ENDIF
   ELSE
      IF(PRESENT(add)) THEN
         CALL exchange_data_3(p_pat, tmp_recv, &
                              add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data_3(p_pat, tmp_recv)
      ENDIF
   ENDIF

   recv(:,:) = tmp_recv(:,1,:)

END SUBROUTINE exchange_data_2

!-------------------------------------------------------------------------

!>
!! Interface for 2D arrays for exchange\_data\_reverse.
!!
!! Just reshapes the arrays and calls exchange\_data\_reverse\_3.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Dec 2009
!!
SUBROUTINE exchange_data_reverse_2(p_pat, recv, send)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   REAL(wp), INTENT(INOUT), TARGET        :: recv(:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)


! !LOCAL VARIABLES:

   REAL(wp) :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))
   REAL(wp) :: send_buf(1,p_pat%n_recv)
   INTEGER :: i

!EOP
!-----------------------------------------------------------------------
!BOC

   IF(use_exchange_delayed) THEN

     CALL add_delayed_request(p_pat, 1, .TRUE.)
     delayed_request(n_delayed_requests)%recv2 => recv

     IF(PRESENT(send)) THEN
       DO i = 1, p_pat%n_pnts
         send_buf(1,p_pat%recv_src(i)) = send(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i))
       ENDDO
     ELSE
       DO i = 1, p_pat%n_pnts
         send_buf(1,p_pat%recv_src(i)) = recv(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i))
       ENDDO
     ENDIF

     CALL buffer_data(p_pat%recv_limits, send_buf)
     RETURN

   ENDIF

   tmp_recv(:,1,:) = recv(:,:)

   IF(PRESENT(send)) THEN
      CALL exchange_data_reverse_3(p_pat, tmp_recv, &
                           send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
   ELSE
      CALL exchange_data_reverse_3(p_pat, tmp_recv)
   ENDIF

   recv(:,:) = tmp_recv(:,1,:)

END SUBROUTINE exchange_data_reverse_2

!--------------------------------------------------------------------------------------------------

!> Starts a delayed exchange, i.e. all following exchange_data_xxx calls are
!! buffered until do_delayed_exchange() is called which actually does the
!! communication calls and distributes the received data.
!!
!! ATTENTION:
!! Since a pointer to the recv and add arrays is stored internally,
!! these arrays must be a real reference (either the target array or a pointer)
!! and not a copy to the real data and  must neither go out of scope nor
!! be deallocated between the call to exchange_data_xxx and the corresponding
!! call to do_delayed_exchange().
!! Same holds for p_pat (but the communication patterns shouldn't be
!! changed after initialization anyways).
!!
!! Initial version by Rainer Johanni, Oct 2010

SUBROUTINE start_delayed_exchange

  INTEGER :: np

  ! Just ignore this call for 1 processor runs
  IF (p_nprocs == 1 .OR. p_pe == p_test_pe) RETURN

  ! If the first time here: initialize send_bufs

  IF(.NOT.ALLOCATED(send_bufs)) THEN

    ALLOCATE(send_bufs(0:p_n_work-1))

    DO np = 0, p_n_work-1
      ALLOCATE(send_bufs(np)%buf(1)) ! will be enlarged later
      send_bufs(np)%nelems = 0 ! nothing filled in yet
    ENDDO

  ENDIF

  use_exchange_delayed = .TRUE.

END SUBROUTINE start_delayed_exchange

!--------------------------------------------------------------------------------------------------
!
!> Adds a delayed request, i.e. increments n_delayed_requests, checks if the limit is not exceeded
!! sets the p_pat, ndim2, reverse members and nullifies the remaining pointers

SUBROUTINE add_delayed_request(p_pat, ndim2, reverse)

  TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
  INTEGER, INTENT(IN) :: ndim2
  LOGICAL, INTENT(IN) :: reverse

  ! Check for limit

  IF(n_delayed_requests >= max_delayed_requests) &
    CALL finish('mo_communication','max_delayed_requests exceeded')

  n_delayed_requests = n_delayed_requests+1

  delayed_request(n_delayed_requests)%p_pat => p_pat
  delayed_request(n_delayed_requests)%ndim2 =  ndim2
  delayed_request(n_delayed_requests)%reverse = reverse
  delayed_request(n_delayed_requests)%recv2 => NULL()
  delayed_request(n_delayed_requests)%recv3 => NULL()
  delayed_request(n_delayed_requests)%add2  => NULL()
  delayed_request(n_delayed_requests)%add3  => NULL()

END SUBROUTINE add_delayed_request

!--------------------------------------------------------------------------------------------------

!> Backend for buffering data for delayed exchanges
!!
!! Initial version by Rainer Johanni, Oct 2010

SUBROUTINE buffer_data(limits, send)

  INTEGER, INTENT(IN) :: limits(0:)
  REAL(wp), INTENT(IN) :: send(:,:)

  REAL(wp), ALLOCATABLE :: buf(:)
  INTEGER :: ndim2, nitems, n, np, iss, ise

  ndim2 = UBOUND(send,1)


  ! Store data in send_bufs

  DO np = 0, p_n_work-1 ! loop over PEs where to send the data

    iss = limits(np)+1
    ise = limits(np+1)

    nitems = ndim2*(ise-iss+1)

    IF(ise<iss) CYCLE

    n = send_bufs(np)%nelems ! Number of elements currently in buffer

    ! Make sure that send_bufs(np)%buf has enough capacity.
    ! There is some deallocating/allocating of the buffers at the beginning
    ! but after the first timestep the final size should be reached.

    IF(n + nitems > SIZE(send_bufs(np)%buf)) THEN
      IF(n>0) THEN
        ALLOCATE(buf(n))
        buf(1:n) = send_bufs(np)%buf(1:n)
      ENDIF
      DEALLOCATE(send_bufs(np)%buf)
      ALLOCATE(send_bufs(np)%buf(n + nitems))
      IF(n>0) THEN
        send_bufs(np)%buf(1:n) = buf(1:n)
        DEALLOCATE(buf)
      ENDIF
    ENDIF

    send_bufs(np)%buf(n+1:n+nitems) = RESHAPE(send(:,iss:ise), (/ nitems /))

    send_bufs(np)%nelems = n+nitems ! new number of elements

  ENDDO

END SUBROUTINE buffer_data

!--------------------------------------------------------------------------------------------------

SUBROUTINE do_delayed_exchange()

!
! !DESCRIPTION:
! Does data exchange according to a communication pattern (in p_pat).
!
! !REVISION HISTORY:
! Initial version by Rainer Johanni, Oct 2010
!

! !INPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

! !LOCAL VARIABLES:

  TYPE(t_buffer), ALLOCATABLE :: recv_bufs(:)
  REAL(wp), ALLOCATABLE :: recv_buf(:,:)
  TYPE(t_comm_pattern), POINTER :: p_pat

  INTEGER :: i, n, np, nr, irs, ire, ndim2
#ifdef __SX__
  INTEGER :: k
#endif

!EOP
!-----------------------------------------------------------------------

  ! Just ignore this call for 1 processor runs
  IF (p_nprocs == 1 .OR. p_pe == p_test_pe) RETURN

  ALLOCATE(recv_bufs(0:p_n_work-1))

  ! Set up irecv's for receive buffers

  DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

    n = 0
    DO nr = 1, n_delayed_requests
      IF(delayed_request(nr)%reverse) THEN
        irs = delayed_request(nr)%p_pat%send_limits(np)+1
        ire = delayed_request(nr)%p_pat%send_limits(np+1)
      ELSE
        irs = delayed_request(nr)%p_pat%recv_limits(np)+1
        ire = delayed_request(nr)%p_pat%recv_limits(np+1)
      ENDIF
      n = n + (ire-irs+1)*delayed_request(nr)%ndim2
    ENDDO

    IF(n>0) THEN
      ALLOCATE(recv_bufs(np)%buf(n))
      CALL p_irecv(recv_bufs(np)%buf, np, 1, p_count=n, comm=p_comm_work)
    ENDIF

  ENDDO

  ! Send our data

  DO np = 0, p_n_work-1 ! loop over PEs where to send the data

     IF(send_bufs(np)%nelems > 0) &
        CALL p_send(send_bufs(np)%buf, np, 1, p_count=send_bufs(np)%nelems, comm=p_comm_work)

  ENDDO

  ! Wait for all outstanding requests to finish

  CALL p_wait

  recv_bufs(:)%nelems = 0 ! counts consumed elements below

  ! Fill in receive buffer

  DO nr = 1, n_delayed_requests

    p_pat => delayed_request(nr)%p_pat
    ndim2 =  delayed_request(nr)%ndim2

    IF(delayed_request(nr)%reverse) THEN
      ALLOCATE(recv_buf(ndim2,p_pat%n_send))
    ELSE
      ALLOCATE(recv_buf(ndim2,p_pat%n_recv))
    ENDIF

    DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

      IF(delayed_request(nr)%reverse) THEN
        irs = p_pat%send_limits(np)+1
        ire = p_pat%send_limits(np+1)
      ELSE
        irs = p_pat%recv_limits(np)+1
        ire = p_pat%recv_limits(np+1)
      ENDIF
      IF(ire<irs) CYCLE

      n = recv_bufs(np)%nelems
      recv_buf(1:ndim2,irs:ire) = RESHAPE(recv_bufs(np)%buf(n+1:n+ndim2*(ire-irs+1)), &
        &                                 (/ndim2,ire-irs+1/))
      recv_bufs(np)%nelems = recv_bufs(np)%nelems + ndim2*(ire-irs+1)
    ENDDO

    IF(delayed_request(nr)%reverse) THEN
      IF(ASSOCIATED(delayed_request(nr)%recv3)) THEN
#ifdef __SX__
!CDIR UNROLL=6
        DO k = 1, ndim2
          DO i = 1, p_pat%n_send
            delayed_request(nr)%recv3(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)) = &
              recv_buf(k,i)
          ENDDO
        ENDDO
#else
        DO i = 1, p_pat%n_send
          delayed_request(nr)%recv3(p_pat%send_src_idx(i),1:ndim2,p_pat%send_src_blk(i)) = &
            recv_buf(1:ndim2,i)
        ENDDO
#endif
      ELSE
        DO i = 1, p_pat%n_send
          delayed_request(nr)%recv2(p_pat%send_src_idx(i),p_pat%send_src_blk(i)) = &
            recv_buf(1,i)
        ENDDO
      ENDIF
    ELSE
      IF(ASSOCIATED(delayed_request(nr)%recv3)) THEN
        IF(ASSOCIATED(delayed_request(nr)%add3)) THEN
#ifdef __SX__
!CDIR UNROLL=6
          DO k = 1, ndim2
            DO i = 1, p_pat%n_pnts
              delayed_request(nr)%recv3(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
                recv_buf(k,p_pat%recv_src(i)) +                                          &
                delayed_request(nr)%add3(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
            ENDDO
          ENDDO
#else
          DO i = 1, p_pat%n_pnts
            delayed_request(nr)%recv3(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i)) = &
              recv_buf(1:ndim2,p_pat%recv_src(i)) +                                          &
              delayed_request(nr)%add3(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i))
          ENDDO
#endif
        ELSE
#ifdef __SX__
!CDIR UNROLL=6
          DO k = 1, ndim2
            DO i = 1, p_pat%n_pnts
              delayed_request(nr)%recv3(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
                recv_buf(k,p_pat%recv_src(i))
            ENDDO
          ENDDO
#else
          DO i = 1, p_pat%n_pnts
            delayed_request(nr)%recv3(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i)) = &
              recv_buf(1:ndim2,p_pat%recv_src(i))
          ENDDO
#endif
        ENDIF
      ELSE
        IF(ASSOCIATED(delayed_request(nr)%add2)) THEN
          DO i = 1, p_pat%n_pnts
            delayed_request(nr)%recv2(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i)) = &
              recv_buf(1,p_pat%recv_src(i)) +                                        &
              delayed_request(nr)%add2(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i))
          ENDDO
        ELSE
          DO i = 1, p_pat%n_pnts
            delayed_request(nr)%recv2(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i)) = &
              recv_buf(1,p_pat%recv_src(i))
          ENDDO
        ENDIF
      ENDIF
    ENDIF

    DEALLOCATE(recv_buf)

  ENDDO

  DO np = 0, p_n_work-1
     IF(ALLOCATED(recv_bufs(np)%buf)) DEALLOCATE(recv_bufs(np)%buf)
  ENDDO

  DEALLOCATE(recv_bufs)

  n_delayed_requests = 0

  DO np = 0, p_n_work-1
     send_bufs(np)%nelems = 0
  ENDDO

  use_exchange_delayed = .FALSE.

END SUBROUTINE do_delayed_exchange

END MODULE mo_communication
