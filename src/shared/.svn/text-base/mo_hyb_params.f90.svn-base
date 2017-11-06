!>
!! module *mo_hyb* - *loop indices and surface-pressure independent
!! variables associated with the vertical finite-difference scheme.
!!
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
MODULE mo_hyb_params

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max, find_next_free_unit
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_impl_constants,     ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_global_variables,   ONLY: nvclev,  nlev, nlevp1, &
    vertical_coordinates_file_name, msg_level
  USE mo_physical_constants, ONLY: g=>grav, rcpd, rd

  IMPLICIT NONE

  PUBLIC

  INTEGER :: nlevm1         ! (number of levels)-1.
  INTEGER :: nplev          ! *number of pressure levels.
  INTEGER :: nplvp1         ! *nplev+1.*
  INTEGER :: nplvp2         ! *nplev+2.*
  INTEGER :: nplvpa         ! *nplvp1,* or 2 if *nplev=0.*
  INTEGER :: nlmsgl         ! *nlev* - (number of sigma levels).
  INTEGER :: nlmslp         ! *nlmsgl+1.*
  INTEGER :: nlmsla         ! *nlmslp,* or 2 if *nlmslp=1.*

  REAL(wp) :: apzero        ! *reference pressure for computation of the
  !                         !  hybrid vertical levels.
  !REAL(wp) :: rpr           ! *reciprocal of reference surface pressure.
  !REAL(wp) :: rdtr          !  rd*(reference temperature).
  REAL(wp) :: t0icao        ! *surface temperatur of reference atmosphere
  REAL(wp) :: tsticao       ! *stratospheric temperature of reference atmosphere
  REAL(wp) :: rdtstic       ! *rd*tsticao
  REAL(wp) :: rdlnp0i       ! *rd*ln(surface pressure) of reference atmosphere
  REAL(wp) :: alrrdic       ! *lapse-rate parameter of reference atmosphere
  REAL(wp) :: rdt0ral       ! *rd*t0icao/alphaic
  REAL(wp) :: ptricao       ! *tropopause pressure of reference atmosphere
  REAL(wp) :: rdlnpti       ! *rd*ln(ptricao)
  REAL(wp) :: gsticao       ! *constant used in geopotential calculation

  REAL(wp), ALLOCATABLE :: ralpha(:) ! rd*alpha at pressure and sigma levels.
  REAL(wp), ALLOCATABLE :: rlnpr(:)  ! rd*ln(p(k+.5)/p(k-.5))
  REAL(wp), ALLOCATABLE :: dela(:)   ! a(k+.5)-a(k-.5).
  REAL(wp), ALLOCATABLE :: delb(:)   ! b(k+.5)-b(k-.5).
  REAL(wp), ALLOCATABLE :: rddelb(:) ! rd*delb.
  REAL(wp), ALLOCATABLE :: cpg(:)    ! a(k+.5)*b(k-.5)-b(k+.5)*a(k-.5).
  REAL(wp), ALLOCATABLE :: delpr(:)  ! p(k+.5)-p(k-.5) of
                                     ! the refrence surface pressure.
  REAL(wp), ALLOCATABLE :: rdelpr(:) ! *reciprocal of *delpr.*
  !REAL(wp), ALLOCATABLE :: ralphr(:) ! *constant array for use by pgrad.
  REAL(wp), ALLOCATABLE :: alpham(:) ! *constant array for use by dyn.
  REAL(wp), ALLOCATABLE :: ardprc(:) ! *constant array for use by dyn.
  !REAL(wp), ALLOCATABLE :: rlnmar(:) ! *constant array for use by pgrad.
  !REAL(wp), ALLOCATABLE :: aktlrd(:) ! *constant array for use by conteq.
  !REAL(wp), ALLOCATABLE :: altrcp(:) ! *constant array for use by conteq.
  REAL(wp), ALLOCATABLE :: ceta(:)   ! *full hybrid vertical levels.
  REAL(wp), ALLOCATABLE :: cetah(:)  ! *half hybrid vertical levels.
  !REAL(wp), ALLOCATABLE :: bb(:,:)   ! *gravity wave matrix

  REAL(wp), ALLOCATABLE :: vct_a(:) ! param. A of the vertical coordinte
  REAL(wp), ALLOCATABLE :: vct_b(:) ! param. B of the vertical coordinate
  REAL(wp), ALLOCATABLE :: vct  (:) ! param. A and B of the vertical coordinate


  PRIVATE :: alloc_hyb, init_hyb_params
  PUBLIC :: read_hyb_params, destruct_vertical_coord

  ! calling tree
  !
  ! (hydro_atmos) -> *init_vertical_coord* |- *read_hyb_params*
  !                                        |- *alloc_hyb*
  !                                        |- *init_hyb_params*

CONTAINS

  !-------------------------------------------------------------------------
  ! !IROUTINE:  init_vertical_coord
  !
  ! !SUBROUTINE INTERFACE:
  SUBROUTINE init_vertical_coord
    ! !DESCRIPTION:
    !  Initialization of the hybrid vertical coordinate
    !
    ! !REVISION HISTORY:
    !  Original version by Hui Wan, MPI-M, 2006-02-09
    !
    INTEGER :: jk
    !EOP
    !-----------------------------------------------------------------------
    !BOC

    ! read the A and B parameters of the vertical coordinate

    CALL read_hyb_params

    IF (msg_level > 10) THEN
      CALL message('mo_hyb_params:init_vertical_coord', '')
      CALL message('', 'Vertical coordinate table')
      CALL message('', '   k     vct_a(k) [Pa]   vct_b(k) []')
      DO jk = 1, SIZE(vct_a)
        WRITE(message_text,'(i4,f18.10,f14.10)')  jk, vct_a(jk), vct_b(jk)
        CALL message('', TRIM(message_text))
      ENDDO
    ENDIF
    ! allocate memory for the auxiliary parameters and arrays
    CALL alloc_hyb

    ! assign values to the the auxiliary parameters and arrays
    CALL init_hyb_params

  END SUBROUTINE init_vertical_coord
  !-------------------------------------------------------------------------

  !>
  !! Read the A and B parameters of the hybrid vertical grid,
  !! which define the half level pressure: ph=A+B*ps [Pa]
  SUBROUTINE  read_hyb_params (dyn_option)

    ! Optional argument to switch into nonhydrostatic mode
    CHARACTER(len=MAX_CHAR_LENGTH), OPTIONAL, INTENT(IN) :: dyn_option

    ! Local variables
    CHARACTER(len=MAX_CHAR_LENGTH),PARAMETER :: routine  = &
         &   'mo_hyb_params:read_hyb_params'
    CHARACTER(len=filename_max), PARAMETER :: hyb_file = 'HYB_PARAMS_'
    CHARACTER(len=filename_max)            :: hyb_file_n, line

    INTEGER :: ist, iunit

    INTEGER :: ik, jk
    !-------------------------------------------------------------------------
    ! allocate memory

    
    ALLOCATE(vct_a(nvclev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of vct_a failed')
    ENDIF

    ALLOCATE(vct_b(nvclev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of vct_b failed')
    ENDIF

    ALLOCATE(vct(nvclev*2), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of vct failed')
    ENDIF

    ! Open file
    iunit = find_next_free_unit(10,20)
    OPEN (UNIT=iunit,file=TRIM(vertical_coordinates_file_name),access='SEQUENTIAL', &
         &form='FORMATTED', IOSTAT=ist)

    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), 'open HYB_PARAMS file failed')
    ENDIF

    ! Skip header line
    READ (iunit,*,IOSTAT=ist) line
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), 'reading header line failed')
    ENDIF

    ! Read A and B
    DO jk=1,nvclev
       READ (iunit,*,IOSTAT=ist) ik, vct_a(jk), vct_b(jk)
       IF(ist/=SUCCESS)THEN
          CALL finish (TRIM(routine), 'reading vct_a and vct_b failed')
       ENDIF
    END DO

    CLOSE(iunit)

    vct(       1:       nvclev) = vct_a(:)
    vct(nvclev+1:nvclev+nvclev) = vct_b(:)

  END SUBROUTINE  read_hyb_params
  !-------------------------------------------------------------------------

  !>
  !!
  SUBROUTINE alloc_hyb

    INTEGER :: ist
    CHARACTER(len=MAX_CHAR_LENGTH),PARAMETER :: routine = &
         & 'mo_hyb_params:alloc_hyb'

    !-----------------------------------------------------------------------
    ALLOCATE (ralpha(nlev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of ralpha failed')
    ENDIF

    ALLOCATE (rlnpr(nlev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of rlnpr failed')
    ENDIF

    ALLOCATE (dela(nlev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of dela failed')
    ENDIF

    ALLOCATE (delb(nlev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of delb failed')
    ENDIF

    ALLOCATE (rddelb(nlev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of rddelb failed')
    ENDIF

    ALLOCATE (cpg(nlev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of cpg failed')
    ENDIF

    ALLOCATE (delpr(nlev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of rdelpr failed')
    ENDIF

    ALLOCATE (rdelpr(nlev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of rdelpr failed')
    ENDIF

    !ALLOCATE (ralphr(nlev), STAT=ist)
    !IF(ist/=SUCCESS)THEN
    !   CALL finish (TRIM(routine), ' allocation of ralphr failed')
    !ENDIF

    ALLOCATE (alpham(nlev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of alpham failed')
    ENDIF

    ALLOCATE (ardprc(nlev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of ardprc failed')
    ENDIF

    !ALLOCATE (rlnmar(nlev), STAT=ist)
    !IF(ist/=SUCCESS)THEN
    !   CALL finish (TRIM(routine), ' allocation of rlnmar failed')
    !ENDIF

    !ALLOCATE (aktlrd(nlev), STAT=ist)
    !IF(ist/=SUCCESS)THEN
    !   CALL finish (TRIM(routine), ' allocation of aktlrd failed')
    !ENDIF

    !ALLOCATE (altrcp(nlev), STAT=ist)
    !IF(ist/=SUCCESS)THEN
    !   CALL finish (TRIM(routine), ' allocation of altrcp failed')
    !ENDIF

    ALLOCATE (ceta(nlev), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of ceta failed')
    ENDIF

    ALLOCATE (cetah(nlevp1), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), ' allocation of cetah failed')
    ENDIF

    !  ALLOCATE (bb(nlev,nlev), STAT=ist)
    !  IF(ist/=SUCCESS)THEN
    !     CALL finish (TRIM(routine), ' allocation of bb failed')
    !  ENDIF

  END SUBROUTINE alloc_hyb
  !-------------------------------------------------------------------------

  !>
  !!  Initializes constants for vertical coordinate calculations.
  !!
  !!  Method:
  !!    Compute loop indices and surface-pressure independent
  !!    variables associated with the vertical finite-difference scheme.
  !!    Output is in module *mo_hyb*
  SUBROUTINE init_hyb_params
    !  Local scalars:
    REAL(wp) :: za, zb, zetam, zetap, zp, zp0icao, zpp, zrd, zs, zsm
    INTEGER  :: ilev, ilevp1, iplev, iplvp1, is, ism, ist, jk, jlev

    !  Intrinsic functions
    INTRINSIC EXP, LOG

    !-----------------------------------------------------------------------
    !BOC

    !  Executable statements

    !-- 1. Initialize variables

!ag    apzero    = 101325._wp ! changed for NCAR summer colloquium!
    apzero    = 100000._wp
    zrd       = rd
    ralpha(1) = zrd*LOG(2._wp)
    rlnpr(1)  = 2._wp*ralpha(1)
    ilev      = nlev
    ilevp1    = ilev + 1
    nlevp1    = ilevp1
    nlevm1    = ilev - 1
    iplev     = 0
    iplvp1    = 1
    is        = nvclev + ilevp1
    ism       = is - 1
    zpp       = vct(1)
    zsm       = vct(is)

    t0icao  = 288._wp
    tsticao = 216.5_wp
    zp0icao = 101320._wp
    rdlnp0i = rd*LOG(zp0icao)
    rdtstic = rd*tsticao
    alrrdic = 0.0065_wp/g
    rdt0ral = t0icao/alrrdic
    rdlnpti = rdlnp0i + (LOG(tsticao/t0icao))/alrrdic
    ptricao = EXP(rdlnpti/rd)
    gsticao = tsticao*(rdlnpti-1._wp/alrrdic)

    zb      = vct(nvclev+iplvp1+1)

    !-- 2. Calculate pressure-level values

    DO WHILE ( zb == 0._wp )

      iplev  = iplvp1
      iplvp1 = iplev + 1
      IF (iplvp1==ilevp1) EXIT    ! if all levels are pressure levels

      zp            = zpp
      zpp           = vct(iplvp1)
      delpr(iplev)  = zpp - zp
      rdelpr(iplev) = 1._wp/delpr(iplev)

      IF ( iplev>1 ) THEN
        rlnpr(iplev)  = zrd*LOG(zpp/zp)
        ralpha(iplev) = zrd - zp*rlnpr(iplev)/delpr(iplev)
      END IF

      alpham(iplev) = ralpha(iplev)*rcpd
      ardprc(iplev) = rlnpr(iplev)*rdelpr(iplev)*rcpd
      zb            = vct(nvclev+iplvp1+1)

    ENDDO


    IF (iplvp1/=ilevp1) THEN   ! All levels are not pressure-levels

      nplev  = iplev
      nplvp1 = iplvp1
      nplvp2 = iplvp1 + 1
      IF (iplev==0) THEN
        nplvpa = 2
      ELSE
        nplvpa = iplvp1
      END IF

      !-- 3. Calculate sigma-level values

      za = vct(ism-nvclev)

      DO WHILE ( za == 0._wp )

        is  = ism
        ism = is - 1
        ist = is - nvclev
        zs  = zsm
        zsm = vct(is)
        IF (ist==1) THEN
          nlmsgl = 0
          nlmslp = 1
          nlmsla = 2
          EXIT
        ELSE
          rlnpr(ist)  = zrd*LOG(zs/zsm)
          ralpha(ist) = zrd - zsm*rlnpr(ist)/(zs-zsm)
        END IF
        za = vct(ism-nvclev)

      END DO

      IF (za>0._wp) THEN
        nlmsgl = ism - nvclev
        nlmslp = nlmsgl + 1
        nlmsla = nlmslp
      END IF

      !-- 4. Calculate dela, delb, rddelb, cpg, and complete alphdb

      DO jk = 1, nlev
        dela(jk)   = vct(jk+1) - vct(jk)
        delb(jk)   = vct(nvclev+jk+1) - vct(nvclev+jk)
        rddelb(jk) = rd*delb(jk)
        cpg(jk)    = vct(nvclev+jk)*vct(jk+1) - vct(nvclev+jk+1)*vct(jk)
      END DO

      DO jk = nlmslp, nlev
        alpham(jk) = ralpha(jk)*delb(jk)
      END DO

    ENDIF  ! If all levels are not pressure-levels

    !-- 5. Compute full level values of the hybrid coordinate

    zetam    = vct(1)/apzero + vct(nvclev+1)
    cetah(1) = zetam

    DO jlev = 1, nlev
      zetap         = vct(jlev+1)/apzero + vct(nvclev+1+jlev)
      ceta(jlev)    = (zetam+zetap)*.5_wp
      cetah(jlev+1) = zetap
      zetam = zetap
    END DO

  END SUBROUTINE init_hyb_params
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE  destruct_vertical_coord
    !-------------------------------------------------------------------------
    DEALLOCATE(vct_a,vct_b,vct)
    DEALLOCATE (ralpha,rlnpr,dela, delb)
    DEALLOCATE (rddelb,cpg,delpr,rdelpr,alpham,ardprc)
    DEALLOCATE (ceta,cetah)

  END SUBROUTINE destruct_vertical_coord
  !-------------------------------------------------------------------------
END MODULE mo_hyb_params

