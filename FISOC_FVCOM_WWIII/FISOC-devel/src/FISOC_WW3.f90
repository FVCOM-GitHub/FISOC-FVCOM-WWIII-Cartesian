!/===========================================================================/
! Copyright (c) 2025 The University of Massachusetts Dartmouth 
! Produced at the School of Marine Science & Technology 
! Marine Ecosystem Dynamics Modeling group
! All rights reserved.
!
! FVCOM has been developed by the joint UMASSD-WHOI research team. For 
! details of authorship and attribution of credit please see the FVCOM
! technical manual or contact the MEDM group.
!
! 
! This file is part of FVCOM. For details, see http://www.fvcom.org. 
! The full copyright notice is contained in the file COPYRIGHT located in the 
! root directory of the FVCOM code. This original header must be maintained
! in all distributed versions.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED.  
!
!/---------------------------------------------------------------------------/
#include "w3macros.h"

# define W3_PDLIB
# define W3_MPI

MODULE FISOC_WW3
  !
  !  1. Readme:
  !
  !     This module was modified from WW3 V7.14 PROGRAM W3SHEL by Jianhua Qi at SMAST,UMASSD
  !
  !  8. Structure :
  !
  !     ----------------------------------------------------------------
  !        0.   Set up data structures.                ( W3NMOD, etc. )
  !        1.   I-O setup.
  !          a  For shell.
  !          b  For WAVEWATCH III.
  !          c  Local parameters.
  !        2.   Define input fields
  !        3.   Set time frame.
  !        4.   Define output
  !          a  Loop over types, do
  !        +--------------------------------------------------------+
  !        | b    Process standard line                             |
  !        | c    If type 1: fields of mean wave parameters         |
  !        | d    If type 2: point output                           |
  !        | e    If type 3: track output                           |
  !        | f    If type 4: restart files                          |
  !        | g    If type 5: boundary output                        |
  !        | h    If type 6: separated wave fields                  |
  !        | i    If type 7: coupling fields                        |
  !        +--------------------------------------------------------+
  !        5.   Initialzations
  !          a  Wave model.                              ( W3INIT )
  !          b  Read homogeneous field data.
  !          c  Prepare input files.                     ( W3FLDO )
  !          d  Set field times.
  !        6.   If no input fields required, run model in a single
  !             sweep and exit.                          ( W3WAVE )
  !        7.   Run model with input
  !             Do until end time is reached
  !        +--------------------------------------------------------+
  !        | a  Determine next time interval and input fields.      |
  !        |   1  Preparation                                       |
  !        |      Loop over input fields                            |
  !        | +------------------------------------------------------|
  !        | | 2  Check if update is needed                         |
  !        | | 3  Update time and fields                 ( W3FLDG ) |
  !        | |                                           ( W3FLDH ) |
  !        | | 4  Update next ending time                           |
  !        | +------------------------------------------------------|
  !        | b  Run wave model.                          ( W3WAVE ) |
  !        | c  If requested, data assimilation.         ( W3WDAS ) |
  !        | d  Final output if needed.                  ( W3WAVE ) |
  !        | e  Check time                                          |
  !        +--------------------------------------------------------+
  !     ----------------------------------------------------------------
  !
  !  9. Switches :
  !
  !       !/SHRD  Switch for shared / distributed memory architecture.
  !       !/DIST  Id.
  !       !/MPI   Id.
  !
  !       !/MGW   Moving grid wind correction.
  !       !/MGP   Moving grid propagation correction.
  !
  !       !/T     Enable test output.
  !       !/O7    Echo input homogeneous fields.
  !
  !       !/NCO   NCEP NCO modifications for operational implementation.
  !
  !
  ! 10. Source code :
  !
  !/ ------------------------------------------------------------------- /

!  use w3servmd, only : print_memcheck
#ifdef W3_PDLIB
  USE CONSTANTS, ONLY: LPDLIB
#endif
  USE W3GDATMD
  USE W3WDATMD, ONLY: TIME, VA, W3NDAT, W3DIMW, W3SETW
#ifdef W3_OASIS
  USE W3WDATMD, ONLY: TIME00, TIMEEND
#endif
#ifdef W3_NL5
  USE W3WDATMD, ONLY: QI5TBEG
#endif
  USE W3ADATMD, ONLY: W3NAUX, W3DIMA, W3SETA
  USE W3IDATMD
#ifdef W3_OASIS
  USE W3ODATMD, ONLY: DTOUT, FLOUT
#endif
  USE W3ODATMD, ONLY: W3NOUT, W3SETO
  USE W3ODATMD, ONLY: NAPROC, IAPROC, NAPOUT, NAPERR, NOGRP,      &
       NGRPP, IDOUT, FNMPRE, IOSTYP, NOTYPE
!  USE W3ODATMD, ONLY: NAPROC, IAPROC, NAPOUT, NAPERR, NOGRP,      &
!       NGRPP, IDOUT, FNMPRE, FNMGRD, FNMPNT, FNMRST, IOSTYP, NOTYPE
  USE W3ODATMD, ONLY: FLOGRR, FLOGR, OFILES
  !/
  USE W3FLDSMD
  USE W3INITMD
  USE W3WAVEMD
  USE W3WDASMD
  !/
  USE W3IOGRMD, ONLY: W3IOGR
  USE W3IOGOMD, ONLY: W3READFLGRD, FLDOUT, W3FLGRDFLAG
  USE W3IORSMD, ONLY: OARST
  USE W3IOPOMD
  USE W3SERVMD, ONLY : NEXTLN, EXTCDE
  USE W3TIMEMD

#ifdef W3_OASIS
  USE W3OACPMD, ONLY: CPL_OASIS_INIT, CPL_OASIS_GRID,            &
       CPL_OASIS_DEFINE, CPL_OASIS_FINALIZE,      &
       ID_OASIS_TIME, CPLT0
#endif
#ifdef W3_OASOCM
  USE W3OGCMMD, ONLY: SND_FIELDS_TO_OCEAN
#endif
#ifdef W3_OASACM
  USE W3AGCMMD, ONLY: SND_FIELDS_TO_ATMOS
#endif
#ifdef W3_OASICM
  USE W3IGCMMD, ONLY: SND_FIELDS_TO_ICE
#endif

#ifdef W3_TIDE
  USE W3TIDEMD
#endif
  !
  USE W3NMLSHELMD

#ifdef W3_OMPG
  USE OMP_LIB
#endif
  IMPLICIT NONE
  !
#ifdef W3_MPI
  INCLUDE "mpif.h"
#endif
  !/
  !/ ------------------------------------------------------------------- /
  !/ Local PARAMETER statements
  !/
  INTEGER, PARAMETER  :: NHMAX =    200
  !/
  !/ ------------------------------------------------------------------- /
  !/ Local parameters
  !/
  TYPE(NML_DOMAIN_T)       :: NML_DOMAIN
  TYPE(NML_INPUT_T)        :: NML_INPUT
  TYPE(NML_OUTPUT_TYPE_T)  :: NML_OUTPUT_TYPE
  TYPE(NML_OUTPUT_DATE_T)  :: NML_OUTPUT_DATE
!  TYPE(NML_OUTPUT_PATH_T)  :: NML_OUTPUT_PATH
  TYPE(NML_HOMOG_COUNT_T)  :: NML_HOMOG_COUNT
  TYPE(NML_HOMOG_INPUT_T), ALLOCATABLE  :: NML_HOMOG_INPUT(:)
  !
  INTEGER             :: NDSI, NDSI2, NDSS, NDSO, NDSE, NDST, NDSL,&
       NDSEN, IERR, J, I, ILOOP, IPTS, NPTS,     &
       NDTNEW, MPI_COMM = -99,                   &
       FLAGTIDE, COUPL_COMM, IH, N_TOT
  INTEGER             :: NDSF(-7:9), NDS(15), NTRACE(2), NDT(7:9), &
       TIME0(2), TIMEN(2), TTIME(2), TTT(2),     &
       NH(-7:10), THO(2,-7:10,NHMAX), RCLD(7:9), &
       NODATA(7:9), ODAT(40), IPRT(6) = 0,       &
       STARTDATE(8), STOPDATE(8), IHH(-7:10)
  !
#ifdef W3_OASIS
  INTEGER             :: OASISED
#endif
#ifdef W3_COU
  INTEGER             :: OFL
#endif
  INTEGER             :: CLKDT1(8), CLKDT2(8), CLKDT3(8)
#ifdef W3_MPI
  INTEGER             :: IERR_MPI
#endif
  !
  REAL                :: FACTOR, DTTST, XX, YY,                    &
       HA(NHMAX,-7:10), HD(NHMAX,-7:10),         &
       HS(NHMAX,-7:10)
  REAL                :: CLKFIN, CLKFEL
  REAL, ALLOCATABLE   :: X(:), Y(:), XXX(:,:), DATA0(:,:),         &
       DATA1(:,:), DATA2(:,:)
  !
  DOUBLE PRECISION    :: STARTJULDAY, STOPJULDAY
  !
  CHARACTER(LEN=1)    :: COMSTR, FLAGTFC(-7:10)
  CHARACTER(LEN=3)    :: IDSTR(-7:10), IDTST
  CHARACTER(LEN=6)    :: YESXNO
  CHARACTER(LEN=40)   :: PN
  CHARACTER(LEN=40),                                               &
       ALLOCATABLE :: PNAMES(:)
  CHARACTER(LEN=13)   :: IDFLDS(-7:10)
  CHARACTER(LEN=20)   :: STRNG
  CHARACTER(LEN=23)   :: DTME21
  CHARACTER(LEN=30)   :: IDOTYP(8)
  CHARACTER(LEN=80)   :: LINE
  CHARACTER(LEN=256)  :: TMPLINE, TEST
  CHARACTER(LEN=1024) :: FLDIN
  CHARACTER(LEN=1024) :: FLDRST=''
  CHARACTER(LEN=80)   :: LINEIN
  CHARACTER(LEN=8)    :: WORDS(7)=''

#ifdef W3_COU
  CHARACTER(LEN=30)   :: OFILE
#endif
  !
  LOGICAL             :: FLLSTL, FLLSTI, FLLSTR, FLFLG, FLHOM,     &
       TFLAGI, PRTFRM, FLAGSCI, FLGNML
  LOGICAL             :: FLGRD(NOGRP,NGRPP), FLGD(NOGRP),          &
       FLGR2(NOGRP,NGRPP), FLG2(NOGRP),          &
       FLAGSTIDE(4), FLH(-7:10), FLGDAS(3),      &
       FLLST_ALL(-7:10)
#ifdef W3_MPI
  LOGICAL             :: FLHYBR = .FALSE.
#endif
#ifdef W3_OMPH
  INTEGER             :: THRLEV
#endif
#ifdef W3_OASIS
  LOGICAL             :: L_MASTER
  LOGICAL             :: FIRST_STEP = .TRUE.
#endif
  character(len=10)   :: jchar
  integer             :: memunit

  LOGICAL                 :: DIR_EXISTS
  INTEGER                 :: DIR_STATUS
  
  DATA IDFLDS / 'ice param. 1 ' , 'ice param. 2 ' ,               &
       'ice param. 3 ' , 'ice param. 4 ' ,               &
       'ice param. 5 ' ,                                 &
       'mud density  ' , 'mud thkness  ' ,               &
       'mud viscos.  ' ,                                 &
       'water levels ' , 'currents     ' ,               &
       'winds        ' , 'ice fields   ' ,               &
       'momentum     ' , 'air density  ' ,               &
       'mean param.  ' , '1D spectra   ' ,               &
       '2D spectra   ' , 'moving grid  ' /
  DATA IDOTYP / 'Fields of mean wave parameters' ,                &
       'Point output                  ' ,                &
       'Track point output            ' ,                &
       'Restart files                 ' ,                &
       'Nesting data                  ' ,                &
       'Partitioned wave field data   ' ,                &
       'Fields for coupling           ' ,                &
       'Restart files second request  '/
  DATA IDSTR  / 'IC1', 'IC2', 'IC3', 'IC4', 'IC5', 'MDN', 'MTH',  &
       'MVS', 'LEV', 'CUR', 'WND', 'ICE', 'TAU', 'RHO',  &
       'DT0', 'DT1', 'DT2', 'MOV' /
  !
  PUBLIC
  CONTAINS
  
  SUBROUTINE FISOC_WW3_INIT(MPI_COMM_opt)
      IMPLICIT NONE
      INTEGER, INTENT(IN),OPTIONAL  :: MPI_COMM_opt
  !
  !/
  !/ ------------------------------------------------------------------- /
  !/
  FLGR2 = .FALSE.
  FLAGSTIDE(:) = .FALSE.
  FLH(:)       = .FALSE.
  !
#ifdef W3_T
  PRTFRM = .TRUE.
#endif
  !
  CALL DATE_AND_TIME ( VALUES=CLKDT1 )
  !
  !--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! 0.  Set up data structures
  !
#ifdef W3_OASIS
  OASISED=1
#endif
#ifdef W3_PDLIB
  LPDLIB = .TRUE.
#endif
  !
  CALL W3NMOD ( 1, 6, 6 )
  CALL W3NDAT (    6, 6 )
  CALL W3NAUX (    6, 6 )
  CALL W3NOUT (    6, 6 )
  CALL W3NINP (    6, 6 )
  !
  CALL W3SETG ( 1, 6, 6 )
  CALL W3SETW ( 1, 6, 6 )
  CALL W3SETA ( 1, 6, 6 )
  CALL W3SETO ( 1, 6, 6 )
  CALL W3SETI ( 1, 6, 6 )

  memunit = 740+IAPROC
!  call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 1')
  !
#ifdef W3_SHRD
  NAPROC = 1
  IAPROC = 1
#endif
  !
#ifdef W3_OMPH
  FLHYBR = .TRUE.
#endif

#ifdef W3_OASIS
!  IF (OASISED.EQ.1) THEN
!    CALL CPL_OASIS_INIT(MPI_COMM)
!  ELSE
#endif
#ifdef W3_OMPH
    ! For hybrid MPI-OpenMP specify required thread level. JGLi06Sep2019
!    IF( FLHYBR ) THEN
!      CALL MPI_INIT_THREAD( MPI_THREAD_FUNNELED, THRLEV, IERR_MPI)
!    ELSE
#endif
#ifdef W3_MPI
!      CALL MPI_INIT      ( IERR_MPI )
#endif
#ifdef W3_OMPH
!    ENDIF
#endif

#ifdef W3_MPI
    MPI_COMM = MPI_COMM_opt
#endif
#ifdef W3_OASIS
!  END IF
#endif
  !
  !
#ifdef W3_MPI
  CALL MPI_COMM_SIZE ( MPI_COMM, NAPROC, IERR_MPI )
#endif
#ifdef W3_MPI
  CALL MPI_COMM_RANK ( MPI_COMM, IAPROC, IERR_MPI )
  IAPROC = IAPROC + 1
#endif
  memunit = 740+IAPROC
  !
#ifdef W3_NCO
  !     IF ( IAPROC .EQ. 1 ) CALL W3TAGB                         &
  !                         ('WAVEFCST',1998,0007,0050,'NP21   ')
#endif
!  call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 2')
  !
  !--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! 1.  IO set-up
  ! 1.a For shell
  !
  NDSI   = 10
  NDSS   = 90
  NDSO   =  6
  NDSE   =  6
  NDST   =  6
  NDSL   = 50
#ifdef W3_COU
  NDSO   =  333
  NDSE   =  333
  NDST   =  333
#endif


  NDSF(-7)  = 1008
  NDSF(-6)  = 1009
  NDSF(-5)  = 1010
  NDSF(-4)  = 1011
  NDSF(-3)  = 1012
  NDSF(-2)  = 1013
  NDSF(-1)  = 1014
  NDSF(0)   = 1015

  NDSF(1)  = 11
  NDSF(2)  = 12
  NDSF(3)  = 13
  NDSF(4)  = 14
  NDSF(5)  = 15
  NDSF(6)  = 16
  NDSF(7)  = 17
  NDSF(8)  = 18
  NDSF(9)  = 19
  !
#ifdef W3_NCO
  !
  ! Redo according to NCO
  !
  NDSI   = 11
  NDSS   = 90
  NDSO   =  6
  NDSE   = NDSO
  NDST   = NDSO
  NDSF(1)  = 12
  NDSF(2)  = 13
  NDSF(3)  = 14
  NDSF(4)  = 15
  NDSF(5)  = 16
  NDSF(6)  = 17
  NDSF(7)  = 18
  NDSF(8)  = 19
  NDSF(9)  = 20
#endif
  !
  NAPOUT = 1
  NAPERR = 1
  !
#ifdef W3_COU
  OFILE  = 'output.ww3'
  OFL    = LEN_TRIM(OFILE)
  J      = LEN_TRIM(FNMPRE)
  IF ( IAPROC .EQ. NAPOUT )             &
       OPEN (333,FILE=FNMPRE(:J)//OFILE(:OFL),ERR=2008,IOSTAT=IERR)
#endif

  IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,900)
  !
  IF ( IAPROC .EQ. NAPERR ) THEN
    NDSEN  = NDSE
  ELSE
    NDSEN  = -1
  END IF
#ifdef W3_OMPH
  IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,905) &
       MPI_THREAD_FUNNELED, THRLEV
#endif
  !
#ifdef W3_OMPG
    IF(IAPROC .EQ. NAPOUT) THEN
      WRITE(NDSO, 906) omp_get_max_threads()
    ENDIF
#endif


  ! 1.b For WAVEWATCH III (See W3INIT)
  !
  NDS( 1) = 20
  NDS( 2) =  6
  NDS( 3) = 21
  NDS( 4) =  6
  NDS( 5) = 30
  NDS( 6) = 30
  NDS( 7) = 31
  NDS( 8) = 32
  NDS( 9) = 33
  NDS(10) = 35
  NDS(11) = 22
  NDS(12) = 23
  NDS(13) = 34
  NDS(14) = 36
  NDS(15) = 37

  !
  NTRACE(1) =  NDS(3)
  NTRACE(2) =  10
  !
#ifdef W3_NCO
  !
  ! Redo according to NCO
  !
  NDS( 1) = 51
  NDS( 2) = NDSO
  NDS( 3) = NDSO
  NDS( 4) = NDSO
  NDS( 5) = 20
  NDS( 6) = 21
  NDS( 7) = 52
  NDS( 8) = 53
  NDS( 9) = 22
  NDS(10) = 71
  NDS(11) = 23
  NDS(12) = 54
  NDS(13) = 55
  NTRACE(1) = NDSO
#endif
  !
#ifdef W3_T
  WRITE (NDST,9000) (NDS(I),I=1,12)
  WRITE (NDST,9001) (NTRACE(I),I=1,2)
#endif
  !
  ! 1.c Local parameters
  !
  ! Default COMSTR to "$" (for when using nml input files)
  COMSTR = "$"
  !
  ! inferred from context: these flags (FL) are to indicate that the last (LST)
  !   field has been read from a file.
  FLLSTL = .FALSE. ! This is associated with J.EQ.1 (wlev)
  FLLSTI = .FALSE. ! This is associated with J.EQ.4 (ice)
  FLLSTR = .FALSE. ! This is associated with J.EQ.6 (rhoa)
  FLLST_ALL = .FALSE. ! For all

  ! If using experimental mud or ice physics, additional lines will
  !  be read in from ww3_shel.inp and applied, so JFIRST is changed from
  !  its initialization setting "JFIRST=1" to some lower value.
#ifdef W3_IC1
  JFIRST=-7
#endif
#ifdef W3_IC2
  JFIRST=-7
#endif
#ifdef W3_IS2
  JFIRST=-7
#endif
#ifdef W3_IC3
  JFIRST=-7
#endif
#ifdef W3_BT8
  JFIRST=-7
#endif
#ifdef W3_BT9
  JFIRST=-7
#endif
#ifdef W3_IC4
  JFIRST=-7
#endif
#ifdef W3_IC5
  JFIRST=-7
#endif

!  call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 2a')
  !
  !--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! 2.  Define input fields
  !

  !
  ! process ww3_prnc namelist
  !
  INQUIRE(FILE=TRIM(FNMPRE)//"ww3_shel.nml", EXIST=FLGNML)
  IF (FLGNML) THEN
    ! Read namelist
    CALL W3NMLSHEL (MPI_COMM, NDSI, TRIM(FNMPRE)//'ww3_shel.nml',  &
         NML_DOMAIN, NML_INPUT, NML_OUTPUT_TYPE,                   &
         NML_OUTPUT_DATE, NML_HOMOG_COUNT,        &
         NML_HOMOG_INPUT, IERR)

    ! 2.1 forcing flags

    FLH(-7:10)=.FALSE.
    FLAGTFC(-7)=TRIM(NML_INPUT%FORCING%ICE_PARAM1)
    FLAGTFC(-6)=TRIM(NML_INPUT%FORCING%ICE_PARAM2)
    FLAGTFC(-5)=TRIM(NML_INPUT%FORCING%ICE_PARAM3)
    FLAGTFC(-4)=TRIM(NML_INPUT%FORCING%ICE_PARAM4)
    FLAGTFC(-3)=TRIM(NML_INPUT%FORCING%ICE_PARAM5)
    FLAGTFC(-2)=TRIM(NML_INPUT%FORCING%MUD_DENSITY)
    FLAGTFC(-1)=TRIM(NML_INPUT%FORCING%MUD_THICKNESS)
    FLAGTFC(0)=TRIM(NML_INPUT%FORCING%MUD_VISCOSITY)
    FLAGTFC(1)=TRIM(NML_INPUT%FORCING%WATER_LEVELS)
    FLAGTFC(2)=TRIM(NML_INPUT%FORCING%CURRENTS)
    FLAGTFC(3)=TRIM(NML_INPUT%FORCING%WINDS)
    FLAGTFC(4)=TRIM(NML_INPUT%FORCING%ICE_CONC)
    FLAGTFC(5)=TRIM(NML_INPUT%FORCING%ATM_MOMENTUM)
    FLAGTFC(6)=TRIM(NML_INPUT%FORCING%AIR_DENSITY)
    FLAGTFC(7)=TRIM(NML_INPUT%ASSIM%MEAN)
    FLAGTFC(8)=TRIM(NML_INPUT%ASSIM%SPEC1D)
    FLAGTFC(9)=TRIM(NML_INPUT%ASSIM%SPEC2D)

    IF (TRIM(NML_INPUT%FORCING%ICE_PARAM1) .EQ. 'H') THEN
      FLAGTFC(-7)='T'
      FLH(-7)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%ICE_PARAM2) .EQ. 'H') THEN
      FLAGTFC(-6)='T'
      FLH(-6)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%ICE_PARAM3) .EQ. 'H') THEN
      FLAGTFC(-5)='T'
      FLH(-5)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%ICE_PARAM4) .EQ. 'H') THEN
      FLAGTFC(-4)='T'
      FLH(-4)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%ICE_PARAM5) .EQ. 'H') THEN
      FLAGTFC(-3)='T'
      FLH(-3)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%MUD_DENSITY) .EQ. 'H') THEN
      FLAGTFC(-2)='T'
      FLH(-2)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%MUD_THICKNESS) .EQ. 'H') THEN
      FLAGTFC(-1)='T'
      FLH(-1)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%MUD_VISCOSITY) .EQ. 'H') THEN
      FLAGTFC(0)='T'
      FLH(0)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%WATER_LEVELS) .EQ. 'H') THEN
      FLAGTFC(1)='T'
      FLH(1)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%CURRENTS) .EQ. 'H') THEN
      FLAGTFC(2)='T'
      FLH(2)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%WINDS) .EQ. 'H') THEN
      FLAGTFC(3)='T'
      FLH(3)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%ICE_CONC) .EQ. 'H') THEN
      FLAGTFC(4)='T'
      FLH(4)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%ATM_MOMENTUM) .EQ. 'H') THEN
      FLAGTFC(5)='T'
      FLH(5)=.TRUE.
    END IF
    IF (TRIM(NML_INPUT%FORCING%AIR_DENSITY) .EQ. 'H') THEN
      FLAGTFC(6)='T'
      FLH(6)=.TRUE.
    END IF
!JQI<
    FLAGTFC(1)='C'
    FLH(1)=.FALSE.       !.TRUE.
    FLAGTFC(2)='C'
    FLH(2)=.FALSE.       !.TRUE.
!JQI>
    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,920)
    DO J=JFIRST, 9
      IF (FLAGTFC(J).EQ.'T') THEN
        INFLAGS1(J)=.TRUE.
        FLAGSC(J)=.FALSE.
      END IF
      IF (FLAGTFC(J).EQ.'F') THEN
        INFLAGS1(J)=.FALSE.
        FLAGSC(J)=.FALSE.
      END IF
      IF (FLAGTFC(J).EQ.'C') THEN
        INFLAGS1(J)=.TRUE.
        FLAGSC(J)=.TRUE.
      END IF
      IF ( J .LE. 6 ) THEN
        FLH(J) = FLH(J) .AND. INFLAGS1(J)
      END IF
      IF ( INFLAGS1(J) ) THEN
        YESXNO = 'YES/--'
      ELSE
        YESXNO = '---/NO'
      END IF
      IF ( FLH(J) ) THEN
        STRNG  = '(homogeneous field) '
      ELSE IF ( FLAGSC(J) ) THEN
        STRNG  = '(coupling field) '
      ELSE
        STRNG  = '                    '
      END IF
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,921) IDFLDS(J), YESXNO, STRNG
    END DO
#ifdef W3_COU
    IF (FLAGSC(1) .AND. INFLAGS1(2) .AND. .NOT. FLAGSC(2)) GOTO 2102
    IF (FLAGSC(2) .AND. INFLAGS1(1) .AND. .NOT. FLAGSC(1)) GOTO 2102
#endif

    INFLAGS1(10) = .FALSE.
#ifdef W3_MGW
    INFLAGS1(10) = .TRUE.
#endif
#ifdef W3_MGP
    INFLAGS1(10) = .TRUE.
#endif
#ifdef W3_MGW
    FLH(10)   = .TRUE.
#endif
#ifdef W3_MGP
    FLH(10)   = .TRUE.
#endif
    IF ( INFLAGS1(10) .AND. IAPROC.EQ.NAPOUT )                         &
         WRITE (NDSO,921) IDFLDS(10), 'YES/--', ' '
    !
    FLFLG  = INFLAGS1(-7) .OR. INFLAGS1(-6) .OR. INFLAGS1(-5) .OR. INFLAGS1(-4) &
         .OR. INFLAGS1(-3) .OR. INFLAGS1(-2) .OR. INFLAGS1(-1)           &
         .OR. INFLAGS1(0)  .OR. INFLAGS1(1)  .OR. INFLAGS1(2)            &
         .OR. INFLAGS1(3)  .OR. INFLAGS1(4)  .OR. INFLAGS1(5)            &
         .OR. INFLAGS1(6)  .OR. INFLAGS1(7)  .OR. INFLAGS1(8)            &
         .OR. INFLAGS1(9)
    FLHOM  = FLH(-7) .OR. FLH(-6) .OR. FLH(-5) .OR. FLH(-4)       &
         .OR. FLH(-3) .OR. FLH(-2) .OR. FLH(-1) .OR. FLH(0)   &
         .OR. FLH(1) .OR. FLH(2) .OR. FLH(3) .OR. FLH(4)      &
         .OR. FLH(5) .OR. FLH(6) .OR. FLH(10)
    !
    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,922)
    !
    !       INFLAGS2 is just "initial value of INFLAGS1", i.e. does *not* get
    !          changed when model reads last record of ice.ww3
    INFLAGS2=INFLAGS1

#ifdef W3_T
    WRITE (NDST,9020) FLFLG, INFLAGS1, FLHOM, FLH
#endif

    ! 2.2 Time setup

    READ(NML_DOMAIN%START,*) TIME0
    CALL T2D(TIME0,STARTDATE,IERR)
    CALL D2J(STARTDATE,STARTJULDAY,IERR)
    READ(NML_DOMAIN%STOP,*) TIMEN
    CALL T2D(TIMEN,STOPDATE,IERR)
    CALL D2J(STOPDATE,STOPJULDAY,IERR)

    ! 2.3 Domain setup

    IOSTYP = NML_DOMAIN%IOSTYP

#ifdef W3_PDLIB
    IF (IOSTYP .gt. 1) THEN
      WRITE(*,*) 'IOSTYP not supported in domain decomposition mode'
      CALL EXTCDE ( 6666 )
    ENDIF
#endif

    CALL W3IOGR ( 'GRID', NDSF(7) )
    IF ( FLAGLL ) THEN
      FACTOR = 1.
    ELSE
      FACTOR = 1.E-3
    END IF

    ! 2.4 Output dates

    READ(NML_OUTPUT_DATE%FIELD%START, *)   ODAT(1), ODAT(2)
    READ(NML_OUTPUT_DATE%FIELD%STRIDE, *)  ODAT(3)
    READ(NML_OUTPUT_DATE%FIELD%STOP, *)    ODAT(4), ODAT(5)

    READ(NML_OUTPUT_DATE%FIELD%OUTFFILE, *)  OFILES(1)
    !        OUTPTS(I)%OUTSTRIDE(1)=ODAT(3,I)

    READ(NML_OUTPUT_DATE%POINT%START, *)   ODAT(6), ODAT(7)
    READ(NML_OUTPUT_DATE%POINT%STRIDE, *)  ODAT(8)
    READ(NML_OUTPUT_DATE%POINT%STOP, *)    ODAT(9), ODAT(10)

    READ(NML_OUTPUT_DATE%POINT%OUTFFILE, *)  OFILES(2)
    !        OUTPTS(I)%OUTSTRIDE(2)=ODAT(8,I)

    READ(NML_OUTPUT_DATE%TRACK%START, *)   ODAT(11), ODAT(12)
    READ(NML_OUTPUT_DATE%TRACK%STRIDE, *)  ODAT(13)
    READ(NML_OUTPUT_DATE%TRACK%STOP, *)    ODAT(14), ODAT(15)
    READ(NML_OUTPUT_DATE%RESTART%START, *)   ODAT(16), ODAT(17)
    READ(NML_OUTPUT_DATE%RESTART%STRIDE, *)  ODAT(18)
    READ(NML_OUTPUT_DATE%RESTART%STOP, *)    ODAT(19), ODAT(20)
    READ(NML_OUTPUT_DATE%RESTART2%START, *)   ODAT(36), ODAT(37)
    READ(NML_OUTPUT_DATE%RESTART2%STRIDE, *)  ODAT(38)
    READ(NML_OUTPUT_DATE%RESTART2%STOP, *)    ODAT(39), ODAT(40)
    READ(NML_OUTPUT_DATE%BOUNDARY%START, *)   ODAT(21), ODAT(22)
    READ(NML_OUTPUT_DATE%BOUNDARY%STRIDE, *)  ODAT(23)
    READ(NML_OUTPUT_DATE%BOUNDARY%STOP, *)    ODAT(24), ODAT(25)
    READ(NML_OUTPUT_DATE%PARTITION%START, *)   ODAT(26), ODAT(27)
    READ(NML_OUTPUT_DATE%PARTITION%STRIDE, *)  ODAT(28)
    READ(NML_OUTPUT_DATE%PARTITION%STOP, *)    ODAT(29), ODAT(30)
    READ(NML_OUTPUT_DATE%COUPLING%START, *)   ODAT(31), ODAT(32)
    READ(NML_OUTPUT_DATE%COUPLING%STRIDE, *)  ODAT(33)
    READ(NML_OUTPUT_DATE%COUPLING%STOP, *)    ODAT(34), ODAT(35)

    ! set the time stride at 0 or more
    ODAT(3) = MAX ( 0 , ODAT(3) )
    ODAT(8) = MAX ( 0 , ODAT(8) )
    ODAT(13) = MAX ( 0 , ODAT(13) )
    ODAT(18) = MAX ( 0 , ODAT(18) )
    ODAT(23) = MAX ( 0 , ODAT(23) )
    ODAT(28) = MAX ( 0 , ODAT(28) )
    ODAT(33) = MAX ( 0 , ODAT(33) )
    ODAT(38) = MAX ( 0 , ODAT(38) )
    !
#ifdef W3_COU
    ! Test the validity of the coupling time step
    IF (ODAT(33) == 0) THEN
      IF ( IAPROC .EQ. NAPOUT ) THEN
        WRITE(NDSO,1010) ODAT(33), INT(DTMAX)
      END IF
      ODAT(33) = INT(DTMAX)
    ELSE IF (MOD(ODAT(33),INT(DTMAX)) .NE. 0) THEN
      GOTO 2009
    END IF
#endif
    !
    ! 2.5 Output types

    NPTS   = 0
    NOTYPE = 6
#ifdef W3_COU
    NOTYPE = 7
#endif
    DO J = 1, NOTYPE
      !          OUTPTS(I)%OFILES(J)=OFILES(J)
      IF ( ODAT(5*(J-1)+3) .NE. 0 ) THEN

        ! Type 1: fields of mean wave parameters
        IF ( J .EQ. 1 ) THEN
          FLDOUT = NML_OUTPUT_TYPE%FIELD%LIST
          CALL W3FLGRDFLAG ( NDSO, NDSO, NDSE, FLDOUT, FLGD,     &
               FLGRD, IAPROC, NAPOUT, IERR )
          IF ( IERR .NE. 0 ) GOTO 2222


          ! Type 2: point output
        ELSE IF ( J .EQ. 2 ) THEN
          OPEN (NDSL, FILE=TRIM(FNMPRE)//TRIM(NML_OUTPUT_TYPE%POINT%FILE), &
               FORM='FORMATTED', STATUS='OLD', ERR=2104, IOSTAT=IERR)

          ! first loop to count the number of points
          ! second loop to allocate the array and store the points
          IPTS = 0
          DO ILOOP=1,2
            REWIND (NDSL)
            !
            IF ( ILOOP.EQ.2) THEN
              NPTS = IPTS
              IF ( NPTS.GT.0 ) THEN
                ALLOCATE ( X(NPTS), Y(NPTS), PNAMES(NPTS) )
                IPTS = 0 ! reset counter to be reused for next do loop
              ELSE
                ALLOCATE ( X(1), Y(1), PNAMES(1) )
                GOTO 2054
              END IF
            END IF
            !
            DO
              READ (NDSL,*,ERR=2004,IOSTAT=IERR) TMPLINE
              ! if end of file or stopstring, then exit
              IF ( IERR.NE.0 .OR. INDEX(TMPLINE,"STOPSTRING").NE.0 ) EXIT
              ! leading blanks removed and placed on the right
              TEST = ADJUSTL ( TMPLINE )
              IF ( TEST(1:1).EQ.COMSTR .OR. LEN_TRIM(TEST).EQ.0 ) THEN
                ! if comment or blank line, then skip
                CYCLE
              ELSE
                ! otherwise, backup to beginning of line
                BACKSPACE ( NDSL, ERR=2004, IOSTAT=IERR)
                READ (NDSL,*,ERR=2004,IOSTAT=IERR) XX, YY, PN
              END IF
              IPTS = IPTS + 1
              IF ( ILOOP .EQ. 1 ) CYCLE
              IF ( ILOOP .EQ. 2 ) THEN
                X(IPTS)      = XX
                Y(IPTS)      = YY
                PNAMES(IPTS) = PN
                IF ( IAPROC .EQ. NAPOUT ) THEN
                  IF ( FLAGLL ) THEN
                    IF ( IPTS .EQ. 1 ) THEN
                      WRITE (NDSO,2945)                     &
                           FACTOR*XX, FACTOR*YY, PN
                    ELSE
                      WRITE (NDSO,2946) IPTS,               &
                           FACTOR*XX, FACTOR*YY, PN
                    END IF
                  ELSE
                    IF ( IPTS .EQ. 1 ) THEN
                      WRITE (NDSO,2955)                     &
                           FACTOR*XX, FACTOR*YY, PN
                    ELSE
                      WRITE (NDSO,2956) IPTS,               &
                           FACTOR*XX, FACTOR*YY, PN
                    END IF
                  END IF
                END IF
              END IF ! ILOOP.EQ.2
            END DO ! end of file
          END DO ! ILOOP
          CLOSE(NDSL)

          ! Type 3: track output
        ELSE IF ( J .EQ. 3 ) THEN
          TFLAGI = NML_OUTPUT_TYPE%TRACK%FORMAT
          IF ( .NOT. TFLAGI ) NDS(11) = -NDS(11)
          IF ( IAPROC .EQ. NAPOUT ) THEN
            IF ( .NOT. TFLAGI ) THEN
              WRITE (NDSO,3945) 'input', 'UNFORMATTED'
            ELSE
              WRITE (NDSO,3945) 'input', 'FORMATTED'
            END IF
          END IF

          ! Type 6: partitioning
        ELSE IF ( J .EQ. 6 ) THEN
          IPRT(1) = NML_OUTPUT_TYPE%PARTITION%X0
          IPRT(2) = NML_OUTPUT_TYPE%PARTITION%XN
          IPRT(3) = NML_OUTPUT_TYPE%PARTITION%NX
          IPRT(4) = NML_OUTPUT_TYPE%PARTITION%Y0
          IPRT(5) = NML_OUTPUT_TYPE%PARTITION%YN
          IPRT(6) = NML_OUTPUT_TYPE%PARTITION%NY
          PRTFRM = NML_OUTPUT_TYPE%PARTITION%FORMAT
          !
          IF ( IAPROC .EQ. NAPOUT ) THEN
            IF ( PRTFRM ) THEN
              YESXNO = 'YES/--'
            ELSE
              YESXNO = '---/NO'
            END IF
            WRITE (NDSO,6945) IPRT, YESXNO
          END IF

#ifdef W3_COU
          ! Type 7: coupling
        ELSE IF ( J .EQ. 7 ) THEN
          FLDOUT = NML_OUTPUT_TYPE%COUPLING%SENT
          CALL W3FLGRDFLAG ( NDSO, NDSO, NDSE, FLDOUT, FLG2,  &
               FLGR2, IAPROC, NAPOUT, IERR )
          IF ( IERR .NE. 0 ) GOTO 2222
          FLDIN = NML_OUTPUT_TYPE%COUPLING%RECEIVED
          CPLT0 = NML_OUTPUT_TYPE%COUPLING%COUPLET0
#endif

        END IF ! J
      END IF ! ODAT
    END DO ! J

    ! Extra fields to be written in the restart
    FLDRST = NML_OUTPUT_TYPE%RESTART%EXTRA
    CALL W3FLGRDFLAG ( NDSO, NDSO, NDSE, FLDRST, FLOGR,  &
         FLOGRR, IAPROC, NAPOUT, IERR )
    IF ( IERR .NE. 0 ) GOTO 2222

    ! force minimal allocation to avoid memory seg fault
    IF ( .NOT.ALLOCATED(X) .AND. NPTS.EQ.0 ) ALLOCATE ( X(1), Y(1), PNAMES(1) )

    ! 2.6 Homogeneous field data

    IF ( FLHOM ) THEN
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,951)                   &
           'Homogeneous field data (and moving grid) ...'

      NH(-7) = NML_HOMOG_COUNT%N_IC1
      NH(-6) = NML_HOMOG_COUNT%N_IC2
      NH(-5) = NML_HOMOG_COUNT%N_IC3
      NH(-4) = NML_HOMOG_COUNT%N_IC4
      NH(-3) = NML_HOMOG_COUNT%N_IC5
      NH(-2) = NML_HOMOG_COUNT%N_MDN
      NH(-1) = NML_HOMOG_COUNT%N_MTH
      NH(0)  = NML_HOMOG_COUNT%N_MVS
      NH(1)  = NML_HOMOG_COUNT%N_LEV
      NH(2)  = NML_HOMOG_COUNT%N_CUR
      NH(3)  = NML_HOMOG_COUNT%N_WND
      NH(4)  = NML_HOMOG_COUNT%N_ICE
      NH(5)  = NML_HOMOG_COUNT%N_TAU
      NH(6)  = NML_HOMOG_COUNT%N_RHO
      NH(10)  = NML_HOMOG_COUNT%N_MOV
      !
      N_TOT = NML_HOMOG_COUNT%N_TOT
      !
      DO J=JFIRST,10
        IF ( NH(J) .GT. NHMAX ) GOTO 2006
      END DO


      ! Store homogeneous fields
      IF ( N_TOT .GT. 0 ) THEN
        IHH(:)=0
        DO IH=1,N_TOT
          READ(NML_HOMOG_INPUT(IH)%NAME,*) IDTST
          SELECT CASE (IDTST)
          CASE ('IC1')
            J=-7
          CASE ('IC2')
            J=-6
          CASE ('IC3')
            J=-5
          CASE ('IC4')
            J=-4
          CASE ('IC5')
            J=-3
          CASE ('MDN')
            J=-2
          CASE ('MTH')
            J=-1
          CASE ('MVS')
            J=0
          CASE ('LEV')
            J=1
          CASE ('CUR')
            J=2
          CASE ('WND')
            J=3
          CASE ('ICE')
            J=4
          CASE ('TAU')
            J=5
          CASE ('RHO')
            J=6
          CASE ('MOV')
            J=10
          CASE DEFAULT
            GOTO 2062
          END SELECT
          IHH(J)=IHH(J)+1
          READ(NML_HOMOG_INPUT(IH)%DATE,*) THO(:,J,IHH(J))
          HA(IHH(J),J) = NML_HOMOG_INPUT(IH)%VALUE1
          HD(IHH(J),J) = NML_HOMOG_INPUT(IH)%VALUE2
          HS(IHH(J),J) = NML_HOMOG_INPUT(IH)%VALUE3
        END DO
      END IF

#ifdef W3_O7
      DO J=JFIRST, 10
        IF ( FLH(J) .AND. IAPROC.EQ.NAPOUT ) THEN
          WRITE (NDSO,952) NH(J), IDFLDS(J)
          DO I=1, NH(J)
            IF ( ( J .LE. 1 ) .OR. ( J .EQ. 4 ) .OR.      &
                 ( J .EQ. 6 ) ) THEN
              WRITE (NDSO,953) I, THO(1,J,I), THO(2,J,I), &
                   HA(I,J)
            ELSE IF ( ( J .EQ. 2 ) .OR. ( J .EQ. 5 ) .OR. &
                 ( J .EQ. 10 ) ) THEN
              WRITE (NDSO,953) I, THO(1,J,I), THO(2,J,I), &
                   HA(I,J), HD(I,J)
            ELSE IF ( J .EQ. 3 ) THEN
              WRITE (NDSO,953) I, THO(1,J,I), THO(2,J,I), &
                   HA(I,J), HD(I,J), HS(I,J)
            END IF
          END DO
        END IF
      END DO
#endif
      !
      IF ( ( FLH(-7) .AND. (NH(-7).EQ.0) ) .OR.                     &
           ( FLH(-6) .AND. (NH(-6).EQ.0) ) .OR.                     &
           ( FLH(-5) .AND. (NH(-5).EQ.0) ) .OR.                     &
           ( FLH(-4) .AND. (NH(-4).EQ.0) ) .OR.                     &
           ( FLH(-3) .AND. (NH(-3).EQ.0) ) .OR.                     &
           ( FLH(-2) .AND. (NH(-2).EQ.0) ) .OR.                     &
           ( FLH(-1) .AND. (NH(-1).EQ.0) ) .OR.                     &
           ( FLH(0)  .AND. (NH(0).EQ.0)  ) .OR.                     &
           ( FLH(1)  .AND. (NH(1).EQ.0)  ) .OR.                     &
           ( FLH(2)  .AND. (NH(2).EQ.0)  ) .OR.                     &
           ( FLH(3)  .AND. (NH(3).EQ.0)  ) .OR.                     &
           ( FLH(4)  .AND. (NH(4).EQ.0)  ) .OR.                     &
           ( FLH(5)  .AND. (NH(5).EQ.0)  ) .OR.                     &
           ( FLH(6)  .AND. (NH(6).EQ.0)  ) .OR.                     &
           ( FLH(10) .AND. (NH(10).EQ.0) ) ) GOTO 2007
      !
    END IF ! FLHOM

    ! USER DEFINED OUTPUT PATH FROM NAMELIST
    ! '/' IS NOT REQUIRED AT THE END OF USER-DEFINED DIRECTORY
!    FNMGRD = TRIM(NML_OUTPUT_PATH%GRD_OUT)
!    IF (FNMGRD(LEN_TRIM(FNMGRD):LEN_TRIM(FNMGRD)) /= '/') THEN
!      FNMGRD = TRIM(FNMGRD) // '/'
!    END IF

!    FNMPNT = TRIM(NML_OUTPUT_PATH%PNT_OUT)
!    IF (FNMPNT(LEN_TRIM(FNMPNT):LEN_TRIM(FNMPNT)) /= '/') THEN
!      FNMPNT = TRIM(FNMPNT) // '/'
!    END IF

!    FNMRST = TRIM(NML_OUTPUT_PATH%RST_OUT)
!    IF (FNMRST(LEN_TRIM(FNMRST):LEN_TRIM(FNMRST)) /= '/') THEN
 !     FNMRST = TRIM(FNMRST) // '/'
 !   END IF

  END IF ! FLGNML



  !
  ! process old ww3_shel.inp format
  !
  IF (.NOT. FLGNML) THEN

    OPEN (NDSI,FILE=TRIM(FNMPRE)//'ww3_shel.inp',STATUS='OLD',IOSTAT=IERR)
    REWIND (NDSI)
    !AR: I changed the error handling for err=2002, see commit message ...
    READ (NDSI,'(A)') COMSTR
    IF (COMSTR.EQ.' ') COMSTR = '$'
    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,901) COMSTR

    ! 2.1 forcing flags

    FLH(-7:10) = .FALSE.
    DO J=JFIRST, 9
      CALL NEXTLN ( COMSTR , NDSI , NDSEN )
      IF ( J .LE. 6 ) THEN
        READ (NDSI,*) FLAGTFC(J), FLH(J)
      ELSE
        READ (NDSI,*) FLAGTFC(J)
      END IF
    END DO

    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,920)
    DO J=JFIRST, 9
      IF (FLAGTFC(J).EQ.'T') THEN
        INFLAGS1(J)=.TRUE.
        FLAGSC(J)=.FALSE.
      END IF
      IF (FLAGTFC(J).EQ.'F') THEN
        INFLAGS1(J)=.FALSE.
        FLAGSC(J)=.FALSE.
      END IF
      IF (FLAGTFC(J).EQ.'C') THEN
        INFLAGS1(J)=.TRUE.
        FLAGSC(J)=.TRUE.
      END IF
      IF ( J .LE. 6 ) THEN
        FLH(J) = FLH(J) .AND. INFLAGS1(J)
      END IF
      IF ( INFLAGS1(J) ) THEN
        YESXNO = 'YES/--'
      ELSE
        YESXNO = '---/NO'
      END IF
      IF ( FLH(J) ) THEN
        STRNG  = '(homogeneous field) '
      ELSE IF ( FLAGSC(J) ) THEN
        STRNG  = '(coupling field) '
      ELSE
        STRNG  = '                    '
      END IF
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,921) IDFLDS(J), YESXNO, STRNG
    END DO
#ifdef W3_COU
    IF (FLAGSC(1) .AND. INFLAGS1(2) .AND. .NOT. FLAGSC(2)) GOTO 2102
    IF (FLAGSC(2) .AND. INFLAGS1(1) .AND. .NOT. FLAGSC(1)) GOTO 2102
#endif
!    call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 2b')
    !
    INFLAGS1(10) = .FALSE.
#ifdef W3_MGW
    INFLAGS1(10) = .TRUE.
#endif
#ifdef W3_MGP
    INFLAGS1(10) = .TRUE.
#endif
#ifdef W3_MGW
    FLH(10)   = .TRUE.
#endif
#ifdef W3_MGP
    FLH(10)   = .TRUE.
#endif
    IF ( INFLAGS1(10) .AND. IAPROC.EQ.NAPOUT )                         &
         WRITE (NDSO,921) IDFLDS(10), 'YES/--', ' '
    !
    FLFLG  = INFLAGS1(-7) .OR. INFLAGS1(-6) .OR. INFLAGS1(-5) .OR. INFLAGS1(-4) &
         .OR. INFLAGS1(-3) .OR. INFLAGS1(-2) .OR. INFLAGS1(-1)           &
         .OR. INFLAGS1(0)  .OR. INFLAGS1(1)  .OR. INFLAGS1(2)            &
         .OR. INFLAGS1(3)  .OR. INFLAGS1(4)  .OR. INFLAGS1(5)            &
         .OR. INFLAGS1(6)  .OR. INFLAGS1(7)  .OR. INFLAGS1(8)            &
         .OR. INFLAGS1(9)
    FLHOM  = FLH(-7) .OR. FLH(-6) .OR. FLH(-5) .OR. FLH(-4)   &
         .OR. FLH(-3) .OR. FLH(-2) .OR. FLH(-1) .OR. FLH(0)   &
         .OR. FLH(1) .OR. FLH(2) .OR. FLH(3) .OR. FLH(4)      &
         .OR. FLH(5) .OR. FLH(6) .OR. FLH(10)
    !
    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,922)
    !
    !       INFLAGS2 is just "initial value of INFLAGS1", i.e. does *not* get
    !          changed when model reads last record of ice.ww3
    INFLAGS2=INFLAGS1

#ifdef W3_T
    WRITE (NDST,9020) FLFLG, INFLAGS1, FLHOM, FLH
#endif


    ! 2.2 Time setup

    CALL NEXTLN ( COMSTR , NDSI , NDSEN )
    READ (NDSI,*) TIME0
!    call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 2c')

    CALL NEXTLN ( COMSTR , NDSI , NDSEN )
    READ (NDSI,*) TIMEN
    !
!    call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 2d')

    ! 2.3 Domain setup

    CALL NEXTLN ( COMSTR , NDSI , NDSEN )
    READ (NDSI,*) IOSTYP
#ifdef W3_PDLIB
    IF (IOSTYP .gt. 1) THEN
      WRITE(*,*) 'IOSTYP not supported in domain decomposition mode'
      CALL EXTCDE ( 6666 )
    ENDIF
#endif
    CALL W3IOGR ( 'GRID', NDSF(7) )
    IF ( FLAGLL ) THEN
      FACTOR = 1.
    ELSE
      FACTOR = 1.E-3
    END IF


    ! 2.4 Output dates

    NPTS   = 0
    NOTYPE = 6
#ifdef W3_COU
    NOTYPE = 7
#endif
    DO J = 1, NOTYPE
      CALL NEXTLN ( COMSTR , NDSI , NDSEN )
      !
      ! CHECKPOINT
      IF(J .EQ. 4) THEN
        ODAT(38)=0
        WORDS(1:7)=''
        READ (NDSI,'(A)') LINEIN
        READ(LINEIN,*,iostat=ierr) WORDS
        READ(WORDS( 1 ), * ) ODAT(16)
        READ(WORDS( 2 ), * ) ODAT(17)
        READ(WORDS( 3 ), * ) ODAT(18)
        READ(WORDS( 4 ), * ) ODAT(19)
        READ(WORDS( 5 ), * ) ODAT(20)
        IF (WORDS(6) .EQ. 'T') THEN
          CALL NEXTLN ( COMSTR , NDSI , NDSEN )
          READ (NDSI,*,END=2001,ERR=2002)(ODAT(I),I=5*(8-1)+1,5*8)
          WRITE(*,*)(ODAT(I),I=5*(8-1)+1,5*8)
        END IF
        IF (WORDS(7) .EQ. 'T') THEN
          CALL NEXTLN ( COMSTR , NDSI , NDSEN )
          READ (NDSI,'(A)',END=2001,ERR=2002) FLDRST
        END IF
        CALL W3FLGRDFLAG ( NDSO, NDSO, NDSE, FLDRST, FLOGR,  &
             FLOGRR, IAPROC, NAPOUT, IERR )
        IF ( IERR .NE. 0 ) GOTO 2222
      ELSE
        !
        !INLINE NEW VARIABLE TO READ IF PRESENT OFILES(J), IF NOT ==0
        !          READ (NDSI,*) (ODAT(I),I=5*(J-1)+1,5*J)
        !          READ (NDSI,*,IOSTAT=IERR) (ODAT(I),I=5*(J-1)+1,5*J),OFILES(J)
        IF(J .LE. 2) THEN
          WORDS(1:6)=''
          !          READ (NDSI,*,END=2001,ERR=2002)(ODAT(I),I=5*(J-1)+1,5*J),OFILES(J)
          READ (NDSI,'(A)') LINEIN
          READ(LINEIN,*,iostat=ierr) WORDS
          !
          IF(J .EQ. 1) THEN
            READ(WORDS( 1 ), * ) ODAT(1)
            READ(WORDS( 2 ), * ) ODAT(2)
            READ(WORDS( 3 ), * ) ODAT(3)
            READ(WORDS( 4 ), * ) ODAT(4)
            READ(WORDS( 5 ), * ) ODAT(5)
          ELSE
            READ(WORDS( 1 ), * ) ODAT(6)
            READ(WORDS( 2 ), * ) ODAT(7)
            READ(WORDS( 3 ), * ) ODAT(8)
            READ(WORDS( 4 ), * ) ODAT(9)
            READ(WORDS( 5 ), * ) ODAT(10)
          END IF

          IF (WORDS(6) .NE. '0' .AND. WORDS(6) .NE. '1') THEN
            OFILES(J)=0
          ELSE
            READ(WORDS( 6 ), * ) OFILES(J)
          END IF


#ifdef W3_COU
        ELSE IF(J .EQ. 7) THEN
          WORDS(1:6)=''
          READ (NDSI,'(A)') LINEIN
          READ(LINEIN,*,iostat=ierr) WORDS

          READ(WORDS( 1 ), * ) ODAT(31)
          READ(WORDS( 2 ), * ) ODAT(32)
          READ(WORDS( 3 ), * ) ODAT(33)
          READ(WORDS( 4 ), * ) ODAT(34)
          READ(WORDS( 5 ), * ) ODAT(35)

          IF (WORDS(6) .EQ. 'T') THEN
            CPLT0 = .TRUE.
          ELSE
            CPLT0 = .FALSE.
          END IF
#endif
        ELSE
          OFILES(J)=0
          READ (NDSI,*,END=2001,ERR=2002)(ODAT(I),I=5*(J-1)+1,5*J)
        END IF
        !          WRITE(*,*) 'OFILES(J)= ', OFILES(J),J
        !
        ODAT(5*(J-1)+3) = MAX ( 0 , ODAT(5*(J-1)+3) )
        !
        write(jchar, '(i0)') j
!        call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL NOTTYPE '//trim(jchar))
        !
        ! 2.5 Output types

        IF ( ODAT(5*(J-1)+3) .NE. 0 ) THEN

          ! Type 1: fields of mean wave parameters
          IF ( J .EQ. 1 ) THEN
            CALL W3READFLGRD ( NDSI, NDSO, 9, NDSEN, COMSTR, FLGD,   &
                 FLGRD, IAPROC, NAPOUT, IERR )
            IF ( IERR .NE. 0 ) GOTO 2222



            ! Type 2: point output
          ELSE IF ( J .EQ. 2 ) THEN
            DO ILOOP=1,2
              IF ( ILOOP .EQ. 1 ) THEN
                NDSI2  = NDSI
                IF ( IAPROC .EQ. 1 ) OPEN                       &
                     (NDSS,FILE=TRIM(FNMPRE)//'ww3_shel.scratch')
              ELSE
                NDSI2  = NDSS
#ifdef W3_MPI
                CALL MPI_BARRIER (MPI_COMM,IERR_MPI)
#endif
                OPEN (NDSS,FILE=TRIM(FNMPRE)//'ww3_shel.scratch')
                REWIND (NDSS)
                !
                IF ( .NOT.ALLOCATED(X) ) THEN
                  IF ( NPTS.GT.0 ) THEN
                    ALLOCATE ( X(NPTS), Y(NPTS), PNAMES(NPTS) )
                  ELSE
                    ALLOCATE ( X(1), Y(1), PNAMES(1) )
                    GOTO 2054
                  END IF
                END IF
              END IF
              !
              NPTS   = 0
              DO
                CALL NEXTLN ( COMSTR , NDSI , NDSEN )
                READ (NDSI2,*) XX, YY, PN
                IF ( ILOOP.EQ.1 .AND. IAPROC.EQ.1 ) THEN
                  BACKSPACE (NDSI)
                  READ (NDSI,'(A)') LINE
                  WRITE (NDSS,'(A)') LINE
                END IF
                IF ( INDEX(PN,"STOPSTRING").NE.0 ) EXIT
                NPTS   = NPTS + 1
                IF ( ILOOP .EQ. 1 ) CYCLE
                X(NPTS)      = XX
                Y(NPTS)      = YY
                PNAMES(NPTS) = PN
                IF ( IAPROC .EQ. NAPOUT ) THEN
                  IF ( FLAGLL ) THEN
                    IF ( NPTS .EQ. 1 ) THEN
                      WRITE (NDSO,2945)                     &
                           FACTOR*XX, FACTOR*YY, PN
                    ELSE
                      WRITE (NDSO,2946) NPTS,               &
                           FACTOR*XX, FACTOR*YY, PN
                    END IF
                  ELSE
                    IF ( NPTS .EQ. 1 ) THEN
                      WRITE (NDSO,2955)                     &
                           FACTOR*XX, FACTOR*YY, PN
                    ELSE
                      WRITE (NDSO,2956) NPTS,               &
                           FACTOR*XX, FACTOR*YY, PN
                    END IF
                  END IF
                END IF
              END DO
              !
              IF ( IAPROC.EQ.1 .AND. ILOOP.EQ.1 ) CLOSE (NDSS)
            END DO
            !
            IF ( NPTS.EQ.0 .AND. IAPROC.EQ.NAPOUT )               &
                 WRITE (NDSO,2947)
            IF ( IAPROC .EQ. 1 ) THEN
#ifdef W3_MPI
              CALL MPI_BARRIER ( MPI_COMM, IERR_MPI )
#endif
              CLOSE (NDSS,STATUS='DELETE')
            ELSE
              CLOSE (NDSS)
#ifdef W3_MPI
              CALL MPI_BARRIER ( MPI_COMM, IERR_MPI )
#endif
            END IF
            !


            ! Type 3: track output
          ELSE IF ( J .EQ. 3 ) THEN
            CALL NEXTLN ( COMSTR , NDSI , NDSEN )
            READ (NDSI,*) TFLAGI
            !
            IF ( .NOT. TFLAGI ) NDS(11) = -NDS(11)
            IF ( IAPROC .EQ. NAPOUT ) THEN
              IF ( .NOT. TFLAGI ) THEN
                WRITE (NDSO,3945) 'input', 'UNFORMATTED'
              ELSE
                WRITE (NDSO,3945) 'input', 'FORMATTED'
              END IF
            END IF


            ! Type 6: partitioning
          ELSE IF ( J .EQ. 6 ) THEN
            !             IPRT: IX0, IXN, IXS, IY0, IYN, IYS
            CALL NEXTLN ( COMSTR , NDSI , NDSEN )
            READ (NDSI,*) IPRT, PRTFRM
            !
            IF ( IAPROC .EQ. NAPOUT ) THEN
              IF ( PRTFRM ) THEN
                YESXNO = 'YES/--'
              ELSE
                YESXNO = '---/NO'
              END IF
              WRITE (NDSO,6945) IPRT, YESXNO
            END IF


#ifdef W3_COU
            ! Type 7: coupling
          ELSE IF ( J .EQ. 7 ) THEN
            CALL W3READFLGRD ( NDSI, NDSO, NDSS, NDSEN, COMSTR, FLG2,     &
                 FLGR2, IAPROC, NAPOUT, IERR )
            IF ( IERR .NE. 0 ) GOTO 2222
            CALL NEXTLN ( COMSTR , NDSI , NDSEN )
            READ (NDSI,'(A)',END=2001,ERR=2002,IOSTAT=IERR) FLDIN
#endif

          END IF ! J
        END IF ! ODAT
      END IF ! IF J=4
    END DO ! J

    ! force minimal allocation to avoid memory seg fault
    IF ( .NOT.ALLOCATED(X) .AND. NPTS.EQ.0 ) ALLOCATE ( X(1), Y(1), PNAMES(1) )

    ! 2.6 Homogeneous field data

    IF ( FLHOM ) THEN
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,951)                  &
           'Homogeneous field data (and moving grid) ...'
      NH     = 0
      !
      ! Start of loop
      DO
        CALL NEXTLN ( COMSTR , NDSI , NDSEN )
        READ (NDSI,*) IDTST


        ! Exit if illegal id
        IF ( IDTST.NE.IDSTR(-7) .AND. IDTST.NE.IDSTR(-6) .AND.   &
             IDTST.NE.IDSTR(-5) .AND. IDTST.NE.IDSTR(-4) .AND.   &
             IDTST.NE.IDSTR(-3) .AND. IDTST.NE.IDSTR(-2) .AND.   &
             IDTST.NE.IDSTR(-1) .AND. IDTST.NE.IDSTR(0)  .AND.   &
             IDTST.NE.IDSTR(1)  .AND. IDTST.NE.IDSTR(2)  .AND.   &
             IDTST.NE.IDSTR(3)  .AND. IDTST.NE.IDSTR(4)  .AND.   &
             IDTST.NE.IDSTR(5)  .AND. IDTST.NE.IDSTR(6)  .AND.   &
             IDTST.NE.IDSTR(10)  .AND. IDTST.NE.'STP' ) GOTO 2005

        ! Stop conditions
        IF ( IDTST .EQ. 'STP' ) THEN
          EXIT
        ELSE
          BACKSPACE ( NDSI )
        END IF

        ! Store data
        DO J=LBOUND(IDSTR,1), 10
          IF ( IDTST .EQ. IDSTR(J) ) THEN
            NH(J)    = NH(J) + 1
            IF ( NH(J) .GT. NHMAX ) GOTO 2006
            IF ( J .LE. 1  ) THEN ! water levels, etc. : get HA
              READ (NDSI,*) IDTST,           &
                   THO(1,J,NH(J)), THO(2,J,NH(J)),            &
                   HA(NH(J),J)
            ELSE IF ( J .EQ. 2 ) THEN ! currents: get HA and HD
              READ (NDSI,*) IDTST,           &
                   THO(1,J,NH(J)), THO(2,J,NH(J)),            &
                   HA(NH(J),J), HD(NH(J),J)
            ELSE IF ( J .EQ. 3 ) THEN ! wind: get HA HD and HS
              READ (NDSI,*) IDTST,           &
                   THO(1,J,NH(J)), THO(2,J,NH(J)),            &
                   HA(NH(J),J), HD(NH(J),J), HS(NH(J),J)
            ELSE IF ( J .EQ. 4 ) THEN ! ice
              READ (NDSI,*) IDTST,           &
                   THO(1,J,NH(J)), THO(2,J,NH(J)),            &
                   HA(NH(J),J)
            ELSE IF ( J .EQ. 5 ) THEN ! atmospheric momentum
              READ (NDSI,*) IDTST,           &
                   THO(1,J,NH(J)), THO(2,J,NH(J)),            &
                   HA(NH(J),J), HD(NH(J),j)
            ELSE IF ( J .EQ. 6 ) THEN ! air density
              READ (NDSI,*) IDTST,           &
                   THO(1,J,NH(J)), THO(2,J,NH(J)),            &
                   HA(NH(J),J)
            ELSE IF ( J .EQ. 10 ) THEN ! mov: HA and HD
              READ (NDSI,*) IDTST,           &
                   THO(1,J,NH(J)), THO(2,J,NH(J)),            &
                   HA(NH(J),J), HD(NH(J),J)
            END IF
          END IF
        END DO
      END DO
!      call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 3')
      !
#ifdef W3_O7
      DO J=JFIRST, 10
        IF ( FLH(J) .AND. IAPROC.EQ.NAPOUT ) THEN
          WRITE (NDSO,952) NH(J), IDFLDS(J)
          DO I=1, NH(J)
            IF ( ( J .LE. 1 ) .OR. ( J .EQ. 4 ) .OR.      &
                 ( J .EQ. 6 ) ) THEN
              WRITE (NDSO,953) I, THO(1,J,I), THO(2,J,I), &
                   HA(I,J)
            ELSE IF ( ( J .EQ. 2 ) .OR. ( J .EQ. 5 ) .OR. &
                 ( J .EQ. 10 ) ) THEN
              WRITE (NDSO,953) I, THO(1,J,I), THO(2,J,I), &
                   HA(I,J), HD(I,J)
            ELSE IF ( J .EQ. 3 ) THEN
              WRITE (NDSO,953) I, THO(1,J,I), THO(2,J,I), &
                   HA(I,J), HD(I,J), HS(I,J)
            END IF
          END DO
        END IF
      END DO
#endif
      !
      !
      IF ( ( FLH(-7) .AND. (NH(-7).EQ.0) ) .OR.                     &
           ( FLH(-6) .AND. (NH(-6).EQ.0) ) .OR.                     &
           ( FLH(-5) .AND. (NH(-5).EQ.0) ) .OR.                     &
           ( FLH(-4) .AND. (NH(-4).EQ.0) ) .OR.                     &
           ( FLH(-3) .AND. (NH(-3).EQ.0) ) .OR.                     &
           ( FLH(-2) .AND. (NH(-2).EQ.0) ) .OR.                     &
           ( FLH(-1) .AND. (NH(-1).EQ.0) ) .OR.                     &
           ( FLH(0)  .AND. (NH(0).EQ.0)  ) .OR.                     &
           ( FLH(1)  .AND. (NH(1).EQ.0)  ) .OR.                     &
           ( FLH(2)  .AND. (NH(2).EQ.0)  ) .OR.                     &
           ( FLH(3)  .AND. (NH(3).EQ.0)  ) .OR.                     &
           ( FLH(4)  .AND. (NH(4).EQ.0)  ) .OR.                     &
           ( FLH(5)  .AND. (NH(5).EQ.0)  ) .OR.                     &
           ( FLH(6)  .AND. (NH(6).EQ.0)  ) .OR.                     &
           ( FLH(10) .AND. (NH(10).EQ.0) ) ) GOTO 2007
      !
    END IF ! FLHOM

  END IF

  !
  ! ----------------
  !

  ! 2.1 input fields

  ! 2.1.a Opening field and data files

  IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,950)
  IF ( FLFLG ) THEN
    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,951)                  &
         'Preparing input files ...'
    !

    DO J=JFIRST, 6
      IF ( INFLAGS1(J) .AND. .NOT. FLAGSC(J)) THEN
        IF ( FLH(J) ) THEN
          IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,954) IDFLDS(J)
        ELSE
          FLAGTIDE = 0
          CALL W3FLDO ('READ', IDSTR(J), NDSF(J), NDST,     &
               NDSEN, NX, NY, GTYPE,               &
               IERR, FPRE=TRIM(FNMPRE), TIDEFLAGIN=FLAGTIDE )
          IF ( IERR .NE. 0 ) GOTO 2222
#ifdef W3_TIDE
          IF (FLAGTIDE.GT.0.AND.J.EQ.1) FLAGSTIDE(1)=.TRUE.
          IF (FLAGTIDE.GT.0.AND.J.EQ.2) FLAGSTIDE(2)=.TRUE.
#endif
          IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,955) IDFLDS(J)
        END IF
      ELSE
        IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,954) IDFLDS(J)
      END IF
    END DO
    !
    DO J=7, 9
      IF ( INFLAGS1(J) .AND. .NOT. FLAGSC(J)) THEN
        CALL W3FLDO ('READ', IDSTR(J), NDSF(J), NDST, NDSEN, &
             RCLD(J), NY, NODATA(J),                 &
             IERR, FPRE=TRIM(FNMPRE) )
        IF ( IERR .NE. 0 ) GOTO 2222
        IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,956) IDFLDS(J),&
             RCLD(J), NODATA(J)
      ELSE
        IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,954) IDFLDS(J)
      END IF
    END DO
    !
  END IF ! FLFLG

!  call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 4')

  ! 2.2 Time setup

  IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,930)
  CALL STME21 ( TIME0 , DTME21 )
  IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,931) DTME21
  TIME = TIME0
  CALL STME21 ( TIMEN , DTME21 )
  IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,932) DTME21
#ifdef W3_OASIS
  TIME00 = TIME0
  TIMEEND = TIMEN
#endif
#ifdef W3_NL5
  QI5TBEG = TIME0
#endif
  !
  DTTST  = DSEC21 ( TIME0 , TIMEN )
  IF ( DTTST .LE. 0. ) GOTO 2003

  ! 2.3 Domain setup

  IOSTYP = MAX ( 0 , MIN ( 3 , IOSTYP ) )
#ifdef W3_PDLIB
  IF (IOSTYP .gt. 1) THEN
    WRITE(*,*) 'IOSTYP not supported in domain decomposition mode'
    CALL EXTCDE ( 6666 )
  ENDIF
#endif

  IF ( IAPROC .EQ. NAPOUT ) THEN
    IF ( IOSTYP .EQ. 0 ) THEN
      WRITE (NDSO,940) 'No dedicated output process, ' //   &
           'parallel file system required.'
    ELSE IF ( IOSTYP .EQ. 1 ) THEN
      WRITE (NDSO,940) 'No dedicated output process, ' //   &
           'any file system.'
    ELSE IF ( IOSTYP .EQ. 2 ) THEN
      WRITE (NDSO,940) 'Single dedicated output process.'
    ELSE IF ( IOSTYP .EQ. 3 ) THEN
      WRITE (NDSO,940) 'Multiple dedicated output processes.'
    ELSE
      WRITE (NDSO,940) 'IOSTYP NOT RECOGNIZED'
    END IF
  END IF


  ! 2.4 Output dates

  DO J = 1, NOTYPE
    !
    IF ( ODAT(5*(J-1)+3) .NE. 0 ) THEN
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,941) J, IDOTYP(J)
      TTIME(1) = ODAT(5*(J-1)+1)
      TTIME(2) = ODAT(5*(J-1)+2)
      CALL STME21 ( TTIME , DTME21 )
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,942) DTME21
      TTIME(1) = ODAT(5*(J-1)+4)
      TTIME(2) = ODAT(5*(J-1)+5)
      CALL STME21 ( TTIME , DTME21 )
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,943) DTME21
      TTIME(1) = 0
      TTIME(2) = 0
      DTTST    = REAL ( ODAT(5*(J-1)+3) )
      CALL TICK21 ( TTIME , DTTST  )
      CALL STME21 ( TTIME , DTME21 )
      IF ( ( ODAT(5*(J-1)+1) .NE. ODAT(5*(J-1)+4) .OR.          &
           ODAT(5*(J-1)+2) .NE. ODAT(5*(J-1)+5) ) .AND.       &
           IAPROC .EQ. NAPOUT ) THEN
        IF ( DTME21(9:9) .NE. '0' ) THEN
          WRITE (NDSO,1944) DTME21( 9:19)
        ELSE IF ( DTME21(10:10) .NE. '0' ) THEN
          WRITE (NDSO,2944) DTME21(10:19)
        ELSE
          WRITE (NDSO,3944) DTME21(12:19)
        END IF
      END IF
    END IF
  END DO
  !
  ! CHECKPOINT
  J=8
  IF (ODAT(38) .NE. 0) THEN
    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,941) J, IDOTYP(J)
    TTIME(1) = ODAT(5*(J-1)+1)
    TTIME(2) = ODAT(5*(J-1)+2)
    CALL STME21 ( TTIME , DTME21 )
    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,942) DTME21
    TTIME(1) = ODAT(5*(J-1)+4)
    TTIME(2) = ODAT(5*(J-1)+5)
    CALL STME21 ( TTIME , DTME21 )
    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,943) DTME21
    TTIME(1) = 0
    TTIME(2) = 0
    DTTST    = REAL ( ODAT(5*(J-1)+3) )
    CALL TICK21 ( TTIME , DTTST  )
    CALL STME21 ( TTIME , DTME21 )
    IF ( ( ODAT(5*(J-1)+1) .NE. ODAT(5*(J-1)+4) .OR.          &
         ODAT(5*(J-1)+2) .NE. ODAT(5*(J-1)+5) ) .AND.       &
         IAPROC .EQ. NAPOUT ) THEN
      IF ( DTME21(9:9) .NE. '0' ) THEN
        WRITE (NDSO,1944) DTME21( 9:19)
      ELSE IF ( DTME21(10:10) .NE. '0' ) THEN
        WRITE (NDSO,2944) DTME21(10:19)
      ELSE
        WRITE (NDSO,3944) DTME21(12:19)
      END IF
    END IF
  END IF

  !
  ! 2.5 Output types

#ifdef W3_T
  WRITE (NDST,9040) ODAT
  WRITE (NDST,9041) FLGRD
  WRITE (NDST,9042) IPRT, PRTFRM
#endif

  !
  ! For outputs with non-zero time step, check dates :
  ! If output ends before run start OR output starts after run end,
  ! deactivate output cleanly with output time step = 0
  ! This is usefull for IOSTYP=3 (Multiple dedicated output processes)
  ! to avoid the definition of dedicated proc. for unused output.
  !
  DO J = 1, NOTYPE
    DTTST  = DSEC21 ( TIME0 , ODAT(5*(J-1)+4:5*(J-1)+5) )
    IF ( DTTST .LT. 0 ) THEN
      ODAT(5*(J-1)+3) = 0
      IF ( IAPROC .EQ. NAPOUT )  WRITE (NDSO,8945) TRIM(IDOTYP(J))
      CONTINUE
    END IF
    DTTST  = DSEC21 ( ODAT(5*(J-1)+1:5*(J-1)+2), TIMEN )
    IF ( DTTST .LT. 0 ) THEN
      ODAT(5*(J-1)+3) = 0
      IF ( IAPROC .EQ. NAPOUT )  WRITE (NDSO,8945) TRIM(IDOTYP(J))
      CONTINUE
    END IF
  END DO
  !
  ! CHECKPOINT
  J = 8
  DTTST  = DSEC21 ( TIME0 , ODAT(5*(J-1)+4:5*(J-1)+5) )
  IF ( DTTST .LT. 0 ) THEN
    ODAT(5*(J-1)+3) = 0
    IF ( IAPROC .EQ. NAPOUT )  WRITE (NDSO,8945) TRIM(IDOTYP(J))
    CONTINUE
  END IF
  DTTST  = DSEC21 ( ODAT(5*(J-1)+1:5*(J-1)+2), TIMEN )
  IF ( DTTST .LT. 0 ) THEN
    ODAT(5*(J-1)+3) = 0
    IF ( IAPROC .EQ. NAPOUT )  WRITE (NDSO,8945) TRIM(IDOTYP(J))
    CONTINUE
  END IF
  !
!  call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 5')
  !
  !--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! 5.  Initializations
  !

  IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,951) 'Wave model ...'
  !
#ifdef W3_TIDE
  IF (FLAGSTIDE(1).OR.FLAGSTIDE(2)) THEN
    CALL VUF_SET_PARAMETERS
    IF (FLAGSTIDE(1)) CALL W3FLDTIDE1 ( 'READ',  NDSF(1), NDST, NDSEN, NX, NY, IDSTR(1), IERR )
    IF (FLAGSTIDE(2)) CALL W3FLDTIDE1 ( 'READ',  NDSF(2), NDST, NDSEN, NX, NY, IDSTR(2), IERR )
  END IF
#endif
  !
#ifdef W3_COU
  ! Sent coupled fields must be written in the restart when coupling at T+0
  IF (CPLT0) THEN
    DO J=1, NOGRP
      FLOGR(J)  = FLOGR(J)  .OR. FLG2(J)
      DO I=1, NGRPP
        FLOGRR(J,I) = FLOGRR(J,I) .OR. FLGR2(J,I)
      END DO
    END DO
  ENDIF
#endif
  !
  OARST = ANY(FLOGR)
  !
  CALL W3INIT ( 1, .FALSE., 'ww3', NDS, NTRACE, ODAT, FLGRD, FLGR2, FLGD,    &
       FLG2, NPTS, X, Y, PNAMES, IPRT, PRTFRM, MPI_COMM,   &
       FLAGSTIDEIN=FLAGSTIDE )
  !
  !      IF (MINVAL(VA) .LT. 0.) THEN
  !        WRITE(740+IAPROC,*) 'NEGATIVE ACTION SHELL 5', MINVAL(VA)
  !        CALL FLUSH(740+IAPROC)
  !        CALL EXTCDE(665)
  !      ENDIF
  !      IF (SUM(VA) .NE. SUM(VA)) THEN
  !        WRITE(740+IAPROC,*) 'NAN in ACTION SHEL1', SUM(VA)
  !        CALL FLUSH(740+IAPROC)
  !        CALL EXTCDE(666)
  !      ENDIF

!  call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 5')
  !
#ifdef W3_TIDE
  IF (FLAGSTIDE(1)) CALL W3FLDTIDE2 ( 'READ',  NDSF(1), NDST, NDSEN, NX, NY, IDSTR(1), 1, IERR )
  IF (FLAGSTIDE(2)) CALL W3FLDTIDE2 ( 'READ',  NDSF(2), NDST, NDSEN, NX, NY, IDSTR(2), 1, IERR )
  ALLOCATE(V_ARG(170,1),F_ARG(170,1),U_ARG(170,1))  ! to be removed later ...
#endif
  !
  ALLOCATE ( XXX(NX,NY) )
  !

  !
#ifdef W3_MPI
  CALL MPI_BARRIER ( MPI_COMM, IERR_MPI )
#endif
  !
  IF ( IAPROC .EQ. NAPOUT ) THEN
    CALL DATE_AND_TIME ( VALUES=CLKDT2 )
  END IF
  !--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
#ifdef W3_OASIS
  ! Initialize L_MASTER, COUPL_COMM
  IF ( IAPROC .EQ. 1) THEN
    L_MASTER = .TRUE.
  ELSE
    L_MASTER = .FALSE.
  ENDIF
  ! Estimate the weights for the spatial interpolation
  IF (DTOUT(7).NE.0) THEN
    CALL CPL_OASIS_GRID(L_MASTER,MPI_COMM)
    CALL CPL_OASIS_DEFINE(NDSO, FLDIN, FLDOUT)
  END IF
#endif

!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     End of shel
!
      GOTO 2222

  !
  ! Error escape locations
  !
2000 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000) IERR
  CALL EXTCDE ( 1000 )
  !
2001 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1001)
  CALL EXTCDE ( 1001 )
  !
2002 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1002) IERR
  CALL EXTCDE ( 1002 )
  !
2102 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1102)
  CALL EXTCDE ( 1102 )
  !
2003 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1003)
  CALL EXTCDE ( 1003 )
  !
2104 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1104) IERR
  CALL EXTCDE ( 1104 )
  !
2004 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1004) IERR
  CALL EXTCDE ( 1004 )
  !
2005 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1005) IDTST
  CALL EXTCDE ( 1005 )
  !
2006 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1006) IDTST, NH(J)
  CALL EXTCDE ( 1006 )
  !
2062 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1062) IDTST
  CALL EXTCDE ( 1062 )
  !
2007 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1007)
  CALL EXTCDE ( 1007 )
  !
2008 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1008) IERR
  CALL EXTCDE ( 1008 )
  !
#ifdef W3_COU
2009 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1009) ODAT(33), NINT(DTMAX)
  CALL EXTCDE ( 1009 )
#endif
  !
2054 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1054)
  CALL EXTCDE ( 1054 )
2222 CONTINUE
  !
#ifdef W3_MPI
  CALL MPI_BARRIER ( MPI_COMM, IERR_MPI )
#endif
  !
  IF ( IAPROC .EQ. NAPOUT ) THEN
    CALL DATE_AND_TIME ( VALUES=CLKDT3 )
    CLKFIN = MAX(TDIFF ( CLKDT1,CLKDT2 ), 0.)
    CLKFEL = MAX(TDIFF ( CLKDT1,CLKDT3 ), 0.)
    WRITE (NDSO,997) CLKFIN
    WRITE (NDSO,998) CLKFEL
    IF ( NDSO .NE. NDS(1) ) THEN
      WRITE (NDS(1),997) CLKFIN
      WRITE (NDS(1),998) CLKFEL
    END IF
    WRITE (NDSO,999)
  END IF
  !
#ifdef W3_NCO
  !     IF ( IAPROC .EQ. 1 ) CALL W3TAGE('WAVEFCST')
#endif
#ifdef W3_OASIS
!  IF (OASISED.EQ.1) THEN
!    CALL CPL_OASIS_FINALIZE
!  ELSE
#endif
#ifdef W3_MPI
!    CALL MPI_FINALIZE  ( IERR_MPI )
#endif
#ifdef W3_OASIS
!  END IF
#endif
  !
  !
  ! Formats
  !
900 FORMAT (/15X,'      *** WAVEWATCH III Program shell ***      '/ &
       15X,'==============================================='/)
901 FORMAT ( '  Comment character is ''',A,''''/)
  !
#ifdef W3_OMPH
905 FORMAT ( '  Hybrid MPI/OMP thread support level:'/        &
       '     Requested: ', I2/                          &
       '      Provided: ', I2/ )
#endif
  !
#ifdef W3_OMPG
906 FORMAT ( '  OMP threading enabled. Number of threads: ', I3 / )
#endif
920 FORMAT (/'  Input fields : '/                                   &
       ' --------------------------------------------------')
921 FORMAT ( '       ',A,2X,A,2X,A)
922 FORMAT ( ' ' )
  !
930 FORMAT (/'  Time interval : '/                                  &
       ' --------------------------------------------------')
931 FORMAT ( '       Starting time : ',A)
932 FORMAT ( '       Ending time   : ',A/)
  !
940 FORMAT (/'  Output requests : '/                                &
       ' --------------------------------------------------'/ &
       '       ',A)
941 FORMAT (/'       Type',I2,' : ',A/                              &
       '      -----------------------------------------')
942 FORMAT ( '            From     : ',A)
943 FORMAT ( '            To       : ',A)
1944 FORMAT ( '            Interval : ', 8X,A11/)
2944 FORMAT ( '            Interval : ', 9X,A10/)
3944 FORMAT ( '            Interval : ',11X,A8/)
2945 FORMAT ( '            Point  1 : ',2F8.2,2X,A)
2955 FORMAT ( '            Point  1 : ',2(F8.1,'E3'),2X,A)
2946 FORMAT ( '              ',I6,' : ',2F8.2,2X,A)
2956 FORMAT ( '              ',I6,' : ',2(F8.1,'E3'),2X,A)
2947 FORMAT ( '            No points defined')
3945 FORMAT ( '            The file with ',A,' data is ',A,'.')
6945 FORMAT ( '            IX first,last,inc :',3I5/                 &
       '            IY first,last,inc :',3I5/                 &
       '            Formatted file    :    ',A)
8945 FORMAT ( '            output dates out of run dates : ', A,     &
       ' deactivated')
  !
950 FORMAT (/'  Initializations :'/                                 &
       ' --------------------------------------------------')
951 FORMAT ( '       ',A)
#ifdef W3_O7
952 FORMAT ( '       ',I6,2X,A)
953 FORMAT ( '          ',I6,I11.8,I7.6,3E12.4)
#endif
954 FORMAT ( '            ',A,': file not needed')
955 FORMAT ( '            ',A,': file OK')
956 FORMAT ( '            ',A,': file OK, recl =',I3,               &
       '  undef = ',E10.3)
  !
960 FORMAT (/'  Running model without input fields'/                &
       ' --------------------------------------------------'/)
  !
970 FORMAT (/'  Running model with input fields'/                   &
       ' --------------------------------------------------')
971 FORMAT (/'  Updating input at ',A)
972 FORMAT ( '     Updating ',A)
973 FORMAT ( '        Past last ',A)
#ifdef W3_TIDE
974 FORMAT ( '     Updating ',A,'using tidal constituents')
#endif
975 FORMAT (/'  Data assimmilation at ',A)
  !
997 FORMAT (/'  Initialization time :',F10.2,' s')
998 FORMAT ( '  Elapsed time        :',F10.2,' s')
  !
999 FORMAT(/'  End of program '/                                    &
       ' ===================================='/               &
       '         WAVEWATCH III Program shell '/)
  !
1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ERROR IN OPENING INPUT FILE'/                    &
       '     IOSTAT =',I5/)
  !
1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     PREMATURE END OF INPUT FILE'/)
  !
1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ERROR IN READING FROM INPUT FILE'/               &
       '     IOSTAT =',I5/)
  !
1102 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     LEVEL AND CURRENT ARE MIXING COUPLED AND FORCED'/&
       '     IT MUST BE FULLY COUPLED OR DISABLED '/)
  !
1003 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ILLEGAL TIME INTERVAL'/)
  !
1104 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ERROR IN OPENING POINT FILE'/                    &
       '     IOSTAT =',I5/)
  !
1004 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ERROR IN READING FROM POINT FILE'/               &
       '     IOSTAT =',I5/)
  !
1005 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ILLEGAL ID STRING HOMOGENEOUS FIELD : ',A/)
  !
1006 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     TOO MANY HOMOGENEOUS FIELDS : ',A,1X,I4/)
  !
1062 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : ***'/             &
       '     HOMOGENEOUS NAME NOT RECOGNIZED : ', A/)
  !
1007 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     INSUFFICIENT DATA FOR HOMOGENEOUS FIELDS'/)
  !
1008 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ERROR IN OPENING OUTPUT FILE'/                   &
       '     IOSTAT =',I5/)
  !
#ifdef W3_COU
1009 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     COUPLING TIME STEP NOT MULTIPLE OF'/             &
       '     MODEL TIME STEP: ',I6, I6/)
#endif
  !
#ifdef W3_COU
1010 FORMAT (/' *** WAVEWATCH III WARNING IN W3SHEL : *** '/         &
       '     COUPLING TIME STEP NOT DEFINED, '/               &
       '     IT WILL BE OVERRIDEN TO DEFAULT VALUE'/          &
       '     FROM ',I6, ' TO ',I6/)
#endif
  !
1054 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     POINT OUTPUT ACTIVATED BUT NO POINTS DEFINED'/)
  !
  !
#ifdef W3_T
9000 FORMAT ( ' TEST W3SHEL : UNIT NUMBERS  :',12I4)
9001 FORMAT ( ' TEST W3SHEL : SUBR. TRACING :',2I4)
#endif
  !
#ifdef W3_T
9020 FORMAT ( ' TEST W3SHEL : FLAGS DEF / HOM  : ',9L2,2X,9L2)
#endif
  !
#ifdef W3_T
9040 FORMAT ( ' TEST W3SHEL : ODAT   : ',I9.8,I7.6,I7,I9.8,I7.6,  &
       4(/24X,I9.8,I7.6,I7,I9.8,I7.6) )
9041 FORMAT ( ' TEST W3SHEL : FLGRD  : ',20L2)
9042 FORMAT ( ' TEST W3SHEL : IPR, PRFRM : ',6I6,1X,L1)
#endif
  !
#ifdef W3_T
9070 FORMAT ( ' TEST W3SHEL : ',A,3X,2(I10.8,I7.6)/                &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,2(I10.8,I7.6)/                &
       '               ',A,L3,2(I10.8,I7.6)/                &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,2(I10.8,I7.6)/                &
       '               ',A,L3,2(I10.8,I7.6)/                &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,2(I10.8,I7.6))
9071 FORMAT ( ' TEST W3SHEL : ',A,', DTTST = ',E10.3)
9072 FORMAT ( ' TEST W3SHEL : ',A,3X,2(I10.8,I7.6)/               &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,2(I10.8,I7.6)/               &
       '               ',A,L3,2(I10.8,I7.6)/               &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,2(I10.8,I7.6)/               &
       '               ',A,L3,2(I10.8,I7.6)/               &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,2(I10.8,I7.6)/           &
       '               ',A,L3,2(I10.8,I7.6))
#endif
  !/
  !/ End of W3SHEL ----------------------------------------------------- /
  !/
END SUBROUTINE FISOC_WW3_INIT

SUBROUTINE FISOC_WW3_RUN(START_DATE_opt,END_DATE_opt)

    
    IMPLICIT NONE

    INTEGER STATUS
    character(len=*),INTENT(IN),OPTIONAL :: START_DATE_opt, END_DATE_opt
    
    character(len=8) DS,DE,TS,TE
    !TYPE(TIME) :: StartTime_sub
    !TYPE(TIME) :: EndTime_sub 
  !--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! 6.  Model without input
  !
  !      IF (MINVAL(VA) .LT. 0.) THEN
  !        WRITE(740+IAPROC,*) 'NEGATIVE ACTION SHELL 6', MINVAL(VA)
  !        CALL FLUSH(740+IAPROC)
  !        CALL EXTCDE(665)
  !      ENDIF
  !      IF (SUM(VA) .NE. SUM(VA)) THEN
  !        WRITE(740+IAPROC,*) 'NAN in ACTION SHEL2', SUM(VA)
  !        CALL FLUSH(740+IAPROC)
  !        CALL EXTCDE(666)
  !      ENDIF
!  call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 6')

   IF(PRESENT(START_DATE_opt).AND.PRESENT(END_DATE_opt)) THEN
    DS(1:4) = START_DATE_opt(1:4)
    DS(5:6) = START_DATE_opt(6:7)
    DS(7:8) = START_DATE_opt(9:10)
    
    TS(1:2) = START_DATE_opt(12:13)
    TS(3:4) = START_DATE_opt(15:16)
    TS(5:6) = START_DATE_opt(18:19)
   
    DE(1:4) = END_DATE_opt(1:4)
    DE(5:6) = END_DATE_opt(6:7)
    DE(7:8) = END_DATE_opt(9:10)
    
    TE(1:2) = END_DATE_opt(12:13)
    TE(3:4) = END_DATE_opt(15:16)
    TE(5:6) = END_DATE_opt(18:19)
    
    READ(DS, '(i8)') TIME0(1) 
    READ(DE, '(i8)') TIMEN(1) 
    READ(TS, '(i6)') TIME0(2) 
    READ(TE, '(i6)') TIMEN(2) 
   END IF 
   
  IF ( .NOT. FLFLG ) THEN
    !
    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,960)
    CALL W3WAVE ( 1, ODAT, TIMEN                      &
        )
    !
    GOTO 2222
    !
  END IF
  !--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! 7.  Model with input
  !
  IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,970)
  !


700 CONTINUE
  !
  !
  ! 7.a Determine next time interval and input fields
  ! 7.a.1 Preparation
  !
  TTIME  = TIMEN
  !
  CALL STME21 ( TIME0 , DTME21 )

  IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,971) DTME21
  !
  !
!  call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 7')

!JQI<
  INFLAGS1(1) = .FALSE.
  INFLAGS1(2) = .FALSE.
!JQI>
  DO J=JFIRST,10
    !
    write(jchar, '(i0)') j
!    call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL UPDATE '//trim(jchar))

    IF ( INFLAGS1(J) ) THEN
      !
      ! 7.a.2 Check if update is needed
      !
      IF (.NOT.FLAGSC(J)) THEN
        TTT(1) = TFN(1,J)
        TTT(2) = TFN(2,J)
        IF ( TTT(1) .EQ. -1 ) THEN
          DTTST  = 0.
        ELSE
          DTTST  = DSEC21 ( TIME0 , TTT )
        END IF
      END IF
      !
      !
      ! 7.a.3 Update time and fields / data
      !
      IF ( DTTST .LE. 0. ) THEN

          IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,972) IDFLDS(J)

! LEV : water levels
            IF ( J .EQ. 1 ) THEN
              IF ( FLH(J) ) THEN
                CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                             TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                             TTT, XXX, XXX, XXX, TLN, XXX, XXX, WLEV, IERR)
              ELSE
                CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                             NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                             TTT, XXX, XXX, XXX, TLN, XXX, XXX, WLEV,   &
                             IERR, FLAGSC(J)                            &
                             )
              END IF
              IF ( IERR .LT. 0 ) FLLSTL = .TRUE.
!could be:    IF ( IERR .LT. 0 ) FLLST_ALL(J) = .TRUE.

! CUR : currents
            ELSE IF ( J .EQ. 2 ) THEN
              IF ( FLH(J) ) THEN
                CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                             TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                             TC0, CX0, CY0, XXX, TCN, CXN, CYN, XXX, IERR)
!
              ELSE
                CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                             NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                             TC0, CX0, CY0, XXX, TCN, CXN, CYN, XXX,    &
                             IERR, FLAGSC(J)                            &
                             )
              END IF

          ! WND : winds
        ELSE IF ( J .EQ. 3 ) THEN
          IF ( FLH(J) ) THEN
            CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                 TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                 TW0, WX0, WY0, DT0, TWN, WXN, WYN, DTN, IERR)
            !
          ELSE
            CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                 NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                 TW0, WX0, WY0, DT0, TWN, WXN, WYN, DTN,    &
                 IERR, FLAGSC(J)                            &
                 )
          END IF

        END IF
        !
        IF ( IERR.GT.0 ) GOTO 2222
        IF ( IERR.LT.0 .AND. IAPROC.EQ.NAPOUT ) WRITE (NDSO,973) IDFLDS(J)


      END IF ! DTTST .LE. 0.
      !
      ! 7.a.4 Update next ending time
      !
      IF ( INFLAGS1(J) ) THEN
        TTT    = TFN(:,J)
        DTTST  = DSEC21 ( TTT , TTIME )
        IF ( DTTST.GT.0. .AND. .NOT.                          &
             ( (FLLSTL .AND. J.EQ.1) .OR.                   &
             (FLLST_ALL(J) .AND. J.EQ.-7) .OR.            &
             (FLLST_ALL(J) .AND. J.EQ.-6) .OR.            &
             (FLLST_ALL(J) .AND. J.EQ.-5) .OR.            &
             (FLLST_ALL(J) .AND. J.EQ.-4) .OR.            &
             (FLLST_ALL(J) .AND. J.EQ.-3) .OR.            &
             (FLLST_ALL(J) .AND. J.EQ.-2) .OR.            &
             (FLLST_ALL(J) .AND. J.EQ.-1) .OR.            &
             (FLLST_ALL(J) .AND. J.EQ.0 ) .OR.            &
             (FLLSTI .AND. J.EQ.4) .OR.                   &
             (FLLSTR .AND. J.EQ.6) ) ) THEN
          TTIME  = TTT
          ! notes: if model has run out beyond field input, then this line should not
          !    be reached.
        END IF
      END IF
      !
    END IF ! INFLAGSC1(J)
    !
  END DO ! J=JFIRST,10

!  call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 8')

  TDN = TTIME
  CALL TICK21 ( TDN, 1. )
  DO J=7, 9
    IF ( INFLAGS1(J) ) THEN
      TTT    = TFN(:,J)
      DTTST  = DSEC21 ( TTT , TDN )
      IF ( DTTST.GT.0. ) TDN = TTT
    END IF
  END DO
  !
  !
  IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,*) ' '
  !
  ! 7.b Run the wave model for the given interval
  !
  TIME0  = TTIME
  !
  CALL W3WAVE ( 1, ODAT, TIME0)
!  call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 9')
  !
  ! The following lines prevents us from trying to read past the end
  ! of the files. This feature existed in v3.14.
  ! "1" is for water levels
  ! "4" is for ice concentration:
  ! "6" is for air density:
  IF ( FLLSTL ) INFLAGS1(1) = .FALSE.
  IF ( FLLSTI ) INFLAGS1(4) = .FALSE.
  IF ( FLLSTR ) INFLAGS1(6) = .FALSE.

  ! We include something like this for mud and ice parameters also:
!  DO J=-7,0
!    IF (FLLST_ALL(J))THEN
!      INFLAGS1(J)=.FALSE.
!    END IF
!  END DO

  !
  ! 7.c Run data assimilation at ending time
  !
!  DTTST  = DSEC21 ( TIME , TDN )
!  IF ( DTTST .EQ. 0 ) THEN
!    CALL STME21 ( TIME0 , DTME21 )
!    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,975) DTME21
    !
!    FLGDAS(1) = DSEC21(TIME,T0N) .EQ. 0.
!    FLGDAS(2) = DSEC21(TIME,T1N) .EQ. 0.
!    FLGDAS(3) = DSEC21(TIME,T2N) .EQ. 0.
    !
!    CALL W3WDAS ( FLGDAS, RCLD, NDT, DATA0, DATA1, DATA2 )
    !
    ! 7.d Call wave model again after data assimilation for output only
    !
!    DTTST  = DSEC21 ( TIME , TIMEN )

!    IF ( DTTST .EQ. 0. ) THEN
!      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,*) ' '
!      CALL W3WAVE ( 1, ODAT, TIME0                                 &
!           )
!    END IF
!  END IF
  !
  ! 7.e Check times
  !
!  call print_memcheck(memunit, 'memcheck_____:'//' WW3_SHEL SECTION 10')

  DTTST  = DSEC21 ( TIME0 , TIMEN )
  
  IF ( DTTST .GT. 0. ) GOTO 700
  !
  !--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !     End of shel
  !
  GOTO 2222
  !
  ! Error escape locations
  !
2000 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000) IERR
  CALL EXTCDE ( 1000 )
  !
2001 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1001)
  CALL EXTCDE ( 1001 )
  !
2002 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1002) IERR
  CALL EXTCDE ( 1002 )
  !
2102 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1102)
  CALL EXTCDE ( 1102 )
  !
2003 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1003)
  CALL EXTCDE ( 1003 )
  !
2104 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1104) IERR
  CALL EXTCDE ( 1104 )
  !
2004 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1004) IERR
  CALL EXTCDE ( 1004 )
  !
2005 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1005) IDTST
  CALL EXTCDE ( 1005 )
  !
2006 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1006) IDTST, NH(J)
  CALL EXTCDE ( 1006 )
  !
2062 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1062) IDTST
  CALL EXTCDE ( 1062 )
  !
2007 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1007)
  CALL EXTCDE ( 1007 )
  !
2008 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1008) IERR
  CALL EXTCDE ( 1008 )
  !
#ifdef W3_COU
2009 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1009) ODAT(33), NINT(DTMAX)
  CALL EXTCDE ( 1009 )
#endif
  !
2054 CONTINUE
  IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1054)
  CALL EXTCDE ( 1054 )
2222 CONTINUE
  !
#ifdef W3_MPI
  CALL MPI_BARRIER ( MPI_COMM, IERR_MPI )
#endif
  !
  IF ( IAPROC .EQ. NAPOUT ) THEN
    CALL DATE_AND_TIME ( VALUES=CLKDT3 )
    CLKFIN = MAX(TDIFF ( CLKDT1,CLKDT2 ), 0.)
    CLKFEL = MAX(TDIFF ( CLKDT1,CLKDT3 ), 0.)
    WRITE (NDSO,997) CLKFIN
    WRITE (NDSO,998) CLKFEL
    IF ( NDSO .NE. NDS(1) ) THEN
      WRITE (NDS(1),997) CLKFIN
      WRITE (NDS(1),998) CLKFEL
    END IF
    WRITE (NDSO,999)
  END IF
  !
#ifdef W3_NCO
  !     IF ( IAPROC .EQ. 1 ) CALL W3TAGE('WAVEFCST')
#endif
#ifdef W3_OASIS
  IF (OASISED.EQ.1) THEN
    CALL CPL_OASIS_FINALIZE
  ELSE
#endif
#ifdef W3_MPI
!JQI    CALL MPI_FINALIZE  ( IERR_MPI )
#endif
#ifdef W3_OASIS
  END IF
#endif
  !
  !
  ! Formats
  !
900 FORMAT (/15X,'      *** WAVEWATCH III Program shell ***      '/ &
       15X,'==============================================='/)
901 FORMAT ( '  Comment character is ''',A,''''/)
  !
#ifdef W3_OMPH
905 FORMAT ( '  Hybrid MPI/OMP thread support level:'/        &
       '     Requested: ', I2/                          &
       '      Provided: ', I2/ )
#endif
  !
#ifdef W3_OMPG
906 FORMAT ( '  OMP threading enabled. Number of threads: ', I3 / )
#endif
920 FORMAT (/'  Input fields : '/                                   &
       ' --------------------------------------------------')
921 FORMAT ( '       ',A,2X,A,2X,A)
922 FORMAT ( ' ' )
  !
930 FORMAT (/'  Time interval : '/                                  &
       ' --------------------------------------------------')
931 FORMAT ( '       Starting time : ',A)
932 FORMAT ( '       Ending time   : ',A/)
  !
940 FORMAT (/'  Output requests : '/                                &
       ' --------------------------------------------------'/ &
       '       ',A)
941 FORMAT (/'       Type',I2,' : ',A/                              &
       '      -----------------------------------------')
942 FORMAT ( '            From     : ',A)
943 FORMAT ( '            To       : ',A)
1944 FORMAT ( '            Interval : ', 8X,A11/)
2944 FORMAT ( '            Interval : ', 9X,A10/)
3944 FORMAT ( '            Interval : ',11X,A8/)
2945 FORMAT ( '            Point  1 : ',2F8.2,2X,A)
2955 FORMAT ( '            Point  1 : ',2(F8.1,'E3'),2X,A)
2946 FORMAT ( '              ',I6,' : ',2F8.2,2X,A)
2956 FORMAT ( '              ',I6,' : ',2(F8.1,'E3'),2X,A)
2947 FORMAT ( '            No points defined')
3945 FORMAT ( '            The file with ',A,' data is ',A,'.')
6945 FORMAT ( '            IX first,last,inc :',3I5/                 &
       '            IY first,last,inc :',3I5/                 &
       '            Formatted file    :    ',A)
8945 FORMAT ( '            output dates out of run dates : ', A,     &
       ' deactivated')
  !
950 FORMAT (/'  Initializations :'/                                 &
       ' --------------------------------------------------')
951 FORMAT ( '       ',A)
#ifdef W3_O7
952 FORMAT ( '       ',I6,2X,A)
953 FORMAT ( '          ',I6,I11.8,I7.6,3E12.4)
#endif
954 FORMAT ( '            ',A,': file not needed')
955 FORMAT ( '            ',A,': file OK')
956 FORMAT ( '            ',A,': file OK, recl =',I3,               &
       '  undef = ',E10.3)
  !
960 FORMAT (/'  Running model without input fields'/                &
       ' --------------------------------------------------'/)
  !
970 FORMAT (/'  Running model with input fields'/                   &
       ' --------------------------------------------------')
971 FORMAT (/'  Updating input at ',A)
972 FORMAT ( '     Updating ',A)
973 FORMAT ( '        Past last ',A)
#ifdef W3_TIDE
974 FORMAT ( '     Updating ',A,'using tidal constituents')
#endif
975 FORMAT (/'  Data assimmilation at ',A)
  !
997 FORMAT (/'  Initialization time :',F10.2,' s')
998 FORMAT ( '  Elapsed time        :',F10.2,' s')
  !
999 FORMAT(/'  End of program '/                                    &
       ' ===================================='/               &
       '         WAVEWATCH III Program shell '/)
  !
1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ERROR IN OPENING INPUT FILE'/                    &
       '     IOSTAT =',I5/)
  !
1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     PREMATURE END OF INPUT FILE'/)
  !
1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ERROR IN READING FROM INPUT FILE'/               &
       '     IOSTAT =',I5/)
  !
1102 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     LEVEL AND CURRENT ARE MIXING COUPLED AND FORCED'/&
       '     IT MUST BE FULLY COUPLED OR DISABLED '/)
  !
1003 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ILLEGAL TIME INTERVAL'/)
  !
1104 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ERROR IN OPENING POINT FILE'/                    &
       '     IOSTAT =',I5/)
  !
1004 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ERROR IN READING FROM POINT FILE'/               &
       '     IOSTAT =',I5/)
  !
1005 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ILLEGAL ID STRING HOMOGENEOUS FIELD : ',A/)
  !
1006 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     TOO MANY HOMOGENEOUS FIELDS : ',A,1X,I4/)
  !
1062 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : ***'/             &
       '     HOMOGENEOUS NAME NOT RECOGNIZED : ', A/)
  !
1007 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     INSUFFICIENT DATA FOR HOMOGENEOUS FIELDS'/)
  !
1008 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     ERROR IN OPENING OUTPUT FILE'/                   &
       '     IOSTAT =',I5/)
  !
#ifdef W3_COU
1009 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     COUPLING TIME STEP NOT MULTIPLE OF'/             &
       '     MODEL TIME STEP: ',I6, I6/)
#endif
  !
#ifdef W3_COU
1010 FORMAT (/' *** WAVEWATCH III WARNING IN W3SHEL : *** '/         &
       '     COUPLING TIME STEP NOT DEFINED, '/               &
       '     IT WILL BE OVERRIDEN TO DEFAULT VALUE'/          &
       '     FROM ',I6, ' TO ',I6/)
#endif
  !
1054 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
       '     POINT OUTPUT ACTIVATED BUT NO POINTS DEFINED'/)
  !
  !
#ifdef W3_T
9000 FORMAT ( ' TEST W3SHEL : UNIT NUMBERS  :',12I4)
9001 FORMAT ( ' TEST W3SHEL : SUBR. TRACING :',2I4)
#endif
  !
#ifdef W3_T
9020 FORMAT ( ' TEST W3SHEL : FLAGS DEF / HOM  : ',9L2,2X,9L2)
#endif
  !
#ifdef W3_T
9040 FORMAT ( ' TEST W3SHEL : ODAT   : ',I9.8,I7.6,I7,I9.8,I7.6,  &
       4(/24X,I9.8,I7.6,I7,I9.8,I7.6) )
9041 FORMAT ( ' TEST W3SHEL : FLGRD  : ',20L2)
9042 FORMAT ( ' TEST W3SHEL : IPR, PRFRM : ',6I6,1X,L1)
#endif
  !
#ifdef W3_T
9070 FORMAT ( ' TEST W3SHEL : ',A,3X,2(I10.8,I7.6)/                &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,2(I10.8,I7.6)/                &
       '               ',A,L3,2(I10.8,I7.6)/                &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,2(I10.8,I7.6)/                &
       '               ',A,L3,2(I10.8,I7.6)/                &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,17X,(I10.8,I7.6)/             &
       '               ',A,L3,2(I10.8,I7.6))
9071 FORMAT ( ' TEST W3SHEL : ',A,', DTTST = ',E10.3)
9072 FORMAT ( ' TEST W3SHEL : ',A,3X,2(I10.8,I7.6)/               &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,2(I10.8,I7.6)/               &
       '               ',A,L3,2(I10.8,I7.6)/               &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,2(I10.8,I7.6)/               &
       '               ',A,L3,2(I10.8,I7.6)/               &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,(I10.8,I7.6)/            &
       '               ',A,L3,17X,2(I10.8,I7.6)/           &
       '               ',A,L3,2(I10.8,I7.6))
#endif
  !/
  !/ End of W3SHEL ----------------------------------------------------- /
  !/
END SUBROUTINE FISOC_WW3_RUN

SUBROUTINE FISOC_WW3_FINAL
END SUBROUTINE FISOC_WW3_FINAL

END MODULE FISOC_WW3
