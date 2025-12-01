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
! This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
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

!
! This is the wave model spcific code for FISOC.  The main purpose is to transfer information 
! between ESMF structures and the WM's internal structures.
!

MODULE FISOC_WM_Wrapper
  
  USE ESMF
  
  USE FISOC_utils_MOD
  USE FISOC_types_MOD
  USE FISOC_WW3, ONLY: FISOC_WW3_INIT, FISOC_WW3_RUN, FISOC_WW3_FINAL
  
  USE MOD_PAR
  USE W3ADATMD, ONLY : HS,FP0,THM,WLM,NSEALM
  USE W3IDATMD, ONLY : WLEV,CXN,CYN
  USE W3GDATMD, ONLY : XGRD,YGRD         !NGRIDS,IGRID,NSEA,NSEAL
  USE W3GDATMD, ONLY: NX,NTRI
  use yowNodepool, only: npa,np,ng,iplg,ipgl
  use yowElementpool, only : INE, NE, ielg
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_WM_Wrapper_Init_Phase1,  FISOC_WM_Wrapper_Init_Phase2,  &
       FISOC_WM_Wrapper_Run, FISOC_WM_Wrapper_Finalize

  
  ! TODO: don't use mpic for anything other than passing to FVCOM.  All mpi operations 
  ! should be through the ESMF VM commands.
  INTEGER, SAVE                :: mpic
  INTEGER, ALLOCATABLE, SAVE   :: ownedNodeIDs(:) ! subset of local nodes owned by current partition

  ! Route handle for mapping from ESMF fields to WW3 fields (one to many mapping)
  TYPE(ESMF_RouteHandle),SAVE  :: RH_ESMF2WW3

  ! nodal distgrid, an ESMF object holding information about the distribution of 
  ! WW3 nodes across partitions, needed for the one to many mapping.
  TYPE(ESMF_distgrid),SAVE     :: distgridWW3

CONTAINS
  
  !--------------------------------------------------------------------------------------
  ! The first phase of initialisation is mainly to initialise the wave model, and access 
  ! grid and variable initial information.
  SUBROUTINE FISOC_WM_Wrapper_Init_Phase1(FISOC_config,vm,WM_ExpFB,WM_mesh,rc)
    
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)              :: vm ! ESMF virtual machine (parallel context)
    TYPE(ESMF_mesh),INTENT(OUT)           :: WM_mesh
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: WM_ExpFB
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                               :: localPet, petCount
    CHARACTER(len=ESMF_MAXSTR)            :: label
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: WM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: WM_configFile, WM_stdoutFile
    LOGICAL                               :: verbose_coupling, first

    first = .TRUE.

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, WM_stdoutFile, label='WM_stdoutFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    WRITE (WM_stdoutFile, "(a,I0)") TRIM(WM_stdoutFile), localPet
    OPEN(unit=WM_outputUnit, file=WM_stdoutFile, STATUS='REPLACE', ERR=101)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "Initialising WW3"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    CALL ESMF_ConfigGetAttribute(FISOC_config, WM_configFile, label='WM_configFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"*******************************************************************************"
       PRINT*,"**********       WM wrapper.  Init phase 1 method.         ********************"
       PRINT*,"********** Creating and registering ESMF WW3 components. ********************"
       PRINT*,"*******************************************************************************"
       PRINT*,""
    END IF

    ! We will pass an MPI communicator to WW3
    CALL FISOC_VM_MPI_Comm_dup(vm,mpic,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (localPet.EQ.0) THEN
       WRITE (WM_outputUnit,*) 'FISOC is about to call WW3 init method.'
       WRITE (WM_outputUnit,*) 'mpic:',mpic
    END IF
    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_WW3_INIT(MPI_COMM_opt=mpic)

    IF (localPet.EQ.0) THEN
      WRITE (WM_outputUnit,*) 'FISOC has just called WW3 init method.'
    END IF
    
    CALL WW32ESMF_mesh(FISOC_config,WM_mesh,vm,rc=rc)

    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    

    ! extract a list of required wave variables from the configuration object
    label = 'FISOC_WM_ReqVars:' ! the FISOC names for the vars
    CALL FISOC_getListFromConfig(FISOC_config, label, WM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_populateFieldBundle(WM_ReqVarList,WM_ExpFB,WM_mesh, &
         init_value=FISOC_missingData,                             &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
 
    CALL getFieldDataFromWM(WM_ExpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    RETURN
    
101 msg = "OM failed to open stdoutFile "//WM_stdoutFile
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  END SUBROUTINE FISOC_WM_Wrapper_Init_Phase1
  
  
  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_WM_Wrapper_Init_Phase2(FISOC_config,vm,WM_ImpFB,WM_ExpFB,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: WM_ImpFB, WM_ExpFB
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    TYPE(ESMF_VM),INTENT(IN)              :: vm
    
    LOGICAL   :: verbose_coupling, WM_initCavityFromOM, OM2WM_init_vars
    INTEGER   :: localpet, urc

    rc = ESMF_FAILURE
    
    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"**********      WM wrapper.  Init phase 2 method.        *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"Here we have access to the initialised OM fields, just in case the WM needs "
       PRINT*,"to know about these in order to complete its initialisation."
       PRINT*,""
    END IF

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM2WM_init_vars, 'OM2WM_init_vars',rc=rc) 
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (OM2WM_init_vars) THEN
       CALL sendFieldDataToWM(WM_ImpFB,FISOC_config,vm,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    
    CALL getFieldDataFromWM(WM_ExpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) & 
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  END SUBROUTINE FISOC_WM_Wrapper_Init_Phase2


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_WM_Wrapper_Run(FISOC_config,vm,WM_ExpFB,WM_ImpFB,rc_local)
    
    TYPE(ESMF_config),INTENT(INOUT)                :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT),OPTIONAL  :: WM_ExpFB, WM_ImpFB 
    TYPE(ESMF_VM),INTENT(IN)                       :: vm
    INTEGER,INTENT(OUT),OPTIONAL                   :: rc_local

    INTEGER                    :: localPet, rc
    LOGICAL                    :: verbose_coupling
    TYPE(ESMF_field)           :: WM_dTdz_l0,WM_z_l0, WM_bmb
    REAL(ESMF_KIND_R8),POINTER :: WM_dTdz_l0_ptr(:,:), WM_z_l0_ptr(:,:), WM_bmb_ptr(:,:)
    INTEGER                    :: WM_dt_sec
    REAL(ESMF_KIND_R8)         :: WM_dt_sec_float
    TYPE(ESMF_TimeInterval)    :: WM_dt
    TYPE(ESMF_time)            :: interval_startTime, interval_endTime
    CHARACTER(len=ESMF_MAXSTR) :: interval_startTime_char, interval_endTime_char

    INTEGER, ALLOCATABLE :: TEND(:,:)
    INTEGER :: I


    rc_local = ESMF_FAILURE
    
    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (PRESENT(WM_ImpFB)) THEN       
       CALL sendFieldDataToWM(WM_ImpFB,FISOC_config,vm,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! Get WM time interval (how long to run WM for) in seconds...
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, WM_dt_sec, label='WM_dt_sec', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    WM_dt_sec_float = REAL(WM_dt_sec,ESMF_KIND_R8)

    ! ...convert to ESMF type...
    CALL ESMF_TimeIntervalSet(WM_dt, s=WM_dt_sec, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... use to create start and end times...
    interval_startTime = FISOC_time
    interval_endTime   = FISOC_time + WM_dt
    

    ! ... and convert these to character strings.
    CALL ESMF_TimeGet(interval_startTime, timeStringISOFrac=interval_startTime_char, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_TimeGet(interval_endTime, timeStringISOFrac=interval_endTime_char, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    IF (localPet.EQ.0) THEN
      WRITE (OM_outputUnit,*) 'FISOC is about to call WW3 run method, period (sec): ',WM_dt_sec_float
      WRITE (OM_outputUnit,*) 'FISOC is about to call WW3 start time: ',trim(interval_startTime_char)
      WRITE (OM_outputUnit,*) 'FISOC is about to call WW3 end   time: ',trim(interval_EndTime_char)
      IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
        msg = "Calling WW3 run method now."
        CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
             line=__LINE__, file=__FILE__, rc=rc)
      END IF
    END IF

    CALL ESMF_VMBarrier(vm, rc=rc)
    CALL FISOC_WW3_RUN(START_DATE_opt=interval_startTime_char,END_DATE_opt=interval_endTime_char)
    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (localPet.EQ.0) THEN
    
!      WRITE (WM_outputUnit,*) 'FISOC has just called WW3 run method.'
  
    END IF
    
    IF (PRESENT(WM_ExpFB)) THEN
      CALL getFieldDataFromWM(WM_ExpFB,FISOC_config,vm,rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************         WM wrapper.  Run method.           **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""

       IF ((PRESENT(WM_ExpFB)).AND.(.NOT.(PRESENT(WM_ImpFB)))) THEN
          PRINT*,"We have no new inputs for the WM from the WM.  We need to call the WM "
          PRINT*,"and record its output in the WM export field bundle."
       END IF

       IF ((PRESENT(WM_ExpFB)).AND.(PRESENT(WM_ImpFB))) THEN
          PRINT*,"We have new inputs for the WM from the WM in the WM import field bundle. "
          PRINT*,"We need to send these inputs to the WM, run one timestep, and record the WM "
          PRINT*,"outputs. "
       END IF
       
       IF ((.NOT.(PRESENT(WM_ExpFB))).AND.(PRESENT(WM_ImpFB))) THEN
          PRINT*,"We have new inputs for the WM from the WM in the WM import field bundle. "
          PRINT*,"We need to send these inputs to the WM and run one timestep. We do not "
          PRINT*,"need to collect WM outputs.  Just run the WM one timestep."
       END IF

       IF ((.NOT.(PRESENT(WM_ExpFB))).AND.(.NOT.(PRESENT(WM_ImpFB)))) THEN
          PRINT*,"We have no new inputs for the WM from the WM, and we do not need to "
          PRINT*,"collect WM outputs.  Just run the WM one timestep."
       END IF

    END IF

    rc_local = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_WM_Wrapper_Run


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_WM_Wrapper_Finalize(FISOC_config,vm,rc)

    TYPE(ESMF_config),INTENT(INOUT)    :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)           :: vm
    INTEGER,INTENT(OUT),OPTIONAL       :: rc

    INTEGER                            :: localPet
    LOGICAL                            :: verbose_coupling

    rc = ESMF_FAILURE


    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************       WM wrapper.  Finalise method.        **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"FISOC has taken care of clearing up ESMF types.  Here we just need to call the "
       PRINT*,"WM finalise method."
    END IF

    IF (localPet.EQ.0) THEN
       WRITE (WM_outputUnit,*) 'FISOC is about to call SWAN finalize.'
    END IF

    CALL FISOC_WW3_FINAL
    IF (localPet.EQ.0) THEN
       WRITE (WM_outputUnit,*) 'FISOC has just called SWAN finalize.'
    END IF

    IF (localPet.EQ.0) THEN
       WRITE (WM_outputUnit,*) 'FISOC has just called SWAN finalize.'
    END IF

!    CLOSE(unit=WM_outputUnit, ERR=102)
    
    rc = ESMF_SUCCESS
    
    RETURN
    
102 msg = "WM failed to close stdoutFile"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  END SUBROUTINE FISOC_WM_Wrapper_Finalize
  

  !--------------------------------------------------------------------------------------
  ! update the fields in the wave export field bundle from the WM
  !--------------------------------------------------------------------------------------
  SUBROUTINE getFieldDataFromWM(WM_ExpFB,FISOC_config,vm,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: WM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    TYPE(ESMF_VM),INTENT(IN)              :: vm

    INTEGER                               :: fieldCount, localPet, petCount
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:)
    INTEGER                               :: IstrR, IendR, JstrR, JendR ! tile start and end coords
    INTEGER                               :: ii, jj, nn
!    INTEGER                               :: LBi, UBi, LBj, UBj ! tile start and end coords including halo

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        
    ! get a list of fields and their names from the WM export field bundle
    fieldCount = 0
    CALL ESMF_FieldBundleGet(WM_ExpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(WM_ExpFB, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    fieldLoop: DO nn = 1,fieldCount
       
       CALL ESMF_FieldGet(fieldList(nn), name=fieldName, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL ESMF_FieldGet(fieldList(nn), farrayPtr=ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       ptr = FISOC_missingData

       SELECT CASE (TRIM(ADJUSTL(fieldName)))
       
       CASE ('WM_HS')
         DO ii = 1,SIZE(ptr)
           ptr(ii) = HS(ownedNodeIds(ii))
         END DO
       CASE ('WM_DIR')
         DO ii = 1,SIZE(ptr)
           ptr(ii) = THM(ownedNodeIds(ii))
         END DO
       CASE ('WM_TPEAK')
         DO ii = 1,SIZE(ptr)
           ptr(ii) = FP0(ownedNodeIds(ii))
         END DO
       CASE ('WM_TMBOT')
         DO ii = 1,SIZE(ptr)
           ptr(ii) = HS(ownedNodeIds(ii))
         END DO
       CASE ('WM_UBOT')
         DO ii = 1,SIZE(ptr)
           ptr(ii) = HS(ownedNodeIds(ii))
         END DO
       CASE ('WM_WLEN')
         DO ii = 1,SIZE(ptr)
           ptr(ii) = WLM(ownedNodeIds(ii))
         END DO
       CASE ('WM_DIRBOT')
         DO ii = 1,SIZE(ptr)
           ptr(ii) = HS(ownedNodeIds(ii))
         END DO
       CASE ('WM_QB1') 
         DO ii = 1,SIZE(ptr)
           ptr(ii) = HS(ownedNodeIds(ii))
         END DO
       CASE DEFAULT
         msg = "ERROR: unknown variable: "//TRIM(ADJUSTL(fieldName))
         CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
              line=__LINE__, file=__FILE__, rc=rc)
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
         
       END SELECT
       
       IF (ASSOCIATED(ptr)) THEN
         NULLIFY(ptr)
       END IF
       
    END DO fieldLoop
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE GetFieldDataFromWM
  
  !--------------------------------------------------------------------------------------
  ! update the fields in the WM from the wave import field bundle
  !--------------------------------------------------------------------------------------
  SUBROUTINE sendFieldDataToWM(WM_ImpFB,FISOC_config,vm,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: WM_ImpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    TYPE(ESMF_VM),INTENT(IN)              :: vm

    INTEGER                               :: fieldCount, localPet, petCount
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:)
    INTEGER                               :: IstrR, IendR, JstrR, JendR ! tile start and end coords
    INTEGER                               :: ii, jj, nn

    REAL, ALLOCATABLE, TARGET :: CXN_TMP(:,:),CYN_TMP(:,:),WLEV_TMP(:,:)
    integer, dimension(:), allocatable :: ivertp   ! vertex index of global grid in own subdomain (without ghost vertices)
    INTEGER                               :: j, k

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        
    ! get a list of fields and their names from the WM export field bundle
    fieldCount = 0
    CALL ESMF_FieldBundleGet(WM_ImpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(WM_ImpFB, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! determine node number per subdomain
    !
    allocate(ivertp(NP))
    k = 0
    do j = 1, NP       !NPA
       if ( iplg(j) > 0 ) then
          k = k + 1
          ivertp(k) = iplg(j) 
       endif
    enddo
    
    fieldLoop: DO nn = 1,fieldCount
      
      CALL ESMF_FieldGet(fieldList(nn), name=fieldName, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      IF (FISOC_OM2WM(fieldName,FISOC_config,rc=rc)) THEN

        CALL FISOC_ArrayRedistFromField(RH_ESMF2WW3,fieldList(nn),distgridWW3,ptr)
	
        SELECT CASE (TRIM(ADJUSTL(fieldName)))
         
       CASE ('OM_UU')
         !FLOW_UU(1:size(ptr)) = ptr
         ALLOCATE(CXN_TMP(NPA,2))
         do ii=1,size(ptr)
         CXN_TMP(ii,1) = ptr(ii)
         CXN_TMP(ii,2) = ptr(ii)
	 end do
	 CXN => CXN_TMP
       DEALLOCATE(CXN_TMP)
       CASE ('OM_VV')
         ALLOCATE(CYN_TMP(NPA,2))
	 do ii=1,size(ptr)
         CYN_TMP(ii,1) = ptr(ii)
         CYN_TMP(ii,2) = ptr(ii)
	 end do
	 CYN => CYN_TMP
       DEALLOCATE(CYN_TMP)
       CASE ('OM_ETA')
         !FLOW_ETA(1:size(ptr)) = ptr
         ALLOCATE(WLEV_TMP(NPA,2))
         do ii=1,size(ptr)
         WLEV_TMP(ii,1) = ptr(ii)
         WLEV_TMP(ii,2) = ptr(ii)
	 end do
	 WLEV => WLEV_TMP
       DEALLOCATE(WLEV_TMP)
       CASE ('OM_HH')
!JQI         ADCIRC_HH2(1:size(ptr)) = ptr
       CASE DEFAULT
         msg = "ERROR: unknown variable: "//TRIM(ADJUSTL(fieldName))
         CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
              line=__LINE__, file=__FILE__, rc=rc)
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
         
       END SELECT
 
       IF (ASSOCIATED(ptr)) THEN
         NULLIFY(ptr)
       END IF
       
     END IF
      
    END DO fieldLoop
    deallocate(ivertp)
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE SendFieldDataToWM
  


  !--------------------------------------------------------------------------------
  ! Extract the WW3 mesh information and use it to create an ESMF_mesh object
  SUBROUTINE WW32ESMF_mesh(FISOC_config,ESMF_WW3mesh,vm,rc)
  
    TYPE(ESMF_config),INTENT(INOUT)  :: FISOC_config
    TYPE(ESMF_mesh),INTENT(INOUT)    :: ESMF_WW3Mesh
    TYPE(ESMF_VM),INTENT(IN)         :: vm
    INTEGER,INTENT(OUT),OPTIONAL     :: rc
    type(ESMF_DistGrid)      :: elementDistgrid, distgrid
    INTEGER                          :: ii, nn, IERR, k, j
    CHARACTER(len=ESMF_MAXSTR)       :: subroutineName = "WW32ESMF_mesh", WM_coords
    INTEGER,ALLOCATABLE              :: elemTypes(:), elemIds(:),elemConn(:), localNodeIDs(:)
    INTEGER,ALLOCATABLE              :: nodeIds(:),nodeOwners(:),nodeOwnersGL(:),nodeOwnersGL_recv(:)
    INTEGER,ALLOCATABLE              :: nodeMask(:),testids(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE   :: nodeCoords(:) 
    INTEGER                          :: localPet, petCount 
    INTEGER                          :: WW3_numNodes, WW3_numElems 
    LOGICAL                          :: verbose_coupling
    integer                          :: nownv    ! number of vertices in own subdomain (without ghost vertices)
    integer, dimension(:), allocatable :: ivertp   ! vertex index of global grid in own subdomain (without ghost vertices)
    real   , dimension(:), allocatable :: xpl      ! user coordinates grid points (without ghost vertices)
    real   , dimension(:), allocatable :: ypl      ! user coordinates grid points (without ghost vertices)

  integer               :: mynp,myne
    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    CALL FISOC_ConfigDerivedAttribute(FISOC_config, WM_coords, label='WM_coords', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       msg = "WW3 mesh: creating ESMF mesh structure"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
  
    !-------------------------------------------------------------------!
    ! Use ww3 grid variables directly, may need to use MODULE ALL_VARS
    ! and MODULE LIMS from WW3 library, see mod_main.F for details
    ! and module par for retreating  EGID NGID
    WW3_numNodes        = NPA    
    WW3_numElems        = NE     
!TODO: just use FVCOM vars directly, no need for intermediate vars, i.e. remove FVCOM_numNodes etc

    nownv = NP 
    !
    allocate(ivertp(nownv))
    allocate(   xpl(nownv))
    allocate(   ypl(nownv))

    ALLOCATE(localNodeIds(nownv))

    ALLOCATE(nodeIds(WW3_numNodes))
    ALLOCATE(nodeCoords(WW3_numNodes*2))
    ALLOCATE(nodeOwners(WW3_numNodes))
    ALLOCATE(nodeMask(WW3_numNodes))  
  
    ALLOCATE(elemIds(WW3_numElems))
    ALLOCATE(elemConn(WW3_numElems*3))
    ALLOCATE(elemTypes(WW3_numElems))

    ALLOCATE(nodeOwnersGL(NX)) 
    ALLOCATE(nodeOwnersGL_recv(NX)) 

    ! Setup mask to avoid regridding ice variables in the open ocean.
    nodeMask = MASK_OPEN_OCEAN !without mask

    ! Construct a global array of node owners in which an arbitrary decision is 
    ! taken about which partition boundary nodes should belong to.
    nodeOwnersGL = -1
    
    ! determine node number per subdomain
    !
    k = 0
    do j = 1, NP       !NPA
       if ( iplg(j) > 0 ) then
          k = k + 1
          ivertp(k) = iplg(j)                           !ipgl(j)                   !ivertg(j)
          xpl   (k) = XGRD(1,iplg(J))                   !xcugrd(j) + XOFFS
          ypl   (k) = YGRD(1,iplg(J))                   !ycugrd(j) + YOFFS
          nodeOwnersGL(iplg(k))=localPET
          localNodeIds(k) = j
       endif
    enddo
    
    CALL MPI_Allreduce(nodeOwnersGL,nodeOwnersGL_recv,NX,MPI_INT,MPI_MAX,mpic,IERR)

    ! sanity check (we initialised owners to -1, but PET count starts at 0, so if it 
    ! works then all nodes should have been assigned an owner .GE. 0)
    IF (MINVAL(nodeOwnersGL_recv).LT.0) THEN
      msg = "ERROR: Some nodes not assigned owners"
      CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
           line=__LINE__, file=__FILE__, rc=rc)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    !populate nodeOwners from nodeOwnersGL
    DO ii = 1, WW3_numNodes
      nodeOwners(ii) = nodeOwnersGL_recv(int(abs(iplg(ii))))
    END DO
    elemTypes         =  ESMF_MESHELEMTYPE_TRI

    IF(localPET.EQ.0)then
      open(1,file='nodeowners_swan.dat',status='replace')
         do ii=1,ww3_numnodes
           write(1,*)ii,int(abs(iplg(ii))),nodeOwners(ii)
         end do
      close(1)
    endif
    
    ! loop over to get nodeIds nodeCoords
    DO ii = 1, WW3_numNodes
       nn = (ii-1)*2
       nodeIds(ii)      = int(abs(iplg(ii)))
       
       SELECT CASE (WM_coords)

       CASE("spherical","Spherical","SPHERICAL")
          !CALL FISOC_MAPLL(nodeCoords(nn+1),nodeCoords(nn+2),LAT(ii),LON(ii))
          nodeCoords(nn+1) = XGRD(1,iplg(ii))
          nodeCoords(nn+2) = YGRD(1,iplg(ii))
       CASE DEFAULT
          nodeCoords(nn+1) = XGRD(1,iplg(ii))
          nodeCoords(nn+2) = YGRD(1,iplg(ii))

       END SELECT

    END DO

    ! loop over to get elemConn
    DO ii = 1, WW3_numElems
       nn = (ii-1)*3
       elemIds(ii)    = ielg(ii)         
       elemConn(nn+1) = INE(1,ii)
       elemConn(nn+2) = INE(2,ii)
       elemConn(nn+3) = INE(3,ii)
    END DO 

    ! get a list od node ids for the subset of nodes that are locally owned
    ! (to be used in mapping SWAN variables to ESMF fields)
    CALL FISOC_locallyOwnedNodes(localPet,localNodeIds,nodeOwners,ownedNodeIDs)

    SELECT CASE (WM_coords)

    CASE("spherical","Spherical","SPHERICAL")
       !----------------------------------------------------------------!       
       ! Create Mesh structure in 1 step
       ESMF_WW3Mesh = ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
            nodeIds=nodeIds, nodeCoords=nodeCoords,                   &
            nodeOwners=nodeOwners,                                    &
            nodeMask=nodeMask, elementIds=elemIds,                    &
            elementTypes=elemTypes, elementConn=elemConn,             &
            coordSys=ESMF_COORDSYS_SPH_DEG,                              &
            rc=rc)
    CASE DEFAULT
       !----------------------------------------------------------------!       
       ! Create Mesh structure in 1 step
       ESMF_WW3Mesh = ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
            nodeIds=nodeIds, nodeCoords=nodeCoords,                   &
            nodeOwners=nodeOwners,                                    &
            nodeMask=nodeMask, elementIds=elemIds,                    &
            elementTypes=elemTypes, elementConn=elemConn,             &
            coordSys=ESMF_COORDSYS_CART,                              &
            rc=rc)

    END SELECT

    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

   ! output mesh information from swan and esmf
    call ESMF_MeshGet(ESMF_WW3Mesh,numOwnedNodes=mynp,numOwnedElements=myne,elementDistgrid=distgrid,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    allocate(testids(myne))
    call ESMF_DistGridGet(distgrid,localDE=0,seqIndexList=testids,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    print*,'localPET=',localPET, 'esmf   owned nodes,elements',mynp,myne
    print*,'localPET=',localPET, 'ww3 owned nodes,elements',nownv,ne
    print*,'localPET=',localPET, 'ww3   all nodes,elements',npa,ne
    !print*,'localPET=',localPET, 'elementIds',elemIds
    !print*,'localPET=',localPET, 'distgridElementIds',testids
    deallocate(testids)	

    CALL FISOC_CreateOneToManyRouteHandle(ESMF_WW3Mesh,nodeIds(1:),RH_ESMF2WW3,distgridWW3,vm)

    DEALLOCATE(nodeIds)
    DEALLOCATE(nodeCoords)
    DEALLOCATE(nodeOwners)
    DEALLOCATE(nodeMask)
    
    DEALLOCATE(elemIds)
    DEALLOCATE(elemConn)
    DEALLOCATE(elemTypes)	

    !TODO: put this in PET logs
    print*,"FINISHED WW3 MESH CREATION"

  END SUBROUTINE WW32ESMF_mesh

END MODULE FISOC_WM_Wrapper
