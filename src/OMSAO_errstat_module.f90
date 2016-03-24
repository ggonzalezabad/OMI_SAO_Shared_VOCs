MODULE OMSAO_errstat_module

  ! =================================================================
  !
  ! This module defines variables associated with error handling. It
  ! also loads/includes all (SDPTK) files that define error messages
  ! and generally deal with error handling.
  !
  ! =================================================================

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen

  IMPLICIT NONE

  ! Toolkit include files (sytem)
  INCLUDE 'PGS_SMF.f'
  INCLUDE 'PGS_PC.f'
  INCLUDE 'PGS_PC_9.f'
  INCLUDE 'PGS_IO.f'

  ! Toolkit include files (PGE specific)
  INCLUDE 'PGS_OMSAO_52500.f'

  ! OMI L1B reader include file
  INCLUDE 'PGS_OMI_1900.f'


  ! ----------------------------------------------------
  ! String carrying a potential error or warning message
  ! ----------------------------------------------------
  CHARACTER (LEN=maxchlen) :: pge_errstat_msg

  ! ----------------
  ! Some error stati
  ! ----------------
  INTEGER (KIND=i4), PARAMETER :: &
       pge_errstat_ok    = 0, pge_errstat_warning = 1, &
       pge_errstat_error = 2, pge_errstat_fatal   = 3
  
  ! -----------------------------------------------------------------
  ! Status variables for READ process; used to identify status of
  ! current READ requests, like "O.K.", "FAILED", "END OF DATA", etc.
  ! -----------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       file_read_ok = 0, file_read_failed = 1, file_read_missing = 3, &
       file_read_eof = -1

  ! --------------------------------------------
  ! Status parameters for HE5 interface routines.
  ! --------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: he5_stat_ok = 0, he5_stat_fail = -1


  ! --------------------------
  ! Verbosity Threshold Levels
  ! --------------------------
  ! NOTE: Maximum accepted level in OMI_smf_setmsg seems to be "2"
  ! --------------------------------------------------------------
  ! vb_lev_default   0  default (one line only, at completion)
  ! vb_lev_omidebug  1  OMI support team debugging
  ! vb_lev_stmdebug  2  STM debugging
  ! vb_lev_develop   3  PGE development
  ! vb_lev_lt1mb     4  LOG files of size <= 1MB
  ! vb_lev_gt1mb     5  LOG files of size >  1MB
  ! vb_lev_screen    6  Screen output (only for PGE developer)
  INTEGER (KIND=i4), PARAMETER :: &
       vb_lev_default = 0, vb_lev_omidebug = 1, vb_lev_stmdebug = 2, &
       vb_lev_develop = 3, vb_lev_1mb      = 4, vb_lev_gt1mb    = 5, vb_lev_screen = 6

  ! ------------------------------------------------
  ! Separator string for Error Message concatination
  ! ------------------------------------------------
  CHARACTER (LEN=5), PARAMETER :: f_sep = " --- "

  ! -------------------
  ! External subroutine
  ! -------------------
  INTEGER (KIND=i4), EXTERNAL :: OMI_SMF_setmsg


CONTAINS

  SUBROUTINE error_check ( &
       iserr, errref, severity, errpoint, addmsg, vblev, errstat )

    IMPLICIT NONE

    ! ---------------------------------------------------------------------
    ! Explanation of subroutine arguments:
    !
    !    iserr ......... error status returned from the last action
    !    errref ........ value to compare locerrstat to
    !    severity ...... severity if locerr is /= "O.K."
    !    errpoint ...... pre-defined error message index from OMSAO_52500.t
    !    addmsg ........ message to add for more details on error
    !    vblev ......... verbosity threshold level
    !    errstat ....... error status returned from the subroutine
    ! ---------------------------------------------------------------------

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER   (KIND=i4), INTENT (IN) :: iserr, errref, severity, vblev, errpoint
    CHARACTER (LEN=*),   INTENT (IN) :: addmsg

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! --------------
    ! Local variable
    ! --------------
    INTEGER   (KIND=i4) :: estat, locerr
    CHARACTER (LEN=300) :: addsev


    ! ---------------------------------------------------
    ! Check the current error value against the reference
    ! ---------------------------------------------------
    IF ( iserr > errref ) THEN   ! Used to be "/="
       locerr = severity
    ELSE
       locerr = pge_errstat_ok
    END IF

    addsev = TRIM(ADJUSTL(addmsg))
    SELECT CASE ( locerr )
    CASE ( pge_errstat_ok )
       ! ------------------------------------------------------------
       ! Nothing to be added here. Otherwise we would get "Severity"
       ! additions to "PROGRESS" messages, which is not what we want.
       ! ------------------------------------------------------------
    CASE ( pge_errstat_warning )
       addsev = TRIM(ADJUSTL(addsev)) !//' Severity is: WARNING'
    CASE ( pge_errstat_error )
       addsev = TRIM(ADJUSTL(addsev))//' Severity is: ERROR'
    CASE ( pge_errstat_fatal )
       addsev = TRIM(ADJUSTL(addsev))//' Severity is: FATAL'
    CASE DEFAULT
       addsev = TRIM(ADJUSTL(addsev))//' Severity is: UNKNOWN'
    END SELECT

    IF( locerr /= pge_errstat_ok ) THEN
       estat = OMI_SMF_setmsg( errpoint, TRIM(ADJUSTL(addsev)), " ", vblev )
       errstat = MAX ( errstat, severity )
    END IF

    RETURN
  END SUBROUTINE error_check


  SUBROUTINE pge_error_status_exit ( pge_error_status, exit_value )

    USE OMSAO_parameters_module, ONLY: maxchlen
    IMPLICIT NONE

    ! ---------------------------------------------------------------------
    ! Explanation of subroutine arguments:
    !
    !    pge_error_status ........ final error status at end of PGE run
    ! ---------------------------------------------------------------------

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_error_status

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4), INTENT (OUT) :: exit_value

    ! --------------
    ! Local variable
    ! --------------
    INTEGER   (KIND=i4)      :: estat
    CHARACTER (LEN=maxchlen) :: errmsg, ecode

    ! -------------------
    ! External subroutine
    ! -------------------
    CHARACTER (LEN=maxchlen), EXTERNAL :: int2string

    exit_value = 0

    SELECT CASE ( pge_error_status )
    CASE ( pge_errstat_ok )
       ! -----------------------------------------------------------------
       ! PGE execution completed successfully. All is well.
       ! -----------------------------------------------------------------
       estat = OMI_SMF_setmsg ( OMSAO_S_ENDOFRUN, " ", " ", vb_lev_default )
       exit_value = 0
    CASE ( pge_errstat_warning )
       ! -----------------------------------------------------------------
       ! PGE execution raised non-terminal warnings. Nothing serious, we
       ! hope, so execution completed but with a non-zero exit status.
       ! -----------------------------------------------------------------
       estat = OMI_SMF_setmsg ( OMSAO_W_ENDOFRUN, " ", " ", vb_lev_default )
       exit_value = 0
    CASE ( pge_errstat_error )
       ! -----------------------------------------------------------------
       ! PGE execution encountered a non-fatal error.
       ! -----------------------------------------------------------------
       estat = OMI_SMF_setmsg ( OMSAO_E_ENDOFRUN, " ", " ", vb_lev_default )
       exit_value = 0
    CASE ( pge_errstat_fatal )
       ! -----------------------------------------------------------------
       ! PGE execution encountered a fatal error that lead to termination
       ! -----------------------------------------------------------------
       estat = OMI_SMF_setmsg ( OMSAO_F_ENDOFRUN, " ", " ", vb_lev_default )
       exit_value = 1
    CASE DEFAULT
       ! -----------------------------------------------------------------
       ! If we ever reach here, then PGE_ERRSTAT has been set to a funny
       ! value. This should never happen, but we buffer this case anyway.
       ! -----------------------------------------------------------------
       ecode = int2string ( estat, 1 )
       errmsg = 'Exit code was "'//TRIM(ADJUSTL(ecode))//'".'
       estat = OMI_SMF_setmsg ( OMSAO_U_ENDOFRUN, errmsg, " ", vb_lev_default )
       exit_value = 1
    END SELECT

    RETURN
  END SUBROUTINE pge_error_status_exit

END MODULE OMSAO_errstat_module
