SUBROUTINE OMSAO_main ( exit_value )

  ! ************************************************************************
  !
  ! This is the main program of the SAO Product Generation Executives (PGEs)
  ! for  the Ozone Monitoring Istrument (OMI).
  !
  ! Authors: Thomas P. Kurosu, Kelly Chance
  !          Smithsonian Astrophysical Observatory
  !          60 Garden Street (MS 50)
  !          Cambridge, MA 02138 (USA)
  !
  !          EMail: tkurosu@cfa.harvard.edu
  !                 kchance@cfa.harvard.edu
  !
  ! ************************************************************************

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: pge_bro_idx
  USE OMSAO_parameters_module
  USE OMSAO_variables_module,  ONLY: pge_idx
  USE OMSAO_errstat_module
  USE OMSAO_he5_datafields_module
  USE OMSAO_slitfunction_module

  IMPLICIT NONE

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (OUT) :: exit_value


  ! -------------------------
  ! Name of module/subroutine
  ! -------------------------
  CHARACTER (LEN=10), PARAMETER :: modulename = 'OMSAO_main'

  ! ------------------------------------------------------------------
  ! The general PGE error status variable. This is a relatively late 
  ! addition and is not used consistently in all routines yet. However,
  ! it is envisaged that it will be used ubiquitously througout the
  ! PGE once the PGE developer gets around to implementing it as such.
  ! ------------------------------------------------------------------
  INTEGER   (KIND=i4)      :: errstat, pge_error_status

  ! --------------------------------------------------------------------
  ! Maximum number of points in any reference spectrum. This is used for
  ! automatic (i.e., SUBROUTINE argument) memory allocation and hence is
  ! defined here instead of in a MODULE.
  ! --------------------------------------------------------------------
  INTEGER (KIND=i4) :: n_max_rspec

  ! ----------------------------
  ! Set PGE_ERROR_STATUS to O.K.
  ! ----------------------------
  pge_error_status = pge_errstat_ok

  ! ----------------------------------------------------------------------------
  CALL unbufferSTDout()                       ! Make PGE write STD/IO unbuffered
  ! ----------------------------------------------------------------------------

  errstat = pge_errstat_ok
  ! ---------------------------------------------------------------------------
  CALL read_pcf_file ( errstat )   ! Read PCF file
  ! ---------------------------------------------------------------------------
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"READ_PCF_FILE.", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) GOTO 666

  errstat = pge_errstat_ok
  ! ------------------------------------------------------------
  CALL init_metadata ( errstat )  ! Initialize MetaData
  ! ------------------------------------------------------------
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"INIT_METADATA.", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666

  errstat = pge_errstat_ok
  ! ----------------------------------------------------------------------------------------
  CALL read_reference_spectra ( pge_idx, n_max_rspec, errstat )     ! Read reference spectra
  ! ----------------------------------------------------------------------------------------
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"READ_REFERENCE_SPECTRA.", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666

  errstat = pge_errstat_ok

  ! ------------------------------------------------------------------------------------
  CALL omi_slitfunc_read ( errstat )                   ! Read OMI slit function
  ! ------------------------------------------------------------------------------------
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"OMI_SLITFUNC_READ.", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666

  ! ---------------------------------------------
  ! Set number of InputPointers and InputVersions
  ! ---------------------------------------------
  CALL set_input_pointer_and_versions ( pge_idx )

  errstat = pge_errstat_ok
  ! --------------------------------------------------------------------------------------------
  CALL omi_pge_fitting_process  ( pge_idx, n_max_rspec, errstat )   ! Where all the work is done
  ! --------------------------------------------------------------------------------------------
  CALL error_check ( errstat, pge_errstat_warning, errstat, OMSAO_A_SUBROUTINE, &
       modulename//f_sep//"OMI_PGE_FITTING_PROCESS.", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666

  ! ------------------------------------
  ! Write END_OF_RUN message to log file
  ! ------------------------------------
666 CALL pge_error_status_exit ( pge_error_status, exit_value )

  RETURN
END SUBROUTINE OMSAO_main
