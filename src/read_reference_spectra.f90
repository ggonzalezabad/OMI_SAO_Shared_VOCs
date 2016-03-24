SUBROUTINE read_reference_spectra ( pge_idx, n_max_rspec, pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: &
       max_rs_idx, wvl_idx, spc_idx, pge_static_input_luns, &
       pge_o3_idx, o3_t1_idx, o3_t2_idx, o3_t3_idx, comm_idx
  USE OMSAO_parameters_module, ONLY: maxchlen, max_spec_pts, zerospec_string, r8_missval
  USE OMSAO_variables_module,  ONLY: &
       verb_thresh_lev, winwav_min, winwav_max, ReferenceSpectrum, refspecs_original, &
       common_mode_spec, yn_solar_comp, solar_comp_typ, solar_comp_orb,           &
       OMSAO_solcomp_filename, l1br_opf_version, l1b_channel, yn_common_iter
  USE OMSAO_he5_datafields_module, ONLY: o3_prefit_he5fields
  USE OMSAO_solcomp_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=22), PARAMETER :: modulename = 'read_reference_spectra'

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (OUT)   :: n_max_rspec
  INTEGER (KIND=i4), INTENT (INOUT) :: pge_error_status

  ! ------------------------
  ! Error handling variables
  ! ------------------------
  INTEGER (KIND=i4) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: i, npts


  errstat = pge_errstat_ok

  ! -------------------------------------------------------------
  ! This variable will hold the maximum number of spectral points
  ! of any reference spectrum read in.
  ! -------------------------------------------------------------
  n_max_rspec = 0

  ! -----------------------------------------------------------
  ! Read spectra one by one. Skip if name of file is ZEROSPEC
  ! -----------------------------------------------------------

  o3_prefit_he5fields(o3_t1_idx:o3_t3_idx,1:2)%SpecTemp = r8_missval
  DO i = 1, max_rs_idx
     IF ( INDEX (TRIM(ADJUSTL(refspecs_original(i)%FileName)), zerospec_string ) == 0 ) THEN
        errstat = pge_errstat_ok
        SELECT CASE ( i )
        CASE ( comm_idx )
           IF ( .NOT. yn_common_iter ) THEN
              ! --------------------------------------------------------------
              ! The common mode spectrum for OMI has 60 across-track positions
              ! --------------------------------------------------------------
              CALL read_commonmode_spec ( &
                   pge_static_input_luns(i), refspecs_original(i)%FileName, &
                   winwav_min, winwav_max, common_mode_spec, errstat )
              
              ! ----------------------------------------------------------------
              ! REFSPECs_ORIGINAL must know the number of spectral points, since
              ! this is checked in PREPARE_DATABASE, where nPoints==0 would lead
              ! to the exclusion of the reference spectrum (for obvious reasons)
              ! ----------------------------------------------------------------
              npts                                     = common_mode_spec%nPoints
              refspecs_original(i)%nPoints             = npts

              ! ------------------------------------------------------------------------
              ! CAREFUL: The Common Mode has 60 potentially different wavelength arrays.
              !          This assignment here can screw up things!
              ! ------------------------------------------------------------------------
              refspecs_original(i)%RefSpecWavs(1:npts) = common_mode_spec%RefSpecWavs(1,1:npts)
           END IF
        CASE DEFAULT
           CALL read_one_refspec ( &
                i, pge_static_input_luns(i), refspecs_original(i)%FileName, &
                winwav_min, winwav_max, refspecs_original(i), npts, errstat )
        END SELECT

        ! --------------------------------------------
        ! Update the maximum number of spectral points
        ! --------------------------------------------
        n_max_rspec = MAX ( npts, n_max_rspec )

        CALL error_check ( &
             errstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
             modulename//f_sep//TRIM(ADJUSTL(refspecs_original(i)%FileName)),       &
             vb_lev_default, pge_error_status )
        IF ( pge_idx == pge_o3_idx .AND. ( i==o3_t1_idx .OR. i==o3_t2_idx .OR. i==o3_t3_idx ) ) &
             o3_prefit_he5fields(i,1:2)%SpecTemp = refspecs_original(i)%Temperature
     END IF
  END DO

  ! -----------------------------
  ! Read composite solar spectrum
  ! -----------------------------------------------------------------------
  ! Practically, this isn't a reference spectrum of the same ilk, but it is
  ! ingested here since it is used in a similar way than all the others. 
  ! -----------------------------------------------------------------------
  IF ( yn_solar_comp ) THEN
     errstat = pge_errstat_ok
     CALL soco_pars_read ( &
          OMSAO_solcomp_filename, solar_comp_typ, l1b_channel, &
          winwav_min, winwav_max, errstat )
     !v002 CALL soco_pars_read ( &
     !v002      OMSAO_solcomp_filename, solar_comp_typ-1, solar_comp_orb, l1b_channel, &
     !v002      l1br_opf_version, winwav_min, winwav_max, errstat )
     CALL error_check ( &
          errstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
          modulename//f_sep//'Composite Solar Spectrum Parameters',       &
          vb_lev_default, pge_error_status )
  END IF

  ! ------------------------------------
  ! Report successful reading of spectra
  ! ------------------------------------
  CALL error_check ( 0, 1, pge_errstat_ok, OMSAO_S_READ_REFSPEC_FILE, &
       modulename, vb_lev_stmdebug, pge_error_status )

  RETURN
END SUBROUTINE read_reference_spectra

SUBROUTINE read_one_refspec ( &
     rs_idx, omi_lun, specname, winwav_min, winwav_max, refspec_orig, nspec, errstat )

  USE OMSAO_precision_module,   ONLY: r8
  USE OMSAO_indices_module,     ONLY: wvl_idx, spc_idx, ring_idx
  USE OMSAO_parameters_module,  ONLY: maxchlen, max_spec_pts, lm_start_of_table
  USE OMSAO_variables_module,   ONLY: verb_thresh_lev, ReferenceSpectrum, l1b_channel
  USE OMSAO_errstat_module
  
  IMPLICIT NONE

  ! ----------------
  ! Input Parameters
  ! ----------------
  INTEGER   (KIND=i4), INTENT (IN) :: rs_idx, omi_lun
  CHARACTER (LEN=*),   INTENT (IN) :: specname
  REAL      (KIND=r8), INTENT (IN) :: winwav_min, winwav_max

  ! -----------------
  ! Output Parameters
  ! -----------------
  INTEGER (KIND=i4),        INTENT (OUT)   :: nspec
  INTEGER (KIND=i4),        INTENT (INOUT) :: errstat
  TYPE (ReferenceSpectrum), INTENT (OUT)   :: refspec_orig
  
  ! ----------------
  ! Local Variables
  ! ----------------
  INTEGER (KIND=i4) :: i, ios, funit, file_read_stat, j1, j2, nskip
  INTEGER (KIND=i4), DIMENSION (max_spec_pts) :: irev
  REAL    (KIND=r8), DIMENSION (max_spec_pts) :: x, y, xtmp, ytmp
  CHARACTER (LEN=maxchlen)                    :: lastline, rs_title, rs_units
  REAL    (KIND=r8)                           :: xdum, rs_temp, ddum, specnorm

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=16), PARAMETER :: modulename = 'read_one_refspec'

  ! ------------------------
  ! Error handling variables
  ! ------------------------
  INTEGER (KIND=i4) :: version, locerrstat

  ! ---------------------------------
  ! External OMI and Toolkit routines
  ! ---------------------------------
  INTEGER (KIND=i4), EXTERNAL :: &
       pgs_smf_teststatuslevel, pgs_io_gen_openf, pgs_io_gen_closef


  locerrstat = pge_errstat_ok

  nspec = 0; rs_temp = 0.0_r8
  x = 0.0_r8 ; y = 0.0_r8

  ! -----------------------
  ! Open reference spectrum
  ! -----------------------
  version = 1
  locerrstat = PGS_IO_GEN_OPENF ( omi_lun, PGSd_IO_Gen_RSeqFrm, 0, funit, version )
  locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
  CALL error_check ( &
       locerrstat, pgs_smf_mask_lev_s, pge_errstat_error, OMSAO_E_OPEN_REFSPEC_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
  IF (  errstat /= pge_errstat_ok ) RETURN

  ! --------------------------------------
  ! Skip comments header to start of table
  ! --------------------------------------
  nskip = 0; file_read_stat = 0
  CALL skip_headerlines ( funit, '#', 'read', nskip, lastline, errstat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
  IF (  errstat /= pge_errstat_ok ) RETURN
  
  ! -----------------------------------------------------------------------
  ! Read dimension, start&end wavelength, and norm of spectrum; note that 
  ! the title of the spectrum has already been read as the last of the 
  ! "skipped" lines.
  ! -----------------------------------------------------------------------
  rs_title = TRIM(ADJUSTL(lastline))
  READ (UNIT=funit, FMT='(A)', IOSTAT=file_read_stat) rs_units
  READ (UNIT=funit, FMT=*, IOSTAT=file_read_stat) rs_temp, specnorm, ddum, ddum, nspec
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
  IF (  errstat /= pge_errstat_ok ) RETURN

  ! ---------------------------------------------------------------------
  ! Find first and last index to read, based on WINWAV_MIN and WINWAV_MAX
  ! ---------------------------------------------------------------------
  j1 = 1  ;  j2 = 0 ; xdum = 0.0_r8
  getidx: DO i = 1, nspec
     READ (UNIT=funit, FMT=*, IOSTAT=ios) xdum
     IF ( i /= nspec .AND. ios /= 0 ) THEN
        CALL error_check ( &
             ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
             modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
        IF (  errstat /= pge_errstat_ok ) RETURN
     END IF
     IF ( xdum <= winwav_min )               j1 = i
     IF ( xdum >= winwav_max .AND. j2 == 0 ) j2 = i
     IF ( j1 > 0 .AND. j2 > 0 ) EXIT getidx
  END DO getidx
  IF ( j2 == 0 ) j2 = nspec

  ! -------------------------------------------
  ! Check for maximum number of spectral points
  ! -------------------------------------------
  IF ( j2-j1+1 > max_spec_pts ) THEN
     CALL error_check ( &
          0, 0, pge_errstat_error, OMSAO_E_REFSPEC_MAXPTS, &
          modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
     IF (  errstat /= pge_errstat_ok ) RETURN
  END IF
  ! --------------------------------------------------------------------------
  ! Rewind file, skip to start of table line, and re-read first entry as dummy
  ! --------------------------------------------------------------------------
  REWIND ( funit )

  ! --------------------------------------
  ! Skip comments header to start of table
  ! --------------------------------------
  CALL skip_headerlines ( funit, '#', 'skip', nskip, lastline, errstat )

  ! -----------------------------------------------
  ! Read dimension, start&end, and norm of spectrum
  ! -----------------------------------------------
  READ (UNIT=funit, FMT='(A)', IOSTAT=file_read_stat) rs_title
  READ (UNIT=funit, FMT='(A)', IOSTAT=file_read_stat) rs_units
  READ (UNIT=funit, FMT=*,     IOSTAT=file_read_stat) rs_temp, specnorm, ddum, ddum, nspec
  CALL error_check ( &
       ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
  IF (  errstat /= pge_errstat_ok ) RETURN

  ! -------------------------------------------------
  ! Skip first J1-1 lines, then read J1 to J2 entries
  ! -------------------------------------------------
  DO i = 1, j1-1
     READ (UNIT=funit, FMT=*, IOSTAT=ios) xdum
     IF ( i /= nspec .AND. ios /= 0 ) THEN
        CALL error_check ( &
             ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
             modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
        IF (  errstat /= pge_errstat_ok ) RETURN
     END IF
  END DO
  DO i = j1, j2
     IF ( rs_idx == ring_idx ) THEN
        SELECT CASE ( l1b_channel )
        CASE ( 'UV1')
           READ (UNIT=funit, FMT=*, IOSTAT=ios) x(i-j1+1), y(i-j1+1)
        CASE ( 'UV2')
           READ (UNIT=funit, FMT=*, IOSTAT=ios) x(i-j1+1), ddum, y(i-j1+1)
        CASE ( 'VIS')
           READ (UNIT=funit, FMT=*, IOSTAT=ios) x(i-j1+1), ddum, ddum, y(i-j1+1)
        END SELECT
     ELSE
        READ (UNIT=funit, FMT=*, IOSTAT=ios) x(i-j1+1), y(i-j1+1)
     END IF
     IF ( i /= nspec .AND. ios /= 0 ) THEN
        CALL error_check ( &
             ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
             modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
        IF (  errstat /= pge_errstat_ok ) RETURN
     END IF
  END DO

  ! -----------------------------------------------
  ! Close fitting control file, report SUCCESS read
  ! -----------------------------------------------
  locerrstat = PGS_IO_GEN_CLOSEF ( funit )
  locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
  CALL error_check ( &
       locerrstat, pgs_smf_mask_lev_s, pge_errstat_warning, OMSAO_W_CLOSE_REFSPEC_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )

  ! ------------------------------------------------------------
  ! Reassign number of spectral points and first/last wavelength
  ! ------------------------------------------------------------
  nspec = j2 - j1 + 1

  ! ---------------------------------------------------
  ! Reorder spectrum so that wavelengths are increasing
  ! ---------------------------------------------------
  IF ( x(nspec) < x(1) ) THEN
     irev = (/ (i, i = nspec, 1, -1) /)
     xtmp(1:nspec) = x   (1:nspec)        ;  ytmp(1:nspec) = y   (1:nspec)
     x   (1:nspec) = xtmp(irev(1:nspec))  ;  y   (1:nspec) = ytmp(irev(1:nspec))
  END IF

  ! ------------------------
  ! Assign output quantities
  ! ------------------------
  refspec_orig%nPoints              = nspec
  refspec_orig%RefSpecWavs(1:nspec) = x(1:nspec)
  refspec_orig%RefSpecData(1:nspec) = y(1:nspec)
  refspec_orig%FirstLastWav         = (/ x(1), x(nspec) /)
  refspec_orig%NormFactor           = specnorm
  refspec_orig%Temperature          = rs_temp
  refspec_orig%Title                = TRIM(ADJUSTL(rs_title))
  refspec_orig%Units                = TRIM(ADJUSTL(rs_units))
  !CHARACTER (LEN=maxchlen)                      :: FittingIdxName

  RETURN
END SUBROUTINE read_one_refspec


SUBROUTINE read_commonmode_spec ( &
     omi_lun, specname, winwav_min, winwav_max, common_orig, errstat )

  USE OMSAO_precision_module,   ONLY: r8
  USE OMSAO_indices_module,     ONLY: wvl_idx, spc_idx, ring_idx
  USE OMSAO_parameters_module,  ONLY: maxchlen, max_spec_pts, lm_start_of_table
  USE OMSAO_variables_module,   ONLY: verb_thresh_lev, CommonModeSpectrum
  USE OMSAO_omidata_module,     ONLY: nxtrack_max
  USE OMSAO_errstat_module
  
  IMPLICIT NONE

  ! ----------------
  ! Input Parameters
  ! ----------------
  INTEGER   (KIND=i4), INTENT (IN) :: omi_lun
  CHARACTER (LEN=*),   INTENT (IN) :: specname
  REAL      (KIND=r8), INTENT (IN) :: winwav_min, winwav_max

  ! -----------------
  ! Output Parameters
  ! -----------------
  INTEGER (KIND=i4),         INTENT (INOUT) :: errstat
  TYPE (CommonModeSpectrum), INTENT (OUT)   :: common_orig
  
  ! ----------------
  ! Local Variables
  ! ----------------
  INTEGER (KIND=i4) :: i, ixt, ios, funit, file_read_stat, j1, j2, nskip, nspec
  INTEGER (KIND=i4), DIMENSION (max_spec_pts)             :: irev
  REAL    (KIND=r8), DIMENSION (max_spec_pts)             :: xtmp, ytmp
  REAL    (KIND=r8), DIMENSION (nxtrack_max,max_spec_pts) :: x, y
  CHARACTER (LEN=maxchlen)                                :: lastline, rs_title, rs_units
  REAL    (KIND=r8)                                       :: xdum, rs_temp, ddum, specnorm

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=20), PARAMETER :: modulename = 'read_commonmode_spec'

  ! ------------------------
  ! Error handling variables
  ! ------------------------
  INTEGER (KIND=i4) :: version, locerrstat

  ! ---------------------------------
  ! External OMI and Toolkit routines
  ! ---------------------------------
  INTEGER (KIND=i4), EXTERNAL :: &
       pgs_smf_teststatuslevel, pgs_io_gen_openf, pgs_io_gen_closef


  locerrstat = pge_errstat_ok

  nspec = 0; rs_temp = 0.0_r8
  x = 0.0_r8 ; y = 0.0_r8

  ! -----------------------
  ! Open reference spectrum
  ! -----------------------
  version = 1
  locerrstat = PGS_IO_GEN_OPENF ( omi_lun, PGSd_IO_Gen_RSeqFrm, 0, funit, version )
  locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
  CALL error_check ( &
       locerrstat, pgs_smf_mask_lev_s, pge_errstat_error, OMSAO_E_OPEN_REFSPEC_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
  IF (  errstat /= pge_errstat_ok ) RETURN

  ! --------------------------------------
  ! Skip comments header to start of table
  ! --------------------------------------
  nskip = 0; file_read_stat = 0
  CALL skip_headerlines ( funit, '#', 'read', nskip, lastline, errstat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
  IF (  errstat /= pge_errstat_ok ) RETURN
  
  ! -----------------------------------------------------------------------
  ! Read dimension, start&end wavelength, and norm of spectrum; note that 
  ! the title of the spectrum has already been read as the last of the 
  ! "skipped" lines.
  ! -----------------------------------------------------------------------
  rs_title = TRIM(ADJUSTL(lastline))
  READ (UNIT=funit, FMT='(A)', IOSTAT=file_read_stat) rs_units
  READ (UNIT=funit, FMT=*, IOSTAT=file_read_stat) rs_temp, specnorm, ddum, ddum, nspec
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
  IF (  errstat /= pge_errstat_ok ) RETURN

  ! ---------------------------------------------------------------------
  ! Find first and last index to read, based on WINWAV_MIN and WINWAV_MAX
  ! ---------------------------------------------------------------------
  j1 = 1  ;  j2 = 0 ; xdum = 0.0_r8
  getidx: DO i = 1, nspec
     READ (UNIT=funit, FMT=*, IOSTAT=ios) xdum
     IF ( i /= nspec .AND. ios /= 0 ) THEN
        CALL error_check ( &
             ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
             modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
        IF (  errstat /= pge_errstat_ok ) RETURN
     END IF
     IF ( xdum <= winwav_min )               j1 = i
     IF ( xdum >= winwav_max .AND. j2 == 0 ) j2 = i
     IF ( j1 > 0 .AND. j2 > 0 ) EXIT getidx
  END DO getidx
  IF ( j2 == 0 ) j2 = nspec

  ! -------------------------------------------
  ! Check for maximum number of spectral points
  ! -------------------------------------------
  IF ( j2-j1+1 > max_spec_pts ) THEN
     CALL error_check ( &
          0, 0, pge_errstat_error, OMSAO_E_REFSPEC_MAXPTS, &
          modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
     IF (  errstat /= pge_errstat_ok ) RETURN
  END IF
  ! --------------------------------------------------------------------------
  ! Rewind file, skip to start of table line, and re-read first entry as dummy
  ! --------------------------------------------------------------------------
  REWIND ( funit )

  ! --------------------------------------
  ! Skip comments header to start of table
  ! --------------------------------------
  CALL skip_headerlines ( funit, '#', 'skip', nskip, lastline, errstat )

  ! -----------------------------------------------
  ! Read dimension, start&end, and norm of spectrum
  ! -----------------------------------------------
  READ (UNIT=funit, FMT='(A)', IOSTAT=file_read_stat) rs_title
  READ (UNIT=funit, FMT='(A)', IOSTAT=file_read_stat) rs_units
  READ (UNIT=funit, FMT=*,     IOSTAT=file_read_stat) rs_temp, specnorm, ddum, ddum, nspec
  CALL error_check ( &
       ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
  IF (  errstat /= pge_errstat_ok ) RETURN

  ! -------------------------------------------------
  ! Skip first J1-1 lines, then read J1 to J2 entries
  ! -------------------------------------------------
  DO i = 1, j1-1
     READ (UNIT=funit, FMT=*, IOSTAT=ios) xdum
     IF ( i /= nspec .AND. ios /= 0 ) THEN
        CALL error_check ( &
             ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
             modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
        IF (  errstat /= pge_errstat_ok ) RETURN
     END IF
  END DO
  DO i = j1, j2
     READ (UNIT=funit, FMT=*, IOSTAT=ios) (x(ixt,i-j1+1), y(ixt,i-j1+1), ixt=1, nxtrack_max)
     IF ( i /= nspec .AND. ios /= 0 ) THEN
        CALL error_check ( &
             ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSPEC_FILE, &
             modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )
        IF (  errstat /= pge_errstat_ok ) RETURN
     END IF
  END DO

  ! -----------------------------------------------
  ! Close fitting control file, report SUCCESS read
  ! -----------------------------------------------
  locerrstat = PGS_IO_GEN_CLOSEF ( funit )
  locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
  CALL error_check ( &
       locerrstat, pgs_smf_mask_lev_s, pge_errstat_warning, OMSAO_W_CLOSE_REFSPEC_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(specname)), vb_lev_default, errstat )

  ! ------------------------------------------------------------
  ! Reassign number of spectral points and first/last wavelength
  ! ------------------------------------------------------------
  nspec = j2 - j1 + 1

  ! ---------------------------------------------------
  ! Reorder spectrum so that wavelengths are increasing
  ! ---------------------------------------------------
  IF ( x(1,nspec) < x(1,1) ) THEN
     irev = (/ (i, i = nspec, 1, -1) /)
     DO i = 1, nxtrack_max
        xtmp(  1:nspec) = x   (i,1:nspec)
        ytmp(  1:nspec) = y   (i,1:nspec)
        x   (i,1:nspec) = xtmp(irev(1:nspec))
        y   (i,1:nspec) = ytmp(irev(1:nspec))
     END DO
  END IF

  ! ------------------------
  ! Assign output quantities
  ! ------------------------
  common_orig%nPoints              = nspec
  common_orig%FirstLastWav         = &
       (/ MINVAL(x(1:nxtrack_max,1:nspec)), MAXVAL(x(1:nxtrack_max,1:nspec)) /)
  common_orig%NormFactor           = specnorm
  common_orig%Temperature          = rs_temp
  common_orig%Title                = TRIM(ADJUSTL(rs_title))
  common_orig%Units                = TRIM(ADJUSTL(rs_units))
  DO i = 1, nxtrack_max
     common_orig%RefSpecWavs(i,1:nspec) = x(i,1:nspec)
     common_orig%RefSpecData(i,1:nspec) = y(i,1:nspec)
  END DO

  RETURN
END SUBROUTINE read_commonmode_spec


SUBROUTINE skip_headerlines ( funit, hstr, read_or_skip, nhead, lastline, errstat )

  USE OMSAO_precision_module, ONLY: i4
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER   (KIND=i4), INTENT (IN) :: funit
  CHARACTER (LEN=1),   INTENT (IN) :: hstr
  CHARACTER (LEN=4),   INTENT (IN) :: read_or_skip

  ! ----------------------
  ! Input/Output variables
  ! ----------------------
  INTEGER   (KIND=i4),      INTENT (INOUT) :: nhead, errstat
  CHARACTER (LEN=maxchlen), INTENT (OUT)   :: lastline

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: i, iskip, file_read_stat

  SELECT CASE ( read_or_skip )
  CASE ( 'read' )
     iskip = 0 ; file_read_stat = 0
     skiplines: DO WHILE ( file_read_stat >= 0 )
        lastline = ""
        READ (UNIT=funit, FMT='(A)', IOSTAT=file_read_stat) lastline
        lastline = TRIM(ADJUSTL(lastline))
        IF (  file_read_stat /= 0 .OR. &
             (LEN_TRIM(ADJUSTL(lastline)) > 0 .AND. lastline(1:1) /= hstr) ) THEN
           EXIT skiplines
        ELSE
           iskip = iskip + 1
        END IF
     END DO skiplines
     nhead = iskip
     errstat = pge_errstat_ok
     IF ( file_read_stat /= 0 ) errstat = pge_errstat_error
  CASE ( 'skip' )
     DO i = 1, nhead
        READ (UNIT=funit, FMT='(A)') lastline
     END DO
  END SELECT

  RETURN
END SUBROUTINE skip_headerlines
