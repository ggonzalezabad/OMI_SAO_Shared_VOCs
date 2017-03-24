SUBROUTINE omi_pge_fitting_process ( pge_idx, n_max_rspec, pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_errstat_module,      ONLY: pge_errstat_ok, pge_errstat_error, pge_errstat_fatal
  USE OMSAO_he5_module,          ONLY: NrofScanLines, NrofCrossTrackPixels
  USE OMSAO_solcomp_module
  USE OMSAO_variables_module,    ONLY: l1b_rad_filename
  USE OMSAO_omidata_module,      ONLY: omi_radiance_swathname
  USE OMSAO_radiance_ref_module, ONLY: yn_radiance_reference, l1b_radref_filename

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, n_max_rspec

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: nTimesRad,   nXtrackRad,   nWvlCCD
  INTEGER (KIND=i4) :: nTimesRadRR, nXtrackRadRR, nWvlCCDrr
  INTEGER (KIND=i4) :: errstat

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_pge_fitting_process'


  pge_error_status = pge_errstat_ok

  ! -------------------------------------------------------------------------------------
  ! Set the swath name of various ESDTs
  ! -------------------------------------------------------------------------------------
  CALL omi_set_fitting_parameters ( pge_idx, errstat )
  ! -------------------------------------------------------------------------------------
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error ) GO TO 666

  ! -----------------------------------------------------------------------------------
  ! Get dimensions the L1B radiance granule
  ! -----------------------------------------------------------------------------------
  errstat = pge_errstat_ok
  CALL omi_read_radiance_paras ( &
       l1b_rad_filename, nTimesRad, nXtrackRad, nWvlCCD, &
       omi_radiance_swathname, errstat )
  NrofScanLines        = nTimesRad
  NrofCrossTrackPixels = nXtrackRad

  ! ---------------------------------------------------------------
  ! Dimensions for Radiance Reference granule
  ! ---------------------------------------------------------------
  IF ( .NOT. yn_radiance_reference ) THEN
     l1b_radref_filename = l1b_rad_filename
     nTimesRadRR = nTimesRad ; nXtrackRadRR = nXtrackRad ; nWvlCCDrr = nWvlCCD
  ELSE
     CALL omi_read_radiance_paras (                                  &
          l1b_radref_filename, nTimesRadRR, nXtrackRadRR, nWvlCCDrr, &
          omi_radiance_swathname, errstat                              )
     ! ----------------------------------------------------------------
     ! Number of cross-track positions must be the same; fold otherwise
     ! ----------------------------------------------------------------
     IF ( nXtrackRad /= nXtrackRadRR ) THEN
        CALL error_check ( 0, 1, pge_errstat_fatal, OMSAO_F_XTRMISRAD, &
             modulename, vb_lev_default, errstat )
        GO TO 666
     END IF
  END IF

  ! -----------------------------------------------------------------------------------
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  ! -----------------------------------------------------------------
  CALL omi_fitting (                                  &
       pge_idx,                                       &
       nTimesRad,   nXtrackRad, nWvlCCD, n_max_rspec, &
       nTimesRadRR, nWvlCCDrr,  pge_error_status       )

  IF ( pge_error_status >= pge_errstat_fatal ) GO TO 666
  
  ! -------------------------------------------------------------
  ! Here is the place to jump to in case some error has occurred.
  ! Naturally, we also reach here when everything executed as it
  ! was supposed to, but that doesn't matter, since we are not
  ! taking any particular action at this point.
  ! -------------------------------------------------------------
666 CONTINUE

  ! -------------------------------------------------
  ! Deallocation of some potentially allocated memory
  ! -------------------------------------------------
  CALL soco_pars_deallocate (errstat)

  IF ( pge_error_status >= pge_errstat_fatal ) RETURN

  RETURN
END SUBROUTINE omi_pge_fitting_process


SUBROUTINE omi_fitting (                                  &
       pge_idx,                                           &
       nTimesRad,   nXtrackRad,   nWvlCCD,   n_max_rspec, &
       nTimesRadRR, nWvlCCDrr, pge_error_status             )

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: sao_molecule_names, pge_hcho_idx, n_max_fitpars
  USE OMSAO_parameters_module, ONLY: i2_missval, r8_missval
  USE OMSAO_variables_module, ONLY: l1b_rad_filename, verb_thresh_lev, hw1e, &
       e_asym, phase, l2_filename, radwavcal_freq, pixnum_lim, radfit_latrange, &
       database, n_fitvar_rad, refspecs_original, fitvar_rad_init, &
       fitvar_rad_saved, yn_solar_comp, yn_diagnostic_run, n_fitres_loop, &
       fitres_range, yn_common_iter, common_latlines, common_latrange, &
       radiance_wavcal_lnums, yn_solmonthave
  USE OMSAO_omidata_module
  USE OMSAO_he5_module,       ONLY:  pge_swath_name
  USE OMSAO_he5_datafields_module
  USE OMSAO_solar_wavcal_module
  USE OMSAO_radiance_ref_module
  USE OMSAO_destriping_module, ONLY: ctr_maxcol
  USE OMSAO_prefitcol_module
  USE OMSAO_errstat_module
  USE OMSAO_solmonthave_module
  USE OMSAO_wfamf_module, ONLY: omi_read_climatology, CmETA

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, nTimesRad, nXtrackRad, nWvlCCD, n_max_rspec
  INTEGER (KIND=i4), INTENT (IN) :: nTimesRadRR, nWvlCCDrr
  
  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: pge_error_status

  ! -------------------------
  ! Local variables (for now)
  ! -------------------------
  INTEGER   (KIND=i4) ::                                                 &
       iline, first_line, last_line, errstat, first_wc_pix, last_wc_pix, &
       first_line_save, last_line_save, first_pix, last_pix
  REAl (KIND=r8), DIMENSION (1) :: zerovec = 0.0_r8

  ! ----------------------------------------------------------------------
  ! Swath dimensions and variables that aren't passed from calling routine
  ! ----------------------------------------------------------------------
  INTEGER   (KIND=i4)      :: nTimesIrr, nXtrackIrr
  CHARACTER (LEN=maxchlen) :: molname

  ! ----------------------------------------------------------
  ! Variables and parameters associated with Spatial Zoom data
  ! and Common Mode spectrum
  ! ----------------------------------------------------------
  INTEGER (KIND=i1), DIMENSION (0:nTimesRad-1)   :: omi_binfac
  INTEGER (KIND=i4), DIMENSION (0:nTimesRad-1,2) :: omi_xtrpix_range
  LOGICAL,           DIMENSION (0:nTimesRad-1)   :: &
       omi_yn_szoom, yn_common_range, yn_radfit_range

  INTEGER (KIND=i1), DIMENSION (0:nTimesRadRR-1)   :: omi_binfac_rr
  INTEGER (KIND=i4), DIMENSION (0:nTimesRadRR-1,2) :: omi_xtrpix_range_rr
  LOGICAL,           DIMENSION (0:nTimesRadRR-1)   :: omi_yn_szoom_rr

  ! ----------------------------------------------------------
  ! OMI L1b latitudes
  ! ----------------------------------------------------------
  REAL (KIND=r4), DIMENSION (1:nXtrackRad, 0:nTimesRad-1) :: l1b_latitudes

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=11), PARAMETER :: modulename = 'omi_fitting'

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER (KIND=i4), EXTERNAL :: &
       he5_init_swath, he5_define_fields, he5_close_output_file, &
       he5_set_field_attributes, he5_write_global_attributes,    &
       he5_write_swath_attributes, he5_open_readwrite

  ! ---------------------------------------------------------------
  ! Some initializations that will save us headaches in cases where 
  ! a proper set-up of those variables failes or is bypassed.
  ! ---------------------------------------------------------------
  first_pix = 1 ; last_pix = 1

  errstat = pge_errstat_ok

  ! --------------------------------
  ! Name of the main output molecule
  ! --------------------------------
  molname = sao_molecule_names(pge_idx)

  ! -------------------------------------------------------------------
  ! Range of cross-track pixels to fit. This is based on the selection
  ! in the fitting control file and whether the granule being processed
  ! is in global or spatial zoom mode, or even a mixture thereof.
  !
  ! NOTE that we set OMI_XTRPIX_RANGE for all swath lines because the
  ! choice of swath lines to process may not contain the radiance 
  ! reference/calibration line.
  ! -------------------------------------------------------------------
  CALL omi_read_binning_factor ( &
       TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
       nTimesRad, omi_binfac(0:nTimesRad-1), omi_yn_szoom(0:nTimesRad-1),   &
       errstat )
  CALL omi_set_xtrpix_range ( &
     nTimesRad, nXtrackRad, pixnum_lim(3:4),                         &
     omi_binfac(0:nTimesRad-1), omi_xtrpix_range(0:nTimesRad-1,1:2), &
     first_wc_pix, last_wc_pix, errstat )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  
  ! --------------------------------------------------------------------
  ! If the radiance reference is obtained from the same L1b file, we can
  ! simply copy the variables we have just read to the corresponding 
  ! "rr" ones (in this case, the dimensions are the same). Otherwise we
  ! have to read them from the radiance reference granule.
  ! --------------------------------------------------------------------
  IF ( TRIM(ADJUSTL(l1b_radref_filename)) /= TRIM(ADJUSTL(l1b_rad_filename)) ) THEN
     CALL omi_read_binning_factor ( &
          TRIM(ADJUSTL(l1b_radref_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
          nTimesRadRR, omi_binfac_rr(0:nTimesRadRR-1), omi_yn_szoom_rr(0:nTimesRadRR-1), &
          errstat )
     CALL omi_set_xtrpix_range ( &
          nTimesRadRR, nXtrackRad, pixnum_lim(3:4),                                 &
          omi_binfac_rr(0:nTimesRadRR-1), omi_xtrpix_range_rr(0:nTimesRadRR-1,1:2), &
          first_wc_pix, last_wc_pix, errstat )
     pge_error_status = MAX ( pge_error_status, errstat )
     IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  ELSE
     omi_binfac_rr      (0:nTimesRad-1)       = omi_binfac      (0:nTimesRad-1)
     omi_yn_szoom_rr    (0:nTimesRadRR-1)     = omi_yn_szoom    (0:nTimesRad-1)
     omi_xtrpix_range_rr(0:nTimesRadRR-1,1:2) = omi_xtrpix_range(0:nTimesRad-1,1:2)
  END IF

  ! --------------------------------------------------------------------
  ! Solar Irradiance Processing: If we don't do a solar composite, we can
  ! use a solar monthly average, if not we have to read the irradiance 
  ! data. 
  ! Otherwise we need to compute them from the solar composite 
  ! parameterization on a equidistant grid.
  ! -------------------------------------------------------------------
  omi_solcal_itnum = i2_missval ; omi_solcal_xflag = i2_missval
  ! --------------------------------------------------------------------------
  ! Check than only one or non of yn_solar_comp are yn_solmonthva are set True
  ! --------------------------------------------------------------------------
  IF ( yn_solar_comp .AND. yn_solmonthave ) THEN
     CALL error_check ( 1, 0, pge_errstat_fatal, OMSAO_F_SOLCOM_VS_SOLAVE, &
          modulename, vb_lev_gt1mb, errstat )
     GO TO 666
  END IF

  IF ( yn_solar_comp ) THEN
     ! -----------------------------------
     ! Compute composite solar irradiances
     ! -----------------------------------
     CALL omi_create_solcomp_irradiance ( nXtrackRad )
  ELSE IF (yn_solmonthave) THEN
     ! -----------------------------------
     ! Read solar monthly mean irradiance
     ! -----------------------------------
     CALL omi_read_monthly_average_irradiance (nTimesIrr, nXtrackIrr, errstat)
     pge_error_status = MAX ( pge_error_status, errstat )
     IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  ELSE
     ! --------------------
     ! Read OMI irradiances
     ! --------------------
     CALL omi_read_irradiance_data ( nTimesIrr, nXtrackIrr, errstat ) 
     pge_error_status = MAX ( pge_error_status, errstat )
     IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  END IF

  ! ---------------------------------------------------------------
  ! Solar wavelength calibration, done even when we use a composite
  ! solar spectrum to avoid un-initialized variables. However, no
  ! actual fitting is performed in the latter case.
  ! ---------------------------------------------------------------
  CALL xtrack_solar_calibration_loop ( first_wc_pix, last_wc_pix, errstat )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

  ! ---------------------------------------------------------------
  ! No matter what, we need a swath line for radiance wavelength 
  ! calibration. This may be a single line or it may be the average
  ! over a block of lines. However, if we are not using a radiance
  ! reference, then we are still doing a radiance calibration and
  ! need to make sure that we are using a radiance from the current
  ! granule.
  ! ---------------------------------------------------------------
  CALL omi_get_radiance_reference (                             &
       l1b_radref_filename, nTimesRadRR, nXtrackRad, nWvlCCDrr, &
       radiance_wavcal_lnums, errstat                             )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

  ! ---------------------------------------------------------------
  ! The Climatology is going to be read here and kept in memory. If
  ! this has a bad impact in the efficiency of the application then
  ! I will find a different way. We are doing this to be able to in
  ! itialize the output he5 with the correct number of levels for
  ! the Scattering weights and Gas_profiles output. We are going to
  ! use the number of levels in the climatology as the number of le
  ! vels of the reported scattering weights.
  ! ---------------------------------------------------------------
  CALL omi_read_climatology ( errstat )

  ! ---------------------------
  ! Initialize omi_ozone_amount
  ! ---------------------------
  ALLOCATE (omi_ozone_amount(1:nXtrackRad,0:nTimesRad-1))
  omi_ozone_amount = r8_missval
  
  ! ----------------------------------------
  ! Initialization of HE5 output data fields
  ! ----------------------------------------
  errstat = HE5_Init_Swath ( l2_filename, pge_swath_name, nTimesRad, nXtrackRad, CmETA )
  CALL he5_initialize_datafields ( )
  errstat = HE5_Define_Fields ( pge_idx, pge_swath_name, nTimesRad, nXtrackRad, CmETA )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

  ! -----------------------------------------------------------------------------------
  ! If we are NOT using a radiance reference, then we need to read the 
  ! swath line for radiance wavelength calibration. In this case, the
  ! value of  RADIANCE_WAVCAL_LNUMS is that of the first line of the radiance reference
  ! selected
  ! -----------------------------------------------------------------------------------
  IF ( .NOT. yn_radiance_reference ) THEN
     ntimes_loop = 1                                           ! The number of scan lines to read
     iline = radiance_wavcal_lnums(1)
     CALL omi_read_radiance_lines (                         &  ! Get NTIMES_LOOP radiance lines
          l1b_rad_filename, iline, nXtrackRad, ntimes_loop, &
          nWvlCCD, errstat )
     pge_error_status = MAX ( pge_error_status, errstat )
     IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  END IF

  ! -----------------------------------------------------
  ! Across-track loop for radiance wavelength calibration
  ! -----------------------------------------------------
  omi_radcal_itnum = i2_missval ; omi_radcal_xflag = i2_missval
  CALL xtrack_radiance_wvl_calibration (                          &
       yn_radiance_reference, yn_solar_comp,                      &
       first_wc_pix, last_wc_pix, n_max_rspec, n_comm_wvl, errstat )
  pge_error_status = MAX ( pge_error_status, errstat )
  IF ( pge_error_status >= pge_errstat_error )  GO TO 666

  ! --------------------------------------------------------------
  ! Terminate on not having any cross-track pixels left to process
  ! --------------------------------------------------------------
  IF ( ALL ( omi_cross_track_skippix ) ) THEN
     CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_NOPIXEL, &
          modulename, vb_lev_default, errstat )
     GO TO 666
  END IF
  
  ! ----------------------------------
  ! CCM - Add lqH2O prefit
  ! ----------------------------------
  SELECT CASE(pge_idx )
  
  CASE( pge_hcho_idx )
  
	  ! -------------------------------------------------------------------
	  ! First access to pre-fitted O3 and BrO columns. At this time we have
	  ! already set up some of the fitting arrays, so any error would lead
	  ! to some major headaches to undo things. Hence we return if access
	  ! fails.
	  ! -------------------------------------------------------------------
	  IF ( ( .NOT. yn_radiance_reference ) .AND. &
	       ( pge_idx == pge_hcho_idx ) .AND. &
	       ANY((/yn_o3_prefit(1),yn_bro_prefit(1)/)) ) THEN
	     CALL init_prefit_files ( pge_idx, nTimesRad, nXtrackRad, errstat )
	     IF ( errstat >= pge_errstat_error ) RETURN
	  END IF
	
	CASE( pge_gly_idx )
		
		! -------------------------------------------------------------------
    ! First access to pre-fitted Liquid Water columns. At this time we
    ! have already set up some of the fitting arrays, so any error would
    ! lead to some major headaches to undo things. Hence we return if
    ! access fails.
    ! -------------------------------------------------------------------
    
    IF ( (.NOT. yn_radiance_reference) .AND. yn_lqh2o_prefit(1) ) THEN
       CALL init_prefit_files ( pge_idx, nTimesRad, nXtrackRad, errstat )
       IF ( errstat >= pge_errstat_error ) RETURN
    END IF
		
	CASE DEFAULT	
		! Nothing here yet
	END SELECT

  ! ---------------------------------------------------------------------
  ! If we are using a radiance reference AND want to remove the target
  ! gas from it (important for BrO, for example), we have to run through
  ! all spectra that go into the reference, compute the average column,
  ! and then remove that from the averaged radiance reference spectrum
  ! (owing to the fact that the average of hundreds of OMI spectra still
  !  doesn't produce a decently fitted column).
  ! ---------------------------------------------------------------------
  IF ( yn_radiance_reference .AND. yn_remove_target ) THEN


     iline = SUM ( radiance_reference_lnums(1:2) ) / 2
     IF ( iline < 0 .OR. iline > nTimesRadRR ) iline = nTimesRadRR / 2

     first_pix = omi_xtrpix_range_rr(iline,1)
     last_pix  = omi_xtrpix_range_rr(iline,2)

     CALL xtrack_radiance_reference_loop (                                   &
          yn_radiance_reference, yn_remove_target,                           &
          nXtrackRad, nWvlCCDrr, first_pix, last_pix, pge_idx, pge_error_status )

     ! -------------------------------------------------------------
     ! Write the output from solar/earthshine wavelength calibration
     ! and radiance reference to file. The latter results will be
     ! overwritten in the call to XTRACK_RADIANCE_REFERENCE_LOOP
     ! below, hence we need to write them out here.
     ! -------------------------------------------------------------
     CALL he5_write_wavcal_output ( nXtrackRad, first_pix, last_pix, errstat )

  END IF

  ! -----------------------------------------------------------------
  ! Before we go any further we need to read the L1b latitude values,
  ! since we base our screening of which swath lines to process on
  ! those values. Both common mode, if used, and the radiance fit
  ! uses the same arrays, so we read this only ones.
  !
  ! We could shave off some fractional minute from the run time by
  ! not reading the latitudes in cases where no radiance reference
  ! is used, i.e., where both radiance granule and radiance reference
  ! granule are the same. The down-side is an increase in virtual
  ! memory program uses, plus some more logic to find out whether to
  ! read the latitudes or not. For now we are going with a second
  ! read, particularly since the current algorithm settings would
  ! require it anyway.
  ! -----------------------------------------------------------------
  CALL read_latitude ( &
       TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
       nTimesRad, nXtrackRad, l1b_latitudes(1:nXtrackRad,0:nTimesRad-1) )

  ! -----------------------------------------------------------------
  ! Now we enter the on-line computation of the common mode spectrum.
  ! -----------------------------------------------------------------
  IF ( yn_common_iter ) THEN
     
     ! ----------------------------------------------------------
     ! Set the logical YN array that determines which swath lines
     ! will be used in the common mode
     ! ----------------------------------------------------------
     yn_common_range(0:nTimesRad-1) = .FALSE.
     CALL find_swathline_range ( &
          TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
          nTimesRad, nXtrackRad, l1b_latitudes(1:nXtrackRad,0:nTimesRad-1),       &
          common_latrange(1:2), yn_common_range(0:nTimesRad-1), errstat             )

     ! -------------------------------------------------------------
     ! First and last swath line number will be overwritten. Hence
     ! we save them for further reference.
     ! -------------------------------------------------------------
     !first_line_save = first_line
     !last_line_save  = last_line

     ! ------------------------------------------------------------------
     ! There is a time saver catch built into the assignment below, 
     ! which we probably want to rethink. If we are only doing a few
     ! lines but the common mode extends over a wide range of a
     ! latitudes, then we would be processing a lot of lines to derive
     ! the common mode. This is currently being excluded, but the better
     ! way may be to remember that we can control everything through
     ! the fitting control file.
     ! ------------------------------------------------------------------
     !IF ( MINVAL(common_latlines(1:2)) >= 0 .AND. &
     !     MAXVAL(common_latlines(1:2)) <= ntimesrad ) THEN
     !   first_line = common_latlines(1) !MAX(first_line, common_latlines(1))
     !   last_line  = common_latlines(2) !MIN(last_line,  common_latlines(2))
     !END IF

     ! ------------------------------------------
     ! Interface to the loop over all swath lines
     ! ------------------------------------------
     CALL omi_pge_swathline_loop (                              &
          pge_idx, nTimesRad, nxtrackRad, nWvlCCD, n_max_rspec, &
          yn_common_range(0:nTimesRad-1),                       &
          omi_xtrpix_range(0:nTimesRad-1,1:2),                  &
          yn_radiance_reference, .FALSE., -1,                   &
          yn_common_iter, pge_error_status )
 
     ! -----------------------------------------------------------------
     ! Reset first and last swath line number to non-common mode values.
     ! -----------------------------------------------------------------
     !first_line = first_line_save
     !last_line  = last_line_save

     ! ---------------------------------------------------
     ! Set the index value of the Common Mode spectrum and
     ! assign values to the fitting parameter arrays
     ! ---------------------------------------------------
     CALL compute_common_mode ( .FALSE., nXtrackRad, 1, zerovec, zerovec, .TRUE. )

     ! -------------------------------------------
     ! Write the just computed common mode to file
     ! -------------------------------------------
     IF ( yn_diagnostic_run ) &
          CALL he5_write_common_mode ( nXtrackRad, n_comm_wvl, pge_error_status )

  END IF
    
  ! ----------------------------------------------------------
  ! Now into the proper fitting, with or without common mode.
  ! ----------------------------------------------------------

  ! -------------------------------------------------------------------
  ! Radiance Reference Fit: Only if we have not selected to remove the 
  ! target gas from the radiance reference (in which case we have done
  ! this fit already), or N_COMM_ITER==0 (which means that we can refit
  ! with a common mode spectrum)
  ! -------------------------------------------------------------------
  IF ( radiance_wavcal_lnums(1) >= 0 ) THEN
  
     ! --------------------------------
     ! The number of scan lines to read
     ! --------------------------------
     ntimes_loop = 1
     iline = radiance_wavcal_lnums(1)

     ! -------------------------------------------------
     ! Only use prefit if not using a Radiance Reference
     ! CCM modify to include lqH2O
     ! -------------------------------------------------
     IF ( ( .NOT. yn_radiance_reference )                             .AND. &
          ( (pge_idx == pge_hcho_idx) .OR. (pge_idx == pge_gly_idx) ) .AND. &
          ANY((/yn_o3_prefit(1),yn_bro_prefit(1),yn_lqh2o_prefit(1)/)) ) THEN
        CALL read_prefit_columns ( pge_idx, nXtrackRad, ntimes_loop, iline, errstat )
        pge_error_status = MAX ( pge_error_status, errstat )
        IF ( errstat >= pge_errstat_error ) GO TO 666
     ELSE
        o3_prefit_col    = 0.0_r8 ; o3_prefit_dcol    = 0.0_r8
        bro_prefit_col   = 0.0_r8 ; bro_prefit_dcol   = 0.0_r8
        lqh2o_prefit_col = 0.0_r8 ; lqh2o_prefit_dcol = 0.0_r8
     END IF

     ! ------------------------------
     ! Get NTIMES_LOOP radiance lines
     ! ------------------------------
     CALL omi_read_radiance_lines ( &
          l1b_rad_filename, iline, nXtrackRad, ntimes_loop, nWvlCCD, errstat )
     pge_error_status = MAX ( pge_error_status, errstat )
     IF ( pge_error_status >= pge_errstat_error )  GO TO 666

     ! -----------------------------------------------
     ! Radiance Reference Fit (or WavCal Radiance Fit)
     ! -----------------------------------------------
     first_pix = omi_xtrpix_range(iline,1)
     last_pix  = omi_xtrpix_range(iline,2)

     ! -------------------------------------
     ! Initialize saved fitting variables
     ! -------------------------------------
     fitvar_rad_saved(1:n_max_fitpars ) = fitvar_rad_init(1:n_max_fitpars)

     IF (yn_radiance_reference) THEN
        CALL xtrack_radiance_reference_loop ( &
             yn_radiance_reference, .FALSE.,  &
             nXtrackRad, nWvlCCDrr, first_pix, last_pix, pge_idx, errstat )
     ENDIF

  END IF

  ! -------------------------------------------------------------------
  ! Output of fit results for solar and radiance wavelength calibration
  ! but ONLY if we haven't done it already (see above). The Radiance
  ! Reference values under YN_REMOVE_TARGET settings carry valuable
  ! information only BEFORE the target has been removed.
  ! -------------------------------------------------------------------
  IF ( .NOT. (yn_radiance_reference .AND. yn_remove_target) ) THEN
     CALL he5_write_wavcal_output ( nXtrackRad, first_pix, last_pix, errstat )
     pge_error_status = MAX ( pge_error_status, errstat )
     IF ( pge_error_status >= pge_errstat_error )  GO TO 666
  END IF

  ! ----------------------------------------------------------
  ! Set the logical YN array that determines which swath lines
  ! will be processed. Unless we have constrained either swath
  ! line numbers or the latitude range, this can be set to
  ! .TRUE. universically.
  ! ----------------------------------------------------------
  ! First, set the range of swath lines to process
  ! ----------------------------------------------
  first_line = 0  ;  last_line = nTimesRad-1
  IF ( pixnum_lim(1) > 0 ) first_line = MIN(pixnum_lim(1), last_line)
  IF ( pixnum_lim(2) > 0 ) last_line  = MAX( MIN(pixnum_lim(2), last_line), first_line )

  yn_radfit_range = .FALSE.
  IF ( first_line         > 0           .OR. &
       last_line          < nTimesRad-1 .OR. &
       radfit_latrange(1) > -90.0_r4    .OR. &
       radfit_latrange(2) < +90.0_r4           ) THEN

     IF ( radfit_latrange(1) > -90.0_r4    .OR. &
          radfit_latrange(2) < +90.0_r4           ) THEN
        CALL find_swathline_range ( &
             TRIM(ADJUSTL(l1b_rad_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
             nTimesRad, nXtrackRad, l1b_latitudes(1:nXtrackRad,0:nTimesRad-1),       &
             radfit_latrange(1:2), yn_radfit_range(0:nTimesRad-1), errstat             )
     ELSE
        yn_radfit_range = .TRUE.
        IF ( first_line > 0           ) yn_radfit_range(0:first_line-1)          = .FALSE.
        IF ( last_line  < nTimesRad-1 ) yn_radfit_range(last_line+1:nTimesRad-1) = .FALSE.
     END IF
  ELSE
     yn_radfit_range = .TRUE.
  END IF
  
  ! ------------------------------------------
  ! Interface to the loop over all swath lines
  ! ------------------------------------------
  CALL omi_pge_swathline_loop (                                  &
       pge_idx, nTimesRad, nXtrackRad, nWvlCCD, n_max_rspec,     &
       yn_radfit_range(0:nTimesRad-1),                           &
       omi_xtrpix_range(0:nTimesRad-1,1:2),                      &
       yn_radiance_reference, .FALSE., -1,                       &
       .FALSE., pge_error_status )

  ! -------------------------------------------------------------
  ! Here is the place to jump to in case some error has occurred.
  ! Naturally, we also reach here when everything executed as it
  ! was supposed to, but that doesn't matter, since we are not
  ! taking any particular action at this point.
  ! -------------------------------------------------------------
666 CONTINUE
  IF ( pge_error_status >= pge_errstat_fatal ) RETURN

  ! ----------------------------------------------
  ! This subroutine completes the following tasks:
  !    (1) Compute pixel geololcation corners
  !    (2) Compute AMFs
  !    (3) Apply cross-track destriping correction
  ! ---------------------------------------
  CALL omi_pge_postprocessing (                                             &
       l1b_rad_filename, pge_idx, nTimesRad, nXtrackRad,                    &
       yn_radfit_range(0:nTimesRad-1), omi_xtrpix_range(0:nTimesRad-1,1:2), &
       omi_yn_szoom(0:nTimesRad-1), n_max_rspec, errstat                 )

  ! ---------------------------
  ! Deallocate omi_ozone_amount
  ! ---------------------------
  IF (ALLOCATED ( omi_ozone_amount ) ) DEALLOCATE ( omi_ozone_amount )

  ! ---------------------
  ! Write some attributes
  ! ---------------------
  errstat = pge_errstat_ok
  errstat = he5_write_global_attributes( )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"HE5_WRITE_GLOBAL_ATTRIBUTES", vb_lev_default, pge_error_status )

  errstat = pge_errstat_ok
  errstat = he5_write_swath_attributes ( pge_idx )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"HE5_WRITE_SWATH_ATTRIBUTES", vb_lev_default, pge_error_status )

  errstat = pge_errstat_ok
  errstat = he5_set_field_attributes   ( pge_idx )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"HE5_SET_FIELD_ATTRIBUTES", vb_lev_default, pge_error_status )
  ! -----------------
  ! Close output file
  ! -----------------
  errstat = pge_errstat_ok
  errstat = he5_close_output_file ( pge_idx)
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"HE5_CLOSE_OUTPUT_FILE", vb_lev_default, pge_error_status )

  ! ----------------------------------------
  ! Get Metadata from MCF, PCF, and L1B file
  ! ----------------------------------------
  errstat = pge_errstat_ok
  CALL check_metadata_consistency ( errstat )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"CHECK_METADATA_CONSISTENCY.", vb_lev_default, pge_error_status )
  errstat = pge_errstat_ok
  CALL set_l2_metadata            ( errstat )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"SET_L2_METADATA.", vb_lev_default, pge_error_status )

  RETURN
END SUBROUTINE omi_fitting
