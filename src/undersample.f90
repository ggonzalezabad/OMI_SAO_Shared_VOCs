SUBROUTINE undersample ( xtrack_pix, n_sensor_pts, curr_wvl, hw1e, e_asym, phase, errstat )

  !     Convolves input spectrum with Gaussian slit function of specified
  !     HW1e, and samples at a particular input phase to give the OMI
  !     undersampling spectrum. This version calculates both phases of the
  !     undersampling spectrum, phase1 - i.e., underspec (1, i) - being the
  !     more common in OMI spectra.

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: solar_idx, us1_idx, us2_idx, wvl_idx, spc_idx
  USE OMSAO_parameters_module, ONLY: max_spec_pts, maxchlen
  USE OMSAO_variables_module,  ONLY: &
       refspecs_original, database, have_undersampling, yn_use_labslitfunc
  USE OMSAO_slitfunction_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                           INTENT (IN) :: n_sensor_pts, xtrack_pix
  REAL    (KIND=r8),                           INTENT (IN) :: hw1e, e_asym, phase
  REAL    (KIND=r8), DIMENSION (n_sensor_pts), INTENT (IN) :: curr_wvl

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  LOGICAL                                     :: yn_full_range
  REAL (KIND=r8), DIMENSION (2,n_sensor_pts)  :: underspec
  REAL (KIND=r8), DIMENSION (max_spec_pts)    :: &
       locwvl, locspec, specmod, tmpwav, over, under, resample

  INTEGER (KIND=i4) :: npts, locerrstat

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=11), PARAMETER :: modulename = 'undersample'

  locerrstat = pge_errstat_ok

  ! ==================================================
  ! Assign solar reference spectrum to local variables
  ! ==================================================
  npts            = refspecs_original(solar_idx)%nPoints
  locwvl (1:npts) = refspecs_original(solar_idx)%RefSpecWavs(1:npts)
  locspec(1:npts) = refspecs_original(solar_idx)%RefSpecData(1:npts)


  IF ( yn_use_labslitfunc ) THEN
     CALL omi_slitfunc_convolve ( &
          xtrack_pix, npts, locwvl(1:npts), locspec(1:npts), specmod(1:npts), locerrstat )
  ELSE
     CALL asymmetric_gaussian_sf ( &
          npts, hw1e, e_asym, locwvl(1:npts), locspec(1:npts), specmod(1:npts))
  END IF
  CALL error_check ( &
       locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename//f_sep//'Convolution', vb_lev_default, errstat )
  IF ( locerrstat >= pge_errstat_error ) RETURN

  ! Phase1 calculation: Calculate spline derivatives for KPNO data
  !                     Calculate solar spectrum at OMI positions

  CALL interpolation (                                                                &
       modulename//f_sep//'Phase 1a',                                                 &
       npts, locwvl(1:npts), specmod(1:npts), n_sensor_pts, curr_wvl(1:n_sensor_pts), &
       resample(1:n_sensor_pts), 'endpoints', 0.0_r8, yn_full_range, locerrstat )
  CALL error_check (                                                    &
       locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename//f_sep//'Phase 1a', vb_lev_default, errstat )
  IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  IF ( .NOT. yn_full_range ) CALL error_check (           &
       0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
       modulename//f_sep//'Phase 1a', vb_lev_develop, errstat )

  ! Calculate solar spectrum at OMI + phase positions, original and resampled.

  ! ------------------------------------------------------------------------------
  ! The original ("modified K.C.) scheme to compute the UNDERSPEC wavelength array
  ! ------------------------------------------------------------------------------
  ! ( assumes ABS(PHASE) < 1.0 )
  ! ----------------------------
  tmpwav(2:n_sensor_pts) = (1.0_r8-phase)*curr_wvl(1:n_sensor_pts-1) + phase*curr_wvl(2:n_sensor_pts)
  tmpwav(1)              = curr_wvl(1)
  tmpwav(n_sensor_pts)   = curr_wvl(n_sensor_pts)
  IF ( tmpwav(2) <= tmpwav(1) ) tmpwav(2) = (tmpwav(1)+tmpwav(3))/2.0_r8

  CALL interpolation (                                                              &
       modulename//f_sep//'Phase 1b',                                               &
       npts, locwvl(1:npts), specmod(1:npts), n_sensor_pts, tmpwav(1:n_sensor_pts), &
       over(1:n_sensor_pts), 'endpoints', 0.0_r8, yn_full_range, locerrstat )
  CALL error_check (                                                    &
       locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename//f_sep//'Phase 1b', vb_lev_default, errstat )
  IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  !IF ( .NOT. yn_full_range ) CALL error_check (                 &
  !     0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE,       &
  !     modulename//f_sep//'Phase 1b', vb_lev_omidebug, errstat )

  CALL interpolation (                                                                 &
       modulename//f_sep//'Phase 1c',                                                  &
       n_sensor_pts, curr_wvl(1:n_sensor_pts), resample(1:n_sensor_pts), n_sensor_pts, &
       tmpwav(1:n_sensor_pts), under(1:n_sensor_pts), 'endpoints', 0.0_r8,             &
       yn_full_range, locerrstat )
  CALL error_check (                                                    &
       locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename//f_sep//'Phase 1c', vb_lev_default, errstat )
  IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  !IF ( .NOT. yn_full_range ) CALL error_check (                 &
  !     0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE,       &
  !     modulename//f_sep//'Phase 1c', vb_lev_omidebug, errstat )

  underspec(1,1:n_sensor_pts) = over(1:n_sensor_pts) - under(1:n_sensor_pts)
  resample (  1:n_sensor_pts) = over(1:n_sensor_pts)

  ! ------------------------------------------------
  ! Save spectra to final arrays for Undersampling 1
  ! ------------------------------------------------
  refspecs_original(us1_idx)%nPoints                     = n_sensor_pts
  refspecs_original(us1_idx)%NormFactor                  = 1.0E+00_R8
  refspecs_original(us1_idx)%RefSpecWavs(1:n_sensor_pts) = tmpwav(1:n_sensor_pts)
  refspecs_original(us1_idx)%RefSpecData(1:n_sensor_pts) = underspec(1,1:n_sensor_pts)

  database(us1_idx,1:n_sensor_pts) = underspec (1,1:n_sensor_pts)

  ! ---------------------------------------------------------
  ! If we haven't selected Undersampling 2 then we return now
  ! ---------------------------------------------------------
  IF ( .NOT. have_undersampling(us2_idx) ) RETURN

  ! --------------------------------------------------------------------------------------
  ! Phase2 calculation: Calculate solar spectrum at OMI positions, original and resampled.
  ! --------------------------------------------------------------------------------------
  CALL interpolation (                                                                &
       modulename//f_sep//'Phase 2a',                                                 &
       npts, locwvl(1:npts), specmod(1:npts), n_sensor_pts, curr_wvl(1:n_sensor_pts), &
       over(1:n_sensor_pts), 'endpoints', 0.0_r8, yn_full_range, locerrstat )
  CALL error_check ( &
       locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename//f_sep//'Phase 2a', vb_lev_default, errstat )

  IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  !IF ( .NOT. yn_full_range ) CALL error_check ( &
  !     0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
  !     modulename//f_sep//'Phase 2a', vb_lev_omidebug, errstat )

  CALL interpolation ( &
       modulename//f_sep//'Phase 2b',                                                &
       n_sensor_pts, tmpwav(1:n_sensor_pts), resample(1:n_sensor_pts), n_sensor_pts, &
       curr_wvl(1:n_sensor_pts), under(1:n_sensor_pts), 'endpoints', 0.0_r8,         &
       yn_full_range, locerrstat )
  CALL error_check ( &
       locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename//f_sep//'Phase 2b', vb_lev_default, errstat )

  IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  !IF ( .NOT. yn_full_range ) CALL error_check ( &
  !     0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
  !     modulename//f_sep//'Phase 2b', vb_lev_omidebug, errstat )

  ! ------------------------------------
  ! Compute final undersampling spectrum
  ! ------------------------------------
  underspec(2,1:n_sensor_pts) = over(1:n_sensor_pts) - under(1:n_sensor_pts)

  ! ------------------------------------------------
  ! Save spectra to final arrays for Undersampling 2
  ! ------------------------------------------------
  refspecs_original(us2_idx)%nPoints                     = n_sensor_pts
  refspecs_original(us2_idx)%NormFactor                  = 1.0E+00_R8
  refspecs_original(us2_idx)%RefSpecWavs(1:n_sensor_pts) = curr_wvl(1:n_sensor_pts)
  refspecs_original(us2_idx)%RefSpecData(1:n_sensor_pts) = underspec(2,1:n_sensor_pts)

  database(us2_idx,1:n_sensor_pts) = underspec (2,1:n_sensor_pts)


  RETURN
END SUBROUTINE undersample


SUBROUTINE undersample_new ( xtrack_pix, n_sensor_pts, curr_wvl, n_solar_pts, solar_wvl, &
                             hw1e, e_asym, errstat )

  !     Convolves input spectrum with Gaussian slit function of specified
  !     HW1e, and samples at a particular input phase to give the OMI
  !     undersampling spectrum. This version calculates both phases of the
  !     undersampling spectrum, phase1 - i.e., underspec (1, i) - being the
  !     more common in OMI spectra.

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: solar_idx, us1_idx, us2_idx, wvl_idx, spc_idx
  USE OMSAO_parameters_module, ONLY: max_spec_pts, maxchlen
  USE OMSAO_variables_module,  ONLY: &
       refspecs_original, database, have_undersampling, yn_use_labslitfunc
  USE OMSAO_slitfunction_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                           INTENT (IN) :: n_sensor_pts, xtrack_pix, &
       n_solar_pts
  REAL    (KIND=r8),                           INTENT (IN) :: hw1e, e_asym
  REAL    (KIND=r8), DIMENSION (n_sensor_pts), INTENT (IN) :: curr_wvl
  REAL    (KIND=r8), DIMENSION (n_solar_pts), INTENT (IN)  :: solar_wvl

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  LOGICAL                                     :: yn_full_range
  REAL (KIND=r8), DIMENSION (2,n_sensor_pts)  :: underspec
  REAL (KIND=r8), DIMENSION (max_spec_pts)    :: &
       locwvl, locspec, specmod, tmpwav, over, under, resample

  INTEGER (KIND=i4) :: npts, locerrstat

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'undersample_new'

  locerrstat = pge_errstat_ok

  ! ==================================================
  ! Assign solar reference spectrum to local variables
  ! ==================================================
  npts            = refspecs_original(solar_idx)%nPoints
  locwvl (1:npts) = refspecs_original(solar_idx)%RefSpecWavs(1:npts)
  locspec(1:npts) = refspecs_original(solar_idx)%RefSpecData(1:npts)

  tmpwav(2:n_solar_pts) = solar_wvl
  tmpwav(1)             = curr_wvl(1)
  tmpwav(n_solar_pts)   = curr_wvl(n_sensor_pts)
  IF ( tmpwav(2) <= tmpwav(1) ) tmpwav(2) = (tmpwav(1)+tmpwav(3))/2.0_r8
  IF ( tmpwav(n_solar_pts-1) >= tmpwav(n_solar_pts) ) tmpwav(n_solar_pts-1) = &
                               (tmpwav(n_solar_pts)+tmpwav(n_solar_pts-2))/2.0_r8



  IF ( yn_use_labslitfunc ) THEN
     CALL omi_slitfunc_convolve ( &
          xtrack_pix, npts, locwvl(1:npts), locspec(1:npts), specmod(1:npts), locerrstat )
  ELSE
     CALL asymmetric_gaussian_sf ( &
          npts, hw1e, e_asym, locwvl(1:npts), locspec(1:npts), specmod(1:npts))
  END IF
  CALL error_check ( &
       locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename//f_sep//'Convolution', vb_lev_default, errstat )
  IF ( locerrstat >= pge_errstat_error ) RETURN

  ! Phase1 calculation: Calculate spline derivatives for KPNO data
  !                     Calculate solar spectrum at OMI radiance positions

  CALL interpolation (                                                                &
       modulename//f_sep//'Phase 1a',                                                 &
       npts, locwvl(1:npts), specmod(1:npts), n_sensor_pts, curr_wvl(1:n_sensor_pts), &
       over(1:n_sensor_pts), 'endpoints', 0.0_r8, yn_full_range, locerrstat )
  CALL error_check (                                                    &
       locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename//f_sep//'Phase 1a', vb_lev_default, errstat )
  IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  IF ( .NOT. yn_full_range ) CALL error_check (           &
       0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
       modulename//f_sep//'Phase 1a', vb_lev_develop, errstat )

  ! Convolved High resolution solar to OMI solar grid
  CALL interpolation (                                                               &
       modulename//f_sep//'Phase 1b',                                                &
       npts, locwvl(1:npts), specmod(1:npts), n_solar_pts, tmpwav(1:n_solar_pts), &
       resample(1:n_solar_pts), 'endpoints', 0.0_r8, yn_full_range, locerrstat )
  CALL error_check (                                                    &
       locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename//f_sep//'Phase 1b', vb_lev_default, errstat )
  IF ( locerrstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------
  ! Issue a warning if we don't have the full interpolation range
  ! -------------------------------------------------------------
  !IF ( .NOT. yn_full_range ) CALL error_check (                 &
  !     0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE,       &
  !     modulename//f_sep//'Phase 1b', vb_lev_omidebug, errstat )

  ! Undersample solar to radiance grid
  CALL interpolation (                                                                 &
       modulename//f_sep//'Phase 1c',                                                  &
       n_solar_pts, tmpwav(1:n_solar_pts), resample(1:n_solar_pts), n_sensor_pts, &
       curr_wvl(1:n_sensor_pts), under(1:n_sensor_pts), 'endpoints', 0.0_r8,             &
       yn_full_range, locerrstat )
  CALL error_check (                                                    &
       locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename//f_sep//'Phase 1c', vb_lev_default, errstat )
  IF ( locerrstat >= pge_errstat_error ) RETURN

  ! Calculate undersample spectrum
  underspec(1,1:n_sensor_pts) = over(1:n_sensor_pts) - under(1:n_sensor_pts)

  ! ------------------------------------------------
  ! Save spectra to final arrays for Undersampling 1
  ! ------------------------------------------------
  refspecs_original(us1_idx)%nPoints                     = n_sensor_pts
  refspecs_original(us1_idx)%NormFactor                  = 1.0E+00_R8
  refspecs_original(us1_idx)%RefSpecWavs(1:n_sensor_pts) = curr_wvl(1:n_sensor_pts)
  refspecs_original(us1_idx)%RefSpecData(1:n_sensor_pts) = underspec(1,1:n_sensor_pts)
  refspecs_original(us2_idx)%nPoints                     = n_sensor_pts
  refspecs_original(us2_idx)%NormFactor                  = 1.0E+00_R8
  refspecs_original(us2_idx)%RefSpecWavs(1:n_sensor_pts) = curr_wvl(1:n_sensor_pts)
  refspecs_original(us2_idx)%RefSpecData(1:n_sensor_pts) = underspec(1,1:n_sensor_pts)

  ! Save undersample spectrum to database, for compliance with previous versions
  ! it is saved to both under sample spectra but only one of them needs to be used.
  database(us1_idx,1:n_sensor_pts) = underspec (1,1:n_sensor_pts)
  database(us2_idx,1:n_sensor_pts) = underspec (1,1:n_sensor_pts)

  RETURN
END SUBROUTINE undersample_new
