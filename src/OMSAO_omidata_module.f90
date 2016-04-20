MODULE OMSAO_omidata_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen, max_spec_pts
  USE OMSAO_indices_module,    ONLY: n_max_fitpars, max_rs_idx, max_calfit_idx, o3_t1_idx, o3_t3_idx
  IMPLICIT NONE

  ! ------------------------------------------------------------
  ! Boundary wavelengths (approximate) for UV-2 and VIS channels
  ! ------------------------------------------------------------
  REAL (KIND=r4), PARAMETER :: uv2_upper_wvl = 385.0_r4, vis_lower_wvl = 350.0_r4

  ! ---------------------------------------
  ! Minimum OMI spectral resolution (in nm)
  ! ---------------------------------------
  REAL (KIND=r8), PARAMETER :: omi_min_specres = 0.5_r8

  ! ---------------------------------
  ! Maximum OMI data/swath dimensions
  ! ---------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
       nxtrack_max    =  60, nwavel_max = 1024, &
       nwavelcoef_max =   5, nlines_max =  100, &
       nUTCdim        =   6

  ! -------------------------------------------------------
  ! Maximum dimension of OMI data fields in SAO L2 products
  ! -------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_field_maxdim = 3

  ! --------------------------------------------------
  ! Parameters defined by the NISE snow cover approach
  ! --------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
       NISE_snowfree =   0, NISE_allsnow = 100, NISE_permice = 101, NISE_drysnow = 103, &
       NISE_ocean    = 104, NISE_suspect = 125, NISE_error   = 127

  ! --------------------------------------------------------------------
  ! Values to go into the diagnostic array that shows how the AMF was
  ! computed, with values that indicate missing cloud products, glint,
  ! and geometric or no AMF
  ! --------------------------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
       omi_cfr_addmiss = 1000, omi_ctp_addmiss = 2000, omi_glint_add = 10000, &
       omi_oob_cfr     =  210, omi_oob_ctp     =  220, &
       omi_geo_amf = -1, omi_oobview_amf = -2, omi_wfmod_amf = -9, omi_bigsza_amf = 500

  ! -----------------------
  ! Arrays for OMI L1b data
  ! -----------------------
  INTEGER (KIND=i2)                                                    :: omi_mflg
  REAL    (KIND=r4), DIMENSION (0:nlines_max-1)                        :: omi_auraalt
  REAL    (KIND=r8), DIMENSION (0:nlines_max-1)                        :: omi_time
  INTEGER (KIND=i4), DIMENSION (0:nlines_max-1)                        :: omi_radiance_errstat
  INTEGER (KIND=i1), DIMENSION (nxtrack_max,0:nlines_max-1)            :: omi_xtrflg_l1b
  INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1)            :: omi_geoflg, omi_xtrflg
  INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1)            :: omi_height, land_water_flg
  REAL    (KIND=r4), DIMENSION (nxtrack_max,0:nlines_max-1)            :: omi_latitude, omi_longitude
  REAL    (KIND=r4), DIMENSION (nxtrack_max,0:nlines_max-1)            :: omi_szenith, omi_sazimuth
  REAL    (KIND=r4), DIMENSION (nxtrack_max,0:nlines_max-1)            :: omi_vzenith, omi_vazimuth
  REAL    (KIND=r8), DIMENSION (nwavel_max,nxtrack_max,0:nlines_max-1) :: omi_radiance_spec
  REAL    (KIND=r8), DIMENSION (nwavel_max,nxtrack_max,0:nlines_max-1) :: omi_radiance_prec
  REAL    (KIND=r8), DIMENSION (nwavel_max,nxtrack_max,0:nlines_max-1) :: omi_radiance_wavl
  INTEGER (KIND=i2), DIMENSION (nwavel_max,nxtrack_max,0:nlines_max-1) :: omi_radiance_qflg
  INTEGER (KIND=i4), DIMENSION (nwavel_max,nxtrack_max,0:nlines_max-1) :: omi_radiance_ccdpix
  INTEGER (KIND=i2), DIMENSION (nUTCdim,0:nlines_max-1)                :: omi_time_utc

  ! ---------------------------------------------------------------
  ! Snow/Ice and Glint flags are used in the AMF computation module 
  ! outside the "nlines_max" loops and hence need to be defined on
  ! the maximum swath dimensions.
  ! ---------------------------------------------------------------
  INTEGER (KIND=i4), DIMENSION (nwavel_max,nxtrack_max) :: omi_irradiance_ccdpix
  INTEGER (KIND=i2), DIMENSION (nwavel_max,nxtrack_max) :: omi_irradiance_qflg, omi_radref_qflg
  REAL    (KIND=r8), DIMENSION (nwavel_max,nxtrack_max) :: &
       omi_irradiance_prec, omi_irradiance_wavl, omi_irradiance_spec, omi_irradiance_wght, &
       omi_radref_spec, omi_radref_wavl, omi_radref_wght
  REAL    (KIND=r4), DIMENSION (nxtrack_max) :: omi_radref_sza, omi_radref_vza

  ! ---------------------------------------------------------------
  ! Scene Albedo, from OMCLDO2, for Wavelength-Modified AMF lookup
  ! ---------------------------------------------------------------
  REAL (KIND=r4), DIMENSION (nxtrack_max) :: omi_scene_albedo


  INTEGER (KIND=i4), DIMENSION (nxtrack_max,4) :: omi_ccdpix_selection
  INTEGER (KIND=i4), DIMENSION (nxtrack_max,2) :: omi_ccdpix_exclusion

  ! ----------------------------------------
  ! Arrays for fitting and/or derived output
  ! ----------------------------------------
  REAL    (KIND=r8), PARAMETER                              :: d_comm_wvl = 0.01_r8
  INTEGER (KIND=i4)                                         :: n_comm_wvl
  INTEGER (KIND=i4), DIMENSION (nxtrack_max)                :: common_cnt
  REAL    (KIND=r8), DIMENSION (nxtrack_max,max_spec_pts)   :: common_spc, common_wvl
  REAL    (KIND=r8), DIMENSION (nxtrack_max,0:nlines_max-1) :: &
       omi_column_amount, omi_column_uncert, omi_fit_rms, omi_radfit_chisq
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION(:,:)            :: omi_ozone_amount
  REAL    (KIND=r4), DIMENSION (nxtrack_max,0:nlines_max-1) :: omi_razimuth
  INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1) :: omi_fitconv_flag
  INTEGER (KIND=i2), DIMENSION (nxtrack_max,0:nlines_max-1) :: omi_itnum_flag

  ! ----------------------------------------------------------------------------
  ! Correlations with main output product. Due to a bug in the HDF-EOS5 routines
  ! (non-TLCF implementation), STRING fields cannot be written to file directly.
  ! A work-around solution is to convert the CHARACTERs to INTEGERs. Thus the
  ! need for the additional array CORRELATION_NAMES_INT.
  ! ----------------------------------------------------------------------------
  CHARACTER (LEN=maxchlen), DIMENSION (n_max_fitpars) :: correlation_names
  CHARACTER (LEN=n_max_fitpars*maxchlen)              :: correlation_names_concat

  ! --------------------------------------------------------
  ! Ozone is a special case: We can have up to 3 temperatues
  ! --------------------------------------------------------
  REAL (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx, nxtrack_max,0:nlines_max-1) :: &
       omi_o3_amount, omi_o3_uncert

  ! ---------------------------------
  ! Dimensions for measurement swaths
  ! ---------------------------------
  !INTEGER (KIND=i4) :: ntimes, ntimessmallpixel, nxtrack, nwavel, ntimes_loop, nwavel_ccd
  INTEGER (KIND=i4) :: ntimessmallpixel, nwavel, ntimes_loop
  INTEGER (KIND=i4), DIMENSION (nxtrack_max)                  :: omi_nwav_irrad, omi_nwav_radref
  INTEGER (KIND=i4), DIMENSION (nxtrack_max,0:nlines_max-1)   :: omi_nwav_rad

  ! ---------------------------------------
  ! Swath attributes for measurement swaths
  ! ---------------------------------------
  INTEGER (KIND=i2) :: ImageBinningFactor, BinnedImageRows, StopColumn
  INTEGER (KIND=i4) :: NumTimes, NumTimesSmallPixel

  INTEGER (KIND=i4), DIMENSION (nxtrack_max)                         :: n_omi_database_wvl
  INTEGER (KIND=i2), DIMENSION (nxtrack_max)                         :: &
       omi_solcal_itnum, omi_radcal_itnum, omi_radref_itnum,            &
       omi_solcal_xflag, omi_radcal_xflag, omi_radref_xflag
  REAL    (KIND=r8), DIMENSION (max_calfit_idx, nxtrack_max)         :: &
       omi_solcal_pars,  omi_radcal_pars,  omi_radref_pars
  REAL    (KIND=r8), DIMENSION (max_rs_idx, nwavel_max, nxtrack_max) :: omi_database
  REAL    (KIND=r8), DIMENSION (            nwavel_max, nxtrack_max) :: omi_database_wvl
  REAL    (KIND=r8), DIMENSION (nxtrack_max)                         :: omi_sol_wav_avg
  REAL    (KIND=r8), DIMENSION (nxtrack_max)                         :: &
       omi_solcal_chisq, omi_radcal_chisq, omi_radref_chisq, &
       omi_radref_col,   omi_radref_dcol,  omi_radref_rms,   &
       omi_radref_xtrcol
  REAL    (KIND=r8), DIMENSION (2,nxtrack_max,0:nlines_max-1)        :: omi_wavwin_rad, omi_fitwin_rad
  REAL    (KIND=r8), DIMENSION (2,nxtrack_max)                       :: omi_wavwin_sol, omi_fitwin_sol

  ! ---------------
  ! OMI swath names
  ! ---------------
  CHARACTER (LEN=maxchlen) :: omi_radiance_swathname, omi_irradiance_swathname
  CHARACTER (LEN=maxchlen) :: l1b_radiance_esdt

  ! ------------------------------
  ! Distance between Earth and Sun
  ! ------------------------------
  REAL (KIND=r4) :: EarthSunDistance

  ! ---------------------------------------------------------
  ! OMI scan line, block line, and across-track pixel numbers
  ! ---------------------------------------------------------
  INTEGER (KIND=i4) :: omi_scanline_no, omi_blockline_no, omi_xtrackpix_no

  ! ---------------------------
  ! OMI L2 output data QA flags
  ! ----------------------------------------------------------------
  ! Per Centages of output column data that are used to classify the
  ! scientific data quality:
  !         Good data >= QA_PERCENT_PASS   : "Passed"
  !         Good data >= QA_PERCENT_SUSPECT: "Suspect"
  !         Good data <  QA_PERCENT_SUSPECT: "Failed"
  ! ----------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: qa_percent_passed = 75, qa_percent_suspect = 50


  ! ------------------------------------------------------
  ! Finally some variables that will be initialized in the
  ! course of the processing.
  ! ------------------------------------------------------
  INTEGER (KIND=i4) :: &
       n_omi_radwvl, n_omi_radrefwvl, n_omi_irradwvl,  &
       nwavelcoef_irrad, nwavelcoef_rad,               &
       ntimes_smapix_irrad, ntimes_smpix_rad, nclenfit

  ! --------------------------------
  ! Current cross-track pixel number
  ! --------------------------------
  INTEGER (KIND=i4) :: curr_xtrack_pixnum

  LOGICAL, DIMENSION (nxtrack_max) ::  omi_cross_track_skippix = .FALSE.

  ! ------------------------------------
  ! Number of significant digits to keep
  ! ------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_roff_dig = 5

  ! ----------------------------------------------------------
  ! Variables and parameters associated with Spatial Zoom data
  ! ----------------------------------------------------------
  CHARACTER (LEN=16), PARAMETER :: uv1_glob_swath = 'Earth UV-1 Swath'
  CHARACTER (LEN=10), PARAMETER ::    &
       uv1_zoom_swath = '(60x159x4)', &
       uv2_zoom_swath = '(60x557x4)', &
       vis_zoom_swath = '(60x751x4)'

  INTEGER (KIND=i4), PARAMETER :: gzoom_spix = 16, gzoom_epix = 45, gzoom_npix = 30
  INTEGER (KIND=i1), PARAMETER :: global_mode = 8_i1, szoom_mode = 4_i1
  INTEGER (KIND=i1)            :: truezoom, fullswath

  INTEGER (KIND=i4) :: omi_l1b_idx

  ! -----------------------------------------
  ! Index to keep track of the line we are in
  ! -----------------------------------------
  INTEGER (KIND=i4) :: omi_iline
END MODULE OMSAO_omidata_module
