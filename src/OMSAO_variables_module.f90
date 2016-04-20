MODULE OMSAO_variables_module

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       max_rs_idx, max_calfit_idx, n_max_fitpars, mxs_idx, sig_idx, icf_idx, &
       o3_t1_idx, o3_t2_idx, o3_t3_idx,                                      &
       us1_idx, us2_idx, n_voc_amf_luns, ccd_idx, radfit_idx
       !n_amftab_ang_max, n_amftab_dim_max, n_amftab_coef_max,                &

  USE OMSAO_parameters_module,   ONLY: maxchlen, max_spec_pts, n_fit_winwav, max_mol_fit
  USE OMSAO_omidata_module,      ONLY: nwavel_max, nxtrack_max, nlines_max
  USE EZspline_obj

  IMPLICIT NONE

  ! -----------------------------
  ! Variables read from PCF file
  ! -----------------------------
  ! * Current PGE name and index
  ! -----------------------------
  INTEGER (KIND=I4) :: pge_idx
  CHARACTER (LEN=6) :: pge_name
  ! ---------------------
  ! * Verbosity threshold
  ! ---------------------
  CHARACTER (LEN=1) :: verb_thresh_char
  INTEGER (KIND=I4) :: verb_thresh_lev
  ! ----------------------------------------------------
  ! * Orbit number, Version ID, L1B radiance OPF version
  ! ----------------------------------------------------
  INTEGER (KIND=I4) :: orbit_number, ecs_version_id, l1br_opf_version

  ! -------------------------------------------------
  ! Variables defined in preamble of original program
  ! -------------------------------------------------
  LOGICAL :: yn_smooth, yn_doas
  LOGICAL :: yn_spectrum_norm

  INTEGER (KIND=I4), DIMENSION (n_max_fitpars)  :: mask_fitvar_rad, mask_fitvar_cal, all_radfit_idx

  INTEGER (KIND=I4)                             :: n_fitvar_rad, n_fitvar_cal

  REAL    (KIND=r8), DIMENSION (max_calfit_idx) :: &
       fitvar_cal, fitvar_cal_saved, fitvar_sol_init, lo_sunbnd, up_sunbnd

  REAL    (KIND=r8), DIMENSION (n_max_fitpars)  :: fitvar_rad, fitvar_rad_init
  REAL    (KIND=r8), DIMENSION (n_max_fitpars)  :: fitvar_rad_saved
  REAL    (KIND=r8), DIMENSION (n_max_fitpars)  :: lo_radbnd, up_radbnd, lobnd, upbnd
  CHARACTER (LEN=6), DIMENSION (n_max_fitpars)  :: fitvar_rad_str, fitvar_sol_str

  REAL    (KIND=r8), DIMENSION (max_rs_idx, nwavel_max) :: database

  ! -------------------------------------
  ! Variables related to Air Mass Factors
  ! -------------------------------------
  !INTEGER (KIND=I4)  :: n_amftab_dim, n_amftab_ang
  LOGICAL            :: have_amftable, yn_o3amf_cor
  !!REAL    (KIND=r8)  :: amf_esza_min, amf_esza_max
  !REAL    (KIND=r8)  :: amf_tab_wvl, amf_tab_alb
  !REAL    (KIND=r8), DIMENSION (n_amftab_ang_max, n_amftab_dim_max) :: amf_table_bro

  !! ---------------------------------------------------------------
  !! Some BrO specific AMF variables: Fitted polynomial coefficients
  !! that are used in lieu of interpolation.
  !! ---------------------------------------------------------------
  !INTEGER (KIND=i4)                                  :: n_amftab_coef
  !REAL    (KIND=r8)                                  :: amftab_coef_norm
  !REAL    (KIND=r8), DIMENSION (0:n_amftab_coef_max) :: amftab_coef

  ! -----------------------------
  ! Previously IMPLICIT variables
  ! -----------------------------
  REAL (KIND=r8) :: phase, szamax, chisq, sol_wav_avg, rad_wav_avg
  REAL (KIND=r8) :: hw1e, e_asym
  REAL (KIND=r4) :: zatmos

  ! -----------------------------------------------------------
  ! Variables related to reference spectra
  !  * number of reference spectra:       N_REfSPEC
  !  * indentification strings:           FITPAR_IDXNAME
  !  * file names with reference spectra: REFSPEC_FNAME
  !  * original (uniterpolated) data:     REFSPEC_ORIG_DATA
  !  * number of spectral points:         N_REFSPEC_PTS
  !  * first and last wavelenghts:        REFSPEC_FIRSTLAST_WAV
  ! -----------------------------------------------------------
  INTEGER   (KIND=I4) :: n_refspec

  ! -------------------------------------
  ! TYPE declaration for Reference Specta
  ! -------------------------------------
  TYPE, PUBLIC :: ReferenceSpectrum
     CHARACTER (LEN=maxchlen)                      :: Title, Units
     CHARACTER (LEN=maxchlen)                      :: FileName
     CHARACTER (LEN=maxchlen)                      :: FittingIdxName
     INTEGER   (KIND=I4)                           :: nPoints
     REAL      (KIND=r8)                           :: NormFactor, Temperature
     REAL      (KIND=r8), DIMENSION (2)            :: FirstLastWav
     REAL      (KIND=r8), DIMENSION (max_spec_pts) :: RefSpecWavs
     REAL      (KIND=r8), DIMENSION (max_spec_pts) :: RefSpecData     
  END TYPE ReferenceSpectrum

  ! -----------------------------------------
  ! TYPE declaration for Common Mode Spectrum
  ! -----------------------------------------
  TYPE, PUBLIC :: CommonModeSpectrum
     CHARACTER (LEN=maxchlen)                                  :: Title, Units
     CHARACTER (LEN=maxchlen)                                  :: FileName
     CHARACTER (LEN=maxchlen)                                  :: FittingIdxName
     INTEGER   (KIND=I4)                                       :: nPoints
     REAL      (KIND=r8)                                       :: NormFactor, Temperature
     REAL      (KIND=r8), DIMENSION (2)                        :: FirstLastWav
     INTEGER   (KIND=I2), DIMENSION (nxtrack_max,2)            :: CCDPixel
     INTEGER   (KIND=I4), DIMENSION (nxtrack_max)              :: RefSpecCount
     REAL      (KIND=r8), DIMENSION (nxtrack_max,max_spec_pts) :: RefSpecWavs
     REAL      (KIND=r8), DIMENSION (nxtrack_max,max_spec_pts) :: RefSpecData     
  END TYPE CommonModeSpectrum

  ! -------------------------------
  ! Array for all Reference Spectra
  ! -------------------------------
  TYPE (ReferenceSpectrum),  DIMENSION (max_rs_idx) :: refspecs_original
  TYPE (CommonModeSpectrum)                         :: common_mode_spec

  ! -------------------------------------------
  ! A special beast: The undersampling spectrum
  ! -------------------------------------------
  LOGICAL, DIMENSION (us1_idx:us2_idx) :: have_undersampling


  CHARACTER (LEN=maxchlen), DIMENSION (icf_idx:max_rs_idx) :: static_input_fnames

  REAL (KIND=r8), DIMENSION (max_spec_pts) :: cubic_x, cubic_y, cubic_w

  ! --------------------------------------
  ! Solar and Earth shine wavlength limits
  ! --------------------------------------
  REAL (KIND=r8), DIMENSION (n_fit_winwav) :: fit_winwav_lim
  REAL (KIND=r8), DIMENSION (2)            :: fit_winexc_lim
  REAL (KIND=r8)                           :: winwav_min, winwav_max

  ! ------------------------------------------------------------------
  ! Indices of fitting window defining wavelengths in current spectrum
  ! ------------------------------------------------------------------
  INTEGER (KIND=I4), DIMENSION (n_fit_winwav) :: fit_winwav_idx

  ! --------------------------------------------------------------------------
  ! The current solar and radiance spectrum, including wavelengths and weights
  ! --------------------------------------------------------------------------
  INTEGER (KIND=i4), DIMENSION (2)                      :: radiance_wavcal_lnums
  REAL    (KIND=r8), DIMENSION (ccd_idx, nwavel_max)    :: curr_rad_spec
  REAL    (KIND=r8), DIMENSION (ccd_idx, nwavel_max)    :: curr_sol_spec
  REAL    (KIND=r8), DIMENSION (nwavel_max)             :: fitwavs, fitweights, currspec
  REAL    (KIND=r8), DIMENSION (nwavel_max,nxtrack_max) :: rad_spec_wavcal, rad_wght_wavcal

  ! ----------------------------------------
  ! Pixel number limits:
  !  * 1,2: First and last scan line number
  !  * 3,4: First and last cross-track pixel
  ! ----------------------------------------
  INTEGER (KIND=I4), DIMENSION (4) :: pixnum_lim

  ! ----------------------------------------------------
  ! Latitude limits: Lower and upper latitude to process
  ! ----------------------------------------------------
  REAL (KIND=r4), DIMENSION (2) :: radfit_latrange

  ! ---------------------------------------------------------------------
  ! Contstraints on the fitting residual: Window and Number of Iterations
  ! ---------------------------------------------------------------------
  INTEGER (KIND=i4), DIMENSION (radfit_idx)  :: n_fitres_loop, fitres_range
  REAL    (KIND=r8), DIMENSION (nxtrack_max) :: xtrack_fitres_limit

  ! --------------------------------------------
  ! Frequency of radiance wavelength calibration
  ! --------------------------------------------
  INTEGER (KIND=I4) :: radwavcal_freq

  ! ------------------------------------------------------------------------
  ! Variables connected with ELSUNC numerical precision/convergence criteria
  ! ------------------------------------------------------------------------
  REAL (KIND=r8) :: tol,  epsrel,  epsabs,  epsx

  ! ----------------------------------------
  ! Variable for +1.0 or -1.0 multiplication
  ! ----------------------------------------
  REAL (KIND=r8) :: pm_one

  ! ----------------------------------------------------------------------
  ! Index for the fitting parameters carrying the fitted column value.
  !
  ! N_MOL_FIT:    Number of "molecules" that carry the final column; this
  !               can be one molecule at different temperatures.
  ! FITCOL_IDX:   The main molecule indices, corresponding to the list of
  !               reference spectra.
  ! FINCOL_IDX:   For the final summation of the fitted column: The total
  !               number is the number of different molecules times the
  !               allowed sub-indices. The second dimension is for the
  !               reference spectrum index - this eases the final sum over
  !               the fitted columns (see RADIANCE_FIT subroutine).
  !                includes
  !               the subindices, hence the dimension.
  ! N_FINCOL_IDX: Number of final column indices.
  ! ----------------------------------------------------------------------
  INTEGER (KIND=I4)                                    :: n_mol_fit, n_fincol_idx
  INTEGER (KIND=I4), DIMENSION (max_mol_fit)           :: fitcol_idx
  INTEGER (KIND=I4), DIMENSION (2,max_mol_fit*mxs_idx) :: fincol_idx

  ! ------------------------------------
  ! Maximum number of fitting iterations
  ! ------------------------------------
  INTEGER (KIND=I4) :: max_itnum_sol, max_itnum_rad
  

  ! ------------------------------------
  ! Maximum good column amount
  ! ------------------------------------
  REAL (KIND=r8) :: max_good_col
  
  ! ----------------------------------------------------------------------------
  ! Number of calls to fitting function. This is counted in the fitting function
  ! itself, in an attempt to catch infinite loops of the ELSUNC routine, which
  ! occur on occasion for reasons not yet known.
  ! ----------------------------------------------------------------------------
  INTEGER (KIND=i4) :: num_fitfunc_calls, num_fitfunc_jacobi, max_fitfunc_calls

  ! ---------------------
  ! L1B and L2 file names
  ! ---------------------
  CHARACTER (LEN=maxchlen) :: l1b_rad_filename, l1b_irrad_filename, l2_filename

  ! -----------------------------------------------------------------
  ! Generic dimension variables (initialized from either GOME or OMI)
  ! -----------------------------------------------------------------
  INTEGER (KIND=I4) :: n_sol_wvl, n_rad_wvl, n_database_wvl
  INTEGER (KIND=I4) :: n_rad_wvl_max

  ! --------------------------------------------
  ! Name of the tabulated OMI slit function data
  ! --------------------------------------------
  CHARACTER (LEN=maxchlen) :: omi_slitfunc_fname
  LOGICAL                  :: yn_use_labslitfunc

  ! ---------------------------------------------------------------------------
  ! And this is the convolved solar spectrum. It is (re)initialized in the
  ! solar fit routines only if the above shift&squeeze parameters have changed.
  ! ---------------------------------------------------------------------------
  REAL (KIND=r8), DIMENSION (max_spec_pts) :: solar_spec_convolved


  ! ---------------------------------------------------------
  ! Filenames specific for the AMF scheme in OMBRO and OMHCHO
  ! ---------------------------------------------------------
  CHARACTER (LEN=maxchlen)                               :: OMBRO_amf_filename
  CHARACTER (LEN=maxchlen), DIMENSION (n_voc_amf_luns)   :: voc_amf_filenames

  ! ---------------------------------------------------------------
  ! Filename, logical and type indices for composite Solar Spectrum
  ! ---------------------------------------------------------------
  CHARACTER (LEN=maxchlen) :: OMSAO_solcomp_filename
  LOGICAL                  :: yn_solar_comp
  INTEGER (KIND=i4)        :: solar_comp_typ, solar_comp_orb

  ! -------------------------------------------------------
  ! Filename and logical for solar monthly average spectrum
  ! -------------------------------------------------------
  CHARACTER (LEN=maxchlen) :: OMSAO_solmonthave_filename
  LOGICAL                  :: yn_solmonthave

  ! -------------------------------
  ! Logicals for Solar I0 correction
  ! -------------------------------
  LOGICAL           :: yn_solar_i0

  ! ------------------------------------------------------
  ! Logical for processing mode (diagnostic or production)
  ! ------------------------------------------------------
  LOGICAL           :: yn_diagnostic_run

  ! --------------------------------------------
  ! * Logical for Common Mode Iteration
  ! * Index position of Common Mode add-on
  ! * Array for initial fitting variables
  ! * Latitude range for common mode computation
  ! --------------------------------------------
  LOGICAL                          :: yn_common_iter
  INTEGER (KIND=i4)                :: common_fitpos
  REAL    (KIND=r8), DIMENSION (3) :: common_fitvar
  REAL    (KIND=r4), DIMENSION (2) :: common_latrange
  INTEGER (KIND=i4), DIMENSION (2) :: common_latlines

  ! ---------------------------------------------------------------------
  ! 3-letter string to identify the OMI channel (UV2 or VIS) we are using
  ! ---------------------------------------------------------------------
  CHARACTER (LEN=3) :: l1b_channel

  ! -------------------------------------------------
  ! Logical for newshift following Xiong comments gga
  ! -------------------------------------------------
  LOGICAL :: yn_newshift

  ! --------------------------------------------------------
  ! Filename and logical for Reference Sector Correction gga
  ! --------------------------------------------------------
  CHARACTER (LEN=maxchlen) :: OMSAO_refseccor_filename
  CHARACTER (LEN=maxchlen) :: OMSAO_refseccor_cld_filename
  LOGICAL :: yn_refseccor

  ! -----------------------------------------------------------------
  ! Logical for Scattering Weights, Gas Profile and Averaging Kernels
  ! Also filename
  ! -----------------------------------------------------------------
  LOGICAL                  :: yn_sw
  CHARACTER (LEN=maxchlen) :: OMSAO_OMLER_filename

  ! ---------------------------------------------
  ! Variables to hold ozone index information for
  ! AMF calculation
  ! ---------------------------------------------
  INTEGER (KIND=i4) :: ozone_idx
  LOGICAL           :: ozone_log = .TRUE.

END MODULE OMSAO_variables_module
