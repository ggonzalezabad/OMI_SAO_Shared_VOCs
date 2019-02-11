MODULE OMSAO_parameters_module

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ================================================
  ! Define array dimensions and other parameters for
  ! SAO Trace Gas Fitting Algorithms
  ! ================================================

  ! -----------------------------------
  ! Maximum length of CHARACTER strings
  ! -----------------------------------
  INTEGER (KIND=i4), PARAMETER :: maxchlen = 256

  ! ----------------------------------
  ! Maximum iteration number for loops
  ! ----------------------------------
  INTEGER (KIND=i4), PARAMETER :: forever = HUGE(1_i4)

  ! -------------------------------------------------------------------------
  ! Maximum numbers for fitting parameters, GOME pixels, spectral points, ...
  ! -------------------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: max_spec_pts = 15001

  ! ==================
  ! Physical constants
  ! ==================
  REAL (KIND=r8), PARAMETER :: pi = 3.14159265358979_r8
  REAL (KIND=r8), PARAMETER :: deg2rad = pi / 180.0_r8, rad2deg = 180.0_r8 / pi
  REAL (KIND=r8), PARAMETER :: smallval = 1.0E-10_r8
  ! ----------------------------------
  ! Edlen 1966 standard air parameters
  ! ----------------------------------
  REAL (KIND=r8), PARAMETER :: &
       c1 = 1.0000834213_r8, c2 = 2.406030E-02_r8, c3 = 130.0_r8, &
       c4 = 1.5997E-04_r8, c5 = 38.9_r8

  ! ------------
  ! Dobson units
  ! ------------
  REAL (KIND=r8), PARAMETER :: dobson_units = 2.68676E+16_r8

  ! -------------------------------------------------------------
  ! "Typical" O3 and O3 loading for cases of I0 Effect correction
  ! -------------------------------------------------------------
  REAL (KIND=r8), PARAMETER :: &
       ozone_300du = 300.0_r8 * dobson_units, &
       no2_1du     =   1.0_r8 * dobson_units

	! ---------------------------------------------------------
	! "Typical" Slant Columns for all species for I0 correction
	! 
	! NB Modified for glyoxal retrieval to work when the 
	!    yn_solar_i0 switch set. May require seperate switch
	!    if others do not want it CCM
	!    
	! ---------------------------------------------------------
	REAL (KIND=r8), DIMENSION(27), PARAMETER :: &
			 solar_i0_scd = (/   1.0_r8 * dobson_units,& ! solar
			 									   1.0_r8 * dobson_units,& ! ring
			 									 300.0_r8 * dobson_units,& ! o3_t1
			 									 300.0_r8 * dobson_units,& ! o3_t2
			 									 300.0_r8 * dobson_units,& ! o3_t3
			 									   2.0d16               ,& ! no2_t1
			 									   2.0d16               ,& ! no2_t2
			 									   5.0d43 * dobson_units,& ! o2o2
			 									   1.0_r8 * dobson_units,& ! so2 
			 									   1.0_r8 * dobson_units,& ! bro
			 									   1.0_r8 * dobson_units,& ! oclo
			 									   1.0_r8 * dobson_units,& ! hcho
				 									 2.0d23               ,& ! h2o
				 									 2.0d01               ,& ! lqh2o
				 									 5.0d14               ,& ! glyox
				 									 1.0_r8 * dobson_units,& ! io
				 									 1.0_r8 * dobson_units,& ! hono
				 									 1.0_r8 * dobson_units,& ! vraman
				 									 1.0_r8 * dobson_units,& ! commod
				 									 1.0_r8 * dobson_units,& ! resid
				 									 1.0_r8 * dobson_units,& ! pseudo
				 									 1.0_r8 * dobson_units,& ! polcor
				 									 1.0_r8 * dobson_units,& ! usamp1
				 									 1.0_r8 * dobson_units,& ! usamp2
				 									 1.0_r8 * dobson_units,& ! bro_tc
				 									 1.0_r8 * dobson_units,& ! o3_tc
			 									   1.0_r8 * dobson_units /)! noname
	
	! Logicals to do correction CCM
	LOGICAL, DIMENSION(27), PARAMETER :: &
			 yn_i0_spc = (/ .FALSE.,& ! solar
			 								.FALSE.,& ! ring
			 								.TRUE. ,& ! o3_t1
			 							        .TRUE. ,& ! o3_t2
			 								.FALSE.,& ! o3_t3
			 								.TRUE., & ! no2_t1
			 								.TRUE., & ! no2_t2
			 								.FALSE.,& ! o2o2
			 								.FALSE.,& ! so2 
			 								.FALSE.,& ! bro
			 								.FALSE.,& ! oclo
			 								.FALSE.,& ! hcho
				 							.TRUE., & ! h2o
				 							.FALSE.,& ! lqh2o
				 							.TRUE., & ! glyox
				 							.TRUE., & ! io
				 							.FALSE.,& ! hono
				 							.FALSE.,& ! vraman
				 							.FALSE.,& ! commod
				 							.FALSE.,& ! resid
				 							.FALSE.,& ! pseudo
				 							.FALSE.,& ! polcor
				 							.FALSE.,& ! usamp1
				 							.FALSE.,& ! usamp2
				 							.FALSE.,& ! bro_tc
				 							.FALSE.,& ! o3_tc
			 								.FALSE. /)! noname

  ! --------------------------------------------
  ! Average Earth Radius and EOS-Aura eleveation
  ! --------------------------------------------
  REAL (KIND=r4), PARAMETER :: earth_radius_avg = 6371.0_r4
  REAL (KIND=r4), PARAMETER :: eos_aura_avgalt  =  705.0_r4

  ! =============================================================
  ! Large weight for wavelengths to be excluded from the fitting;
  ! small weight for wavelengths to be included in the fitting
  ! =============================================================
  REAL (KIND=r8), PARAMETER :: downweight = 1.0E-30_r8, normweight = 1.0_r8


  ! =========================================
  ! Order of polynomial for DOAS baseline fit
  ! =========================================
  INTEGER (KIND=i4), PARAMETER :: doas_npol = 4

  ! ===================================================
  ! Dimension parameters for ELSUNC auxiliary variables
  ! ===================================================
  INTEGER (KIND=i4), PARAMETER :: elsunc_np = 11, elsunc_nw = 6

  
  ! ---------------------------------------
  ! ELSUNC Abnormal Termination Exit Values
  ! ---------------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
       elsunc_nostart_eval =    -1, &  ! Wrong dimensions or wrong starting point
       elsunc_maxiter_eval =    -2, &  ! Interation exceeded maximum allowed number
       elsunc_hessian_eval =    -3, &  ! Interation exceeded maximum allowed number
       elsunc_no2deri_eval =    -4, &  ! No Second Derivative
       elsunc_newtstp_eval =    -5, &  ! Undamped Newton Step is a failure
       elsunc_nodecst_eval =    -6, &  ! Last step was not descending
       elsunc_onesolu_eval =    -7, &  ! Only one feasible point
       elsunc_parsoob_eval =   -11, &  ! User-defined: Fitting parameters out of bounds
       elsunc_infloop_eval =   -12, &  ! User-defined: Hit "infinite loop" snag
       elsunc_highest_eval = 12344, &  ! Largest possible exit value
       elsunc_usrstop_eval = elsunc_infloop_eval
  REAL    (KIND=r8), PARAMETER ::   &    ! R8 versions for Valid entries
       elsunc_usrstop_eval_r8 =   -10_r8, &
       elsunc_highest_eval_r8 = 12344_r8

  ! -------------------------------------------------------------
  ! Elsunc exit value below which we are computing at noise level
  ! -------------------------------------------------------------
  INTEGER (KIND=i2), PARAMETER :: elsunc_less_is_noise = 300_i2

  ! ----------------------------------
  ! Entries for main quality data flag
  ! ----------------------------------
  INTEGER (KIND=i2), PARAMETER :: &
       main_qa_missing = -1_i2, main_qa_good = 0_i2, main_qa_suspect = 1_i2, main_qa_bad = 2_i2
  INTEGER (KIND=i2), PARAMETER :: &
       main_qa_min_flag = main_qa_missing, main_qa_max_flag = main_qa_bad
  REAL    (KIND=r8), PARAMETER :: &
       main_qa_min_flag_r8 = REAL(main_qa_missing, KIND=r8), &
       main_qa_max_flag_r8 = REAL(main_qa_bad,     KIND=r8)

  ! ===============================================================
  ! Number of wavelengths defining the Solar and Earthshine fitting
  ! windows. N windows require 2*N+2 bounding wavelenghts: Two each
  ! for each window, plus first and last wavlength.
  !
  ! ===============================================================
  INTEGER (KIND=i4), PARAMETER :: n_fit_winwav = 4

  ! ----------------
  ! ZERO_SPEC string
  ! ----------------
  CHARACTER (LEN=9), PARAMETER :: zerospec_string = 'Zero_Spec'

  ! -----------------------------------------------------------------
  ! Maximum number of molecules to fit. This number can be greater
  ! than ONE: We may include multiple reference spectra for the same
  ! molecule in the fit, e.g., at different temperatures. In this
  ! case we need to collect the fitted columns from multiple indices.
  ! -----------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: max_mol_fit = 3

  ! -------------------------
  ! "start of table" landmark
  ! -------------------------
  CHARACTER (LEN=39), PARAMETER :: &
       lm_start_of_table = "start of table (don't delete this line)"


  ! -------------------------------------------------------------------------------
  ! Missing values: We define two sets of the same value, but with different names.
  ! --------------- one set follows the naming convention of the OMI-GDPS-IODS
  ! document (Table 4-14), the other one follows more closely the KIND definition
  ! used in the present PGEs.
  !
  ! NOTE that the values for R4_MISSVAL and R8_MISSVAL are identical. This is for
  ! the cases of quantities that are defined as R8 in the PGE but are written as
  ! R4 to the output file. If R8_MISSVAL was a truly R8 value (e.g., -1.0E+300),
  ! then the conversionto R4 and the HE5 write would fail.
  ! -------------------------------------------------------------------------------
  CHARACTER (LEN=9),   PARAMETER :: str_missval     = "undefined"
  INTEGER   (KIND=i1), PARAMETER :: int8_missval    = -100         !-127
  INTEGER   (KIND=i2), PARAMETER :: int16_missval   = -30000       !-32767
  INTEGER   (KIND=i2), PARAMETER :: int16_missval_l1= -32767
  INTEGER   (KIND=i4), PARAMETER :: int32_missval   = -2000000000  !-2147483647
  REAL      (KIND=r4), PARAMETER :: float32_missval = -1.0E+30_r4  !-1.0_r4*(2.0_r4**100)
  REAL      (KIND=r8), PARAMETER :: float64_missval = -1.0E+30_r8  !-HUGE(1.0_r8) !-1.0_r8*(2.0_r8**100)

  INTEGER   (KIND=i1), PARAMETER :: i1_missval    = int8_missval
  INTEGER   (KIND=i2), PARAMETER :: i2_missval    = int16_missval
  INTEGER   (KIND=i2), PARAMETER :: i2_missval_l1 = int16_missval_l1
  INTEGER   (KIND=i4), PARAMETER :: i4_missval    = int32_missval
  REAL      (KIND=r4), PARAMETER :: r4_missval    = float32_missval
  REAL      (KIND=r8), PARAMETER :: r8_missval    = float64_missval


  ! --------------------------------------------------
  ! Some generic Valid Ranges of the output quantities
  ! --------------------------------------------------
  REAL (KIND=r8), PARAMETER :: &
       valid_min_r8 = -1.0E+30_r8, valid_max_r8 = +1.0E+30_r8, zero_r8 = 0.0_r8, one_r8 = 1.0_r8, &
       valid_max_i2 = REAL(HUGE(1_i2), KIND=r8), valid_max_i4 = REAL(HUGE(1_i4), KIND=r8)
  

  ! -----------------------------------------------------------------
  ! Blank strings of various lengths.
  !
  ! Some compilers don't allow to define CHARACTER PARAMETER arrays
  ! which are initialized with field of unequal length. The following
  ! are a few padding strings that we attach to shorter entries. Not
  ! very stylish, but effective. Note that what matters is the LEN
  ! declaration - we don't need to initialize the strings with the
  ! appropriate number of blanks. Using "" is perfectly fine.
  ! -----------------------------------------------------------------
  CHARACTER (LEN=13), PARAMETER :: blank13 = ""
  CHARACTER (LEN=21), PARAMETER :: blank21 = ""
  CHARACTER (LEN=23), PARAMETER :: blank23 = ""
  CHARACTER (LEN=24), PARAMETER :: blank24 = ""
  CHARACTER (LEN=25), PARAMETER :: blank25 = ""
  CHARACTER (LEN=27), PARAMETER :: blank27 = ""
  CHARACTER (LEN=30), PARAMETER :: blank30 = ""

  ! ----------------------------------------------------------
  ! Some parameters connected with acceptable viewing geometry
  ! ----------------------------------------------------------
  REAL (KIND=r4), PARAMETER :: &
       min_zenith    =     0.0_r4, max_zenith    =  90.0_r4, &  ! "non-inclusive"
       min_azimuth   =  -360.0_r4, max_azimuth   = 360.0_r4, &  ! "inclusive"
       min_latitude  =   -90.0_r4, max_latitude  =  90.0_r4, &  ! "inclusive"
       min_longitude =  -180.0_r4, max_longitude = 180.0_r4     ! "inclusive"

END MODULE OMSAO_parameters_module
