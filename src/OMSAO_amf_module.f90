MODULE OMSAO_AMF_module

  ! ===========================================
  ! Module for the OMBRO and OMHCHO AMF schemes
  ! ===========================================

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY:                 &
       maxchlen, i2_missval, r4_missval, r8_missval, &
       min_zenith, max_zenith, lm_start_of_table
  USE OMSAO_indices_module,    ONLY:                               &
       pge_bro_idx, pge_hcho_idx, pge_gly_idx,                     &
       n_voc_amf_luns, voc_amf_idx, voc_isccp_idx, voc_omicld_idx, &
       voc_amf_luns, voc_amf_lun_str, OMBRO_amf_lun
  USE OMSAO_omidata_module,    ONLY: &
       nxtrack_max, nlines_max,   &
       NISE_snowfree, NISE_ocean, &
       omi_latitude, omi_longitude, omi_szenith, omi_vzenith,               &  ! r4 arrays
       gzoom_spix, gzoom_epix, gzoom_npix,                                  &  ! i4 integers
       omi_cfr_addmiss, omi_ctp_addmiss, omi_glint_add, omi_geo_amf,        &  ! I2 integers
       omi_oobview_amf                                                         ! I2 integers
  USE OMSAO_variables_module,  ONLY:                            &
       voc_amf_filenames, OMBRO_amf_filename, have_amftable, &
       static_input_fnames, verb_thresh_lev
  USE OMSAO_he5_module
  USE OMSAO_errstat_module

  IMPLICIT NONE


  ! -------------------------------------------
  ! Parameters and Variables relevant for OMBRO
  ! -------------------------------------------
  INTEGER (KIND=i4), PARAMETER, PRIVATE :: n_amftab_sza_max = 30, n_amftab_vza_max = 30
  INTEGER (KIND=I4),                                                 PRIVATE :: n_amftab_sza, n_amftab_vza
  REAL    (KIND=r4),                                                 PRIVATE :: amf_tab_wvl, amf_tab_alb
  REAL    (KIND=r8), DIMENSION (n_amftab_sza_max),                   PRIVATE :: amf_table_sza
  REAL    (KIND=r8), DIMENSION (n_amftab_vza_max),                   PRIVATE :: amf_table_vza
  REAL    (KIND=r8), DIMENSION (n_amftab_sza_max, n_amftab_vza_max), PRIVATE :: amf_table_bro

  !INTEGER (KIND=i4)                                  :: n_amftab_coef
  !REAL    (KIND=r8)                                  :: amftab_coef_norm
  !REAL    (KIND=r8), DIMENSION (0:n_amftab_coef_max) :: amftab_coef


  ! --------------------------------------------------------------------
  ! Logical flag for when we are doing Raman clouds, in which case
  ! we compute an additional AMF based on the "forO3" cloud fields.
  ! --------------------------------------------------------------------
  LOGICAL, PRIVATE :: yn_raman_clouds

  ! -------------------------------
  ! TYPE declaration for OMI clouds
  ! -------------------------------
  TYPE, PUBLIC :: OMI_CloudBlock
     REAL (KIND=r8)                                         :: &
          CFRmissing, CFRscale, CFRoffset, CTPmissing, CTPscale, CTPoffset
     REAL (KIND=r8), DIMENSION (nxtrack_max,0:nlines_max-1) :: CFR
     REAL (KIND=r8), DIMENSION (nxtrack_max,0:nlines_max-1) :: CTP
     CHARACTER (LEN=maxchlen)                               :: PGEversion
  END TYPE OMI_CloudBlock


  ! ----------------------------------------------------------------------
  ! AMF tables are computed on the GEOS-Chem lat/lon grid, which has it's
  ! own idiosyncrasies. The grid values correspond to the CENTER points of
  ! the model and have been taken directly from the  the GEOS-Chem manual,
  ! http://www.as.harvard.edu/chemistry/trop/geos/doc/man/index.html
  !
  ! The first grid box (centered around -180deg or similar, i.e., western
  ! longitudes) wraps around to positve (i.e., eastern) longitudes. Hence
  ! any LON >= lonvals(mlon)+dlon/2, where dlon is the distance between
  ! two grid points, is covered by lonvals(1).
  ! ----------------------------------------------------------------------
  
  ! --------------------------------------------------------
  ! All arrays of the HCHO/CHOCHO AMF table are ALLOCATABLE.
  ! --------------------------------------------------------
  REAL    (KIND=r8)                                     :: scale_ncns, scale_ncys, scale_yc
  REAL    (KIND=r4)                                     :: delta_lat, delta_lon, lon_max_geos
  REAL    (KIND=r4), DIMENSION (:),         ALLOCATABLE :: latvals
  REAL    (KIND=r4), DIMENSION (:),         ALLOCATABLE :: lonvals
  REAL    (KIND=r4), DIMENSION (:),         ALLOCATABLE :: szavals
  REAL    (KIND=r4), DIMENSION (:),         ALLOCATABLE :: vzavals
  REAL    (KIND=r4), DIMENSION (:),         ALLOCATABLE :: ctpvals
  INTEGER (KIND=i4), DIMENSION (:,:,:,:),   ALLOCATABLE :: amf_ncns, amf_ncys
  INTEGER (KIND=i4), DIMENSION (:,:,:,:,:), ALLOCATABLE :: amf_ycld

  INTEGER   (KIND=i4),                    PARAMETER, PRIVATE :: nmonths = 12
  CHARACTER (LEN=9), DIMENSION (nmonths), PARAMETER, PRIVATE :: &
  months = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December '    /)


  ! --------------------------------------------
  ! TYPE declaration for ISCCP cloud climatology
  ! --------------------------------------------
  INTEGER (KIND=i4), PARAMETER, PRIVATE :: nlat_isccp=72, nlon_isccp=6596
  TYPE, PUBLIC :: CloudClimatology
     REAL    (KIND=r8)                         :: &
          scale_ctp, scale_cfr, delta_lat, missval_cfr, missval_ctp
     REAL    (KIND=r4), DIMENSION (nlat_isccp) :: latvals, delta_lon
     INTEGER (KIND=i4), DIMENSION (nlat_isccp) :: n_lonvals
     REAL    (KIND=r4), DIMENSION (nlon_isccp) :: lonvals
     REAL    (KIND=r8), DIMENSION (nlon_isccp) :: cfr, ctp
  END TYPE CloudClimatology

  ! ----------------------------------------------
  ! Composite variable for ISCCP Cloud Climatology
  ! ----------------------------------------------
  TYPE (CloudClimatology), PRIVATE :: ISCCP_CloudClim


  CHARACTER (LEN=maxchlen), DIMENSION (n_voc_amf_luns) :: amf_swath_names    = 'undefined'
  INTEGER   (KIND=i4),      DIMENSION (n_voc_amf_luns) :: amf_swath_ids      = -1
  INTEGER   (KIND=i4),      DIMENSION (n_voc_amf_luns) :: amf_swath_file_ids = -1

  ! ----------------------------------------------
  ! Names of various HE5 fields to read from files
  ! ----------------------------------------------
  !(1) OMI Cloud Fields
  ! -------------------
  CHARACTER (LEN=13), PARAMETER, PRIVATE :: omicld_cfrac_field      = 'CloudFraction'
  CHARACTER (LEN=22), PARAMETER, PRIVATE :: omicld_cfracd_field     = 'CloudFractionPrecision'
  CHARACTER (LEN=13), PARAMETER, PRIVATE :: omicld_cpres_field      = 'CloudPressure'
  CHARACTER (LEN=22), PARAMETER, PRIVATE :: omicld_cpresd_field     = 'CloudPressurePrecision'
  ! ----------------------------
  !(2) OMHCHO/OMCHOCHO AMF Table
  ! ---------------------------------------------------------------------------
  ! Note that there is only a "Cloud" AMF but no "Cloud_NoSnow" and "Cloud_Snow"
  ! AMFs. This is because the cloud top is assumed to be a Lambertian surface,
  ! hence whatever is beneath it is of no relevance. Obviously, this is not the
  ! case for a scattering cloud. However, the AMF code used for the computation
  ! of the tables has used a scattering cloud of optical thickness 75 (cloud
  ! top albedo ~0.8) with a black surface underneath to simulate a Lambertian
  ! cloud top.
  ! ---------------------------------------------------------------------------
  ! NOTE: "HCHO_" or "CHOCHO_" must be pre-pended to the field names,
  !       depending on which PGE we are actually running. This is a new
  !       development necessitated by the availability of CHOCHO AMFs.
  !       (tpk, 8 May 2006)
  ! ---------------------------------------------------------------------------
  CHARACTER (LEN= 8), PARAMETER, PRIVATE :: amf_lat_field  = 'AMF_Lats'
  CHARACTER (LEN= 8), PARAMETER, PRIVATE :: amf_lon_field  = 'AMF_Lons'
  CHARACTER (LEN=10), PARAMETER, PRIVATE :: amf_sza_field  = 'AMF_SolZen'
  CHARACTER (LEN=10), PARAMETER, PRIVATE :: amf_vza_field  = 'AMF_LosZen'
  CHARACTER (LEN=10), PARAMETER, PRIVATE :: amf_ctp_field  = 'AMF_CTPres'
  CHARACTER (LEN=18), PARAMETER, PRIVATE :: amf_ncns_field = 'AMF_NoCloud_NoSnow'
  CHARACTER (LEN=16), PARAMETER, PRIVATE :: amf_ncys_field = 'AMF_NoCloud_Snow'
  CHARACTER (LEN= 9), PARAMETER, PRIVATE :: amf_ycld_field = 'AMF_Cloud'

  CHARACTER (LEN= 5), PARAMETER, PRIVATE :: hcho_field   = "HCHO_"
  CHARACTER (LEN= 7), PARAMETER, PRIVATE :: chocho_field = "CHOCHO_"

  ! --------------------------
  !(3) ISCCP Cloud Climatology
  ! --------------------------
  CHARACTER (LEN=10), PARAMETER, PRIVATE :: isccp_lat_field  = 'ISCCP_Lats'
  CHARACTER (LEN=15), PARAMETER, PRIVATE :: isccp_dlon_field = 'ISCCP_DeltaLons'
  CHARACTER (LEN=13), PARAMETER, PRIVATE :: isccp_nlon_field = 'ISCCP_NumLons'
  CHARACTER (LEN=10), PARAMETER, PRIVATE :: isccp_lon_field  = 'ISCCP_Lons'
  CHARACTER (LEN=30), PARAMETER, PRIVATE :: isccp_mcfr_field = 'ISCCP_MonthlyAVG_CloudFraction'
  CHARACTER (LEN=30), PARAMETER, PRIVATE :: isccp_mctp_field = 'ISCCP_MonthlyAVG_CloudPressure'


  ! ------------------------------------------------------------------------
  ! Special parameters for HCHO AMF: Limits for cloud fraction and cloud top
  ! pressure beyond which we flag the HCHO column as "suspect".
  ! ------------------------------------------------------------------------
  !REAL (KIND=r8), PARAMETER :: HCHO_max_cfr = 0.9_r8, HCHO_max_ctp = 900.0_r8

  ! ---------------------------
  ! 32bit/64bit C_LONG integers
  ! ---------------------------
  INTEGER (KIND=C_LONG), PARAMETER, PRIVATE :: zerocl = 0, onecl = 1

CONTAINS

  SUBROUTINE amf_calculation (              &
       pge_idx, nt, nx, lat, lon, sza, vza, &
       snow, glint, xtrange, yn_szoom,      &
       saocol, saodco, saoamf,              &
       errstat                                )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx, pge_idx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon, sza, vza
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: snow, glint
    LOGICAL,           DIMENSION (     0:nt-1), INTENT (IN) :: yn_szoom
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! -----------------------------
    ! Output and modified variables
    ! -----------------------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: saocol, saodco, saoamf
    INTEGER (KIND=i4),                          INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                          :: locerrstat
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1) :: amfdiag
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: amfgeo
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: l2cfr, l2ctp


    locerrstat = pge_errstat_ok

    ! -------------------------
    ! Compute the geometric AMF
    ! -------------------------
    CALL compute_geometric_amf ( nt, nx, sza, vza, xtrange, amfgeo, amfdiag )

    ! -------------------------------------------------------
    ! Initialize molecular AMF with geometric AMF. Subsequent
    ! subroutines will replace any entries where the true
    ! molecular AMF can be computed.
    ! -------------------------------------------------------
    saoamf = amfgeo

    ! -----------------------------
    ! Compute molecule-specific AMF
    ! -----------------------------
    SELECT CASE ( pge_idx )
    CASE ( pge_bro_idx ) 
       CALL ombro_amf_calculation ( &
            nt, nx, sza, vza, snow, glint, xtrange, saoamf, amfdiag, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       WHERE ( saoamf > 0.0_r8 .AND. saocol > r8_missval .AND. saodco > r8_missval ) 
          saocol = saocol / saoamf
          saodco = saodco / saoamf
       END WHERE
    CASE ( pge_hcho_idx ) 
       CALL voc_amf_calculation (                             &
            pge_hcho_idx,                                     &
            nt, nx, lat, lon, sza, vza, snow, glint, xtrange, &
            yn_szoom, saoamf, amfdiag, l2cfr, l2ctp, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       WHERE ( saoamf > 0.0_r8 .AND. saocol > r8_missval .AND. saodco > r8_missval ) 
          saocol = saocol / saoamf
          saodco = saodco / saoamf
       END WHERE
    CASE ( pge_gly_idx ) 
       CALL voc_amf_calculation (                             &
            pge_gly_idx,                                      &
            nt, nx, lat, lon, sza, vza, snow, glint, xtrange, &
            yn_szoom, saoamf, amfdiag, l2cfr, l2ctp, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       WHERE ( saoamf > 0.0_r8 .AND. saocol > r8_missval .AND. saodco > r8_missval ) 
          saocol = saocol / saoamf
          saodco = saodco / saoamf
       END WHERE
    CASE DEFAULT
       ! Nothing to be done here
    END SELECT

    ! -----------------------------------------------
    ! Write AMFs, AMF diagnosting, and AMF-adjusted
    ! columns and column uncertainties to output file
    ! -----------------------------------------------
    CALL he5_write_amf ( &
         pge_idx, nx, nt, saocol, saodco, saoamf, amfgeo, amfdiag, &
         l2cfr, l2ctp, errstat )
    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE amf_calculation


  SUBROUTINE voc_amf_calculation (                       &
       pge_idx,                                          &
       nt, nx, lat, lon, sza, vza, snow, glint, xtrange, &
       yn_szoom, saoamf, amfdiag, l2cfr, l2ctp, errstat    )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: pge_idx, nt, nx
    REAL    (KIND=r4), DIMENSION (nx,0:nt-1),   INTENT (IN) :: lat, lon, sza, vza
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: snow, glint
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange
    LOGICAL,           DIMENSION (   0:nt-1),   INTENT (IN) :: yn_szoom

    ! ----------------------------
    ! Output and Modified varables
    ! ----------------------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (OUT)   :: l2cfr, l2ctp
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: saoamf
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: amfdiag

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, current_month
    INTEGER (KIND=i4) :: mlat, mlon, msza, mvza, mctp

    ! ----------------------
    ! Name of the subroutine
    ! ----------------------
    CHARACTER (LEN=22), PARAMETER :: modulename = 'voc_amf_calculation'


    ! -------------------------
    ! Initialize error variable
    ! -------------------------
    locerrstat = pge_errstat_ok

    ! ------------------------------------------------
    ! GRANULE_MONTH is passed through OMSAO_he5_module
    ! ------------------------------------------------
    current_month = granule_month

    ! -----------------------------
    ! Read the OMI L2 cloud product
    ! -----------------------------
    locerrstat = pge_errstat_ok
    CALL voc_amf_read_omiclouds ( nt, nx, yn_szoom, l2cfr, l2ctp, locerrstat )
    errstat = MAX ( errstat, locerrstat )
    IF ( locerrstat >= pge_errstat_error ) THEN
       l2cfr = r8_missval
       l2ctp = r8_missval
    END IF

    ! ----------------------------
    ! Read ISCCP cloud climatology
    ! ----------------------------
    locerrstat = pge_errstat_ok
    CALL voc_amf_readisccp  ( current_month, locerrstat )
    errstat = MAX ( errstat, locerrstat )
    IF ( locerrstat >= pge_errstat_error ) THEN
       ISCCP_CloudClim%cfr = r8_missval
       ISCCP_CloudClim%ctp = r8_missval
    END IF

    ! --------------------------------------------------------------------
    ! Read AMF table(s). This call will also allocate all table arrays.
    ! Note that all arrays are defined in the header of this module and
    ! are thus known to all subroutines.
    ! --------------------------------------------------------------------
    locerrstat = pge_errstat_ok
    CALL voc_amftable_read ( &
         pge_idx, current_month, mlat, mlon, msza, mvza, mctp, locerrstat )

    ! ---------------------------------------------------------------------
    ! Here we adjust any missing values of the OMI clouds with ISCCP
    ! climatology. We also updateAMFDIAG, which tells us exactly how
    ! the AMF was calculated, i.e., whether anything had to be
    ! initialized from climatology.
    ! ---------------------------------------------------------------------
    CALL voc_amf_diagnostic (                                    &
         nt, nx, lat, lon, sza, snow, glint, xtrange,            &
         MINVAL(ctpvals), MAXVAL(ctpvals), l2cfr, l2ctp, amfdiag  )

    ! ------------------------------------------------------------------------
    ! The AMF calucation routine. We have done all the reading and have set 
    ! the diagnostic flag. The rest is the (tedious) computation of the AMF.
    ! ------------------------------------------------------------------------
    locerrstat = pge_errstat_ok
    CALL voc_amf_compute (                                   &
         nt, nx, lat, lon, sza, vza, snow, xtrange, amfdiag, &
         MINVAL(ctpvals), MAXVAL(ctpvals), l2cfr, l2ctp,     &
         mlat, mlon, msza, mvza, mctp, saoamf, locerrstat      )
    errstat = MAX ( errstat, locerrstat )

    ! -------------------------------
    ! De-allocate AMF table arrays
    ! -------------------------------
    CALL voc_amftable_allocate ( "d", mlat, mlon, msza, mvza, mctp, locerrstat )

    RETURN
  END SUBROUTINE voc_amf_calculation

  SUBROUTINE voc_amf_diagnostic ( &
       nt, nx, lat, lon, sza, snow, glint, xtrange, ctpmin, ctpmax, l2cfr, l2ctp, amfdiag )
    
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    REAL    (KIND=r4),                          INTENT (IN) :: ctpmin, ctpmax
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon, sza
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: snow, glint
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! ----------------
    ! Modified variabe
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: amfdiag
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: l2cfr, l2ctp

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: j1, j2, ix, it, ilat, ilon, spix, epix
    REAL    (KIND=r8) :: latdp, londp


    ! -------------------------------------------------------------------
    ! AMFDIAG has already been set to "geometric" AMF where SZA and VZA
    ! information was available, and "missing"/"out of bounds" otherwise.
    ! Here we need to check for cloud information as well as for snow 
    ! and glint.
    ! -------------------------------------------------------------------


    ! ----------------------------------------------------------------
    ! Checking is done within a loop over NT to assure that we have
    ! "missing" values in all the right places. A single comprehensive
    ! WHERE statement over "1:nx,0:nt-1" would be more efficient but
    ! would also overwrite missing values with real diagnostic flags.
    ! ----------------------------------------------------------------

    DO it = 0, nt-1
       spix = xtrange(it,1) ; epix = xtrange(it,2)

       ! ----------------------------------------------------
       ! Set AMFDIAG to 0 whereever we have "good" clouds and
       ! an SZA that is within the table boundaries.
       ! ----------------------------------------------------
       WHERE (                                           &
            l2cfr(spix:epix,it) >= 0.0_r8          .AND. &
            l2ctp(spix:epix,it) >= 0.0_r8          .AND. &
            sza  (spix:epix,it) >= MINVAL(szavals) .AND. &
            sza  (spix:epix,it) <= MAXVAL(szavals) )
          amfdiag(spix:epix,it) = 0_i2
       END WHERE

       ! -----------------------
       ! Start with the ice flag
       ! -----------------------
       WHERE (                                       &
            amfdiag     (spix:epix,it) >= 0_i2 .AND. &
            snow(spix:epix,it) >= 0_i2         )
          amfdiag(spix:epix,it) = snow(spix:epix,it)
       END WHERE

       !! ----------------------------------------------------
       !! Adjust L2CTP to fall within the AMF table boundaries
       !! ----------------------------------------------------
       !WHERE (                                  &
       !     amfdiag(spix:epix,it) >= 0_i2 .AND. &
       !     l2ctp  (spix:epix,it) > ctpmax         )
       !   l2ctp(spix:epix,it) = REAL ( ctpmax, KIND=r8 )
       !END WHERE
       !WHERE (                                  &
       !     amfdiag(spix:epix,it) >= 0_i2 .AND. &
       !     l2ctp  (spix:epix,it) < ctpmin         )
       !   l2ctp(spix:epix,it) = REAL ( ctpmin, KIND=r8 )
       !END WHERE


       ! ---------------------------------------------------------------------
       ! Where the OMI cloud information is missing, select ISCCP climatology,
       ! but skip pixels where the SZA values are out of the AMF table bounds.
       ! ---------------------------------------------------------------------
       IF ( ( ANY( l2cfr(spix:epix,it) < 0.0_r8 ) ) .OR. &
            ( ANY( l2ctp(spix:epix,it) < 0.0_r8 ) )       ) THEN
          DO ix = spix, epix

             IF ( (l2cfr(ix,it) >= 0.0_r8   .AND. l2ctp(ix,it) >= 0.0_r8   ) .OR. &
                  (sza(ix,it) <  MINVAL(szavals) .OR. sza(ix,it) > MAXVAL(szavals) )  ) CYCLE

             latdp = REAL ( lat(ix,it), KIND=r8 ) ; londp = REAL ( lon(ix,it), KIND=r8 ) ;
             ilat = MAXVAL(MINLOC( ABS(ISCCP_CloudClim%latvals-latdp) ))
             j1   = SUM(ISCCP_CloudClim%n_lonvals(1:ilat-1)) + 1
             j2   = ISCCP_CloudClim%n_lonvals(ilat) + j1
             ilon = MAXVAL(MINLOC( ABS(ISCCP_CloudClim%lonvals(j1:j2)-londp) ))

             IF ( l2ctp(ix,it) < 0.0_r8 .AND. ISCCP_CloudClim%ctp(ilon) >= 0.0_r8 ) THEN
                amfdiag(ix,it) = amfdiag(ix,it) + omi_ctp_addmiss
                l2ctp  (ix,it) = ISCCP_CloudClim%ctp(ilon)
             END IF
             IF ( l2cfr(ix,it) < 0.0_r8 .AND. ISCCP_CloudClim%cfr(ilon) >= 0.0_r8 ) THEN
                amfdiag(ix,it) = amfdiag(ix,it) + omi_cfr_addmiss
                l2cfr  (ix,it) = ISCCP_CloudClim%cfr(ilon)
             END IF

          END DO
       END IF

       ! ---------------------------
       ! Check for glint possibility
       ! ---------------------------
       WHERE (                                    &
            amfdiag  (spix:epix,it) >= 0_i2 .AND. &
            glint(spix:epix,it) > 0_i2           )
          amfdiag(spix:epix,it) = amfdiag(spix:epix,it) + omi_glint_add
       END WHERE
 
    END DO

    RETURN
  END SUBROUTINE voc_amf_diagnostic

  SUBROUTINE voc_amf_compute (                                 &
       nt, nx, lat, lon, sza, vza, snow, xtrange, amfdiag,     &
       ctpmin, ctpmax, l2cfr, l2ctp, mlat, mlon, msza, mvza,   &
       mctp, amf, errstat  )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    INTEGER (KIND=i4),                          INTENT (IN) :: mlat, mlon, msza, mvza, mctp
    REAL    (KIND=r4),                          INTENT (IN) :: ctpmin, ctpmax
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon, sza, vza
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfdiag
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: l2cfr, l2ctp
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: snow
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: amf
    INTEGER (KIND=i4),                          INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                             :: j1, j2, j3, k, k1, k2, ix, it, spix, epix
    REAL    (KIND=r8)                             :: szadp, vzadp, l2cfrloc, l2ctploc
    INTEGER (KIND=i4)                             :: locerrstat
    INTEGER (KIND=i4), DIMENSION (1:nx,0:nt-1)    :: idx_lat, idx_lon, idx_vza
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1)    :: sfr
    REAL    (KIND=r8), DIMENSION (msza)           :: szavdp
    REAL    (KIND=r8), DIMENSION (mvza)           :: vzavdp
    REAL    (KIND=r8), DIMENSION (mctp)           :: ctpvdp
    REAL    (KIND=r8), DIMENSION (msza,mvza)      :: amf_tmp_nc
    REAL    (KIND=r8), DIMENSION (msza,mvza,mctp) :: amf_tmp_yc


    locerrstat = pge_errstat_ok

    ! -----------------------------------------------------------------------
    ! Determine the minimum and maximum indices for the AMF table slice to be
    ! read from file.
    ! -----------------------------------------------------------------------
    CALL voc_amf_findpositions (                                    &
         nt, nx, lat, lon, vza, xtrange, amfdiag, mlat, mlon, mvza, &
         idx_lat, idx_lon, idx_vza                                    )

    ! ------------------------------------------------------------------------------
    ! Now determine the ice fraction; we only consider "non error" and "non suspect"
    ! values for snow and ice cover. The rest is assumed snow free.
    ! ------------------------------------------------------------------------------
    sfr = 0.0_r8
    DO it = 0, nt-1
       spix = xtrange(it,1) ; epix = xtrange(it,2)
       WHERE ( &
            snow(spix:epix,it) > NISE_snowfree .AND. &
            snow(spix:epix,it) < NISE_ocean            )
          sfr(spix:epix,it) = REAL(snow(spix:epix,it), KIND=r8) / 100.0_r8
       ENDWHERE
       ! ----------------------------------------------------------------------
       ! Some of the NISE snow values are > 100, but ice cover has to be <= 1.0
       ! Also, glint conditions are treated as 100% ice cover.
       ! ----------------------------------------------------------------------
       WHERE ( sfr(spix:epix,it) > 1.0_r8 )
          sfr(spix:epix,it) = 1.0_r8
       ENDWHERE
    END DO

    ! -------------------------------------------------------------------------------
    ! Convert AMF table boundaries for SZA and CTP to R8 (required for interpolation)
    ! -------------------------------------------------------------------------------
    ctpvdp = REAL ( ctpvals, KIND=r8 )
    szavdp = REAL ( szavals, KIND=r8 )
    vzavdp = REAL ( vzavals, KIND=r8 )

    ! -------------------------------------------------------------------
    ! Go through the OMI granule pixel by pixel and compute and AMF for
    ! all those that have AMFDIAG >= 0. For all others the AMF has either
    ! been set to MISSING or to the geometric AMF.
    ! -------------------------------------------------------------------

    DO it = 0, nt-1
       spix = xtrange(it,1) ; epix = xtrange(it,2)
       DO ix = spix, epix
          IF ( amfdiag(ix,it) < 0_i2 ) CYCLE

          ! --------------------------------------------
          ! The single (i.e., non-interpolating) indices
          ! --------------------------------------------
          j1 = idx_lat(ix,it) ; j2 = idx_lon(ix,it) ; j3 = idx_vza(ix,it)


          ! --------------------------------------------------
          ! Find slice of 5 VZA around the current J3 location
          ! --------------------------------------------------
          IF      ( j3 <= 3 ) THEN
             k1 = 1      ; k2 = 5
          ELSE IF ( j3 >= mvza-2 ) THEN
             k1 = mvza-4 ; k2 = mvza
          ELSE
             k1 = j3-2 ; k2 = j3+2
          END IF

          ! ----------------------------------------------------------------
          ! Variables cannot be outside the AMF table bounds. For SZA this
          ! has been taken into account through a "missing" AMFDIAG, and the
          ! cloud fraction has already been checked. Below we adjust a local
          ! cloud top variable to fall within the table range.
          ! -----------------------------------------------------------------

          ! ----------------------------------------------------
          ! Adjust L2CTP to fall within the AMF table boundaries
          ! ----------------------------------------------------
          l2ctploc = MAX ( MIN ( l2ctp(ix,it), REAL(ctpmax,KIND=r8) ), REAL(ctpmin,KIND=r8) )

          ! ---------------------
          ! Convert OMI SZA to R8
          ! ---------------------
          szadp = REAL ( sza(ix,it), KIND=r8 )
          vzadp = REAL ( vza(ix,it), KIND=r8 )

          ! ---------------------------------------------------
          ! The cloud-free AMF part (weighted by snow fraction)
          ! ---------------------------------------------------
          amf_tmp_nc(1:msza,k1:k2) = &
               (1.0_r8-sfr(ix,it)) * REAL(amf_ncns(j1,j2,1:msza,k1:k2), KIND=r8) * scale_ncns + &
               (       sfr(ix,it)) * REAL(amf_ncys(j1,j2,1:msza,k1:k2), KIND=r8) * scale_ncys

          ! ---------------------------------------------------------------------
          ! Further action depends on cloud fraction: If the pixel is cloud free,
          ! all we need to do is interpolate a two-dimensional array to the OMI
          ! SZA and VZA. Else we compute the cloudy AMF part and interpolate to
          ! SZA, VZA, and L2CTPLOC
          ! ---------------------------------------------------------------------
          IF ( l2cfr(ix,it) <= 0.0_r8 ) THEN
             CALL ezspline_2d_interpolation ( &
                  msza, k2-k1+1, szavdp(1:msza), vzavdp(k1:k2), amf_tmp_nc(1:msza,k1:k2), &
                  1, 1, szadp, vzadp, amf(ix,it), locerrstat )
          ELSE
             amf_tmp_yc(1:msza,k1:k2,1:mctp) = REAL(amf_ycld(j1,j2,1:msza,k1:k2,1:mctp), KIND=r8) * scale_yc
             DO k = 1, mctp
                amf_tmp_yc (1:msza,k1:k2,k) = &
                     (       l2cfr(ix,it)) * amf_tmp_yc(1:msza,k1:k2,k) + &
                     (1.0_r8-l2cfr(ix,it)) * amf_tmp_nc(1:msza,k1:k2  )
             END DO

             CALL ezspline_3d_interpolation ( &
                  msza, k2-k1+1, mctp, szavdp(1:msza), vzavdp(k1:k2), ctpvdp(1:mctp), &
                  amf_tmp_yc(1:msza,k1:k2,1:mctp), &
                  1, 1, 1, szadp, vzadp, l2ctploc, amf(ix,it), locerrstat )
          END IF

       END DO
    END DO

    errstat = MAX ( errstat, locerrstat )


    RETURN
  END SUBROUTINE voc_amf_compute


  SUBROUTINE voc_amf_read_omiclouds ( nt, nx, yn_szoom, l2cfr, l2ctp, errstat )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),           INTENT (IN) :: nt, nx
    LOGICAL, DIMENSION (0:nt-1), INTENT (IN) :: yn_szoom
    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: l2cfr, l2ctp

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)        :: it, nt_loc, nx_loc, locerrstat, swath_id
    REAL    (KIND=r4)        :: scale_cfr, offset_cfr, missval_cfr, scale_ctp, offset_ctp
    INTEGER (KIND=i2)        :: missval_ctp
    CHARACTER (LEN=5)        :: addstr
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1) :: o4ctp
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1) :: cfr, ctp

    ! ----------------------
    ! Name of the subroutine
    ! ----------------------
    CHARACTER (LEN=22), PARAMETER :: modulename = 'voc_amf_read_omiclouds'

    locerrstat = pge_errstat_ok

    ! -------------------------
    ! Attach to OMI Cloud Swath
    ! -------------------------
    CALL he5_init_input_file ( &
         voc_amf_filenames(voc_omicld_idx), amf_swath_names(voc_omicld_idx), &
         swath_id, amf_swath_file_ids(voc_omicld_idx), &
         nt_loc, nx_loc, locerrstat )
    amf_swath_ids(voc_omicld_idx) = swath_id

    ! ------------------------------------------------------------------------
    ! Usually we would compare nx and nt here, but OMIL2 contains a
    ! different value for nt (=1). Thus we skip the test.
    ! ------------------------------------------------------------------------
    IF ( locerrstat /= pge_errstat_ok ) THEN
       CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
            modulename//f_sep//'OMIL2 access failed.', vb_lev_default, errstat )
       RETURN
    END IF

    IF ( INDEX( voc_amf_filenames(voc_omicld_idx), 'CLDRR' ) /= 0 ) THEN
       yn_raman_clouds = .TRUE.
       addstr          = ""
       !addstr          = "forO3" (this is no longer required; tpk, 8 May 2008)
    ELSE
       yn_raman_clouds = .FALSE.
       addstr          = ""
    END IF

    ! -----------------------
    ! This doesn't work. Why?
    ! -----------------------
    !locerrstat = HE5_SWrdattr  ( swath_id, 'PGEVersion', pgeversion )

    ! --------------------------------------------
    ! Read scaling of cloud data fields (working?)
    ! --------------------------------------------
    scale_cfr = 1.0_r4 ; offset_cfr = 0.0_r4 ; missval_cfr = 0.0_r4
    locerrstat = HE5_SWrdlattr ( swath_id, omicld_cfrac_field//TRIM(ADJUSTL(addstr)), 'MissingValue', missval_cfr )
    !locerrstat = HE5_SWrdlattr ( swath_id, omicld_cfrac_field, 'ScaleFactor',  scale_cfr   )
    !locerrstat = HE5_SWrdlattr ( swath_id, omicld_cfrac_field, 'Offset',       offset_cfr  )
    !IF ( locerrstat /= pge_errstat_ok ) &
    !     CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
    !     modulename//f_sep//'OMIL2 CFR attributes access failed.', vb_lev_default, errstat )

    scale_ctp = 1.0_r4 ; offset_ctp = 0.0_r4 ; missval_ctp = 0
    locerrstat = HE5_SWrdlattr ( swath_id, omicld_cpres_field//TRIM(ADJUSTL(addstr)), 'MissingValue', missval_ctp )
    !locerrstat = HE5_SWrdlattr ( swath_id, omicld_cpres_field, 'ScaleFactor',  scale_ctp  )
    !locerrstat = HE5_SWrdlattr ( swath_id, omicld_cpres_field, 'Offset',       offset_ctp )
    !IF ( locerrstat /= pge_errstat_ok ) &
    !     CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
    !     modulename//f_sep//'OMIL2 CTP attributes access failed.', vb_lev_default, errstat )
        
    swath_id = amf_swath_ids(voc_omicld_idx)

    ! -----------------------
    ! Read current data block
    ! -----------------------
    he5_start_2d = (/ 0, 0 /) ; he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nx, nt /)

    ! ----------------------------------------------------------------
    ! Read cloud fraction and cloud top pressure, and check for error.
    ! Eventually we may read the cloud uncertainties also, but for the
    ! first version we stick with just the basic cloud products.
    ! ----------------------------------------------------------------

    ! ---------------------------------------------------------------------------
    ! (1) Cloud Fraction is of type REAL*4 in both Raman and O2-O2 cloud products
    ! ---------------------------------------------------------------------------
    locerrstat = HE5_SWrdfld ( &
         amf_swath_ids(voc_omicld_idx), omicld_cfrac_field//TRIM(ADJUSTL(addstr)),   &
         he5_start_2d, he5_stride_2d, he5_edge_2d, cfr(1:nx,0:nt-1) )
    IF ( locerrstat /= pge_errstat_ok ) &
         CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'OMIL2 CFR access failed.', vb_lev_default, errstat )

    ! ---------------------------------------------------------------
    ! Check for rebinned zoom data swath storage ("1-30" vs. "16-45")
    ! ---------------------------------------------------------------
    DO it = 0, nt-1
       IF ( yn_szoom(it) .AND. &
            ALL ( cfr(gzoom_epix:nx,it) <= missval_cfr ) ) THEN
          cfr(gzoom_spix:gzoom_epix,it) = cfr(1:gzoom_npix,it)
          cfr(1:gzoom_spix-1,       it) = missval_cfr
       END IF
    END DO

    ! -----------------------------------------------------------
    ! Assign the cloud fraction array used in the AMF calculation
    ! -----------------------------------------------------------
    l2cfr = REAL(cfr, KIND=r8)
    WHERE ( cfr > r4_missval )
       l2cfr = l2cfr * scale_cfr + offset_cfr
    END WHERE

    ! ---------------------------------------------------------------------------
    ! (2) Cloud Pressure of type REAL*4 in Raman but INT*2 in O2-O2
    ! ---------------------------------------------------------------------------
    IF ( yn_raman_clouds ) THEN
       locerrstat = HE5_SWrdfld ( &
            amf_swath_ids(voc_omicld_idx), omicld_cpres_field//TRIM(ADJUSTL(addstr)),   &
            he5_start_2d, he5_stride_2d, he5_edge_2d, ctp(1:nx,0:nt-1) )
    ELSE
       locerrstat = HE5_SWrdfld ( &
            amf_swath_ids(voc_omicld_idx), omicld_cpres_field,   &
            he5_start_2d, he5_stride_2d, he5_edge_2d, o4ctp(1:nx,0:nt-1) )

       ! -------------------------------------------
       ! Temporary copy of O4CTP to an R4 type array
       ! -------------------------------------------
       ctp(1:nx,0:nt-1) = REAL ( o4ctp(1:nx,0:nt-1), KIND=r4 )
    END IF

    ! ---------------------------------------------------------------
    ! Check for rebinned zoom data swath storage ("1-30" vs. "16-45")
    ! ---------------------------------------------------------------
    DO it = 0, nt-1
       IF ( yn_szoom(it) .AND. &
            ALL ( ctp(gzoom_epix:nx,it) <= REAL(missval_ctp, KIND=r4) ) ) THEN
          ctp(gzoom_spix:gzoom_epix,it)  = ctp(1:gzoom_npix,it)
          ctp(1:gzoom_spix-1,       it)  = REAL(missval_ctp, KIND=r4)
       END IF
    END DO

    ! ---------------------------------------------------------------
    ! Assign the cloud top pressure array used in the AMF calculation
    ! ---------------------------------------------------------------
    l2ctp  = REAL(ctp, KIND=r8)
    WHERE ( ctp > REAL(missval_ctp, KIND=r4) )
       l2ctp = l2ctp * scale_ctp + offset_ctp
    END WHERE

    IF ( locerrstat /= pge_errstat_ok ) &
         CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'OMIL2 CTP access failed.', vb_lev_default, errstat )

    ! ------------------------------------------------
    ! Force the cloud parameters into physical bounds.
    ! But make sure not to remove MissingValues.
    ! ------------------------------------------------
    WHERE ( l2cfr > REAL(missval_cfr, KIND=r8) .AND. l2cfr < 0.0_r8 )
       l2cfr = 0.0_r8
    ENDWHERE
    WHERE ( l2cfr > 1.0_r8 )
       l2cfr = 1.0_r8
    ENDWHERE
    WHERE ( l2ctp > REAL(missval_ctp, KIND=r8) .AND. l2ctp < 0.0_r8 )
       l2ctp = 0.0_r8
    ENDWHERE

    ! ------------------------------------------------------
    ! Replace cloud missing values by SAO PGE missing values
    ! ------------------------------------------------------
    WHERE ( l2cfr <= REAL(missval_cfr, KIND=r8) )
       l2cfr = r8_missval
    ENDWHERE
    WHERE ( l2ctp <= REAL(missval_ctp, KIND=r8) )
       l2ctp = r8_missval
    ENDWHERE
    
    RETURN
  END SUBROUTINE voc_amf_read_omiclouds


  SUBROUTINE voc_amf_readisccp ( i0, errstat )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: i0

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)                         :: locerrstat, swath_id, swath_file_id, swlen, he5stat
    CHARACTER (LEN=maxchlen)                    :: swath_file, swath_name
    INTEGER   (KIND=i4), DIMENSION (nlon_isccp) :: tmparr

    ! -----------------------
    ! Name of this subroutine
    ! -----------------------
    CHARACTER (LEN=18), PARAMETER :: modulename = 'voc_amf_read_isccp'

    locerrstat = pge_errstat_ok

    swath_file = TRIM(ADJUSTL(voc_amf_filenames(voc_isccp_idx)))
    
    ! -----------------------------------------------------------
    ! Open HE5 output file and check SWATH_FILE_ID ( -1 if error)
    ! -----------------------------------------------------------
    swath_file_id = HE5_SWopen ( swath_file, he5f_acc_rdonly )
    amf_swath_file_ids(voc_isccp_idx) = swath_file_id
    IF ( swath_file_id == he5_stat_fail ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_error, OMSAO_E_HE5SWOPEN, modulename, vb_lev_default, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       RETURN
    END IF
    
    ! ---------------------------------------------
    ! Check for existing HE5 swath and attach to it
    ! ---------------------------------------------
    locerrstat  = HE5_SWinqswath  ( TRIM(ADJUSTL(swath_file)), swath_name, swlen )
    amf_swath_names(voc_isccp_idx) = TRIM(ADJUSTL(swath_name))
    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
    amf_swath_ids(voc_isccp_idx) = swath_id
    IF ( swath_id == he5_stat_fail ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_error, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       RETURN
    END IF

    ! ----------------------------
    ! Read ISCCP Cloud Climatoloty
    ! ----------------------------
    ! * Latitudes
    he5_start_1d = 0 ; he5_stride_1d = 1 ; he5_edge_1d = nlat_isccp
    he5stat = HE5_SWrdfld ( swath_id, isccp_lat_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, ISCCP_CloudClim%latvals(1:nlat_isccp) )
    ! * Number of longitudes per latitude
    he5stat = HE5_SWrdfld ( swath_id, isccp_nlon_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, ISCCP_CloudClim%n_lonvals(1:nlat_isccp) )
    ! * Delta-Longitudes    
    he5stat = HE5_SWrdfld ( swath_id, isccp_dlon_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, ISCCP_CloudClim%delta_lon(1:nlat_isccp) )
    ! * Longitudes
    he5_start_1d = 0 ; he5_stride_1d = 1 ; he5_edge_1d = nlon_isccp
    he5stat = HE5_SWrdfld ( swath_id, isccp_lon_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, ISCCP_CloudClim%lonvals(1:nlon_isccp) )

    ! * Cloud fraction
    he5_start_2d  = (/ i0-1, 0 /) ; he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ 1, nlon_isccp /)
    he5stat = HE5_SWrdfld ( swath_id, isccp_mcfr_field,                 &
         he5_start_2d, he5_stride_2d, he5_edge_2d, tmparr(1:nlon_isccp) )
    ISCCP_CloudClim%cfr(1:nlon_isccp) = REAL (tmparr(1:nlon_isccp), KIND=r8)
    ! * Cloud top pressure
    he5stat = HE5_SWrdfld ( swath_id, isccp_mctp_field,                 &
         he5_start_2d, he5_stride_2d, he5_edge_2d, tmparr(1:nlon_isccp) )
    ISCCP_CloudClim%ctp(1:nlon_isccp) = REAL (tmparr(1:nlon_isccp), KIND=r8)

    ! --------------------
    ! Read some attributes
    ! --------------------
    he5stat = HE5_SWrdlattr ( swath_id, isccp_lat_field,  "DeltaGrid",    ISCCP_CloudClim%delta_lat   )
    he5stat = HE5_SWrdlattr ( swath_id, isccp_mcfr_field, "ScaleFactor",  ISCCP_CloudClim%scale_cfr   )
    he5stat = HE5_SWrdlattr ( swath_id, isccp_mctp_field, "ScaleFactor",  ISCCP_CloudClim%scale_ctp   )
    he5stat = HE5_SWrdlattr ( swath_id, isccp_mcfr_field, "MissingValue", ISCCP_CloudClim%missval_cfr )
    he5stat = HE5_SWrdlattr ( swath_id, isccp_mctp_field, "MissingValue", ISCCP_CloudClim%missval_ctp )

    ! -------------------------------------------------------
    ! Scale the ISCCP cloud values with their scaling factors
    ! -------------------------------------------------------
    WHERE ( ISCCP_CloudClim%cfr /= ISCCP_CloudClim%missval_ctp )
       ISCCP_CloudClim%cfr = ISCCP_CloudClim%cfr * ISCCP_CloudClim%scale_cfr
    END WHERE
    WHERE ( ISCCP_CloudClim%ctp /= ISCCP_CloudClim%missval_ctp )
       ISCCP_CloudClim%ctp = ISCCP_CloudClim%ctp * ISCCP_CloudClim%scale_ctp
    END WHERE

    ! ---------------------
    ! The final error check
    ! ---------------------
    IF ( he5stat /= pge_errstat_ok ) &
         CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'OMHCHO AMF attributes access failed.', vb_lev_default, errstat )

    ! -----------------------------------------------
    ! Detach from HE5 swath and close HE5 output file
    ! -----------------------------------------------
    he5stat = HE5_SWdetach ( swath_id )
    he5stat = HE5_SWclose  ( swath_file_id )

    RETURN
  END SUBROUTINE voc_amf_readisccp

  SUBROUTINE voc_amftable_read ( pge_idx, ismonth, mlat, mlon, msza, mvza, mctp, errstat )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, ismonth

    ! ----------------
    ! Output variables
    ! ----------------
    ! (all others are defined in the header of this module)
    INTEGER (KIND=i4), INTENT (OUT) :: mlat, mlon, msza, mvza, mctp

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)         :: nswath, locerrstat, swath_id, swath_file_id, swlen, he5stat
    INTEGER   (KIND=C_LONG)     :: nswathcl, swlencl
    CHARACTER (LEN=   maxchlen) :: swath_file, locswathname, voc_str
    CHARACTER (LEN=10*maxchlen) :: swath_name

    ! -------------------------
    ! KIND=4 or KIND=8 integers
    ! -------------------------
    INTEGER (KIND=C_LONG) :: mlatcl, mloncl, mszacl, mvzacl, mctpcl

    ! -----------------------
    ! Name of this subroutine
    ! -----------------------
    CHARACTER (LEN=17), PARAMETER :: modulename = 'voc_amftable_read'


    locerrstat = pge_errstat_ok

    swath_file = TRIM(ADJUSTL(voc_amf_filenames(voc_amf_idx)))
    

    ! ---------------------------------------------------------
    ! Check for the type of PGE we are running; this determines
    ! which string we need to prepend to the field names.
    ! ---------------------------------------------------------
    SELECT CASE ( pge_idx )
    CASE ( pge_hcho_idx )
       voc_str = hcho_field
    CASE ( pge_gly_idx )
       voc_str = chocho_field
    CASE DEFAULT
       voc_str = ""
    END SELECT

    ! -----------------------------------------------------------
    ! Open HE5 output file and check SWATH_FILE_ID ( -1 if error)
    ! -----------------------------------------------------------
    swath_file_id = HE5_SWopen ( swath_file, he5f_acc_rdonly )
    amf_swath_file_ids(voc_amf_idx) = swath_file_id
    IF ( swath_file_id == he5_stat_fail ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_error, OMSAO_E_HE5SWOPEN, modulename, vb_lev_default, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       RETURN
    END IF
    
    ! ---------------------------------------------
    ! Check for existing HE5 swath and attach to it
    ! ---------------------------------------------
    nswathcl = HE5_SWinqswath  ( TRIM(ADJUSTL(swath_file)), swath_name, swlen )
    nswath   = INT ( nswathcl, KIND=i4 )
    ! -------------------------------------------------------------------------
    ! If there is only one swath in the file, we can attach to it right away.
    ! But if there are more (NSWATH > 1), then we must find the swath that
    ! corresponds to the current month we are processing.
    !
    ! Note that we are not peforming a "month" check for a single-swath file.
    ! There is always the possibility that we want to experiment with different
    ! swaths, and forcing name compliance here potentially makes this kind of
    ! experimentation more difficult. Not likely to happen, but who knows.
    ! -------------------------------------------------------------------------
    IF ( nswath > 1 ) THEN
       CALL voc_extract_swathname ( &
            nswath, TRIM(ADJUSTL(swath_name)), TRIM(ADJUSTL(months(ismonth))), &
            locswathname )

       ! ---------------------------------------------------------------------------
       ! Check if we found the correct swath name. If not, report an error and exit.
       ! ---------------------------------------------------------------------------
       IF ( INDEX (TRIM(ADJUSTL(locswathname)),TRIM(ADJUSTL(months(ismonth)))) == 0 ) THEN
          CALL error_check ( &
               0, 1, pge_errstat_error, OMSAO_E_HE5SWLOCATE, modulename, &
               vb_lev_default, locerrstat )
          errstat = MAX ( errstat, locerrstat )
          RETURN
       END IF
    ELSE
       locswathname = TRIM(ADJUSTL(swath_name))
    END IF

    ! -----------------------------------------------------
    ! Proceed with attaching to swath and reading the table
    ! -----------------------------------------------------

    amf_swath_names(voc_amf_idx) = TRIM(ADJUSTL(locswathname))
    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(locswathname)) )
    amf_swath_ids(voc_amf_idx) = swath_id
    IF ( swath_id == he5_stat_fail ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_error, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       RETURN
    END IF

    ! -------------------------------
    ! Read dimensions of AMF table(s)
    ! -------------------------------
    locerrstat = pge_errstat_ok
    CALL voc_amftable_getdim ( swath_id, mlat, mlon, msza, mvza, mctp, locerrstat )
    IF ( ismonth < 1 .OR. ismonth > nmonths .OR. locerrstat /= pge_errstat_ok ) THEN
       errstat = MAX ( errstat, locerrstat ) 
       RETURN
    END IF

    ! ----------------------------------------------------------------------------
    ! Create KIND=4/KIND=8 variables. We have to use those dimensions a few times
    ! in this subroutine, so it saves some typing if we do the conversion once and
    ! save them in new variables.
    ! ----------------------------------------------------------------------------
    mlatcl = INT ( mlat, KIND=C_LONG )
    mloncl = INT ( mlon, KIND=C_LONG )
    mszacl = INT ( msza, KIND=C_LONG )
    mvzacl = INT ( mvza, KIND=C_LONG )
    mctpcl = INT ( mctp, KIND=C_LONG )

    ! -------------------------------
    ! Allocate AMF table arrays
    ! -------------------------------
    locerrstat = pge_errstat_ok
    CALL voc_amftable_allocate ( "a", mlat, mlon, msza, mvza, mctp, locerrstat )
    IF ( locerrstat /= pge_errstat_ok ) THEN
       errstat = MAX ( errstat, locerrstat ) 
       CALL voc_amftable_allocate ( "d", mlat, mlon, msza, mvza, mctp, locerrstat )
       RETURN
    END IF


    ! -------------------------------
    ! Read dimension-defining arrays
    ! -------------------------------
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = mlatcl
    he5stat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(voc_str))//amf_lat_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, latvals(1:mlat) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = mloncl
    he5stat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(voc_str))//amf_lon_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, lonvals(1:mlon) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = mszacl
    he5stat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(voc_str))//amf_sza_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, szavals(1:msza) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = mvzacl
    he5stat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(voc_str))//amf_vza_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, vzavals(1:mvza) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = mctpcl
    he5stat = HE5_SWrdfld ( swath_id, TRIM(ADJUSTL(voc_str))//amf_ctp_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, ctpvals(1:mctp) )
    IF ( he5stat /= pge_errstat_ok ) &
         CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'AMF arrays access failed.', vb_lev_default, errstat )
    ! --------------------
    ! Read some attributes
    ! --------------------
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(voc_str))//amf_lat_field,  "DeltaGrid",   delta_lat )
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(voc_str))//amf_lon_field,  "DeltaGrid",   delta_lon )
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(voc_str))//amf_ncns_field, "ScaleFactor", scale_ncns )
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(voc_str))//amf_ncys_field, "ScaleFactor", scale_ncys )
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(voc_str))//amf_ycld_field, "ScaleFactor", scale_yc )
    IF ( he5stat /= pge_errstat_ok ) &
         CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'AMF attributes access failed.', vb_lev_default, errstat )


    ! --------------------------------------------------------------------
    ! Compute the maximum longitude (note the GEOS-Chem grid idiosyncrasy)
    ! --------------------------------------------------------------------
    lon_max_geos = lonvals(mlon) + delta_lon / 2.0_r8


    ! ---------------
    ! Read the tables
    ! ---------------
    he5_start_4d  = (/   zerocl,  zerocl,  zerocl,  zerocl /)
    he5_stride_4d = (/   onecl,   onecl,   onecl,   onecl /)
    he5_edge_4d   = (/ mlatcl, mloncl, mszacl, mvzacl /)
    he5stat = HE5_SWrdfld (                                &
         swath_id, TRIM(ADJUSTL(voc_str))//amf_ncns_field, &
         he5_start_4d, he5_stride_4d, he5_edge_4d,         &
         amf_ncns(1:mlat,1:mlon,1:msza,1:mvza)               )
    he5stat = HE5_SWrdfld (                                &
         swath_id, TRIM(ADJUSTL(voc_str))//amf_ncys_field, &
         he5_start_4d, he5_stride_4d, he5_edge_4d,         &
         amf_ncys(1:mlat,1:mlon,1:msza,1:mvza)               )
    IF ( he5stat /= pge_errstat_ok ) &
         CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'AMF NoCloud access failed.', vb_lev_default, errstat )

    he5_start_5d  = (/  zerocl, zerocl, zerocl, zerocl, zerocl /)
    he5_stride_5d = (/  onecl,  onecl,  onecl,  onecl,  onecl /)
    he5_edge_5d   = (/ mlatcl, mloncl, mszacl, mvzacl, mctpcl /)

    he5stat = HE5_SWrdfld (                                &
         swath_id, TRIM(ADJUSTL(voc_str))//amf_ycld_field, &
         he5_start_5d, he5_stride_5d, he5_edge_5d,         &
         amf_ycld(1:mlat,1:mlon,1:msza,1:mvza,1:mctp)        )
    IF ( he5stat /= pge_errstat_ok ) &
         CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'AMF Cloud access failed.', vb_lev_default, errstat )

    RETURN
  END SUBROUTINE voc_amftable_read


  SUBROUTINE voc_extract_swathname ( nswath, multi_swath, swathstr, single_swath )

    ! ---------------------------------------------------------------------
    ! Extracts SINGLE_SWATH from MULTI_SWATH, based on presence of SWATHSTR
    ! ---------------------------------------------------------------------
 
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER   (KIND=i4), INTENT (IN) :: nswath
    CHARACTER (LEN=*),   INTENT (IN) :: multi_swath, swathstr

    ! ----------------
    ! Output variables
    ! ----------------
    CHARACTER (LEN=*), INTENT (OUT) :: single_swath

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4), DIMENSION(0:nswath) :: swsep
    INTEGER   (KIND=i4)                      :: mslen, k, j1, j2, nsep
    CHARACTER (LEN=LEN(multi_swath))         :: tmpstr


    ! --------------------------
    ! Initialize output variable
    ! --------------------------
    single_swath = ''

    ! ---------------------------------
    ! Find length of MULTI_SWATH string
    ! ---------------------------------
    tmpstr = TRIM(ADJUSTL(multi_swath))
    mslen  = LEN_TRIM(ADJUSTL(tmpstr))

    ! --------------------------------------
    ! First find the number of "," in TMPSTR
    ! --------------------------------------
    nsep = 0 ; swsep(0:nswath) = 0
    DO k = 1, mslen
       IF ( tmpstr(k:k) == ',' ) THEN
          nsep = nsep + 1
          swsep(nsep) = k
       END IF
    END DO
    IF ( nsep == nswath-1 ) THEN
       nsep = nsep + 1 ; swsep(nsep) = mslen + 1
    END IF


    ! ----------------------------------------------------------------------
    ! Hangle along the positions of separators (commas, ",") between the
    ! concatinated Swath Name entries. The first Swath Name to contain 
    ! SWATHSTR is taken as the match - not a perfect rationale but simple
    ! enough if we have set up the AMF table file correctly.
    ! ----------------------------------------------------------------------    
    getswath: DO k = 1, nswath
       j1 = swsep(k-1)+1  ;  j2 = swsep(k)-1
       IF ( INDEX ( tmpstr(j1:j2), TRIM(ADJUSTL(swathstr)) ) > 0 ) THEN
          single_swath = TRIM(ADJUSTL(tmpstr(j1:j2)))
          EXIT getswath
       END IF
    END DO getswath

    RETURN
  END SUBROUTINE voc_extract_swathname


  SUBROUTINE voc_amftable_getdim ( &
       swath_id, mlat, mlon, msza, mvza, mctp, errstat )

    ! ---------------------------------
    ! Return dimensions of AMF table(s)
    ! ---------------------------------
 
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: swath_id

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    INTEGER (KIND=i4), INTENT (OUT)   :: mlat, mlon, msza, mvza, mctp

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4), PARAMETER             :: maxdim = 100
    INTEGER   (KIND=i4)                        :: ndim, nsep
    INTEGER   (KIND=C_LONG)                       :: ndimcl
    INTEGER   (KIND=i4)                        :: i, j, swlen, iend, istart
    INTEGER   (KIND=i4),  DIMENSION(0:maxdim)  :: dim_array, dim_seps
    INTEGER   (KIND=C_LONG), DIMENSION(0:maxdim)  :: dim_arraycl
    CHARACTER (LEN=10*maxdim)                  :: dim_chars


    ! ---------------------------
    ! Initialize output variables
    ! ---------------------------
    mlat = -1 ; mlon = -1 ; msza = -1 ; mvza = -1 ; mctp = -1

    ! ------------------------------
    ! Inquire about swath dimensions
    ! ------------------------------
    ndimcl = HE5_SWinqdims  ( swath_id, dim_chars, dim_arraycl(0:maxdim) )
    ndim   = INT ( ndimcl, KIND=i4 )
    IF ( ndim <= 0 ) THEN
       errstat = MAX ( errstat, pge_errstat_error )
       RETURN
    END IF
    dim_array(0:maxdim) = INT ( dim_arraycl(0:maxdim), KIND=i4 )

    dim_chars =     TRIM(ADJUSTL(dim_chars))
    swlen     = LEN_TRIM(ADJUSTL(dim_chars))

    ! ----------------------------------------------------------------------
    ! Find the positions of separators (commas, ",") between the dimensions.
    ! Add a "pseudo separator" at the end to fully automate the consecutive
    ! check for nTimes and nXtrack.
    ! ----------------------------------------------------------------------
    nsep = 0 ; dim_seps(0:ndim) = 0
    getseps: DO i = 1, swlen 
       IF ( dim_chars(i:i) == ',' ) THEN
          nsep           = nsep + 1
          dim_seps(nsep) = i
       END IF
    END DO getseps
    nsep = nsep + 1 ; dim_seps(nsep) = swlen+1

    ! --------------------------------------------------------------------
    ! Hangle along the NSEP indices until we have found the two dimensions
    ! we are interested in.
    ! --------------------------------------------------------------------
    getdims:DO j = 0, nsep-1
       istart = dim_seps(j)+1 ; iend = dim_seps(j+1)-1

       SELECT CASE ( dim_chars(istart:iend) )
       CASE ( "n_AMF_lat" )
          mlat = dim_array(j)
       CASE ( "n_AMF_lon" )
          mlon = dim_array(j)
       CASE ( "n_AMF_sza" )
          msza = dim_array(j)
       CASE ( "n_AMF_vza")
          mvza = dim_array(j)
       CASE ( "n_AMF_ctp" )
          mctp = dim_array(j)
       CASE DEFAULT
          ! Whatever. Nothing to be done here.
       END SELECT

    END DO getdims

    RETURN
  END SUBROUTINE voc_amftable_getdim


  SUBROUTINE voc_amftable_allocate ( ad, mlat, mlon, msza, mvza, mctp, errstat )

    USE OMSAO_casestring_module, ONLY: lower_case
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=1), INTENT (IN) :: ad
    INTEGER (KIND=i4), INTENT (IN) :: mlat, mlon, msza, mvza, mctp

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4) :: estat
    CHARACTER (LEN=1)   :: adlow


    estat = pge_errstat_ok

    ! ------------------------------------------
    ! Make sure AD ("a" or "d") is in lower case
    ! ------------------------------------------
    adlow = lower_case ( ad )

    SELECT CASE ( adlow )
    CASE ('a')
       ALLOCATE (latvals (mlat),                     STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (lonvals (mlon),                     STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (szavals (msza),                     STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (vzavals (mvza),                     STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (ctpvals (mctp),                     STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (amf_ncns(mlat,mlon,msza,mvza),      STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (amf_ncys(mlat,mlon,msza,mvza),      STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (amf_ycld(mlat,mlon,msza,mvza,mctp), STAT=estat ) ; errstat = MAX ( errstat, estat )
    CASE ('d')
       IF ( ALLOCATED ( latvals  ) )  DEALLOCATE ( latvals  )
       IF ( ALLOCATED ( lonvals  ) )  DEALLOCATE ( lonvals  )
       IF ( ALLOCATED ( szavals  ) )  DEALLOCATE ( szavals  )
       IF ( ALLOCATED ( vzavals  ) )  DEALLOCATE ( vzavals  )
       IF ( ALLOCATED ( ctpvals  ) )  DEALLOCATE ( ctpvals  )
       IF ( ALLOCATED ( amf_ncns ) )  DEALLOCATE ( amf_ncns )
       IF ( ALLOCATED ( amf_ncys ) )  DEALLOCATE ( amf_ncys )
       IF ( ALLOCATED ( amf_ycld ) )  DEALLOCATE ( amf_ycld )
    CASE DEFAULT
       ! Whatever. Nothing to be done here
    END SELECT

    RETURN
  END SUBROUTINE voc_amftable_allocate


  SUBROUTINE voc_amf_findpositions (                              &
       nt, nx, lat, lon, vza, xtrange, amfdiag, mlat, mlon, mvza, &
       idx_lat, idx_lon, idx_vza                                    )
    
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                         INTENT (IN) :: nt, nx, mlat, mlon, mvza
    INTEGER (KIND=i2), DIMENSION(1:nx,0:nt-1), INTENT (IN) :: amfdiag
    REAL    (KIND=r4), DIMENSION(1:nx,0:nt-1), INTENT (IN) :: lat, lon, vza
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2), INTENT (IN) :: xtrange

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), DIMENSION(1:nx,0:nt-1), INTENT (OUT) :: idx_lat, idx_lon, idx_vza

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)  :: ix, it, spix, epix

    idx_lat = -1 ; idx_lon = -1 ; idx_vza = -1

    DO it = 0, nt-1
       spix = xtrange(it,1) ; epix = xtrange(it,2)

       ! ----------------------------------------------------------
       ! Set all longitude indices for LON > lon_max_geos to 1
       ! ----------------------------------------------------------
       WHERE (                                         &
            amfdiag(spix:epix,it) >= 0_2         .AND. &
            lon    (spix:epix,it) > lon_max_geos         )
          idx_lon(spix:epix,it) = 1
       END WHERE

       DO ix = spix, epix
          ! -----------------------------------------------
          ! Skip "bad" AMF diagnostics (missing, geometric)
          ! -----------------------------------------------
          IF ( amfdiag(ix,it) < 0_i2 ) CYCLE

          ! --------------------------------------------------------
          ! LAT/LON/VZA is "drop in the box", i.e., nearest neighbor
          ! --------------------------------------------------------
          idx_vza(ix,it) = MINVAL( MINLOC( ABS(vzavals(1:mvza)-vza(ix,it)) ) )
          idx_lat(ix,it) = MINVAL( MINLOC( ABS(latvals(1:mlat)-lat(ix,it)) ) )
          IF ( idx_lon(ix,it) /= 1 ) &
               idx_lon(ix,it) = MINVAL( MINLOC( ABS(lonvals(1:mlon)-lon(ix,it)) ) ) 

       END DO
    END DO

    RETURN
  END SUBROUTINE voc_amf_findpositions


  SUBROUTINE compute_geometric_amf ( nt, nx, sza, vza, xtrange, amfgeo, amfdiag )

    USE OMSAO_parameters_module, ONLY: deg2rad, rad2deg
    USE OMSAO_variables_module,  ONLY: szamax

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                         INTENT (IN) :: nx, nt
    REAL    (KIND=r4), DIMENSION (nx,0:nt-1),  INTENT (IN) :: sza, vza
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2), INTENT (IN) :: xtrange

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: amfgeo
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: amfdiag

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: it, spix, epix


    ! ---------------------------------------------
    ! Compute geometric AMF and set diagnostic flag
    ! ---------------------------------------------
    amfgeo  = r8_missval
    amfdiag = i2_missval

    ! ----------------------------------------------------------------
    ! Checking is done within a loop over NT to assure that we have
    ! "missing" values in all the right places. A single comprehensive
    ! WHERE statement over "1:nx,0:nt-1" would be more efficient but
    ! would also overwrite missing values with OMI_OOBVIEW_AMF.
    ! ----------------------------------------------------------------
    DO it = 0, nt-1
       spix = xtrange(it,1) ; epix = xtrange(it,2)
       WHERE ( &
            sza(spix:epix,it) /= r4_missval .AND. &
            sza(spix:epix,it) >=     0.0_r4 .AND. &
            sza(spix:epix,it) <     90.0_r4 .AND. &
            vza(spix:epix,it) /= r4_missval .AND. &
            vza(spix:epix,it) >=     0.0_r4 .AND. &
            vza(spix:epix,it) <     90.0_r4         )
          amfgeo(spix:epix,it) = &
               1.0_r8 / COS ( REAL(sza(spix:epix,it),KIND=r8)*deg2rad ) + &
               1.0_r8 / COS ( REAL(vza(spix:epix,it),KIND=r8)*deg2rad )
          amfdiag(spix:epix,it) = omi_geo_amf
       ELSEWHERE
          amfdiag(spix:epix,it) = omi_oobview_amf
       ENDWHERE
    END DO

    RETURN
  END SUBROUTINE compute_geometric_amf


  SUBROUTINE ombro_amf_read_table ( errstat )


    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat


    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=20), PARAMETER :: modulename = 'OMBRO_AMF_read_table'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)      :: locerrstat, version, file_read_stat
    INTEGER   (KIND=i4)      :: i, funit
    CHARACTER (LEN=maxchlen) :: lastline

    ! ---------------------------------
    ! External OMI and Toolkit routines
    ! ---------------------------------
    INTEGER (KIND=i4) :: pgs_smf_teststatuslevel, pgs_io_gen_openf, pgs_io_gen_closef

    ! ----------------------------------------------------------------------------------
    ! First check whether HAVE_AMFTABLE=.TRUE. from READ_PCF_FILE. HAVE_AMFTABLE=.FALSE.
    ! means that we failed to get the AMF table file from the PCF and thus can't 
    ! determine the AMFs. HAVE_AMFTABLE will also be set .FALSE. if any error occurs in
    ! this routine, in which case no slant-to-vertical conversion of the column amount 
    ! will be performed.
    ! ----------------------------------------------------------------------------------
    IF ( .NOT. have_amftable ) RETURN

    locerrstat = pge_errstat_ok

    ! -----------------------
    ! Open reference spectrum
    ! -----------------------
    version = 1
    locerrstat = PGS_IO_GEN_OPENF ( OMBRO_amf_lun, PGSd_IO_Gen_RSeqFrm, 0, funit, version )
    locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)

    CALL error_check ( &
         locerrstat, PGS_SMF_MASK_LEV_S, pge_errstat_warning, OMSAO_W_OPEN_AMFTABLE_FILE, &
         modulename//f_sep//TRIM(ADJUSTL(OMBRO_amf_filename)), vb_lev_default, &
         errstat )
    IF ( errstat /= pge_errstat_ok ) THEN
       have_amftable = .FALSE. ; locerrstat = PGS_IO_GEN_CLOSEF ( funit ) ; RETURN
    END IF

    ! --------------------------------------------------
    ! Re-initialize the local error value, which is most
    ! likely set to 512, the TOOLKIT success value.
    ! --------------------------------------------------
    locerrstat = pge_errstat_ok
    
    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    n_amftab_sza = 0      ; n_amftab_vza   = 0
    amf_tab_wvl   = 0.0_r4 ; amf_tab_alb   = 0.0_r4
    amf_table_sza = 0.0_r8 ; amf_table_vza = 0.0_r8
    amf_table_bro = 0.0_r8

    ! --------------------------------------
    ! Skip comments header to start of table
    ! --------------------------------------
    CALL skip_to_filemark ( funit, lm_start_of_table, lastline, file_read_stat )

    CALL error_check ( &
         file_read_stat, file_read_ok, pge_errstat_warning, OMSAO_W_READ_AMFTABLE_FILE, &
         modulename//f_sep//TRIM(ADJUSTL(OMBRO_amf_filename)), vb_lev_default, &
         locerrstat )
    IF ( locerrstat /= pge_errstat_ok ) THEN
       errstat = MAX ( errstat, locerrstat )
       have_amftable = .FALSE. ;   locerrstat = PGS_IO_GEN_CLOSEF ( funit ) ; RETURN
    END IF

    ! ---------------------------------------------------------
    ! Read table dimensions, table wavelength, and table albedo
    ! ---------------------------------------------------------
    READ (UNIT=funit, FMT=*, IOSTAT=locerrstat) n_amftab_sza, n_amftab_vza, amf_tab_wvl, amf_tab_alb

    ! ------------------------------------------
    ! Check whether we can proceed to read table
    ! ------------------------------------------
    IF ( (locerrstat   /= file_read_ok)     .OR. &   ! READ successful
         (n_amftab_sza <= 0)                .OR. &   ! Number of SZAs in table >  0
         (n_amftab_sza >  n_amftab_sza_max) .OR. &   ! Table dimensions <= maximum table dim
         (n_amftab_vza <= 0)                .OR. &   ! Number of VZAs in table >  0
         (n_amftab_vza >  n_amftab_vza_max)      &   ! Table dimensions <= maximum table dim
         ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_warning, OMSAO_W_READ_AMFTABLE_FILE, &
            modulename//f_sep//TRIM(ADJUSTL(OMBRO_amf_filename)), vb_lev_default, &
            locerrstat )
       errstat = MAX ( errstat, locerrstat )
       have_amftable = .FALSE. ; locerrstat = PGS_IO_GEN_CLOSEF ( funit ) ; RETURN
    ELSE
       READ (UNIT=funit, FMT=*, IOSTAT=locerrstat) amf_table_vza(1:n_amftab_vza)       
       DO i = 1, n_amftab_sza
          READ (UNIT=funit, FMT=*, IOSTAT=locerrstat) amf_table_sza(i), amf_table_bro(i,1:n_amftab_vza)
          IF ( i /= n_amftab_sza .AND. locerrstat /= 0 ) THEN
             CALL error_check ( &
                  0, 1, pge_errstat_warning, OMSAO_W_READ_AMFTABLE_FILE, &
                  modulename//f_sep//TRIM(ADJUSTL(OMBRO_amf_filename)), vb_lev_default, &
                  locerrstat )
             errstat = MAX ( errstat, locerrstat )
             have_amftable = .FALSE.  ;  RETURN
          END IF
       END DO
    END IF

    ! --------------------------------------
    ! Close BrO AMF file
    ! --------------------------------------
    locerrstat = PGS_IO_GEN_CLOSEF ( funit )
    locerrstat = PGS_SMF_TESTSTATUSLEVEL ( errstat )

    ! ------------------------------------
    ! Report successful reading of spectra
    ! ------------------------------------
    CALL error_check ( &
         0, 1, pge_errstat_ok, OMSAO_S_READ_AMFTABLE_FILE, &
         modulename//f_sep//TRIM(ADJUSTL(OMBRO_amf_filename)), vb_lev_stmdebug, &
         locerrstat )

    RETURN
  END SUBROUTINE OMBRO_AMF_read_table

  SUBROUTINE ombro_amf_calculation ( &
       nt, nx, sza, vza, snow, glint, xtrange, amf, amfdiag, errstat )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    REAL (KIND=r4),    DIMENSION (1:nx,0:nt-1), INTENT (IN) :: sza, vza
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: snow, glint
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! -------------------------
    ! Output/modified variables
    ! -------------------------
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: amfdiag
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: amf
    INTEGER (KIND=i4),                          INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat

    ! ----------------------
    ! Name of the subroutine
    ! ----------------------
    CHARACTER (LEN=24), PARAMETER :: modulename = 'ombro_amf_calculation'

    ! -------------------------
    ! Initialize error variable
    ! -------------------------
    locerrstat = pge_errstat_ok

    ! --------------------------------------------------------------------
    ! Read the table - it's short enough.
    ! --------------------------------------------------------------------
    CALL ombro_amf_read_table ( locerrstat )
    errstat = MAX ( errstat, locerrstat )
    IF ( .NOT. have_amftable ) RETURN

    ! ---------------------------------------------------------------------
    ! Here we adjust any missing values of the OMI clouds with climatology.
    ! We also update AMFDIAG, which tells us exactly how the AMF was
    ! calculated, i.e., whether anything had to be initialized from
    ! climatology.
    ! ---------------------------------------------------------------------
    CALL ombro_amf_diagnostic ( nt, nx, sza, vza, snow, glint, xtrange, amfdiag )

    ! --------------------------
    ! The AMF calucation routine
    ! --------------------------
    CALL ombro_amf_compute ( nt, nx, amfdiag, sza, vza, xtrange, amf )


    RETURN
  END SUBROUTINE ombro_amf_calculation

  SUBROUTINE ombro_amf_diagnostic ( nt, nx, sza, vza, snow, glint, xtrange, amfdiag )
 
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: sza, vza
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: snow, glint
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: amfdiag

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: it, spix, epix


    ! -------------------------------------------------------------------
    ! AMFDIAG has already been set to "geometric" AMF where SZA and VZA
    ! information was available, and "missing" otherwise. Here we need
    ! to check for SZA and VZA being within/out range. For all good angle
    ! values (i.e., those within the table bounds) we check for snow and
    ! glint. Otherwise we set values to "out of bounds" or "missing".
    ! -------------------------------------------------------------------

    DO it = 0, nt-1

       spix = xtrange(it,1) ; epix = xtrange(it,2)

       ! ----------------
       ! Missing SZA, VZA
       ! ----------------
       WHERE (                                   &
            sza(spix:epix,it) <= r8_missval .OR. &
            vza(spix:epix,it) <= r8_missval        )
          amfdiag(spix:epix,it) = i2_missval
       END WHERE

       ! ----------------------------------------
       ! Out-of-Bound SZA, VZA (but not missing!)
       ! ----------------------------------------
       WHERE ( &
            ( sza(spix:epix,it) < MINVAL(amf_table_sza(1:n_amftab_sza)) ) .OR. &
            ( sza(spix:epix,it) > MAXVAL(amf_table_sza(1:n_amftab_sza)) ) .OR. &
            ( vza(spix:epix,it) < MINVAL(amf_table_vza(1:n_amftab_vza)) ) .OR. &
            ( vza(spix:epix,it) > MAXVAL(amf_table_vza(1:n_amftab_vza)) )      )
          amfdiag(spix:epix,it) = omi_oobview_amf
       END WHERE

       ! ------------------------------------------------------
       ! And AMFDIAG values > OOB must be good and are set to 0
       ! ------------------------------------------------------
       WHERE ( amfdiag(spix:epix,it) > omi_oobview_amf )
          amfdiag(spix:epix,it) = 0_i2
       END WHERE
        
       ! ---------------------------------------------------------------
       ! Now we can check the "good" values for ice and glint conditions
       ! ---------------------------------------------------------------
       WHERE (                                       &
            amfdiag     (spix:epix,it) == 0_i2 .AND. &
            snow(spix:epix,it) >= 0_i2         )
          amfdiag(spix:epix,it) = snow(spix:epix,it)
       END WHERE
       WHERE (                                    &
            amfdiag  (spix:epix,it) >= 0_i2 .AND. &
            glint(spix:epix,it) >  0_2           )
          amfdiag(spix:epix,it) = amfdiag(spix:epix,it) + omi_glint_add
       END WHERE
 
    END DO

    RETURN
  END SUBROUTINE ombro_amf_diagnostic


  SUBROUTINE ombro_amf_compute ( nt, nx, amfdiag, sza, vza, xtrange, amf )
    
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfdiag
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: sza, vza
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: amf

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: ix, it, ks, ks1, ks2, kv, kv1, kv2, locerrstat
    INTEGER (KIND=i4) :: spix, epix
    REAL    (KIND=r8) :: szadp, vzadp


    locerrstat = pge_errstat_ok

    ! -----------------------
    ! Compute AMFs one by one
    ! -----------------------

    DO it = 0, nt-1

       spix = xtrange(it,1) ; epix = xtrange(it,2)
       DO ix = spix, epix

          ! --------------------------
          ! Skip and "bad" AMF entries
          ! --------------------------
          IF ( amfdiag(ix,it) < 0_i2 ) CYCLE

          ! -----------------------------
          ! Convert OMI SZA and VZA to R8
          ! -----------------------------
          szadp = REAL ( sza(ix,it), KIND=r8 )
          vzadp = REAL ( vza(ix,it), KIND=r8 )

          ! -----------------------------------------------------
          ! Find postions of current SZA and VZA and table arrays
          ! -----------------------------------------------------
          ks = MINVAL( MINLOC( ABS(amf_table_sza(1:n_amftab_sza)-szadp) ) )
          kv = MINVAL( MINLOC( ABS(amf_table_vza(1:n_amftab_vza)-vzadp) ) )

          ! ----------------------------------------------------------------
          ! Find slice of 5 SZAs/VZAs around the current KS and KV locations
          ! ----------------------------------------------------------------
          ks1 = MIN ( MAX(1, ks-3), n_amftab_sza-6 ) ; ks2 = MIN(ks1+6, n_amftab_sza)
          kv1 = MIN ( MAX(1, kv-3), n_amftab_vza-6 ) ; kv2 = MIN(kv1+6, n_amftab_vza)


          CALL ezspline_2d_interpolation ( &
               ks2-ks1+1, kv2-kv1+1, amf_table_sza(ks1:ks2), amf_table_vza(kv1:kv2), &
               amf_table_bro(ks1:ks2,kv1:kv2), &
               1, 1, szadp, vzadp, amf(ix,it), locerrstat )


       END DO
    END DO

    RETURN
  END SUBROUTINE ombro_amf_compute

END MODULE OMSAO_AMF_module
