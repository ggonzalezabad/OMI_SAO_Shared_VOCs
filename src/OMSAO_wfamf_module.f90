MODULE OMSAO_wfamf_module

  ! ====================================================================
  ! This module defines variables associated with the wavelength depende
  ! nt AMF calculations and contains necessary subroutines to read files
  ! and calculate them
  ! ====================================================================
  USE OMSAO_precision_module, ONLY: i4, r8, C_LONG, r4
  USE OMSAO_parameters_module
  USE OMSAO_errstat_module
  USE OMSAO_he5_module
  USE OMSAO_indices_module, ONLY: n_voc_amf_luns
  USE HDF5

  IMPLICIT NONE

  ! ====================================================================
  ! Wavelength dependent AMF factor specific variables
  ! ====================================================================
  LOGICAL                                      :: yn_amf_wfmod
  INTEGER                                      :: amf_wfmod_idx
  REAL(KIND=r8)                                :: amf_alb_lnd, amf_alb_sno, amf_wvl, &
       amf_wvl2, amf_alb_cld, amf_max_sza

  ! ---------------------------------
  ! GMAO GEOS-5 hybrid grid Ap and Bp
  ! ---------------------------------
  REAL(KIND=r4), DIMENSION(48), PARAMETER :: Ap=(/0.000000E+00, 4.804826E-02, 6.593752E+00, 1.313480E+01, &
       1.961311E+01, 2.609201E+01, 3.257081E+01, 3.898201E+01, &
       4.533901E+01, 5.169611E+01, 5.805321E+01, 6.436264E+01, &
       7.062198E+01, 7.883422E+01, 8.909992E+01, 9.936521E+01, &
       1.091817E+02, 1.189586E+02, 1.286959E+02, 1.429100E+02, &
       1.562600E+02, 1.696090E+02, 1.816190E+02, 1.930970E+02, &
       2.032590E+02, 2.121500E+02, 2.187760E+02, 2.238980E+02, &
       2.243630E+02, 2.168650E+02, 2.011920E+02, 1.769300E+02, &
       1.503930E+02, 1.278370E+02, 1.086630E+02, 9.236572E+01, &
       7.851231E+01, 5.638791E+01, 4.017541E+01, 2.836781E+01, &
       1.979160E+01, 9.292942E+00, 4.076571E+00, 1.650790E+00, &
       6.167791E-01, 2.113490E-01, 6.600001E-02, 1.000000E-02/)
  REAL(KIND=r4), DIMENSION(48), PARAMETER :: Bp=(/1.000000E+00, 9.849520E-01, 9.634060E-01, 9.418650E-01, &
       9.203870E-01, 8.989080E-01, 8.774290E-01, 8.560180E-01, &
       8.346609E-01, 8.133039E-01, 7.919469E-01, 7.706375E-01, &
       7.493782E-01, 7.211660E-01, 6.858999E-01, 6.506349E-01, &
       6.158184E-01, 5.810415E-01, 5.463042E-01, 4.945902E-01, &
       4.437402E-01, 3.928911E-01, 3.433811E-01, 2.944031E-01, &
       2.467411E-01, 2.003501E-01, 1.562241E-01, 1.136021E-01, &
       6.372006E-02, 2.801004E-02, 6.960025E-03, 8.175413E-09, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00/)
  INTEGER(KIND=i4), PARAMETER :: ngeos5 = 48

  ! ---------------------------------------
  ! Data obtained from the climatology file
  ! ---------------------------------------
  REAL(KIND=r4), DIMENSION(:),     ALLOCATABLE :: latvals, lonvals
  REAL(KIND=r4), DIMENSION(:,:),   ALLOCATABLE :: Psurface
  REAL(KIND=r4), DIMENSION(:,:,:), ALLOCATABLE :: Temperature, Gas_profiles, H2O_profiles

  ! ---------------------------------------
  ! Data obtained from Vlidort lookup table
  ! Look up table variables
  ! ---------------------------------------
  ! To hold data:
  REAL(KIND=r4),    DIMENSION(:),           ALLOCATABLE :: lut_alb, lut_clp, lut_sza, &
       lut_vza, lut_wavelength, lut_srf
  CHARACTER(LEN=5), DIMENSION(:),           ALLOCATABLE :: lut_toms
  REAL(KIND=r4),    DIMENSION(:,:),         ALLOCATABLE :: lut_alt_lay, lut_alt_lev, lut_pre_lay, lut_pre_lev, &
       lut_pre_lev_cld, lut_pre_lay_cld
  REAL(KIND=r4),    DIMENSION(:,:,:),       ALLOCATABLE :: lut_air, lut_ozo, lut_tem
  REAL(KIND=r4),    DIMENSION(:,:),         ALLOCATABLE :: lut_Sb_clr, lut_Sb_cld
  REAL(KIND=r4),    DIMENSION(:,:,:,:),     ALLOCATABLE :: lut_I0_clr, lut_I1_clr, lut_I2_clr, lut_Ir_clr, &
       lut_I0_cld, lut_I1_cld, lut_I2_cld, lut_Ir_cld
  REAL(KIND=r4),    DIMENSION(:,:,:,:,:),   ALLOCATABLE :: lut_dI0_cld, lut_dI1_cld, lut_dI2_cld
  REAL(KIND=r4),    DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: lut_dI0_clr, lut_dI1_clr, lut_dI2_clr
  ! To read the data:
  INTEGER(HSIZE_T), DIMENSION(1) :: alb_dim, alb_maxdim, &
       clp_dim, clp_maxdim, &
       srf_dim, srf_maxdim, &
       sza_dim, sza_maxdim, &
       toz_dim, toz_maxdim, &
       vza_dim, vza_maxdim, &
       wav_dim, wav_maxdim
  INTEGER(HSIZE_T), DIMENSION(2) :: Sb_clr_dim, Sb_clr_maxdim, &
       alt_lev_dim, alt_lev_maxdim, &
       alt_lay_dim, alt_lay_maxdim, &
       alt_lev_cld_dim, alt_lev_cld_maxdim, &
       alt_lay_cld_dim, alt_lay_cld_maxdim, Sb_cld_dim, Sb_cld_maxdim
  INTEGER(HSIZE_T), DIMENSION(3) :: air_dim, air_maxdim, tem_dim, tem_maxdim
  INTEGER(HSIZE_T), DIMENSION(4) :: I0_clr_dim,  I0_clr_maxdim, I0_cld_dim, I0_cld_maxdim
  INTEGER(HSIZE_T), DIMENSION(5) :: dI0_cld_dim, dI0_cld_maxdim
  INTEGER(HSIZE_T), DIMENSION(6) :: dI0_clr_dim,  dI0_clr_maxdim

  ! ---------
  ! PCF stuff
  ! ---------
  INTEGER(KIND=i4), PARAMETER :: wfamf_table_lun = 700250
  INTEGER(KIND=i4), PARAMETER :: climatology_lun = 700270
  CHARACTER(LEN=maxchlen)     :: OMSAO_wfamf_table_filename
  CHARACTER(LEN=maxchlen)     :: OMSAO_climatology_filename

  ! -----------------------------
  ! Dimensions of the climatology
  ! -----------------------------
  INTEGER (KIND=i4) :: Cmlat, Cmlon, CmETA, CmEp1

  ! -----------------------
  ! To find the right swath
  ! -----------------------
  INTEGER   (KIND=i4),                    PARAMETER, PRIVATE :: nmonths = 12
  CHARACTER (LEN=9), DIMENSION (nmonths), PARAMETER, PRIVATE :: &
  months = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December '    /)

  ! ---------------------------
  ! 32bit/64bit C_LONG integers
  ! ---------------------------
  INTEGER (KIND=C_LONG), PARAMETER, PRIVATE :: zerocl = 0, onecl = 1

  ! -----------
  ! ISCCP stuff
  ! -----------
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

  ! --------------------------
  !(3) ISCCP Cloud Climatology
  ! --------------------------
  CHARACTER (LEN=10), PARAMETER, PRIVATE :: isccp_lat_field  = 'ISCCP_Lats'
  CHARACTER (LEN=15), PARAMETER, PRIVATE :: isccp_dlon_field = 'ISCCP_DeltaLons'
  CHARACTER (LEN=13), PARAMETER, PRIVATE :: isccp_nlon_field = 'ISCCP_NumLons'
  CHARACTER (LEN=10), PARAMETER, PRIVATE :: isccp_lon_field  = 'ISCCP_Lons'
  CHARACTER (LEN=30), PARAMETER, PRIVATE :: isccp_mcfr_field = 'ISCCP_MonthlyAVG_CloudFraction'
  CHARACTER (LEN=30), PARAMETER, PRIVATE :: isccp_mctp_field = 'ISCCP_MonthlyAVG_CloudPressure'


CONTAINS

  SUBROUTINE amf_calculation (                   &
       pge_idx, nt, nx, lat, lon, sza, vza,      &
       saa, vaa, snow, glint, xtrange, yn_szoom, &
       saocol, saodco, saoamf, terrain_height,   &
       yn_write, errstat                                )

    ! =================================================================
    ! This subroutine computes the AMF factor using the following eleme
    ! nts:
    !     - Kleipool OMLER database
    !     - GEOS Chem climatology
    !     - VLIDORT calculated scattering weights
    ! =================================================================
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx, pge_idx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon, sza, vza, saa, vaa, &
         terrain_height
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: snow, glint
    LOGICAL,           DIMENSION (     0:nt-1), INTENT (IN) :: yn_szoom
    LOGICAL                                                 :: yn_write
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! -----------------------------
    ! Output and modified variables
    ! -----------------------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: saocol, saodco, saoamf
    INTEGER (KIND=i4),                          INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                :: locerrstat, itt, spixx, epixx, current_month
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1)       :: amfdiag
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1)       :: amfgeo
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1)       :: l2cfr, l2ctp
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1)       :: albedo, cli_psurface
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA) :: climatology, cli_temperature, cli_heights
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA) :: scattw !, akernels


    locerrstat  = pge_errstat_ok

    ! ------------------------------------
    ! Initialize variables that are output
    ! ------------------------------------
    albedo       = r8_missval
    climatology  = r8_missval
    cli_heights  = r8_missval
    cli_psurface = r8_missval
    scattw       = r8_missval
    saoamf       = r8_missval
    amfgeo       = r8_missval
    amfdiag      = i2_missval

    ! -----------------------------------------
    ! If amf_wvl < 0.0 then the slant column is
    ! reported and AMFs equal to 1 with scattw, 
    ! akernels and climatology set to missval
    ! -----------------------------------------
    IF (amf_wvl .LT. 0.0) THEN

       DO itt = 0, nt-1
          spixx = xtrange(itt,1) ; epixx = xtrange(itt,2)
          saoamf(spixx:epixx,itt) = 1.0_r8
       END DO

    ELSE

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
       
       ! ----------------------------------------------------
       ! Read OMLER albedo database stored in variable albedo
       ! ----------------------------------------------------
       CALL omi_omler_albedo ( lat, lon, albedo, nt, nx, xtrange, locerrstat)
       
       ! ---------------------------------------
       ! Write the albedo to the output file he5
       ! ---------------------------------------
       IF (yn_write) CALL write_albedo_he5 ( albedo, nt, nx, locerrstat)
       
       ! -----------------------------
       ! Read the OMI L2 cloud product
       ! -----------------------------
       locerrstat = pge_errstat_ok
       CALL amf_read_omiclouds ( nt, nx, yn_szoom, l2cfr, l2ctp, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       IF ( locerrstat >= pge_errstat_error ) THEN
          l2cfr = r8_missval
          l2ctp = r8_missval
       END IF
       
       ! ----------------------------
       ! Read ISCCP cloud climatology
       ! ----------------------------
       locerrstat = pge_errstat_ok
       current_month = granule_month
       CALL voc_amf_readisccp  ( current_month, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       IF ( locerrstat >= pge_errstat_error ) THEN
          ISCCP_CloudClim%cfr = r8_missval
          ISCCP_CloudClim%ctp = r8_missval
       END IF

       ! ---------------------------------------------------------------------
       ! The climatology has already been read, inside omi_pge_fitting_process
       ! Now it is only needed to interpolate to the pixels of the granule.
       ! It was read there to obtain the dimensions of the number of levels.
       ! ---------------------------------------------------------------------
       CALL omi_climatology (climatology, cli_heights, cli_psurface, cli_temperature, terrain_height, &
            lat, lon, nt, nx, xtrange, locerrstat)
       
       ! -------------------------------------
       ! Write the climatology to the he5 file
       ! -------------------------------------
       IF (yn_write) CALL write_climatology_he5 (climatology, cli_psurface, nt, nx, CmETA, locerrstat)

       ! ------------------------------------------------------------------
       ! Read VLIDORT look up table. Variables are declared at module level
       ! ------------------------------------------------------------------
       CALL read_lookup_table (locerrstat)

       ! ----------------------------------------------------------------------
       ! amfdiag is used to keep track of the pixels were enough information is
       ! available to carry on the AMFs calculation.
       ! ----------------------------------------------------------------------
       CALL amf_diagnostic (                                                   &
            nt, nx, lat, lon, sza, vza, snow, glint, xtrange,                  &
            MINVAL(lut_clp), MAXVAL(lut_clp), l2cfr, l2ctp,&
            amfdiag  )

       ! --------------------------------------------------------
       ! Compute Scattering weights in the look up table grid but
       ! with the correct albedo. amfdiag is used to skip pixel
       ! ---------------------------------------------------------
       CALL compute_scatt ( nt, nx, albedo, lat, sza, vza, saa, vaa, l2ctp, l2cfr, terrain_height, amfgeo, amfdiag, &
            scattw)

       ! ----------------------------
       ! Deallocate Vlidort variables
       ! ----------------------------
       CALL vlidort_allocate ("d", INT(toz_dim(1)),INT(srf_dim(1)),INT(clp_dim(1)),INT(alb_dim(1)), &
            INT(sza_dim(1)),INT(vza_dim(1)),INT(wav_dim(1)), &
            INT(alt_lay_dim(2)),INT(alt_lev_dim(2)),errstat)

       ! -----------------------------------------------------------------
       ! Work out the AMF using the scattering weights and the climatology
       ! Work out Averaging Kernels
       ! -----------------------------------------------------------------
       CALL compute_amf ( nt, nx, CmETA, climatology, cli_heights, cli_temperature, cli_psurface, &
            scattw, saoamf, amfdiag, locerrstat)

       ! -----------------------------------------------------------------
       ! Write out scattering weights, altitude grid and averaging kernels
       ! -----------------------------------------------------------------
       IF (yn_write) CALL write_scatt_he5 (scattw, nt, nx, CmETA, locerrstat)
       
    END IF
    
    ! --------------------------
    ! Apply the air mass factors
    ! --------------------------
    WHERE ( saoamf > 0.0_r8 .AND. saocol > r8_missval .AND. saodco > r8_missval ) 
       saocol = saocol / saoamf
       saodco = saodco / saoamf
    END WHERE
    
    ! -----------------------------------------------
    ! Write AMFs, AMF diagnosting, and AMF-adjusted
    ! columns and column uncertainties to output file
    ! -----------------------------------------------
    IF (yn_write) CALL he5_amf_write ( pge_idx, nx, nt, saocol, saodco, saoamf, &
         amfgeo, amfdiag, l2cfr, l2ctp, locerrstat )
    
    errstat = MAX ( errstat, locerrstat )
    
  END SUBROUTINE amf_calculation
  
  SUBROUTINE omi_climatology (climatology, local_heights, local_psurf, local_temperature, terrain_height, &
                              lat, lon, nt, nx, xtrange, locerrstat)

    ! =========================================
    ! Extract Gas climatology to granule pixels
    ! No interpolation or something like that,
    ! Just pick the closest model grid
    ! =========================================
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon, terrain_height
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4),                                INTENT (INOUT) :: locerrstat
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1, CmETA), INTENT (INOUT) :: climatology, local_temperature, local_heights
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1),        INTENT (INOUT) :: local_psurf  

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: itimes, ixtrack, spix, epix, idx_lat, idx_lon, ilevel, n, n1
    REAL    (KIND=r8)                      :: rho, lhgt, aircolumn
    REAL    (KIND=r8), DIMENSION (0:CmETA) :: lpre, level_press
    REAL    (KIND=r8), DIMENSION (1:CmETA) :: ltmp, layer_press
    REAL    (KIND=r8) :: thish2omxr, Mwet, Rwet, detlnp

    ! -----------------------
    ! Some physical constants
    ! -----------------------
    REAL (KIND=r8), PARAMETER ::  &
         Mdry = 0.02896,   &
         Mh2o = 0.018,     &
         Rstar = 8.314,    &
         Navogadro = 6.02214e+23, &
         gplanet = 9.806    
 
    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=30), PARAMETER :: modulename = 'omi_climatology' 

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    locerrstat = pge_errstat_ok
    
    ! ------------------------------------------------------------
    ! Find the Climatology corresponding to each lat and lon pixel
    ! ------------------------------------------------------------
    DO itimes = 0, nt-1

       spix = xtrange(itimes,1); epix = xtrange(itimes,2)
       DO ixtrack = spix, epix

          IF (lon(ixtrack,itimes) .LT. -180.0_r4 .OR. &
              lat(ixtrack,itimes) .LT.  -90.0_r4 .OR. &
              lon(ixtrack,itimes) .GT.  180.0_r4 .OR. &
              lat(ixtrack,itimes) .GT.  90.0_r4) CYCLE

          ! ----------------------------------------------------
          ! Just selecting the closest location to lat lon pixel
          ! ----------------------------------------------------
          idx_lat = MINVAL(MINLOC(ABS(latvals(1:Cmlat) - lat(ixtrack,itimes) )))
          idx_lon = MINVAL(MINLOC(ABS(lonvals(1:Cmlon) - lon(ixtrack,itimes) )))

          local_temperature(ixtrack,itimes,1:CmETA) = Temperature(idx_lon,idx_lat,1:CmETA)
          ltmp(1:CmETA)                             = Temperature(idx_lon,idx_lat,1:CmETA)

          ! ----------------------------------------------
          ! Convert pixel terrain height to pressure using
          ! Xiong suggested to use pressure altitude:
          !  Z = -16 alog10 (P / Po) Z in km and P in hPa.
          ! ----------------------------------------------
          local_psurf(ixtrack,itimes) = 1013.0_r8 * &
               (10.0_r8 ** (REAL(terrain_height(ixtrack,itimes),KIND=r8) / 1000.0_r8 / (-16.0_r8)))

          ! lpre(0:CmETA) is pressure at layer boundaries, same unit as
          ! local_surf, convert from hPa to Pa, it is a CmETA+1 = CmEp1 array
          ! Ap(1:CmEp1) & Bp(1:CmEp1) are coeff at layer boundary for pressure
          ! calculation
          DO ilevel = 0, CmETA
             lpre(ilevel) = (Ap(ilevel+1) + local_psurf(ixtrack,itimes) * Bp(ilevel+1))*1.E2 ! Pa
          END DO


          ! layer_press(1:CmETA) are pressure at layer centers in Pa
          ! local_heights are pressures in hPa at layer centers (I should get rid of this variable)
          DO ilevel = 1, CmETA ! for layer center
             layer_press(ilevel) = (lpre(ilevel-1) + lpre(ilevel))/2.0 ! Pa
             local_heights(ixtrack,itimes,ilevel)   = layer_press(ilevel)/1.E2 ! hPa
          END DO

          ! -------------------
          ! Develop air density
          ! -------------------
          
          DO n = 1, CmETA

             rho       = 0.0_r8
             aircolumn = 0.0_r8
             lhgt      = 0.0_r8
             
             ! Convert input water vapor mixing ratio from PPB to unitless
             thish2omxr = H2O_profiles(idx_lon,idx_lat,n) / 1.0E9

             ! Calculate mean molecular weight of wet air 
             Mwet = (1.0_r8 - thish2omxr)*Mdry + thish2omxr*Mh2o
             
             ! Calculate gas constant for wet air
             Rwet = Rstar / Mwet
             
             ! Calculate layer thickness using the following
             ! dz = -(R*T/g) * dlnP
             n1                   = n - 1
             detlnp = log(lpre(n1)) - log(lpre(n))
             lhgt = Rwet * ltmp(n) * detlnp / gplanet ! meter
             
             rho                = layer_press(n) / ltmp(n) / Rstar 
             aircolumn          = rho*lhgt*Navogadro*1.0E-4 ! # air/cm^2
             
             ! -------------------------------------------------------------
             climatology(ixtrack,itimes,n) = aircolumn * Gas_profiles(idx_lon,idx_lat,n) /1.0E9 ! [GAS]/cm^2
             
          END DO
          
          !  Set non-physical entries to zero.
          WHERE ( climatology(ixtrack,itimes,1:CmETA) < 0.0_r8 )
             climatology(ixtrack,itimes,1:CmETA) = 0.0_r8
          END WHERE
          
       END DO
    END DO

  END SUBROUTINE omi_climatology

  SUBROUTINE omi_read_climatology(errstat)
    
    ! ==========================================================
    ! This subroutine reads in the climatology from GEOS-Chem or
    ! other source. The climatology needs to be produced keeping
    ! the format so the reader will be flexible enough to move
    ! from one climatology to another one
    ! ==========================================================
    
    USE OMSAO_indices_module,   ONLY: sao_molecule_names
    USE OMSAO_variables_module, ONLY: pge_idx, pge_name
    IMPLICIT NONE    

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)         :: nswath, locerrstat, swath_id, swath_file_id, swlen, he5stat, &
                                   ismonth, ndatafields, h2o_cli_idx
    INTEGER   (KIND=i4), DIMENSION(10) :: datafield_rank, datafield_type
    INTEGER   (KIND=C_LONG)     :: nswathcl, Cmlatcl, Cmloncl, CmETAcl, CmEp1cl
    CHARACTER (LEN=   maxchlen) :: swath_file, locswathname, gasdatafieldname, datafield_name
    CHARACTER (LEN=10*maxchlen) :: swath_name

    CHARACTER (LEN= 9), PARAMETER :: cli_lat_field         = 'Latitudes'
    CHARACTER (LEN=10), PARAMETER :: cli_lon_field         = 'Longitudes'
    CHARACTER (LEN=15), PARAMETER :: cli_Psurf_field       = 'SurfacePressure'
    CHARACTER (LEN=18), PARAMETER :: cli_Temperature_field = 'TemperatureProfile'

    REAL      (KIND=r4) :: scale_lat, scale_lon, scale_gas, scale_Psurf, &
                           scale_temperature, scale_H2O

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER (KIND=i4) :: version

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=35), PARAMETER :: modulename = 'omi_read_climatology' 

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    locerrstat = pge_errstat_ok
    
    swath_file  = TRIM(ADJUSTL(OMSAO_climatology_filename))
    ismonth     = granule_month

    ! Index to read H2O water vapor using sao_molecule_names
    h2o_cli_idx = 18

    ! --------------------------------------------------------------
    ! Open he5 OMI climatology and check SWATH_FILE_ID (-1 if error)
    ! --------------------------------------------------------------
    swath_file_id = HE5_SWOPEN (swath_file, he5f_acc_rdonly)
    IF (swath_file_id == he5_stat_fail) THEN
       CALL error_check (0, 1, pge_errstat_error, OMSAO_E_HE5SWOPEN, modulename, &
            vb_lev_default, locerrstat)
       errstat = MAX (errstat, locerrstat)
       RETURN
    END IF
    
    ! -----------------------------------------------------------
    ! Check for existing HE5 swathw and attach to the one we need
    ! -----------------------------------------------------------
    nswathcl = HE5_SWinqswath(TRIM(ADJUSTL(swath_file)), swath_name, swlen )
    nswath   = INT(nswathcl, KIND=i4 )

    ! ----------------------------------------------------------------
    ! If there is only one swath in the file, we can attach to it but
    ! if there are more (NSWATH > 1), then we must find the swath that
    ! corresponds to the current month.
    ! ----------------------------------------------------------------
    IF (nswath > 1) THEN
       CALL extract_swathname(nswath, TRIM(ADJUSTL(swath_name)), &
            TRIM(ADJUSTL(months(granule_month))), locswathname)
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

    ! -----------------------------
    ! Attach to current month swath
    ! -----------------------------
    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(locswathname)) )
    IF ( swath_id == he5_stat_fail ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_error, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, locerrstat )
       errstat = MAX ( errstat, locerrstat )
       RETURN
    END IF

    ! ------------------------------------
    ! Read dimensions of Climatology swath
    ! ------------------------------------
    locerrstat = pge_errstat_ok
    CALL climatology_getdim ( swath_id, Cmlat, Cmlon, CmETA, CmEp1, locerrstat )
    IF ( ismonth < 1 .OR. ismonth > nmonths .OR. locerrstat /= pge_errstat_ok ) THEN
       errstat = MAX ( errstat, locerrstat ) 
       RETURN
    END IF

    ! ----------------------------------------------------------------------------
    ! Create KIND=4/KIND=8 variables. We have to use those dimensions a few times
    ! in this subroutine, so it saves some typing if we do the conversion once and
    ! save them in new variables.
    ! ----------------------------------------------------------------------------
    Cmlatcl = INT ( Cmlat, KIND=C_LONG )
    Cmloncl = INT ( Cmlon, KIND=C_LONG )
    CmETAcl = INT ( CmETA, KIND=C_LONG )
    CmEp1cl = INT ( CmEp1, KIND=C_LONG )

    ! ---------------------------
    ! Allocate Climatology arrays
    ! ---------------------------
    locerrstat = pge_errstat_ok
    CALL climatology_allocate ( "a", Cmlat, Cmlon, CmETA, CmEp1, locerrstat )
    IF ( locerrstat /= pge_errstat_ok ) THEN
       errstat = MAX ( errstat, locerrstat ) 
       CALL climatology_allocate ( "d", Cmlat, Cmlon, CmETA, CmEp1, locerrstat )
       RETURN
    END IF

    ! -------------------------------
    ! Read dimension-defining arrays
    ! -------------------------------
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = Cmlatcl
    he5stat = HE5_SWrdfld ( swath_id, cli_lat_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, latvals(1:Cmlat) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = Cmloncl
    he5stat = HE5_SWrdfld ( swath_id, cli_lon_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, lonvals(1:Cmlon) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = CmEp1cl
    IF ( he5stat /= pge_errstat_ok ) &
         CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'Climatology arrays access failed.', vb_lev_default, errstat )

    ! -----------------------------------------------
    ! Read dimension-defining scale factor attributes
    ! -----------------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, cli_lat_field, "ScaleFactor", scale_lat )
    he5stat = HE5_SWrdlattr ( swath_id, cli_lon_field, "ScaleFactor", scale_lon )
    IF ( he5stat /= pge_errstat_ok ) &
         CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'Climatology attributes access failed.', vb_lev_default, errstat )

    ! -----------------------------------
    ! Apply scaling factors to geo fields
    ! -----------------------------------
    lonvals = lonvals * scale_lon
    latvals = latvals * scale_lat

    ! -----------------------------------------------
    ! Read the tables: Psurface, Heights, Temperature
    ! -----------------------------------------------
    he5_start_2d  = (/ zerocl, zerocl /)
    he5_stride_2d = (/  onecl,  onecl /)
    he5_edge_2d   = (/ Cmloncl, Cmlatcl /)

    he5_start_3d  = (/ zerocl, zerocl, zerocl/)
    he5_stride_3d = (/  onecl,  onecl,  onecl /)
    he5_edge_3d   = (/ Cmloncl, Cmlatcl, CmETAcl /)
    he5stat = HE5_SWrdfld (                                &
         swath_id, cli_Psurf_field,                        &
         he5_start_2d, he5_stride_2d, he5_edge_2d,         &
         Psurface(1:Cmlon,1:Cmlat) )
    he5stat = HE5_SWrdfld (                                &
         swath_id, cli_Temperature_field,                  &
         he5_start_3d, he5_stride_3d, he5_edge_3d,         &
         Temperature(1:Cmlon,1:Cmlat,1:CmETA) )

    ! --------------------------------------
    ! Read datafields scale factor attribute
    ! --------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, cli_Psurf_field,       "ScaleFactor", scale_Psurf       )
    he5stat = HE5_SWrdlattr ( swath_id, cli_Temperature_field, "ScaleFactor", scale_Temperature )

    IF ( he5stat /= pge_errstat_ok ) &
         CALL error_check ( he5stat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
         modulename//f_sep//'Climatology data fields access failed.', vb_lev_default, errstat )

    ! -----------------------------------------------------------------------
    ! Finding out the data field for the gas of interest (.eq. to target gas)
    ! -----------------------------------------------------------------------
    ndatafields = HE5_swinqdflds(swath_id, datafield_name, datafield_rank, datafield_type)
    CALL extract_swathname(nswath, TRIM(ADJUSTL(datafield_name)), &
         TRIM(ADJUSTL(sao_molecule_names(pge_idx))), gasdatafieldname)

    ! ---------------------------------------------------------------------------
    ! Check if we found the correct swath name. If not, report an error and exit.
    ! ---------------------------------------------------------------------------
    IF ( INDEX (TRIM(ADJUSTL(gasdatafieldname)),TRIM(ADJUSTL(sao_molecule_names(pge_idx)))) == 0 ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_error, OMSAO_E_HE5SWLOCATE, modulename, &
            vb_lev_default, locerrstat )
       WRITE(*,*) "Climatology file does not contain data for ", sao_molecule_names(pge_idx)
       errstat = MAX ( errstat, locerrstat )
       RETURN
    END IF

    ! -----------------------------
    ! Read data from this datafield
    ! -----------------------------
    he5stat = HE5_SWrdfld (                                &
         swath_id, TRIM(ADJUSTL(gasdatafieldname)),        &
         he5_start_3d, he5_stride_3d, he5_edge_3d,         &
         Gas_profiles(1:Cmlon,1:Cmlat,1:CmETA) )

    ! -----------------------------------------
    ! Read gas datafield scale factor attribute
    ! -----------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(gasdatafieldname)),&
              "ScaleFactor", scale_gas       )

    ! -------------------------------------------------------------------
    ! Finding out the data field for water vapor (to compute air density)
    ! -------------------------------------------------------------------
    ndatafields = HE5_swinqdflds(swath_id, datafield_name, datafield_rank, datafield_type)
    CALL extract_swathname(nswath, TRIM(ADJUSTL(datafield_name)), &
         TRIM(ADJUSTL(sao_molecule_names(h2o_cli_idx))), gasdatafieldname)

    ! ---------------------------------------------------------------------------
    ! Check if we found the correct swath name. If not, report an error and exit.
    ! ---------------------------------------------------------------------------
    IF ( INDEX (TRIM(ADJUSTL(gasdatafieldname)),TRIM(ADJUSTL(sao_molecule_names(pge_idx)))) == 0 ) THEN
       CALL error_check ( &
            0, 1, pge_errstat_error, OMSAO_E_HE5SWLOCATE, modulename, &
            vb_lev_default, locerrstat )
       WRITE(*,*) "Climatology file does not contain data for ", sao_molecule_names(h2o_cli_idx)
       errstat = MAX ( errstat, locerrstat )
       RETURN
    END IF

    ! -----------------------------
    ! Read data from this datafield
    ! -----------------------------
    he5stat = HE5_SWrdfld (                                &
         swath_id, TRIM(ADJUSTL(gasdatafieldname)),        &
         he5_start_3d, he5_stride_3d, he5_edge_3d,         &
         H2O_profiles(1:Cmlon,1:Cmlat,1:CmETA) )

    ! -----------------------------------------
    ! Read gas datafield scale factor attribute
    ! -----------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(gasdatafieldname)),&
              "ScaleFactor", scale_H2O       )

    ! ------------------------------------
    ! Apply scaling factors to data fields
    ! ------------------------------------
    Temperature  = Temperature  * scale_Temperature
    Psurface     = Psurface     * scale_Psurf
    Gas_profiles = Gas_profiles * scale_gas
    H2O_profiles = H2O_profiles * scale_H2O

  END SUBROUTINE omi_read_climatology

  SUBROUTINE omi_omler_albedo( lat, lon, albedo, nt, nx, xtrange, &
                               errstat)

    ! ==================================================================
    ! This subroutine reads the OMLER albedo data base for the month of
    ! the orbit to processed. Then it interpolates the values for each
    ! one of the pixels of the orbit to be analyzed
    ! ==================================================================

    USE OMSAO_variables_module, ONLY: OMSAO_OMLER_filename, &
                                      winwav_min, winwav_max

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4),                         INTENT (INOUT) :: errstat
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1), INTENT (INOUT) :: albedo

    ! ------------------------------------------------------------------
    ! Local variables, the variables to hold the OMLER data are going to
    ! be allocated and deallocated within this subroutine.
    ! ------------------------------------------------------------------
    REAL    (KIND=r4), ALLOCATABLE, DIMENSION(:) :: OMLER_longitude,   &
                                                    OMLER_latitude,    &
                                                    OMLER_wvl
    INTEGER (KIND=i2), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
                                    OMLER_monthly_albedo
    REAL    (KIND=r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
                                    OMLER_albedo, OMLER_wvl_albedo

    ! --------------------
    ! More Local variables
    ! --------------------
    CHARACTER (LEN=34),  PARAMETER :: grid_name = &
                                   'EarthSurfaceReflectanceClimatology'
    CHARACTER (LEN=maxchlen)       :: grid_file

    INTEGER   (KIND=i4), DIMENSION(1)   :: minwvl, maxwvl, minlon,    &
                                           maxlon, minlat, maxlat
    INTEGER   (KIND=i4), PARAMETER      :: OMLER_n_latitudes   = 360, &
                                           OMLER_n_longitudes  = 720, &
                                           OMLER_n_wavelenghts =  23, &
                                           one                 =   1
    INTEGER   (KIND=i4)                 :: itimes, ixtrack, spix,     &
                                           epix, ilon, ilat, nlon,    &
                                           nlat, OMnwvl
    INTEGER   (KIND=i4)                 :: grid_id, grid_file_id, month


    REAL      (KIND=r8), DIMENSION(1)   :: plon, plat, midwvl
    REAL      (KIND=r4)                 :: scale_factor, offset
    REAL      (KIND=r8)                 :: lonp, latp

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER (KIND=i4) :: version, locerrstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=16), PARAMETER :: modulename = 'omi_omler_albedo' 

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    locerrstat = pge_errstat_ok

    grid_file = TRIM(ADJUSTL(OMSAO_OMLER_filename))

    ! -------------------------------------------------------------------------------
    ! Open he5 OMI OMLER grid file and check GRID_FILE_ID (-1 if error)
    ! -------------------------------------------------------------------------------
    grid_file_id = HE5_GDOPEN (grid_file, he5f_acc_rdonly)
    IF (grid_file_id == he5_stat_fail) THEN
       CALL error_check (0, 1, pge_errstat_error, OMSAO_E_HE5GDOPEN, modulename, &
            vb_lev_default, locerrstat)
       errstat = MAX (errstat, locerrstat)
       RETURN
    END IF

    ! ----------------------------------------------
    ! Attach to grid and check GRID_ID (-1 if error)
    ! ----------------------------------------------
    grid_id = HE5_GDattach (grid_file_id, grid_name)
    IF (grid_id == he5_stat_fail) THEN
       CALL error_check (0, 1, pge_errstat_error, OMSAO_E_HE5GDATTACH, modulename, &
            vb_lev_default, locerrstat)
       errstat = MAX( errstat, locerrstat)
    END IF

    ! -------------------------
    ! Read longitude data field
    ! -------------------------
    ALLOCATE (OMLER_longitude(OMLER_n_longitudes))
    he5_start_1d = 0; he5_stride_1d = 1; he5_edge_1d = OMLER_n_longitudes
    locerrstat = HE5_GDRDFLD(grid_id, "Longitude", he5_start_1d, he5_stride_1d, &
                             he5_edge_1d, OMLER_longitude)
    errstat = MAX (errstat, locerrstat)

    ! ------------------------
    ! Read latitude data field
    ! ------------------------
    ALLOCATE (OMLER_latitude(OMLER_n_latitudes))
    he5_start_1d = 0; he5_stride_1d = 1; he5_edge_1d = OMLER_n_latitudes
    locerrstat = HE5_GDRDFLD(grid_id, "Latitude", he5_start_1d, he5_stride_1d, &
                             he5_edge_1d, OMLER_latitude)
    errstat = MAX (errstat, locerrstat)

    ! ---------------------------
    ! Read wavelenghts data field
    ! ---------------------------
    ALLOCATE (OMLER_wvl(OMLER_n_wavelenghts))
    he5_start_1d = 0; he5_stride_1d = 1; he5_edge_1d = OMLER_n_wavelenghts
    locerrstat = HE5_GDRDFLD(grid_id, "Wavelength", he5_start_1d, he5_stride_1d, &
                             he5_edge_1d, OMLER_wvl)
    errstat = MAX (errstat, locerrstat)

    ! -----------------------------------------------
    ! Select the wavelenghts; finding array positions
    ! -----------------------------------------------
    midwvl = (amf_wvl + amf_wvl2) / 2.0_r8
    minwvl = MINLOC(OMLER_wvl, OMLER_wvl .GE. REAL(winwav_min,KIND=r4))
    maxwvl = MAXLOC(OMLER_wvl, OMLER_wvl .LE. REAL(winwav_max,KIND=r4))
    OMnwvl = maxwvl(1)-minwvl(1)+1

    ! ---------------------------
    ! Read the albedo data field:
    ! -Month
    ! -Selected wavelenghts
    ! ---------------------------
    ALLOCATE (OMLER_monthly_albedo(OMLER_n_longitudes, &
              OMLER_n_latitudes, OMnwvl,1))
    ALLOCATE (OMLER_wvl_albedo(OMLER_n_longitudes,     &
              OMLER_n_latitudes, OMnwvl,1))
    ALLOCATE (OMLER_albedo(OMLER_n_longitudes,         &
              OMLER_n_latitudes, 1,1))

    month = granule_month - 1

    he5_start_4d  = (/  0,  0, minwvl-1, month/)
    he5_stride_4d = (/  1,  1,        1,     1/)
    he5_edge_4d   = (/OMLER_n_longitudes,OMLER_n_latitudes, &
                      OMnwvl, 1/)
    locerrstat = HE5_GDRDFLD(grid_id, "MonthlySurfaceReflectance",  &
      he5_start_4d, he5_stride_4d, he5_edge_4d, OMLER_monthly_albedo)
    errstat = MAX (errstat, locerrstat)

    ! ----------------------------
    ! Read the albedo scale factor
    ! ----------------------------
    locerrstat = HE5_GDRDLATTR(grid_id, "MonthlySurfaceReflectance", &
         "ScaleFactor", scale_factor)
    locerrstat = HE5_GDRDLATTR(grid_id, "MonthlySurfaceReflectance", &
         "Offset", offset)

    ! --------------------
    ! Deattached from grid
    ! --------------------
    locerrstat = HE5_GDDETACH(grid_id)

    ! -------------------------------------------
    ! Close he5 OMI OMLER grid file (-1 if error)
    ! -------------------------------------------
    locerrstat = HE5_GDclose ( grid_file_id)
    IF ( locerrstat == he5_stat_fail) THEN
       CALL error_check (0, 1, pge_errstat_error, OMSAO_W_HE5GDCLOSE, modulename, &
            vb_lev_default, locerrstat)
       errstat = MAX (errstat, locerrstat)
       RETURN
    END IF
  
    OMLER_wvl_albedo = REAL(offset, KIND = r8) +          &
                      REAL(scale_factor, KIND = r8)*      &
                      REAL(OMLER_monthly_albedo, KIND=r8)
 
    ! ----------------------------------------------------------------
    ! Interpolate for each pixel to one single wavelenght, just at the
    ! mid point of the fitting window: midwvl.
    ! ----------------------------------------------------------------
    DO ilon = 1, OMLER_n_longitudes
       DO ilat = 1, OMLER_n_latitudes
          CALL ezspline_1d_interpolation (                            &
               OMnwvl, REAL(OMLER_wvl(minwvl(1):maxwvl(1)), KIND=r8), &
               OMLER_wvl_albedo(ilon,ilat,1:OMnwvl,1),                &
               one, midwvl, OMLER_albedo(ilon,ilat,1,1), errstat )
       END DO
    END DO

    ! --------------------------------------------------
    ! Interpolate to the lat and longitude of each pixel
    ! --------------------------------------------------
    DO itimes = 0, nt-1

       spix = xtrange(itimes,1); epix = xtrange(itimes,2)
       DO ixtrack = spix, epix
          
          plon(1) = r8_missval
          plat(1) = r8_missval
          plon(1) = REAL(lon(ixtrack,itimes), KIND=r8)
          plat(1) = REAL(lat(ixtrack,itimes), KIND=r8)
          lonp = plon(1)
          latp = plat(1)

          IF (lon(ixtrack,itimes) .LT. -180.0_r4 .OR. &
              lat(ixtrack,itimes) .LT.  -90.0_r4 .OR. &
              lon(ixtrack,itimes) .GT.  180.0_r4 .OR. &
              lat(ixtrack,itimes) .GT.  90.0_r4) CYCLE

          ! -------------------------------------------------
          ! To speed up the interpolation a field of +- 2 deg
          ! rees on latitude and longitude will be used from
          ! OMLER centered around the OMI pixel.
          ! ----------------------------------------------------------
          ! Finding which OMLER pixels go in to the interpolation
          ! If we are too close to the boundaries I just move them to
          ! the interior. Good for the latitude, no so good for the 
          ! longitude, but over the Pacific it should make no much dif
          ! ference
          ! ----------------------------------------------------------
          IF (latp .LT. -88.0) latp = -89.0
          IF (latp .GT.  88.0) latp =  89.0
          IF (lonp .LT. -178.0) lonp = -179.0
          IF (lonp .GT.  178.0) lonp =  179.0

          minlon = MINLOC(OMLER_longitude, OMLER_longitude .GE. lonp-1.0)
          maxlon = MAXLOC(OMLER_longitude, OMLER_longitude .LE. lonp+1.0)
          minlat = MINLOC(OMLER_latitude,  OMLER_latitude  .GE. latp-1.0)
          maxlat = MAXLOC(OMLER_latitude,  OMLER_latitude  .LE. latp+1.0)
          nlon   = maxlon(1)-minlon(1)+1
          nlat   = maxlat(1)-minlat(1)+1                

          CALL ezspline_2d_interpolation ( nlon, nlat,             &
               REAL(OMLER_longitude(minlon(1):maxlon(1)),KIND=r8), &
               REAL(OMLER_latitude (minlat(1):maxlat(1)),KIND=r8), &
               OMLER_albedo(minlon(1):maxlon(1),                   &
                            minlat(1):maxlat(1), 1, 1),            &
               one, one, plon, plat, albedo(ixtrack,itimes),       &
               locerrstat )
       END DO
    END DO
    
    ! --------------------
    ! Deallocate variables
    ! --------------------
    DEALLOCATE (OMLER_monthly_albedo)
    DEALLOCATE (OMLER_wvl_albedo)
    DEALLOCATE (OMLER_albedo)
    DEALLOCATE (OMLER_longitude)
    DEALLOCATE (OMLER_latitude)
    DEALLOCATE (OMLER_wvl)
   
    errstat = MAX(errstat, locerrstat)
    
  END SUBROUTINE omi_omler_albedo
    
  SUBROUTINE extract_swathname ( nswath, multi_swath, swathstr, single_swath )

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
  END SUBROUTINE extract_swathname

  SUBROUTINE climatology_getdim ( &
       swath_id, Cmlat, Cmlon, CmETA, CmEp1, errstat )

    ! --------------------------------
    ! Return dimensions of Climatology
    ! --------------------------------
 
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: swath_id

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    INTEGER (KIND=i4), INTENT (OUT)   :: Cmlat, Cmlon, CmETA, CmEp1

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4), PARAMETER               :: maxdim = 100
    INTEGER   (KIND=i4)                          :: ndim, nsep
    INTEGER   (KIND=C_LONG)                      :: ndimcl
    INTEGER   (KIND=i4)                          :: i, j, swlen, iend, istart
    INTEGER   (KIND=i4),  DIMENSION(0:maxdim)    :: dim_array, dim_seps
    INTEGER   (KIND=C_LONG), DIMENSION(0:maxdim) :: dim_arraycl
    CHARACTER (LEN=10*maxdim)                    :: dim_chars


    ! ---------------------------
    ! Initialize output variables
    ! ---------------------------
    Cmlat = -1 ; Cmlon = -1 ; CmETA = -1; CmEp1 = -1

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
       CASE ( "nLat" )
          Cmlat = dim_array(j)
       CASE ( "nLon" )
          Cmlon = dim_array(j)
       CASE ( "nETA" )
          CmETA = dim_array(j)
       CASE ( "nEp1")
          CmEp1 = dim_array(j)
       CASE DEFAULT
          ! Whatever. Nothing to be done here.
       END SELECT

    END DO getdims

    RETURN
  END SUBROUTINE climatology_getdim

  SUBROUTINE climatology_allocate ( ad, Cmlat, Cmlon, CmETA, CmEp1, errstat )

    USE OMSAO_casestring_module, ONLY: lower_case
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=1), INTENT (IN) :: ad
    INTEGER (KIND=i4), INTENT (IN) :: Cmlat, Cmlon, CmETA, CmEp1

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
       ALLOCATE (latvals (Cmlat),                 STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (lonvals (Cmlon),                 STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (Temperature (Cmlon, Cmlat, CmETA), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (Gas_profiles(Cmlon, Cmlat, CmETA), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (H2O_profiles(Cmlon, Cmlat, CmETA), STAT=estat ) ; errstat = MAX ( errstat, estat )
       ALLOCATE (Psurface(Cmlon,Cmlat),            STAT=estat ) ; errstat = MAX ( errstat, estat )
    CASE ('d')
       IF ( ALLOCATED ( latvals      ) )  DEALLOCATE ( latvals      )
       IF ( ALLOCATED ( lonvals      ) )  DEALLOCATE ( lonvals      )
       IF ( ALLOCATED ( Temperature  ) )  DEALLOCATE ( Temperature  )
       IF ( ALLOCATED ( Gas_profiles ) )  DEALLOCATE ( Gas_profiles )
       IF ( ALLOCATED ( H2O_profiles ) )  DEALLOCATE ( H2O_profiles )
       IF ( ALLOCATED ( Psurface     ) )  DEALLOCATE ( Psurface     )
    CASE DEFAULT
       ! Whatever. Nothing to be done here
    END SELECT

    RETURN
  END SUBROUTINE climatology_allocate

  SUBROUTINE vlidort_allocate ( ad, anozo, ansrf, ancld, analb, ansza, anvza, anwav, anlay, anlev, errstat )

    USE OMSAO_casestring_module, ONLY: lower_case
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=1), INTENT (IN) :: ad
    INTEGER (KIND=i4), INTENT (IN) :: anozo, ansrf, ancld, analb, ansza, anvza, anwav, anlay, anlev

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

       ALLOCATE (lut_alb(1:analb),     STAT=estat)
       ALLOCATE (lut_srf(1:ansrf),        STAT=estat)
       ALLOCATE (lut_clp(1:ancld),        STAT=estat)
       ALLOCATE (lut_sza(1:ansza),        STAT=estat)
       ALLOCATE (lut_vza(1:anvza),        STAT=estat)
       ALLOCATE (lut_wavelength(1:anwav), STAT=estat)
       ALLOCATE (lut_toms(1:anozo),       STAT=estat)
       
       ALLOCATE (lut_air(1:anozo,1:ansrf,1:anlay), STAT=estat)
       ALLOCATE (lut_alt_lay(1:ansrf,1:anlay),   STAT=estat)
       ALLOCATE (lut_alt_lev(1:ansrf,1:anlev),   STAT=estat)
       ALLOCATE (lut_ozo(1:anozo,1:ansrf,1:anlay), STAT=estat)
       ALLOCATE (lut_pre_lay(1:ansrf,1:anlay),   STAT=estat)
       ALLOCATE (lut_pre_lev_cld(1:ancld,1:anlev),   STAT=estat)
       ALLOCATE (lut_pre_lay_cld(1:ancld,1:anlay),   STAT=estat)
       ALLOCATE (lut_pre_lev(1:ansrf,1:anlev),   STAT=estat)

       ALLOCATE (lut_tem(1:anozo,1:ansrf,1:anlev), STAT=estat)
       
       ALLOCATE (lut_I0_clr(1:anozo,1:ansrf,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_I1_clr(1:anozo,1:ansrf,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_I2_clr(1:anozo,1:ansrf,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_Ir_clr(1:anozo,1:ansrf,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_Sb_clr(1:anozo,1:ansrf),                 STAT=estat)
       
       ALLOCATE (lut_I0_cld(1:anozo,1:ancld,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_I1_cld(1:anozo,1:ancld,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_I2_cld(1:anozo,1:ancld,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_Ir_cld(1:anozo,1:ancld,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_Sb_cld(1:anozo,1:ancld),                 STAT=estat)
       
       ALLOCATE (lut_dI0_clr(1:anozo,1:ansrf,1:anlay,1:analb,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_dI1_clr(1:anozo,1:ansrf,1:anlay,1:analb,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_dI2_clr(1:anozo,1:ansrf,1:anlay,1:analb,1:anvza,1:ansza), STAT=estat)
       
       ALLOCATE (lut_dI0_cld(1:anozo,1:anlay,1:ancld,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_dI1_cld(1:anozo,1:anlay,1:ancld,1:anvza,1:ansza), STAT=estat)
       ALLOCATE (lut_dI2_cld(1:anozo,1:anlay,1:ancld,1:anvza,1:ansza), STAT=estat)    

    CASE ('d')

       IF ( ALLOCATED ( lut_alb ) ) DEALLOCATE ( lut_alb )
       IF ( ALLOCATED ( lut_srf ) ) DEALLOCATE ( lut_srf )
       IF ( ALLOCATED ( lut_clp ) ) DEALLOCATE ( lut_clp )
       IF ( ALLOCATED ( lut_sza ) ) DEALLOCATE ( lut_sza )
       IF ( ALLOCATED ( lut_vza ) ) DEALLOCATE ( lut_vza )
       IF ( ALLOCATED ( lut_wavelength ) ) DEALLOCATE ( lut_wavelength )
       IF ( ALLOCATED ( lut_toms ) ) DEALLOCATE ( lut_toms )
       
       IF ( ALLOCATED ( lut_air ) ) DEALLOCATE ( lut_air )
       IF ( ALLOCATED ( lut_alt_lev ) ) DEALLOCATE ( lut_alt_lev )
       IF ( ALLOCATED ( lut_alt_lay ) ) DEALLOCATE ( lut_alt_lay )
       IF ( ALLOCATED ( lut_ozo ) ) DEALLOCATE ( lut_ozo )
       IF ( ALLOCATED ( lut_pre_lev ) ) DEALLOCATE ( lut_pre_lev )
       IF ( ALLOCATED ( lut_pre_lay ) ) DEALLOCATE ( lut_pre_lay )
       IF ( ALLOCATED ( lut_pre_lay_cld ) ) DEALLOCATE ( lut_pre_lay_cld )
       IF ( ALLOCATED ( lut_pre_lev_cld ) ) DEALLOCATE ( lut_pre_lev_cld )

       IF ( ALLOCATED ( lut_tem ) ) DEALLOCATE ( lut_tem )
       
       IF ( ALLOCATED ( lut_I0_clr ) ) DEALLOCATE ( lut_I0_clr )
       IF ( ALLOCATED ( lut_I1_clr ) ) DEALLOCATE ( lut_I1_clr )
       IF ( ALLOCATED ( lut_I2_clr ) ) DEALLOCATE ( lut_I2_clr )
       IF ( ALLOCATED ( lut_Ir_clr ) ) DEALLOCATE ( lut_Ir_clr )
       IF ( ALLOCATED ( lut_Sb_clr ) ) DEALLOCATE ( lut_Sb_clr )
       
       IF ( ALLOCATED ( lut_I0_cld ) ) DEALLOCATE ( lut_I0_cld )
       IF ( ALLOCATED ( lut_I1_cld ) ) DEALLOCATE ( lut_I1_cld )
       IF ( ALLOCATED ( lut_I2_cld ) ) DEALLOCATE ( lut_I2_cld )
       IF ( ALLOCATED ( lut_Ir_cld ) ) DEALLOCATE ( lut_Ir_cld )
       IF ( ALLOCATED ( lut_Sb_cld ) ) DEALLOCATE ( lut_Sb_cld )
       
       IF ( ALLOCATED ( lut_dI0_clr ) ) DEALLOCATE ( lut_dI0_clr ) 
       IF ( ALLOCATED ( lut_dI1_clr ) ) DEALLOCATE ( lut_dI1_clr ) 
       IF ( ALLOCATED ( lut_dI2_clr ) ) DEALLOCATE ( lut_dI2_clr ) 
       
       IF ( ALLOCATED ( lut_dI0_cld ) ) DEALLOCATE ( lut_dI0_cld ) 
       IF ( ALLOCATED ( lut_dI1_cld ) ) DEALLOCATE ( lut_dI1_cld ) 
       IF ( ALLOCATED ( lut_dI2_cld ) ) DEALLOCATE ( lut_dI2_cld ) 

    CASE DEFAULT
       ! Whatever. Nothing to be done here
    END SELECT

    RETURN
  END SUBROUTINE vlidort_allocate

  SUBROUTINE compute_geometric_amf ( nt, nx, sza, vza, xtrange, amfgeo, amfdiag )

    USE OMSAO_parameters_module, ONLY: deg2rad, rad2deg
    USE OMSAO_variables_module,  ONLY: szamax
    USE OMSAO_omidata_module,    ONLY: omi_geo_amf, omi_oobview_amf

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

  SUBROUTINE amf_read_omiclouds ( nt, nx, yn_szoom, l2cfr, l2ctp, errstat )

    USE OMSAO_variables_module,  ONLY: voc_amf_filenames
    USE OMSAO_indices_module,    ONLY: voc_omicld_idx
    USE OMSAO_omidata_module,    ONLY: gzoom_spix, gzoom_epix, gzoom_npix

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
    LOGICAL                  :: yn_raman_clouds

    ! ---------------------------------------
    ! For accesing the file (local variables)
    ! ---------------------------------------
    CHARACTER (LEN=maxchlen) :: omicloud_swath_name    = 'undefined'
    INTEGER   (KIND=i4)      :: omicloud_swath_id      = -1
    INTEGER   (KIND=i4)      :: omicloud_swath_file_id = -1

    ! ----------------------------------------------
    ! Names of various HE5 fields to read from files
    ! ----------------------------------------------
    !(1) OMI Cloud Fields
    ! -------------------
    CHARACTER (LEN=13), PARAMETER :: omicld_cfrac_field      = 'CloudFraction'
    CHARACTER (LEN=13), PARAMETER :: omicld_cpres_field      = 'CloudPressure'

    ! ----------------------
    ! Name of the subroutine
    ! ----------------------
    CHARACTER (LEN=22), PARAMETER :: modulename = 'amf_read_omiclouds'

    locerrstat = pge_errstat_ok

    ! -------------------------
    ! Attach to OMI Cloud Swath
    ! -------------------------
    CALL he5_init_input_file ( &
         voc_amf_filenames(voc_omicld_idx), omicloud_swath_name, &
         omicloud_swath_id, omicloud_swath_file_id,              &
         nt_loc, nx_loc, locerrstat )

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
    ELSE
       yn_raman_clouds = .FALSE.
       addstr          = ""
    END IF

    ! --------------------------------------------
    ! Read scaling of cloud data fields (working?)
    ! --------------------------------------------
    scale_cfr = 1.0_r4 ; offset_cfr = 0.0_r4 ; missval_cfr = 0.0_r4
    locerrstat = HE5_SWrdlattr ( omicloud_swath_id, omicld_cfrac_field//TRIM(ADJUSTL(addstr)), 'MissingValue', missval_cfr )

    scale_ctp = 1.0_r4 ; offset_ctp = 0.0_r4 ; missval_ctp = 0
    locerrstat = HE5_SWrdlattr ( omicloud_swath_id, omicld_cpres_field//TRIM(ADJUSTL(addstr)), 'MissingValue', missval_ctp )

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
         omicloud_swath_id, omicld_cfrac_field//TRIM(ADJUSTL(addstr)),   &
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
            omicloud_swath_id, omicld_cpres_field//TRIM(ADJUSTL(addstr)),   &
            he5_start_2d, he5_stride_2d, he5_edge_2d, ctp(1:nx,0:nt-1) )
    ELSE
       locerrstat = HE5_SWrdfld ( &
            omicloud_swath_id, omicld_cpres_field,   &
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
  END SUBROUTINE amf_read_omiclouds

  SUBROUTINE voc_amf_readisccp ( i0, errstat )

    USE OMSAO_variables_module,  ONLY: voc_amf_filenames
    USE OMSAO_indices_module,    ONLY: voc_isccp_idx

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
  
  SUBROUTINE amf_diagnostic ( &
       nt, nx, lat, lon, sza, vza, snow, glint, xtrange, ctpmin, ctpmax, l2cfr, l2ctp, amfdiag )

    USE OMSAO_omidata_module,   ONLY: omi_oobview_amf, omi_glint_add, omi_height, &
         omi_geo_amf, omi_bigsza_amf, omi_cfr_addmiss, omi_ctp_addmiss, &
         omi_oob_cfr, omi_oob_ctp
    USE OMSAO_variables_module, ONLY: winwav_min, winwav_max
    
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    REAL    (KIND=r4),                          INTENT (IN) :: ctpmin, ctpmax
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon, sza, vza
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
    ! -------------------------------------------------------------------
    ! If the file provided for Scattering weights has no information with
    ! in the fitting window, then nothing is to be done here.
    ! The geometric AMFs flag will remain in amfdiag and no further calcu
    ! lation will be performed inside the calculate_amf and calculate_sca
    ! subroutines.
    ! -------------------------------------------------------------------
    ! ----------------------------------------------------------
    ! Check that the value in amf_wvl is within the range of the
    ! suplied scattering weights file.
    ! ----------------------------------------------------------
    IF ( &
         (amf_wvl  .LT. MINVAL(lut_wavelength) ) .OR. &
         (amf_wvl2 .GT. MAXVAL(lut_wavelength) ) ) THEN
       RETURN
    END IF

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
            ( sza(spix:epix,it) < MINVAL(lut_sza) ) .OR. &
            ( sza(spix:epix,it) > MAXVAL(lut_sza) ) .OR. &
            ( vza(spix:epix,it) < MINVAL(lut_vza) ) .OR. &
            ( vza(spix:epix,it) > MAXVAL(lut_vza) )      )
          amfdiag(spix:epix,it) = omi_oobview_amf
       END WHERE

       ! ---------------------------------------------------------------------
       ! Where the OMI cloud information is missing, select ISCCP climatology,
       ! but skip pixels where the SZA values are out of the AMF table bounds.
       ! ---------------------------------------------------------------------
       IF ( ( ANY( l2cfr(spix:epix,it)   < 0.0_r8          ) ) .OR. &
            ( ANY( l2ctp(spix:epix,it)   < 0.0_r8          ) ) ) THEN

          DO ix = spix, epix
             IF (amfdiag(ix,it) > omi_oobview_amf) THEN
                latdp = REAL ( lat(ix,it), KIND=r8 ) ; londp = REAL ( lon(ix,it), KIND=r8 ) ;
                ilat = MAXVAL(MINLOC( ABS(ISCCP_CloudClim%latvals-latdp) ))
                j1   = SUM(ISCCP_CloudClim%n_lonvals(1:ilat-1)) + 1
                j2   = ISCCP_CloudClim%n_lonvals(ilat) + j1
                ilon = MAXVAL(MINLOC( ABS(ISCCP_CloudClim%lonvals(j1:j2)-londp) ))
                
                amfdiag(ix,it) = 0_i2
                IF ( l2ctp(ix,it) < 0.0_r8 .AND. ISCCP_CloudClim%ctp(ilon) >= 0.0_r8 ) THEN
                   amfdiag(ix,it) = omi_ctp_addmiss + amfdiag(ix,it)
                   l2ctp  (ix,it) = ISCCP_CloudClim%ctp(ilon)
                END IF
                IF ( l2cfr(ix,it) < 0.0_r8 .AND. ISCCP_CloudClim%cfr(ilon) >= 0.0_r8 ) THEN
                   amfdiag(ix,it) = omi_cfr_addmiss + amfdiag(ix,it)
                   l2cfr  (ix,it) = ISCCP_CloudClim%cfr(ilon)
                END IF
             ENDIF
          END DO
       END IF
     
       ! ---------------------------------------------------
       ! Out of bounds clouds (too high or too low), make it
       ! the highest possible value in the look up table.
       ! -------------------------------------------------
       WHERE ( &
            ( l2ctp(spix:epix,it)   <  MINVAL(lut_clp) ) .AND. &
            ( amfdiag(spix:epix,it) >  omi_oobview_amf ) )
          l2ctp(spix:epix,it) = MINVAL(lut_clp)
          amfdiag(spix:epix,it) = amfdiag(spix:epix,it) + omi_oob_ctp
       END WHERE

       WHERE ( &
            ( l2ctp(spix:epix,it)   >  MAXVAL(lut_clp) ) .AND. &
            ( amfdiag(spix:epix,it) >  omi_oobview_amf ) )
          l2ctp(spix:epix,it) = MAXVAL(lut_clp)
          amfdiag(spix:epix,it) = amfdiag(spix:epix,it) + omi_oob_ctp
       END WHERE

       ! -------------------------------------------------------
       ! For pixel without cloud information set flag to oobview
       ! -------------------------------------------------------
       WHERE ( &
            ( l2cfr(spix:epix,it) .EQ. r8_missval ) .OR. &
            ( l2ctp(spix:epix,it) .EQ. r8_missval)       )
              amfdiag(spix:epix,it) = omi_oobview_amf
       END WHERE

       ! ------------------------------------------------------
       ! And AMFDIAG values > OOB must be good and are set to 0
       ! if we have "good" clouds
       ! ------------------------------------------------------
       WHERE ( &
            ( amfdiag(spix:epix,it) > omi_oobview_amf ) .AND. &
            ( amfdiag(spix:epix,it) < omi_cfr_addmiss )  )
             amfdiag(spix:epix,it) = 0_i2
       END WHERE

       ! --------------------------------------------------
       ! Angles above the top value set on the control file
       ! are calculated "using this maximum value".
       ! --------------------------------------------------
       WHERE ( &
            ( sza(spix:epix,it)     .GE. amf_max_sza       ) .AND. &
            ( amfdiag(spix:epix,it) .GT. omi_oobview_amf ) )
              amfdiag(spix:epix,it) = omi_bigsza_amf + amfdiag(spix:epix,it)
       END WHERE

       ! -----------------------
       ! Start with the ice flag
       ! -----------------------
       WHERE (                                       &
            amfdiag     (spix:epix,it) >= 0_i2 .AND. &
            snow(spix:epix,it) >= 0_i2         )
          amfdiag(spix:epix,it) = snow(spix:epix,it) + amfdiag(spix:epix,it)
       END WHERE

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
  END SUBROUTINE amf_diagnostic

  SUBROUTINE compute_scatt ( nt, nx, albedo, lat, sza, vza, saa, vaa, l2ctp, l2cfr, terrain_height, amfgeo, &
       amfdiag, scattw)

    USE OMSAO_lininterpolation_module
    USE OMSAO_variables_module, ONLY: verb_thresh_lev
    USE OMSAO_omidata_module, ONLY: omi_ozone_amount

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                INTENT (IN) :: nt, nx
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1),       INTENT (IN) :: amfdiag
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1),       INTENT (IN) :: lat, sza, vza, saa, vaa, terrain_height
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1),       INTENT (IN) :: albedo, l2ctp, l2cfr, amfgeo

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,ngeos5-1), INTENT (INOUT) :: scattw

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: itime, ixtrack, status, one, &
         iozo, iwav, isza, ivza, icld, isrf, ipre, ialt, ialb, ilev
    INTEGER (KIND=i4), DIMENSION(1) :: iwavs, iwavf
    REAL    (KIND=r8) :: grad, raa, tmp_saa, tmp_vaa
    REAL    (KIND=r8) :: crf, nwavs
    REAL    (KIND=r8), DIMENSION(1), PARAMETER :: one_r8 =(/1.0_r8/)
    REAL    (KIND=r8), DIMENSION(srf_dim(1), vza_dim(1), sza_dim(1)) :: Inte_clear_3D
    REAL    (KIND=r8), DIMENSION(vza_dim(1), sza_dim(1))             :: Inte_clear_2D
    REAL    (KIND=r8), DIMENSION(sza_dim(1))                         :: Inte_clear_1D
    REAL    (KIND=r8), DIMENSION(1)                                  :: Radiance_clr
    REAL    (KIND=r8), DIMENSION(clp_dim(1), vza_dim(1), sza_dim(1)) :: Inte_cloud_3D
    REAL    (KIND=r8), DIMENSION(vza_dim(1), sza_dim(1))             :: Inte_cloud_2D
    REAL    (KIND=r8), DIMENSION(sza_dim(1))                         :: Inte_cloud_1D
    REAL    (KIND=r8), DIMENSION(1)                                  :: Radiance_cld
    REAL    (KIND=r8), DIMENSION(srf_dim(1),alt_lay_dim(2),alb_dim(1),vza_dim(1),sza_dim(1)) :: sw_clear_5D
    REAL    (KIND=r8), DIMENSION(srf_dim(1),alt_lay_dim(2),vza_dim(1),sza_dim(1))            :: sw_clear_4D
    REAL    (KIND=r8), DIMENSION(alt_lay_dim(2),vza_dim(1),sza_dim(1))                       :: sw_clear_3D
    REAL    (KIND=r8), DIMENSION(alt_lay_dim(2),sza_dim(1))                                  :: sw_clear_2D
    REAL    (KIND=r8), DIMENSION(alt_lay_dim(2))                                             :: sw_clear_1D
    REAL    (KIND=r8), DIMENSION(ngeos5-1)                                                   :: sw_clear
    REAL    (KIND=r8), DIMENSION(alt_lay_dim(2),clp_dim(1),vza_dim(1),sza_dim(1)) :: sw_cloud_4D
    REAL    (KIND=r8), DIMENSION(alt_lay_dim(2),vza_dim(1),sza_dim(1))            :: sw_cloud_3D
    REAL    (KIND=r8), DIMENSION(alt_lay_dim(2),sza_dim(1))                       :: sw_cloud_2D
    REAL    (KIND=r8), DIMENSION(alt_lay_dim(2))                                  :: sw_cloud_1D
    REAL    (KIND=r8), DIMENSION(ngeos5-1)                                        :: sw_cloud
    REAL    (KIND=r8), DIMENSION(alt_lay_dim(2)) :: re_alt
    REAL    (KIND=r8), DIMENSION(srf_dim(1))     :: re_srf
    REAL    (KIND=r8), DIMENSION(alb_dim(1))     :: re_alb
    REAL    (KIND=r8), DIMENSION(clp_dim(1))     :: re_cld
    REAL    (KIND=r8), DIMENSION(sza_dim(1))     :: re_sza
    REAL    (KIND=r8), DIMENSION(vza_dim(1))     :: re_vza
    REAL    (KIND=r8), DIMENSION(1)              :: local_alb, local_sza, local_vza, local_raa, &
         local_srf, local_cld, local_cfr
    REAL    (KIND=r8), DIMENSION(ngeos5)         :: local_lev_pre
    REAL    (KIND=r8), DIMENSION(ngeos5-1)       :: local_lay_pre
    REAL    (KIND=r8), DIMENSION(srf_dim(1),alt_lay_dim(2)) :: re_pre_2D, re_pre_cld_2D
    REAL    (KIND=r8), DIMENSION(alt_lay_dim(2))            :: lay_pre_1D, lay_pre_cld_1D
    REAL    (KIND=r8), PARAMETER :: d2r = 3.141592653589793d0/180.0  !! JED fix

    ! To select TOMS profile
    REAL    (KIND=r8), PARAMETER :: du  = 2.69e16 ! molecules/cm^2
    INTEGER (KIND=i4)  :: toms_idx
    REAL    (KIND=r8), DIMENSION(10), PARAMETER :: hxxx = (/125,175,225,275,325,375,425,475,523,575/)
    REAL    (KIND=r8), DIMENSION(10), PARAMETER :: mxxx = (/125,175,225,275,325,375,425,475,523,575/)
    REAL    (KIND=r8), DIMENSION(6),  PARAMETER :: lxxx = (/225,275,325,375,425,475/)

    ! -----------------------------------
    ! Find look up table wavelength index
    ! No interpolation, closest available
    ! is selected.
    ! -----------------------------------
    one = 1_i4
    iwavs = MINLOC(ABS(lut_wavelength - REAL(amf_wvl,  KIND = r4) ))
    iwavf = MINLOC(ABS(lut_wavelength - REAL(amf_wvl2, KIND = r4) ))
    nwavs = REAL(iwavf(1), KIND=r8) - REAL(iwavs(1), KIND=r8) + 1.0_r8

    ! --------------------------------------------------------
    ! Re-order the pressure, altitude, sza, and vza dimensions
    ! from (look up table). Needed for interpolation. 
    ! This should be moved to the program that creates the
    ! look up tables
    ! ------------------------------- ------------------------
    DO isza = 1, sza_dim(1)
       re_sza(isza) = cos(d2r*REAL(lut_sza(sza_dim(1)+1-isza), KIND = r8))  ! JED fix
    END DO
    DO ivza = 1, vza_dim(1)
       re_vza(ivza) = cos(d2r*REAL(lut_vza(vza_dim(1)+1-ivza), KIND = r8)) ! JED fix
    END DO
    DO isrf = 1, srf_dim(1)
       re_srf(isrf) = REAL(lut_srf(srf_dim(1)+1-isrf), KIND = r8)
       re_pre_2D(isrf,1:alt_lay_dim(2)) = lut_pre_lay(srf_dim(1)+1-isrf,1:alt_lay_dim(2))
    END DO
    DO icld = 1, clp_dim(1)
       re_cld(icld) = REAL(lut_clp(clp_dim(1)+1-icld), KIND = r8)
       re_pre_cld_2D(icld,1:alt_lay_dim(2)) = lut_pre_lay_cld(clp_dim(1)+1-icld,1:alt_lay_dim(2))
    END DO
    re_alb(1:alb_dim(1)) = REAL(lut_alb(1:alb_dim(1)), KIND = r8)

    ! ---------------
    ! Loop over lines
    ! ---------------
    DO itime = 0, nt-1
       ! --------------------------
       ! Loop over xtrack positions
       ! --------------------------
       DO ixtrack = 1, nx

          toms_idx = -1
          IF (amfdiag(ixtrack,itime) .LT. 0) CYCLE

          ! -----------------------------------------
          ! Find out which TOMS profile we should use
          ! in function of latitude and ozone column.
          ! -----------------------------------------        
          IF ( ABS(lat(ixtrack,itime)) .GT. 60.0 ) THEN
             toms_idx = MINLOC(ABS(hxxx-omi_ozone_amount(ixtrack,itime)/du/amfgeo(ixtrack,itime)), 1) + 0
          ELSE IF ( ABS(lat(ixtrack,itime)) .GT. 30.0 .AND. ABS(lat(ixtrack,itime)) .LE. 60.0 ) THEN
             toms_idx = MINLOC(ABS(mxxx-omi_ozone_amount(ixtrack,itime)/du/amfgeo(ixtrack,itime)), 1) + 17
          ELSE IF ( ABS(lat(ixtrack,itime)) .LE. 30.0 ) THEN
             toms_idx = MINLOC(ABS(lxxx-omi_ozone_amount(ixtrack,itime)/du)/amfgeo(ixtrack,itime), 1) + 11
          ENDIF

          ! ----------------------------------------------
          ! Work out relative azimuth angle for this pixel
          ! ----------------------------------------------
          tmp_saa = REAL(saa(ixtrack,itime),KIND=r8); tmp_vaa = REAL(vaa(ixtrack,itime),KIND=r8)
          ! (1) Map [-180, +180] to [0, 360]
          IF ( tmp_saa .NE. r8_missval .AND. tmp_saa .LT. 0.0_r8 .AND. &
               tmp_saa .GE. -180.0_r8 ) THEN
             tmp_saa = 360.0_r4 + tmp_saa
          ENDIF
          IF ( tmp_vaa .NE. r8_missval .AND. tmp_vaa .LT. 0.0_r8 .AND. &
               tmp_vaa .GE. -180.0_r8 ) THEN
             tmp_vaa = 360.0_r8 + tmp_vaa
          ENDIF
          ! (2) Compute relative azimuth angle (RELATIVE means absolute value)
          IF ( tmp_saa .GE. 0.0_r8 .AND. tmp_vaa .GE. 0.0_r8 ) THEN
             raa = ABS( tmp_saa - tmp_vaa )
          ENDIF
          IF ( raa .GE. 180.0_r8 ) THEN
             raa = 360.0_r8 - raa
          ENDIF

          ! ----------------------------------------------
          ! If this point is reached then scattw should be
          ! different from r8_missval and it needs to be
          ! initialized to 0.0 to work out the average
          ! ----------------------------------------------
          ! Initialize variables
          ! --------------------
          scattw(ixtrack,itime,:) = 0.0_r8

          ! ----------------------------------------------
          ! If sza > amf_max_sza set it for calculation to
          ! amf_max_sza
          ! ----------------------------------------------
          local_sza(1) = REAL(sza(ixtrack,itime), KIND = r8)
          IF (local_sza(1) .GT. amf_max_sza) local_sza(1) = amf_max_sza

          local_alb(1) = REAL(albedo(ixtrack,itime), KIND=r8)
          local_cld(1) = REAL(l2ctp(ixtrack,itime), KIND=r8)
          local_cfr(1) = REAL(l2cfr(ixtrack,itime), KIND=r8)
          local_sza(1) = cos(d2r*local_sza(1))  ! JED fix
          local_vza(1) = cos(d2r*REAL(vza(ixtrack,itime), KIND = r8))  ! JED fix
          local_raa(1) = d2r*REAL(raa,KIND = r8)
          local_srf(1) = REAL(terrain_height(ixtrack,itime), KIND = r8)

          ! ----------------------------------------------
          ! Convert pixel terrain height to pressure using
          ! Xiong suggested to use pressure altitude:
          !  Z = -16 alog10 (P / Po) Z in km and P in hPa.
          ! ----------------------------------------------
          local_srf(1) = 1013.0_r8 * (10.0_r8 ** (local_srf(1) / 1000.0_r8 / (-16.0_r8))) 

          !Bringing surface pressure to highest available in lookup table if needed
          !Current highest surface pressure is 1030 hPa.
          IF (local_srf(1) .GT. MAXVAL(re_srf)) local_srf(1) = MAXVAL(re_srf)
          !Bringing clouds heights to highest available pressure if needed. We avoid extrapolation
          !Current highest cloud pressure is 1030 hPa.
          IF ( local_cld(1) .GT. MAXVAL(re_cld) ) local_cld(1) = MAXVAL(re_cld)
          !Be sure that clouds are below lower cloud pressure. (Ap, Bp calculation
          !of pressures above 350 hPa breaks down). Current lowest pressure 350 hPa ~8 km
          IF ( local_cld(1) .LT. MINVAL(re_cld) ) local_cld(1) = MINVAL(re_cld)
          !Be sure that clouds are above or at the surface.
          IF ( local_cld(1) .GT. local_srf(1) ) local_cld(1) = local_srf(1)

          ! ------------------------------------------------
          ! Now that we have the surface pressure I work out
          ! the pressure of each layer and layer edges using
          ! GEOS5 a_p and b_p values.
          ! ------------------------------------------------
          ! Level
          DO ilev = 1, ngeos5
             local_lev_pre(ngeos5+1-ilev) = (Ap(ilev) + local_srf(1) * Bp(ilev)) !hPa
          END DO
          ! Layer
          DO ialt = 1, ngeos5-1 ! for layer center
             local_lay_pre(ialt) = (local_lev_pre(ialt+1) + local_lev_pre(ialt))/2.0 ! hPa
          END DO      
          ! -----------------------
          ! Interpolate lut_pre_lay
          ! to local_srf(1)
          ! -----------------------
          DO ialt = 1, alt_lay_dim(2)
             CALL ezspline_1d_interpolation (INT(srf_dim(1),KIND=i4), re_srf, &
                  re_pre_2D(1:srf_dim(1),ialt), &
                  one, local_srf(1), &
                  lay_pre_1D(ialt), status)
          END DO
          ! ---------------------------
          ! Interpolate lut_pre_lay_cld
          ! to local_cld(1)
          ! ---------------------------
          DO ialt = 1, alt_lay_dim(2)
             CALL ezspline_1d_interpolation (INT(clp_dim(1),KIND=i4), re_cld, &
                  re_pre_cld_2D(1:clp_dim(1),ialt), &
                  one, local_cld(1), &
                  lay_pre_cld_1D(ialt), status)
          END DO

          ! --------------------
          ! Compute sky radiance
          ! --------------------
          DO isza = 1, sza_dim(1)
             DO ivza = 1, vza_dim(1)
                ! ---------
                ! Clear sky
                ! ---------
                DO isrf = 1, srf_dim(1)
                   ! For the intensity the TOMRAD formula works perfect so we need no albedo loop
                   Inte_clear_3D(srf_dim(1)+1-isrf,vza_dim(1)+1-ivza,sza_dim(1)+1-isza) = &
                        REAL(lut_I0_clr(toms_idx,isrf,ivza,isza), KIND=r8)  +  &
                        REAL(lut_I1_clr(toms_idx,isrf,ivza,isza), KIND=r8) * cos(local_raa(1))        + &
                        REAL(lut_I2_clr(toms_idx,isrf,ivza,isza), KIND=r8) * cos(2.0_r8*local_raa(1))
                   Inte_clear_3D(srf_dim(1)+1-isrf,vza_dim(1)+1-ivza,sza_dim(1)+1-isza) = &
                        Inte_clear_3D(srf_dim(1)+1-isrf,vza_dim(1)+1-ivza,sza_dim(1)+1-isza) + &
                        REAL(lut_Ir_clr(toms_idx,isrf,ivza,isza), KIND=r8) * local_alb(1) / &
                        ( one_r8(1) - local_alb(1) * REAL(lut_Sb_clr(toms_idx,isrf), KIND=r8) )
                ENDDO
                ! ---------
                ! Cloud sky
                ! ---------
                DO icld = 1, clp_dim(1)
                   Inte_cloud_3D(clp_dim(1)+1-icld,vza_dim(1)+1-ivza,sza_dim(1)+1-isza) = &
                        REAL(lut_I0_cld(toms_idx,icld,ivza,isza), KIND=r8)  +  &
                        REAL(lut_I1_cld(toms_idx,icld,ivza,isza), KIND=r8) * cos(local_raa(1))      + &
                        REAL(lut_I2_cld(toms_idx,icld,ivza,isza), KIND=r8) * cos(2.0_r8*local_raa(1))
                ENDDO ! End cloud loop
                ! ----------------------------------------------
                ! Interpolate clear sky part to surface pressure
                ! ----------------------------------------------
                CALL ezspline_1d_interpolation (INT(srf_dim(1),KIND=i4), re_srf, &
                     Inte_clear_3D(1:srf_dim(1),vza_dim(1)+1-ivza,sza_dim(1)+1-isza), &
                     one, local_srf(1), &
                     Inte_clear_2D(vza_dim(1)+1-ivza,sza_dim(1)+1-isza), status)
                ! --------------------------------------------
                ! Interpolate cloud sky part to cloud pressure
                ! --------------------------------------------
                CALL ezspline_1d_interpolation (INT(clp_dim(1),KIND=i4), re_cld, &
                     Inte_cloud_3D(1:clp_dim(1),vza_dim(1)+1-ivza,sza_dim(1)+1-isza), &
                     one, local_cld(1), &
                     Inte_cloud_2D(vza_dim(1)+1-ivza,sza_dim(1)+1-isza), status)
             ENDDO ! End VZA loop
             ! ------------------
             ! Interpolate to VZA
             ! ------------------
             CALL ezspline_1d_interpolation (INT(vza_dim(1),KIND=i4), re_vza, &
                  Inte_clear_2D(1:vza_dim(1),sza_dim(1)+1-isza), &
                  one, local_vza(1), &
                  Inte_clear_1D(sza_dim(1)+1-isza), status)
             CALL ezspline_1d_interpolation (INT(vza_dim(1),KIND=i4), re_vza, &
                  Inte_cloud_2D(1:vza_dim(1),sza_dim(1)+1-isza), &
                  one, local_vza(1), &
                  Inte_cloud_1D(sza_dim(1)+1-isza), status)
          ENDDO ! End SZA loop
          Radiance_clr(1) = 0.0_r8
          Radiance_cld(1) = 0.0_r8
          ! ------------------
          ! Interpolate to SZA
          ! ------------------
          CALL ezspline_1d_interpolation (INT(sza_dim(1),KIND=i4), re_sza, &
               Inte_clear_1D(1:sza_dim(1)), &
               one, local_sza(1), &
               Radiance_clr(1), status)
          CALL ezspline_1d_interpolation (INT(sza_dim(1),KIND=i4), re_sza, &
               Inte_cloud_1D(1:sza_dim(1)), &
               one, local_sza(1), &
               Radiance_cld(1), status)
          ! ---------------------------------------
          ! Working out the cloud radiance fraction
          ! See note below, Boersma et al. 2011 and
          ! Martin et al. 2003 (crf used below)
          ! ---------------------------------------                 
          crf = 0.0_r8
          crf = local_cfr(1) * Radiance_cld(1) / &
               ( local_cfr(1) * Radiance_cld(1) + &
               (1.0_r8 - local_cfr(1)) * Radiance_clr(1) )

          ! ---------------------------
          ! Work out scattering weights
          ! ---------------------------
          DO ialt = 1, alt_lay_dim(2)
             DO isza = 1, sza_dim(1)
                DO ivza = 1, vza_dim(1)
                   ! ---------
                   ! Clear sky
                   ! ---------
                   DO isrf = 1, srf_dim(1)
                      ! ------------------------------
                      ! Albedo loop only for clear sky
                      ! ------------------------------
                      DO ialb = 1, alb_dim(1)
                         sw_clear_5D(srf_dim(1)+1-isrf,ialt,ialb,vza_dim(1)+1-ivza,sza_dim(1)+1-isza) = &
                              REAL( lut_dI0_clr(toms_idx,isrf,ialt,ialb,ivza,isza), KIND=r8 ) + &
                              REAL( lut_dI1_clr(toms_idx,isrf,ialt,ialb,ivza,isza), KIND=r8 ) * cos( local_raa(1) ) + &
                              REAL( lut_dI2_clr(toms_idx,isrf,ialt,ialb,ivza,isza), KIND=r8 ) * cos(2.0_r8 * local_raa(1) )
                      END DO
                      ! ---------------------------
                      ! Interpolate to local albedo
                      ! ---------------------------
                      sw_clear_4D(srf_dim(1)+1-isrf,ialt,vza_dim(1)+1-ivza,sza_dim(1)+1-isza) = &
                           linInterpol(INT(alb_dim(1),KIND=i4), re_alb, &
                           sw_clear_5D(srf_dim(1)+1-isrf,ialt,1:alb_dim(1),vza_dim(1)+1-ivza,sza_dim(1)+1-isza), &
                           local_alb(1), status=status)
                   END DO
                   ! --------------------------------------
                   ! Cloud sky
                   ! --------------------------------------
                   DO icld = 1, clp_dim(1)
                      sw_cloud_4D(ialt,clp_dim(1)+1-icld,vza_dim(1)+1-ivza,sza_dim(1)+1-isza) = &
                           REAL( lut_dI0_cld(toms_idx,ialt,icld,ivza,isza), KIND=r8 ) + &
                           REAL( lut_dI1_cld(toms_idx,ialt,icld,ivza,isza), KIND=r8 ) * cos( local_raa(1) ) + &
                           REAL( lut_dI2_cld(toms_idx,ialt,icld,ivza,isza), KIND=r8 ) * cos(2.0_r8 * local_raa(1) )
                   END DO
                   ! -------------------------------------
                   ! Interpolate to local surface pressure
                   ! If cloud pressure is greater than
                   ! surface pressure then cloud SW are to
                   ! be cero.
                   ! -------------------------------------
                CALL ezspline_1d_interpolation (INT(srf_dim(1),KIND=i4), re_srf, &
                     sw_clear_4D(1:srf_dim(1),ialt,vza_dim(1)+1-ivza,sza_dim(1)+1-isza), &
                     one, local_srf(1), &
                     sw_clear_3D(ialt,vza_dim(1)+1-ivza,sza_dim(1)+1-isza), status)
                CALL ezspline_1d_interpolation (INT(clp_dim(1),KIND=i4), re_cld, &
                     sw_cloud_4D(ialt,1:clp_dim(1),vza_dim(1)+1-ivza,sza_dim(1)+1-isza), &
                     one, local_cld(1), &
                     sw_cloud_3D(ialt,vza_dim(1)+1-ivza,sza_dim(1)+1-isza), status)
                END DO ! End VZA loop
                ! ------------------------
                ! Interpolate to local VZA
                ! ------------------------
                CALL ezspline_1d_interpolation (INT(vza_dim(1),KIND=i4), re_vza, &
                     sw_clear_3D(ialt,1:vza_dim(1),sza_dim(1)+1-isza), &
                     one, local_vza(1), &
                     sw_clear_2D(ialt,sza_dim(1)+1-isza), status)
                CALL ezspline_1d_interpolation (INT(vza_dim(1),KIND=i4), re_vza, &
                     sw_cloud_3D(ialt,1:vza_dim(1),sza_dim(1)+1-isza), &
                     one, local_vza(1), &
                     sw_cloud_2D(ialt,sza_dim(1)+1-isza), status)
             END DO ! End SZA angle 
             ! ------------------------
             ! Interpolate to local SZA
             ! ------------------------
             CALL ezspline_1d_interpolation (INT(sza_dim(1),KIND=i4), re_sza, &
                  sw_clear_2D(ialt,1:sza_dim(1)), &
                  one, local_sza(1), &
                  sw_clear_1D(ialt), status)
             CALL ezspline_1d_interpolation (INT(sza_dim(1),KIND=i4), re_sza, &
                  sw_cloud_2D(ialt,1:sza_dim(1)), &
                  one, local_sza(1), &
                  sw_cloud_1D(ialt), status)
          END DO ! End alt layer loop

          ! --------------------------------
          ! Finally interpolate to local_pre
          ! work out from local_srf(1) and
          ! GEOS5 Ap and Bp (or just from
          ! interpolating re_pre to
          ! local_srf the cloud SW
          ! --------------------------------
          DO ialt = 1, ngeos5-1
             IF ( local_lay_pre(ialt) .GT. local_cld(1) ) THEN
                sw_cloud(ialt) = 0.0_r8
             ELSE
                CALL ezspline_1d_interpolation (INT(alt_lay_dim(2),KIND=i4), lay_pre_cld_1D(1:alt_lay_dim(2)), &
                     sw_cloud_1D(1:alt_lay_dim(2)), &
                     one, local_lay_pre(ialt), &
                     sw_cloud(ialt), status)
             END IF
          END DO

          ! ---------------------------------------------------------------------------------
          ! Boersma et al. 2011 AMT, 4, 2011
          ! Cloud radiance fraction: Crf= Cfr * Icl / Ir
          !  We define Icl = From the Vlidort calculation, see above
          !            Icr = From the Vlidort calculation, see above
          !             Ir = Cfr * Icl + (1 - Cfr) * Icr (Total pixel radiance)
          ! 
          ! Now the scattering weights become w = crf * scatt_cloud + (1 - crf) * scatt_clear
          !  We add the scattweights calculated in the previous wavelengths.
          ! ---------------------------------------------------------------------------------
          DO ialt = 1, ngeos5-1 
             scattw(ixtrack,itime,ialt) = crf * sw_cloud(ialt) + (1.0_r8 - crf) * sw_clear_1d(ialt)
          END DO
          !  Set non-physical entries to zero.
          WHERE ( scattw(ixtrack,itime,1:ngeos5-1) < 0.0_r8 )
             scattw(ixtrack,itime,1:ngeos5-1) = 0.0_r8
          END WHERE
       END DO ! End loop xtrack
       IF ( verb_thresh_lev .GE. vb_lev_screen ) WRITE(*,*) 'Scattering weights line', itime

    END DO ! End loop lines
    
  END SUBROUTINE COMPUTE_SCATT

  SUBROUTINE compute_amf ( nt, nx, CmETA, climatology, cli_heights, cli_temperature, cli_psurface, &
                           scattw, saoamf, amfdiag, errstat)

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                INTENT(IN) :: nt, nx, CmETA
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA), INTENT(IN) :: climatology, cli_heights, &
                                                                    cli_temperature, scattw
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1),       INTENT(IN) :: cli_psurface
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1),       INTENT(IN) :: amfdiag

    ! -----------------------------
    ! Output and modified variables
    ! -----------------------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1),       INTENT (INOUT) :: saoamf
    INTEGER (KIND=i4),                                INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                      :: locerrstat, n, n1, ixtrack, itimes, itest
   
    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=11), PARAMETER :: modulename = 'compute_amf'

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    DO itimes  = 0, nt-1 ! Swath lines loop
       DO ixtrack = 1, nx ! Xtrack pixel loop

          ! ------------------------------------------------
          ! Set all amf with amfdiag < 0 to r8_missval
          ! ------------------------------------------------
          IF ( amfdiag(ixtrack,itimes) .LT. 0 ) THEN
             saoamf(ixtrack,itimes) = r8_missval
             CYCLE
          ENDIF

          ! -------------------------
          ! Finally work out the AMFs
          ! -------------------------
          saoamf(ixtrack,itimes) = SUM(scattw(ixtrack, itimes, 1:CmETA) * &
                                        climatology(ixtrack,itimes,1:CmETA))     / &
                                   SUM(climatology(ixtrack,itimes,1:CmETA))
       END DO ! Finish xtrack pixel loop
    END DO ! Finish 
    
  END SUBROUTINE compute_amf

  SUBROUTINE write_albedo_he5(albedo, nt, nx, errstat)

    ! ==================================================================
    ! This routines writes the albedos obtained from the OMLER climatolo
    ! gy to the output file.
    ! ==================================================================
    USE OMSAO_omidata_module,   ONLY: n_roff_dig

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                         INTENT (IN) :: nt, nx
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1), INTENT (IN) :: albedo

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4),                         INTENT (INOUT) :: errstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=16), PARAMETER :: modulename = 'write_albedo_he5'
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                          :: locerrstat
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: colloc

    locerrstat = pge_errstat_ok

    he5_start_2d  = (/ 0, 0 /)
    he5_stride_2d = (/ 1, 1 /)
    he5_edge_2d   = (/ nx, nt /)      


    colloc = albedo
    CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                               TRIM(ADJUSTL(albedo_field)),             &
                               he5_start_2d, he5_stride_2d, he5_edge_2d,&
                               colloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )
    
  END SUBROUTINE write_albedo_he5

  SUBROUTINE write_climatology_he5(climatology, cli_psurface, nt, nx, nl, errstat)

    ! ===============================================================
    ! This routines writes the Target Gas Profiles from the GEOS-Chem
    ! climatology to the output file.
    ! ===============================================================
    USE OMSAO_omidata_module,   ONLY: n_roff_dig

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                              INTENT (IN) :: nt, nx, nl
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1)     , INTENT (IN) :: cli_psurface
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1,1:nl), INTENT (IN) :: climatology

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4),                         INTENT (INOUT) :: errstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=24), PARAMETER :: modulename = 'write_climatology_he5' ! JED fix
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                               :: locerrstat
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1)      :: colloc_2d
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,1:nl) :: colloc_3d

    locerrstat = pge_errstat_ok

    he5_start_2d  = (/ 0, 0   /)
    he5_stride_2d = (/ 1, 1   /)
    he5_edge_2d   = (/ nx, nt /)      

    colloc_2d = cli_psurface
    CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc_2d(1:nx,0:nt-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                               TRIM(ADJUSTL(surfacepre_field)),         &
                               he5_start_2d, he5_stride_2d, he5_edge_2d,&
                               colloc_2d(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

    he5_start_3d  = (/ 0, 0, 0 /)
    he5_stride_3d = (/ 1, 1, 1 /)
    he5_edge_3d   = (/ nx, nt, nl /)      

    colloc_3d = climatology
    CALL roundoff_3darr_r8 ( n_roff_dig, nx, nt, nl, colloc_3d(1:nx,0:nt-1,1:nl) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                               TRIM(ADJUSTL(gasprofile_field)),         &
                               he5_start_3d, he5_stride_3d, he5_edge_3d,&
                               colloc_3d(1:nx,0:nt-1,1:nl) )
    errstat = MAX ( errstat, locerrstat )  
    
  END SUBROUTINE write_climatology_he5

  SUBROUTINE write_scatt_he5(scattw, nt, nx, nl, errstat)

    ! ===============================================================
    ! This routines writes the scattering weigths to the output file.
    ! ===============================================================
    USE OMSAO_omidata_module,   ONLY: n_roff_dig

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                              INTENT (IN) :: nt, nx, nl
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1,1:nl), INTENT (IN) :: scattw !, akernels

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4),                         INTENT (INOUT) :: errstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=16), PARAMETER :: modulename = 'write_scatt_he5'
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                               :: locerrstat
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,1:nl) :: colloc

    locerrstat = pge_errstat_ok

    he5_start_3d  = (/ 0, 0, 0 /)
    he5_stride_3d = (/ 1, 1, 1 /)
    he5_edge_3d   = (/ nx, nt, nl /)      


    colloc = scattw
    CALL roundoff_3darr_r8 ( n_roff_dig, nx, nt, nl, colloc(1:nx,0:nt-1,1:nl) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                               TRIM(ADJUSTL(scaweights_field)),         &
                               he5_start_3d, he5_stride_3d, he5_edge_3d,&
                               colloc(1:nx,0:nt-1,1:nl) )

    errstat = MAX ( errstat, locerrstat )
    
  END SUBROUTINE write_scatt_he5


SUBROUTINE he5_amf_write ( &
         pge_idx, nx, nt, saocol, saodco, amfmol, amfgeo, amfdiag, &
         amfcfr, amfctp, errstat )

  USE OMSAO_precision_module, ONLY: i2, i4, r8
  USE OMSAO_he5_module
  USE OMSAO_errstat_module
  USE OMSAO_omidata_module,   ONLY: n_roff_dig
  USE OMSAO_indices_module,   ONLY: pge_hcho_idx, pge_gly_idx, pge_bro_idx, pge_h2o_idx

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=13), PARAMETER :: modulename = 'he5_write_amf'


  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                          INTENT (IN) :: pge_idx, nx, nt
  REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: saocol, saodco
  REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfmol, amfgeo
  REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfcfr, amfctp
  INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfdiag

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4)                          :: locerrstat
  REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1) :: amfloc
  REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: colloc

  
  ! -------------------------------------------------------------
  ! Air mass factor plus diagnostic.
  ! -------------------------------------------------------------
  ! As yet, only OMBRO and OMHCHO have true, non-geometric AMFs.
  ! But we try to have symmetric data fields as much as possible,
  ! hence the presence of the "molecule specific" AMF and its
  ! diagnostic for all PGEs. Non-OMBRO and -OMHCHO PGEs carry a
  ! geometric AMF here.
  ! 
  ! For completeness, the geometric AMF is added.
  ! -------------------------------------------------------------

  locerrstat = pge_errstat_ok

  he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nx, nt /)

  ! ----------------------------------------
  ! (1) AMF diagnostic. No rounding required
  ! ----------------------------------------
  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfdiag_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, amfdiag(1:nx,0:nt-1) )
  errstat = MAX ( errstat, locerrstat )

  ! -----------------
  ! (2) Geometric AMF
  ! -----------------
  amfloc = REAL ( amfgeo, KIND=r4 )
  CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfgeo_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
  errstat = MAX ( errstat, locerrstat )

  ! -----------------
  ! (3) Molecular AMF
  ! -----------------
  amfloc = REAL ( amfmol, KIND=r4 )
  CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfmol_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
  errstat = MAX ( errstat, locerrstat )

  ! ----------------------------------------------------------
  ! (4) OMHCHO, OMCHOCHO only: AMF cloud fraction and pressure
  ! ----------------------------------------------------------
  IF ( pge_idx == pge_hcho_idx .OR. pge_idx == pge_gly_idx .OR. pge_idx==pge_h2o_idx) THEN
     amfloc = REAL ( amfcfr, KIND=r4 )
     CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
     locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfcfr_field)), &
          he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
     errstat = MAX ( errstat, locerrstat )
  
     amfloc = REAL ( amfctp, KIND=r4 )
     CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
     locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfctp_field)), &
          he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
     errstat = MAX ( errstat, locerrstat )
  END IF

  ! -----------------------------------------------------------------------
  ! (5) All PGEs: Output of columns and column uncertainties. For some PGEs
  !     (e.g., OMBRO, OMHCHO, OMCHOCHO) those have been adjusted by the AMF,
  !     but we have as yet to perform the rounding for any of them.
  ! -----------------------------------------------------------------------
  colloc = saocol
  CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(col_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
  errstat = MAX ( errstat, locerrstat )
  
  colloc = saodco
  CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(dcol_field)), &
       he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE he5_amf_write

SUBROUTINE read_lookup_table (errstat)

    ! ====================================================
    ! This subroutine reads in the VLIDORT calculations to
    ! compute the Scattering Weights.
    ! ====================================================

    IMPLICIT NONE    

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: hdferr

    INTEGER(HID_T) :: input_file_id                                  ! File identifier
    INTEGER(HID_T) :: alb_did, srf_did, clp_did, sza_did, toz_did, vza_did, wav_did, & ! Dataset identifiers
         I0_clr_did, I1_clr_did, I2_clr_did, Ir_clr_did, Sb_clr_did,        &
         I0_cld_did, I1_cld_did, I2_cld_did, Ir_cld_did, Sb_cld_did,        &
         dI0_clr_did, dI1_clr_did, dI2_clr_did, & 
         dI0_cld_did, dI1_cld_did, dI2_cld_did, & 
         air_did, alt_lev_did, alt_lay_did, pre_lev_did, pre_lay_did, ozo_did, tem_did, &
         dspace, toms_datatype_id, pre_lev_cld_did, pre_lay_cld_did

    INTEGER(SIZE_T)                :: size
    LOGICAL, SAVE :: h5inited = .FALSE.
  
    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=26), PARAMETER :: modulename = 'read_lookup_table' 

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    errstat = pge_errstat_ok

    ! ---------------------------------
    ! Initialize hdf5 FORTRAN Interface
    ! ---------------------------------
    if (.NOT.h5inited) then
      CALL h5open_f(hdferr)
      h5inited = .TRUE.
    endif 

    ! ------------------
    ! Dataset data types
    ! ------------------
    size = 5
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER, toms_datatype_id, hdferr)
    CALL h5tset_size_f(toms_datatype_id, size, hdferr)

    ! ******************************************
    ! Find out the dimensions of the input file:
    !  # of pressure levels
    !  # of SZA
    !  # of Ozone profiles
    !  # of VZA
    !  # of wavelenghts
    !  # of altitude levels  
    ! ******************************************
    ! -------------------
    ! Opening input TABLE
    ! -------------------
    CALL h5fopen_f(TRIM(ADJUSTL(OMSAO_wfamf_table_filename)), H5F_ACC_RDONLY_F, &
                   input_file_id, hdferr)
    IF (hdferr .eq. -1) THEN
       WRITE(*,100) 'ERROR: Opening '//TRIM(ADJUSTL(OMSAO_wfamf_table_filename))
    END IF

    ! --------------------------------------------------------------------------
    ! Open grid, intensity, scattering weights and profile datasets
    ! --------------------------------------------------------------------------
    CALL h5dopen_f(input_file_id,'/Grid/Albedo', alb_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/Surface Pressure', srf_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/Cloud Pressure', clp_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/SZA', sza_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/TOMS', toz_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/VZA', vza_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/Wavelength', wav_did,hdferr)
    
    CALL h5dopen_f(input_file_id,'/Intensity/Clear Sky/I0', I0_clr_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Clear Sky/I1', I1_clr_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Clear Sky/I2', I2_clr_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Clear Sky/Ir', Ir_clr_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Clear Sky/Sb', Sb_clr_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Cloud Sky/I0', I0_cld_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Cloud Sky/I1', I1_cld_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Cloud Sky/I2', I2_cld_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Cloud Sky/Ir', Ir_cld_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Cloud Sky/Sb', Sb_cld_did, hdferr)
    
    
    CALL h5dopen_f(input_file_id,'/Profiles/Air Column Layer', air_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Altitude Level', alt_lev_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Altitude Layer', alt_lay_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Pressure Level', pre_lev_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Pressure Layer', pre_lay_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Cloud pressure Level', pre_lev_cld_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Cloud pressure Layer', pre_lay_cld_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Ozone Column Layer', ozo_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Temperature Level', tem_did, hdferr)

    CALL h5dopen_f(input_file_id,'/Scattering Weights/Clear Sky/dI0', dI0_clr_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Scattering Weights/Clear Sky/dI1', dI1_clr_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Scattering Weights/Clear Sky/dI2', dI2_clr_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Scattering Weights/Cloud Sky/dI0', dI0_cld_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Scattering Weights/Cloud Sky/dI1', dI1_cld_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Scattering Weights/Cloud Sky/dI2', dI2_cld_did, hdferr)
    
    ! -----------------------
    ! Find out the dimensions
    ! -----------------------
    CALL h5dget_space_f(alb_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, alb_dim, alb_maxdim, hdferr)
    CALL h5dget_space_f(srf_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, srf_dim, srf_maxdim, hdferr)
    CALL h5dget_space_f(clp_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, clp_dim, clp_maxdim, hdferr)
    CALL h5dget_space_f(sza_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, sza_dim, sza_maxdim, hdferr)
    CALL h5dget_space_f(toz_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, toz_dim, toz_maxdim, hdferr)
    CALL h5dget_space_f(vza_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, vza_dim, vza_maxdim, hdferr)
    CALL h5dget_space_f(wav_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, wav_dim, wav_maxdim, hdferr)

    CALL h5dget_space_f(I0_clr_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, I0_clr_dim, I0_clr_maxdim, hdferr)
    CALL h5dget_space_f(I0_cld_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, I0_cld_dim, I0_cld_maxdim, hdferr)

    CALL h5dget_space_f(Sb_clr_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, Sb_clr_dim, Sb_clr_maxdim, hdferr)
    CALL h5dget_space_f(Sb_cld_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, Sb_cld_dim, Sb_cld_maxdim, hdferr)

    CALL h5dget_space_f(dI0_clr_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, dI0_clr_dim, dI0_clr_maxdim, hdferr)
    CALL h5dget_space_f(dI0_cld_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, dI0_cld_dim, dI0_cld_maxdim, hdferr)

    CALL h5dget_space_f(alt_lev_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, alt_lev_dim, alt_lev_maxdim, hdferr)
    CALL h5dget_space_f(alt_lay_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, alt_lay_dim, alt_lay_maxdim, hdferr)

    CALL h5dget_space_f(pre_lev_cld_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, alt_lev_cld_dim, alt_lev_cld_maxdim, hdferr)
    CALL h5dget_space_f(pre_lay_cld_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, alt_lay_cld_dim, alt_lay_cld_maxdim, hdferr)

    CALL h5dget_space_f(air_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, air_dim, air_maxdim, hdferr)
    CALL h5dget_space_f(tem_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, tem_dim, tem_maxdim, hdferr)
    
    ! ---------------------------------------------------------------
    ! Allocate & initialize variables now that we have the dimensions
    ! ---------------------------------------------------------------
    CALL vlidort_allocate('a', INT(toz_dim(1)),INT(srf_dim(1)),INT(clp_dim(1)),INT(alb_dim(1)), &
         INT(sza_dim(1)),INT(vza_dim(1)),INT(wav_dim(1)), &
         INT(alt_lay_dim(2)),INT(alt_lev_dim(2)),errstat)
    
    ! ----------------------------------------------------
    ! Read from the h5 file all these small size variables
    ! ----------------------------------------------------
    CALL h5dread_f(alb_did, H5T_NATIVE_REAL,  lut_alb(1:alb_dim(1)),  alb_dim, hdferr)
    CALL h5dread_f(srf_did, H5T_NATIVE_REAL,  lut_srf(1:srf_dim(1)),  srf_dim, hdferr)
    CALL h5dread_f(clp_did, H5T_NATIVE_REAL,  lut_clp(1:clp_dim(1)),  clp_dim, hdferr)
    CALL h5dread_f(sza_did, H5T_NATIVE_REAL,  lut_sza(1:sza_dim(1)),  sza_dim, hdferr)
    CALL h5dread_f(toz_did, toms_datatype_id, lut_toms(1:toz_dim(1)), toz_dim, hdferr)
    CALL h5dread_f(vza_did, H5T_NATIVE_REAL,  lut_vza(1:vza_dim(1)),  vza_dim, hdferr)
    CALL h5dread_f(wav_did, H5T_NATIVE_REAL,  lut_wavelength(1:wav_dim(1)),  wav_dim, hdferr)
    
    CALL h5dread_f(I0_clr_did, H5T_NATIVE_REAL, &
         lut_I0_clr(1:toz_dim(1),1:srf_dim(1),1:vza_dim(1),1:sza_dim(1)), I0_clr_dim, hdferr)
    CALL h5dread_f(I1_clr_did, H5T_NATIVE_REAL, &
         lut_I1_clr(1:toz_dim(1),1:srf_dim(1),1:vza_dim(1),1:sza_dim(1)), I0_clr_dim, hdferr)
    CALL h5dread_f(I2_clr_did, H5T_NATIVE_REAL, &
         lut_I2_clr(1:toz_dim(1),1:srf_dim(1),1:vza_dim(1),1:sza_dim(1)), I0_clr_dim, hdferr)
    CALL h5dread_f(Ir_clr_did, H5T_NATIVE_REAL, &
         lut_Ir_clr(1:toz_dim(1),1:srf_dim(1),1:vza_dim(1),1:sza_dim(1)), I0_clr_dim, hdferr)
    CALL h5dread_f(Sb_clr_did, H5T_NATIVE_REAL, &
         lut_Sb_clr(1:toz_dim(1),1:srf_dim(1)), Sb_clr_dim, hdferr)

    CALL h5dread_f(I0_cld_did, H5T_NATIVE_REAL, &
         lut_I0_cld(1:toz_dim(1),1:clp_dim(1),1:vza_dim(1),1:sza_dim(1)), I0_cld_dim, hdferr)
    CALL h5dread_f(I1_cld_did, H5T_NATIVE_REAL, &
         lut_I1_cld(1:toz_dim(1),1:clp_dim(1),1:vza_dim(1),1:sza_dim(1)), I0_cld_dim, hdferr)
    CALL h5dread_f(I2_cld_did, H5T_NATIVE_REAL, &
         lut_I2_cld(1:toz_dim(1),1:clp_dim(1),1:vza_dim(1),1:sza_dim(1)), I0_cld_dim, hdferr)
    CALL h5dread_f(Ir_cld_did, H5T_NATIVE_REAL, &
         lut_Ir_cld(1:toz_dim(1),1:clp_dim(1),1:vza_dim(1),1:sza_dim(1)), I0_cld_dim, hdferr)
    CALL h5dread_f(Sb_cld_did, H5T_NATIVE_REAL, &
         lut_Sb_cld(1:toz_dim(1),1:clp_dim(1)), Sb_cld_dim, hdferr)

    CALL h5dread_f(dI0_clr_did, H5T_NATIVE_REAL,  &
         lut_dI0_clr(1:toz_dim(1),1:srf_dim(1),1:alt_lay_dim(2),1:alb_dim(1),1:vza_dim(1),1:sza_dim(1)), &
         dI0_clr_dim, hdferr)
    CALL h5dread_f(dI1_clr_did, H5T_NATIVE_REAL,  &
         lut_dI1_clr(1:toz_dim(1),1:srf_dim(1),1:alt_lay_dim(2),1:alb_dim(1),1:vza_dim(1),1:sza_dim(1)), &
         dI0_clr_dim, hdferr)
    CALL h5dread_f(dI2_clr_did, H5T_NATIVE_REAL,  &
         lut_dI2_clr(1:toz_dim(1),1:srf_dim(1),1:alt_lay_dim(2),1:alb_dim(1),1:vza_dim(1),1:sza_dim(1)), &
         dI0_clr_dim, hdferr)

    CALL h5dread_f(dI0_cld_did, H5T_NATIVE_REAL,  &
         lut_dI0_cld(1:toz_dim(1),1:alt_lay_dim(2),1:clp_dim(1),1:vza_dim(1),1:sza_dim(1)), &
         dI0_cld_dim, hdferr)
    CALL h5dread_f(dI1_cld_did, H5T_NATIVE_REAL,  &
         lut_dI1_cld(1:toz_dim(1),1:alt_lay_dim(2),1:clp_dim(1),1:vza_dim(1),1:sza_dim(1)), &
         dI0_cld_dim, hdferr)
    CALL h5dread_f(dI2_cld_did, H5T_NATIVE_REAL,  &
         lut_dI2_cld(1:toz_dim(1),1:alt_lay_dim(2),1:clp_dim(1),1:vza_dim(1),1:sza_dim(1)), &
         dI0_cld_dim, hdferr)
    
    CALL h5dread_f(air_did,     H5T_NATIVE_REAL, lut_air(1:toz_dim(1),1:srf_dim(1),1:alt_lay_dim(2)), &
         air_dim, hdferr)
    CALL h5dread_f(alt_lev_did, H5T_NATIVE_REAL, lut_alt_lev(1:srf_dim(1),1:alt_lev_dim(2)), &
         alt_lev_dim, hdferr)
    CALL h5dread_f(alt_lay_did, H5T_NATIVE_REAL, lut_alt_lay(1:srf_dim(1),1:alt_lay_dim(2)), &
         alt_lay_dim, hdferr)
    CALL h5dread_f(pre_lev_did, H5T_NATIVE_REAL, lut_pre_lev(1:srf_dim(1),1:alt_lev_dim(2)), &
         alt_lev_dim, hdferr)
    CALL h5dread_f(pre_lay_did, H5T_NATIVE_REAL, lut_pre_lay(1:srf_dim(1),1:alt_lay_dim(2)), &
         alt_lay_dim, hdferr)
    CALL h5dread_f(pre_lev_cld_did, H5T_NATIVE_REAL, lut_pre_lev_cld(1:clp_dim(1),1:alt_lev_dim(2)), &
         alt_lev_cld_dim, hdferr)
    CALL h5dread_f(pre_lay_cld_did, H5T_NATIVE_REAL, lut_pre_lay_cld(1:clp_dim(1),1:alt_lay_dim(2)), &
         alt_lay_cld_dim, hdferr)
    CALL h5dread_f(ozo_did,     H5T_NATIVE_REAL, lut_ozo(1:toz_dim(1),1:srf_dim(1),1:alt_lay_dim(2)), &
         air_dim, hdferr)
    CALL h5dread_f(tem_did,     H5T_NATIVE_REAL, lut_tem(1:toz_dim(1),1:srf_dim(1),1:alt_lev_dim(2)), &
         tem_dim, hdferr)
    
    ! --------------
    ! Close datasets
    ! --------------    
    CALL h5dclose_f(alb_did, hdferr)
    CALL h5dclose_f(srf_did, hdferr)
    CALL h5dclose_f(clp_did, hdferr)
    CALL h5dclose_f(sza_did, hdferr)
    CALL h5dclose_f(toz_did, hdferr)
    CALL h5dclose_f(vza_did, hdferr)
    CALL h5dclose_f(wav_did, hdferr)
    
    CALL h5dclose_f(I0_clr_did, hdferr); CALL h5dclose_f(I0_cld_did, hdferr)
    CALL h5dclose_f(I1_clr_did, hdferr); CALL h5dclose_f(I1_cld_did, hdferr)
    CALL h5dclose_f(I2_clr_did, hdferr); CALL h5dclose_f(I2_cld_did, hdferr)
    CALL h5dclose_f(Ir_clr_did, hdferr); CALL h5dclose_f(Ir_cld_did, hdferr)
    CALL h5dclose_f(Sb_clr_did, hdferr); CALL h5dclose_f(Sb_cld_did, hdferr)
    
    CALL h5dclose_f(dI0_clr_did, hdferr); CALL h5dclose_f(dI0_cld_did, hdferr)
    CALL h5dclose_f(dI1_clr_did, hdferr); CALL h5dclose_f(dI1_cld_did, hdferr)
    CALL h5dclose_f(dI2_clr_did, hdferr); CALL h5dclose_f(dI2_cld_did, hdferr)
    
    CALL h5dclose_f(air_did, hdferr)
    CALL h5dclose_f(alt_lay_did, hdferr); CALL h5dclose_f(alt_lev_did, hdferr)
    CALL h5dclose_f(pre_lay_did, hdferr); CALL h5dclose_f(pre_lev_did, hdferr)
    CALL h5dclose_f(pre_lay_cld_did, hdferr); CALL h5dclose_f(pre_lev_cld_did, hdferr)
    CALL h5dclose_f(ozo_did, hdferr)
    CALL h5dclose_f(tem_did, hdferr)
    
    ! ----------
    ! Close file
    ! ----------
    CALL h5fclose_f(input_file_id, hdferr)

    errstat = hdferr

100 FORMAT (A)

END SUBROUTINE read_lookup_table

END MODULE OMSAO_wfamf_module

