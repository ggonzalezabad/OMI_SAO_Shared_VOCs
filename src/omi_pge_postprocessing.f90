SUBROUTINE omi_pge_postprocessing ( &
     l1bfile, pge_idx, ntimes, nxtrack, yn_process, xtrange, yn_szoom, n_max_rspec, errstat )

  ! ---------------------------------------------------------
  ! In this subroutine we collect all those computations that
  ! are done "post fitting". These include:
  !
  ! (1) AMF calculation
  ! (2) Fitting statistics
  ! (3) Cross-track destriping
  ! (4) Ground-pixel corner computation
  ! (5) Reference Sector Background Correction for HCHO
  ! ---------------------------------------------------------

  USE OMSAO_precision_module
  USE OMSAO_he5_module
  USE OMSAO_pixelcorner_module
  USE OMSAO_destriping_module
  USE OMSAO_errstat_module
  USE OMSAO_variables_module, ONLY: yn_refseccor
  USE OMSAO_indices_module, ONLY: pge_hcho_idx
  USE OMSAO_Reference_sector_module
  USE OMSAO_radiance_ref_module, ONLY: yn_radiance_reference
  USE OMSAO_wfamf_module, ONLY: amf_calculation_bis, climatology_allocate, Cmlat, Cmlon, CmETA, CmEp1

  IMPLICIT NONE


  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=*),                              INTENT (IN) :: l1bfile  
  INTEGER (KIND=i4),                              INTENT (IN) :: ntimes, nxtrack, n_max_rspec, pge_idx
  INTEGER (KIND=i4), DIMENSION (0:ntimes-1,1:2),  INTENT (IN) :: xtrange
  LOGICAL,           DIMENSION (0:ntimes-1),      INTENT (IN) :: yn_process, yn_szoom

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ----------------
  ! Local variables
  ! ----------------
  ! (1) OMI data
  ! ----------------
  REAL    (KIND=r4), DIMENSION (1:nxtrack,0:ntimes-1) :: lat, lon, sza, vza, thg, extr
  REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1) :: saocol, saodco, saorms, saoamf
  INTEGER (KIND=i2), DIMENSION (1:nxtrack,0:ntimes-1) :: saofcf, saomqf
  INTEGER (KIND=i2), DIMENSION (1:nxtrack,0:ntimes-1) :: glint_flg, snow_ice_flg
  LOGICAL                                             :: yn_write

  ! --------------
  ! Error handling
  ! --------------
  INTEGER (KIND=i4) :: locerrstat

  ! -------------------------
  ! Initialize error variable
  ! -------------------------
  locerrstat = pge_errstat_ok

  ! ----------------------------------------
  ! Read geolocation fields (Lat/Lon/SZA/VZA
  ! ----------------------------------------
  CALL  saopge_geofield_read ( ntimes, nxtrack, lat_field,  lat,  locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, lon_field,  lon,  locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, sza_field,  sza,  locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, vza_field,  vza,  locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, thgt_field, thg,  locerrstat )
  CALL  saopge_geofield_read ( ntimes, nxtrack, extr_field, extr, locerrstat )
 
  ! ----------------------------------------------------
  ! Compute ground pixel corner latitudes and longitudes
  ! ----------------------------------------------------
  CALL compute_pixel_corners ( ntimes, nXtrack, lat, lon, yn_szoom, locerrstat )

  ! ----------------------------------------
  ! Read geolocation fields (Lat/Lon/SZA/VZA
  ! ----------------------------------------
  CALL saopge_columninfo_read (                 &
       ntimes, nxtrack, saocol, saodco, saorms, &
       saoamf, saofcf, locerrstat                 )

  ! ----------------------------------
  ! Read L1b glint and snow/ice flags
  ! ----------------------------------
  CALL omi_read_glint_ice_flags ( &
       l1bfile, nxtrack, ntimes, snow_ice_flg, glint_flg, errstat )

  ! ----------------
  ! Comnpute AMF bis
  ! ----------------
  yn_write = .TRUE.
  CALL amf_calculation_bis (                             &
       pge_idx, ntimes, nxtrack, lat, lon, sza, vza,     &
       snow_ice_flg, glint_flg, xtrange, yn_szoom,       &
       saocol, saodco, saoamf, thg, yn_write, &
       locerrstat              )
 
  ! ----------------------------------
  ! Compute average fitting statistics
  ! ----------------------------------
  CALL compute_fitting_statistics (                      &
       pge_idx, ntimes, nxtrack, xtrange,                &
       saocol, saodco, saorms, saofcf, saomqf, locerrstat )

  ! ---------------------------------------
  ! Apply cross-track destriping correction
  ! ---------------------------------------
  CALL xtrack_destriping (                                    &
       pge_idx, ntimes, nxtrack, yn_process, xtrange,         &
       lat, saocol, saodco, saoamf, saofcf, saomqf, locerrstat )

  ! ---------------------------------------------------------------
  ! Apply Reference Sector Correction; Only for HCHO retrieval !gga
  ! ---------------------------------------------------------------
  IF ((yn_refseccor) .AND. ( pge_idx == pge_hcho_idx ) .AND.  &
      (yn_radiance_reference)) THEN
     CALL Reference_Sector_correction (ntimes, nxtrack, xtrange,     &
          lat, saocol, saodco, saoamf, saomqf, saorms, saofcf,       &
          extr, pge_idx, n_max_rspec, locerrstat)
  ENDIF

  ! --------------------------------
  ! Deallocate Climatology variables
  ! --------------------------------
  CALL climatology_allocate ( "d", Cmlat, Cmlon, CmETA, CmEp1, locerrstat )


  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE omi_pge_postprocessing
