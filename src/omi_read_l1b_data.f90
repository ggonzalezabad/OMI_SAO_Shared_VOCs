SUBROUTINE omi_read_irradiance_data ( ntimes, nxtrack, errstat ) 

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_variables_module,  ONLY: &
       l1b_irrad_filename, l1b_channel, verb_thresh_lev, fit_winwav_lim, &
       fit_winexc_lim
  USE OMSAO_omidata_module,    ONLY: &
       nwavel_max, nxtrack_max, omi_irradiance_swathname, omi_irradiance_spec,        &
       omi_irradiance_qflg, omi_irradiance_prec, omi_irradiance_wavl, omi_nwav_irrad, &
       omi_irradiance_ccdpix, omi_ccdpix_selection, omi_ccdpix_exclusion,             &
       omi_sol_wav_avg
  USE hdfeos4_parameters
  USE L1B_Reader_class
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat
  INTEGER (KIND=i4), INTENT (OUT)   :: ntimes, nxtrack

  ! ---------------
  ! Local variables
  ! ---------------
  TYPE      (L1B_block_type)   :: omi_data_block
  INTEGER   (KIND=i4)          :: &
       nwl, nwavel, imin, imax, i, j, locerrstat, ix, icnt, nwavelcoef

  REAL      (KIND=r4), DIMENSION (nwavel_max)             :: tmp_www, tmp_sss, tmp_ppp
  INTEGER   (KIND=i2), DIMENSION (nwavel_max)             :: tmp_fff
  REAL      (KIND=r4), DIMENSION (nwavel_max,nxtrack_max) :: tmp_spc, tmp_wvl, tmp_prc
  INTEGER   (KIND=i2), DIMENSION (nwavel_max,nxtrack_max) :: tmp_flg

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=24), PARAMETER :: modulename = 'omi_read_irradiance_data' 

  locerrstat = pge_errstat_ok

  ntimes = 0 ; nxtrack = 0 ; nwavel = 0 ; nwavelcoef = 0


  ! ------------------------------------------------------
  ! Open data block structure with default size of 1 lines
  ! ------------------------------------------------------
  locerrstat = L1Br_open ( &
       omi_data_block, l1b_irrad_filename, TRIM(ADJUSTL(omi_irradiance_swathname)) )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_E_READ_L1B_FILE, &
       modulename//f_sep//"L1Br_open failed.", vb_lev_default, errstat )
  IF ( errstat >= pge_errstat_error ) RETURN

  ! ----------------------------------
  ! Obtain irradiance swath dimensions
  ! ----------------------------------
  locerrstat = L1Br_getSWdims ( omi_data_block, &
       NumTimes_k=ntimes, nXtrack_k=nxtrack, nWavel_k=nwavel, nWavelCoef_k=nwavelcoef )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_E_READ_L1B_FILE, &
       modulename//f_sep//"L1Br_getSWdims failed.", vb_lev_default, errstat )
  IF ( errstat >= pge_errstat_error ) RETURN

  ! ----------------------------
  ! Initialize irradiance arrays
  ! ----------------------------
  omi_irradiance_spec (1:nwavel,1:nxtrack) = 0.0_r8
  omi_irradiance_prec (1:nwavel,1:nxtrack) = 0.0_r8
  omi_irradiance_qflg (1:nwavel,1:nxtrack) = 0_i2
  omi_irradiance_wavl (1:nwavel,1:nxtrack) = 0.0_r8

  ! -----------------------------------------------------------------
  ! Read Irradiances from L1b file. Only limit the upper end of the
  ! spectrum because we want to save the pixel numbers and hence need
  ! the wavelengths from the first detector pixel on.
  ! -----------------------------------------------------------------
  locerrstat = L1Br_getSIGline ( omi_data_block, 0, &
       Signal_k            = tmp_spc,   &
       SignalPrecision_k   = tmp_prc,   &
       PixelQualityFlags_k = tmp_flg,   &
       Wavelength_k        = tmp_wvl,   &
       Nwl_k               = nwl                       )

  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_E_READ_L1B_FILE, &
       modulename//f_sep//"L1Br_getSIGline failed.", vb_lev_default, errstat )

  IF ( errstat >= pge_errstat_error ) RETURN

  ! -------------------------------
  ! Reverse arrays for UV-1 channel
  ! -------------------------------
  IF ( l1b_channel == 'UV1' ) THEN
     DO ix = 1, nxtrack
        tmp_www(1:nwl) = tmp_wvl(1:nwl,ix)
        tmp_sss(1:nwl) = tmp_spc(1:nwl,ix)
        tmp_ppp(1:nwl) = tmp_prc(1:nwl,ix)
        tmp_fff(1:nwl) = tmp_flg(1:nwl,ix)
        DO i = 1, nwl
           j = nwl - i + 1
           tmp_wvl(j,ix) = tmp_www(i)
           tmp_spc(j,ix) = tmp_sss(i)
           tmp_prc(j,ix) = tmp_ppp(i)
           tmp_flg(j,ix) = tmp_fff(i)
        END DO
     END DO
  END IF

  ! ----------------------------------------------------
  ! Limit irradiance arrays to fitting window. Check for
  ! strictly ascending wavelengths in the process.
  ! ----------------------------------------------------
  DO ix = 1, nxtrack

     ! -------------------------------------------------------------------------------
     ! Determine the CCD pixel numbers based on the selected wavelength fitting window
     ! -------------------------------------------------------------------------------
     DO j = 1, 3, 2
        CALL array_locate_r4 ( &
             nwl, tmp_wvl(1:nwl,ix), REAL(fit_winwav_lim(j  ),KIND=r4), 'LE', &
             omi_ccdpix_selection(ix,j  ) )
        CALL array_locate_r4 ( &
             nwl, tmp_wvl(1:nwl,ix), REAL(fit_winwav_lim(j+1),KIND=r4), 'GE', &
             omi_ccdpix_selection(ix,j+1) )
     END DO


     imin = omi_ccdpix_selection(ix,1)
     imax = omi_ccdpix_selection(ix,4)

     icnt = imax - imin + 1
     omi_irradiance_wavl(1:icnt,ix) = REAL(tmp_wvl(imin:imax,ix), KIND=r8)
     omi_irradiance_spec(1:icnt,ix) = REAL(tmp_spc(imin:imax,ix), KIND=r8)
     omi_irradiance_prec(1:icnt,ix) = REAL(tmp_prc(imin:imax,ix), KIND=r8)
     omi_irradiance_qflg(1:icnt,ix) = tmp_flg(imin:imax,ix)
     omi_nwav_irrad (ix) = icnt
     omi_sol_wav_avg(ix) = SUM( tmp_wvl(imin:imax,ix) ) / REAL(icnt, KIND=r8)

     ! ------------------------------------------------------------------------------
     ! If any window is excluded, find the corresponding indices. This has to be done
     ! after the array assignements above because we need to know which indices to
     ! exclude from the final arrays, not the complete ones read from the HE4 file.
     ! ------------------------------------------------------------------------------
     omi_ccdpix_exclusion(ix,1:2) = -1
     IF ( MINVAL(fit_winexc_lim(1:2)) > 0.0_r8 ) THEN
        CALL array_locate_r4 ( &
             nwl, tmp_wvl(1:nwl,ix), REAL(fit_winexc_lim(1),KIND=r4), 'GE', omi_ccdpix_exclusion(ix,1) )
        CALL array_locate_r4 ( &
             nwl, tmp_wvl(1:nwl,ix), REAL(fit_winexc_lim(2),KIND=r4), 'LE', omi_ccdpix_exclusion(ix,2) )
     END IF

  END DO

  ! --------------------------
  ! Close data block structure
  ! --------------------------
  locerrstat = L1Br_close ( omi_data_block )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_warning, OMSAO_W_CLOS_L1B_FILE, &
       modulename//f_sep//"L1Br_close failed.", vb_lev_omidebug, errstat )

  RETURN
END SUBROUTINE omi_read_irradiance_data


SUBROUTINE omi_read_radiance_paras ( &
     l1bfile, ntimes, nxtrack, nwavel_ccd, l1bswath, errstat )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen, r4_missval, r8_missval
  USE OMSAO_variables_module,  ONLY : verb_thresh_lev, l1b_channel
  USE OMSAO_omidata_module,    ONLY : EarthSunDistance
  USE OMSAO_errstat_module
  USE hdfeos4_parameters
  USE L1B_Reader_class

  IMPLICIT NONE

  ! --------------
  ! Input Variable
  ! --------------
  CHARACTER (LEN=*), INTENT (IN) :: l1bfile

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4),                             INTENT (INOUT) :: errstat
  INTEGER (KIND=i4),                             INTENT (OUT)   :: ntimes, nxtrack, nwavel_ccd
  CHARACTER (LEN=*),                             INTENT (OUT)   :: l1bswath

  ! ---------------
  ! Local variables
  ! ---------------
  TYPE (L1B_block_type) :: omi_data_block
  INTEGER   (KIND=i4)   :: locerrstat, ntimes_small

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_read_radiance_paras'

  ! ------------------
  ! External Functions
  ! ------------------
  REAL    (KIND=r4) :: L1Bga_EarthSunDistance

  ! --------------------------
  ! Initialize OUTPUT variable
  ! --------------------------
  locerrstat = pge_errstat_ok
  ntimes = 0 ; nxtrack = 0 ; nwavel_ccd = 0

  ! ---------------------------------------------------
  ! Retrieve the swath in the current L1b radiance file
  ! ---------------------------------------------------
  CALL omi_xtract_swathname ( l1bfile, l1b_channel, l1bswath )

  ! ----------------------------------------------------------------------
  ! Open data block called 'omi_data_block' with default size of 100 lines
  ! ----------------------------------------------------------------------
  locerrstat = L1Br_open ( omi_data_block, l1bfile, l1bswath, nL=1 )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_F_OPEN_L1B_FILE, &
       modulename//f_sep//"L1Br_open failed.", vb_lev_default, errstat )
  IF ( errstat >= pge_errstat_error ) RETURN

  ! -------------------------------
  ! Get dimensions of current Swath
  ! -------------------------------
  locerrstat = L1Br_getSWdims ( omi_data_block, &
       NumTimes_k=ntimes, nXtrack_k=nxtrack, NumTimesSmallPixel_k=ntimes_small, &
       nWavel_k=nwavel_ccd )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_E_READ_L1B_FILE, &
       modulename//f_sep//"L1Br_getSWdims failed.", vb_lev_default, errstat )
  IF ( errstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------------
  ! Close data block structure and reopen with larger line number block
  ! -------------------------------------------------------------------
  locerrstat = L1Br_close ( omi_data_block )
  locerrstat = L1Br_open ( omi_data_block, l1bfile, l1bswath, nL=ntimes )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_F_OPEN_L1B_FILE, &
       modulename//f_sep//"L1Br_open failed.", vb_lev_default, errstat )
  IF ( errstat >= pge_errstat_error ) RETURN


  ! --------------------------
  ! Close data block structure
  ! --------------------------
  locerrstat = L1Br_close ( omi_data_block )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_warning, OMSAO_W_CLOS_L1B_FILE, &
       modulename//f_sep//"L1Br_close failed.", vb_lev_omidebug, errstat )

  ! -------------------------------
  ! Get EarthSunDistance
  ! -------------------------------
  EarthSunDistance = L1Bga_EarthSunDistance( l1bfile, l1bswath )


  RETURN
END SUBROUTINE omi_read_radiance_paras


SUBROUTINE omi_read_binning_factor ( &
     l1bfile, l1bswath, ntimes, binfac, yn_szoom, errstat )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen, r4_missval, r8_missval
  USE OMSAO_variables_module,  ONLY : verb_thresh_lev, l1b_channel
  USE OMSAO_omidata_module,    ONLY : global_mode, szoom_mode
  USE OMSAO_errstat_module
  USE hdfeos4_parameters
  USE L1B_Reader_class

  IMPLICIT NONE

  ! --------------
  ! Input Variable
  ! --------------
  CHARACTER (LEN=*), INTENT (IN) :: l1bfile, l1bswath
  INTEGER (KIND=i4), INTENT (IN) :: ntimes

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4),                         INTENT (INOUT) :: errstat
  INTEGER (KIND=i1), DIMENSION (0:ntimes-1), INTENT (OUT)   :: binfac
  LOGICAL,           DIMENSION (0:ntimes-1), INTENT (OUT)   :: yn_szoom

  ! ---------------
  ! Local variables
  ! ---------------
  TYPE (L1B_block_type) :: omi_data_block
  INTEGER   (KIND=i4)   :: locerrstat, iline

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_read_binning_factor'


  ! --------------------------
  ! Initialize OUTPUT variable
  ! --------------------------
  locerrstat = pge_errstat_ok

  ! ----------------------------------------------------------------------
  ! Open data block called 'omi_data_block' with size of nTimes lines
  ! ----------------------------------------------------------------------
  locerrstat = L1Br_open ( omi_data_block, l1bfile, l1bswath, nL=ntimes )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_fatal, OMSAO_F_OPEN_L1B_FILE, &
       modulename//f_sep//"L1Br_open failed.", vb_lev_default, errstat )

  IF ( errstat >= pge_errstat_error ) RETURN

  ! --------------------
  ! Read binning factors
  ! --------------------
  DO iline = 0, ntimes-1
     locerrstat = L1Br_getDATA ( omi_data_block, iline, ImageBinningFactor_k=binfac(iline) )
  END DO

  ! ----------------------------------------------------------------------
  ! Check whether we have a Spatial Zoom granule, in which case we need to
  ! process 60 cross-track positions rather than 30.
  ! ----------------------------------------------------------------------
  IF ( ( INDEX(l1bfile, 'OML1BRUZ') > 0 ) .OR. &
       ( INDEX(l1bfile, 'OML1BRVZ') > 0 ) )    &
       binfac(0:ntimes-1) = global_mode

  ! --------------------------
  ! Close data block structure
  ! --------------------------
  locerrstat = L1Br_close ( omi_data_block )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_warning, OMSAO_W_CLOS_L1B_FILE, &
       modulename//f_sep//"L1Br_close failed.", vb_lev_omidebug, errstat )

  ! ------------------------------------------------------------------------------
  ! Check for GLOBAL and SPATIAL ZOOM mode and set up arrays for index adjustment.
  ! ------------------------------------------------------------------------------
  WHERE ( binfac(0:ntimes-1) == szoom_mode )
     yn_szoom (0:nTimes-1) = .TRUE.
  ELSEWHERE
     yn_szoom (0:nTimes-1) = .FALSE.
  END WHERE

  RETURN
END SUBROUTINE omi_read_binning_factor


FUNCTION L1Bga_EarthSunDistance( he4filename, swathname ) RESULT( ESdistance )

  USE OMSAO_precision_module
  IMPLICIT NONE

  INCLUDE 'hdf.f90'

  ! ---------------
  ! Input Variables
  ! ---------------
  CHARACTER (LEN=*), INTENT(IN) :: he4filename, swathname

  ! ---------------
  ! Result Variable
  ! ---------------
  REAL (KIND=r4) :: ESdistance

  ! ---------------
  ! Local Variables
  ! ---------------
  INTEGER (KIND=i4) :: swfid, swid, status

  ! ------------------
  ! External Functions
  ! ------------------
  INTEGER (KIND=i4) :: swopen, swattach, swrdattr, swdetach, swclose


  ESdistance = -1.0_r4

  swfid = swopen( he4filename, DFACC_READ )

  IF( swfid /=  -1) THEN
     swid = swattach( swfid, swathname )
     IF( swid /= -1 ) THEN
        status = swrdattr( swid, "EarthSunDistance", ESdistance )
        !IF( status == -1 ) THEN
        !   WRITE(msg,'(A)') "Get swath attribute EarthSunDistance"// &
        !        "failed from "//TRIM(swathname)//","  // &
        !        TRIM( he4filename )  
        !   ierr = OMI_SMF_setmsg( OZT_E_INPUT,  msg, &
        !        "L1Bga_EarthSunDistance", zero )
        !ENDIF
        status = swdetach(swid)
     ENDIF
     status = swclose(swfid)
  ENDIF

END FUNCTION L1Bga_EarthSunDistance


SUBROUTINE omi_read_radiance_lines ( &
     l1bfile, iline, nxtrack, nloop, nwavel_ccd, errstat )

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: pge_o3_idx
  USE OMSAO_parameters_module, ONLY: &
       i1_missval, i2_missval, i2_missval_l1, r4_missval, &
       min_zenith, max_zenith, min_azimuth, max_azimuth, &  ! "non-inclusive"
       min_latitude, max_latitude, min_longitude, max_longitude, &
       earth_radius_avg, eos_aura_avgalt
  USE OMSAO_variables_module,  ONLY: verb_thresh_lev, zatmos, pge_idx
  USE OMSAO_omidata_module,  ONLY: &
       nxtrack_max, omi_radiance_swathname, omi_radiance_spec, omi_radiance_prec,  &
       omi_radiance_wavl, omi_radiance_qflg, omi_height, omi_geoflg, omi_latitude,             &
       omi_longitude, omi_szenith, omi_sazimuth, omi_vzenith, omi_vazimuth,                    &
       omi_razimuth, omi_auraalt, omi_time, omi_nwav_rad, omi_radiance_errstat,                &
       omi_min_specres,                                                                        &
       omi_irradiance_ccdpix, omi_radiance_ccdpix, omi_nwav_irrad, omi_irradiance_wavl,        &
       omi_ccdpix_selection,                                                                   &
       omi_xtrflg_l1b, omi_xtrflg
  USE OMSAO_errstat_module
  USE hdfeos4_parameters
  USE L1B_Reader_class

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=*),                                 INTENT (IN) :: l1bfile
  INTEGER (KIND=i4),                                 INTENT (IN) :: iline, nloop, nxtrack, nwavel_ccd

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)                           :: &
       iloop, blockline, nwl, imin, imax, locerrstat, ix, icnt
  REAL      (KIND=r4), DIMENSION (nxtrack)      :: tmp_sazm, tmp_vazm
  REAL      (KIND=r4), DIMENSION (nwavel_ccd,nxtrack) :: tmp_wvl, tmp_spc, tmp_prc
  INTEGER   (KIND=i2), DIMENSION (nwavel_ccd,nxtrack) :: tmp_flg

  TYPE (L1B_block_type)  :: omi_data_block

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_read_radiance_lines'

  locerrstat = pge_errstat_ok

  ! -----------------------------------------------------------
  ! Open data block called 'omi_data_block' with default size of 100 lines
  ! -----------------------------------------------------------
  locerrstat = L1Br_open ( omi_data_block, l1bfile, omi_radiance_swathname, nL=nloop )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_F_OPEN_L1B_FILE, &
       modulename//f_sep//"L1Br_open failed.", vb_lev_default, errstat )
  IF ( errstat >= pge_errstat_error ) RETURN

  lineloop: DO iloop = 0, nloop-1
     ! -------------------------------------------
     ! The current scan line number we are reading
     ! -------------------------------------------
     blockline = iline + iloop


     ! --------------------------------
     ! Initialize all local data arrays
     ! --------------------------------
     omi_auraalt       (iloop)                        = r4_missval
     omi_latitude      (1:nxtrack,iloop)              = r4_missval
     omi_longitude     (1:nxtrack,iloop)              = r4_missval
     omi_szenith       (1:nxtrack,iloop)              = r4_missval
     omi_sazimuth      (1:nxtrack,iloop)              = r4_missval
     omi_vzenith       (1:nxtrack,iloop)              = r4_missval
     omi_vazimuth      (1:nxtrack,iloop)              = r4_missval
     omi_razimuth      (1:nxtrack,iloop)              = r4_missval
     omi_height        (1:nxtrack,iloop)              = i2_missval
     omi_geoflg        (1:nxtrack,iloop)              = i2_missval
     omi_xtrflg        (1:nxtrack,iloop)              = i2_missval
     omi_xtrflg_l1b    (1:nxtrack,iloop)              = i1_missval
     omi_radiance_spec (1:nwavel_ccd,1:nxtrack,iloop) = REAL ( r4_missval, KIND=r8 )
     omi_radiance_prec (1:nwavel_ccd,1:nxtrack,iloop) = REAL ( r4_missval, KIND=r8 )
     omi_radiance_qflg (1:nwavel_ccd,1:nxtrack,iloop) = i2_missval
     omi_radiance_wavl (1:nwavel_ccd,1:nxtrack,iloop) = REAL ( r4_missval, KIND=r8 )

     ! ---------------------------------------------------
     ! Set error status of current line to O.K. by default
     ! ---------------------------------------------------
     omi_radiance_errstat (iloop) = pge_errstat_ok
     locerrstat = L1Br_getGEOline ( omi_data_block, blockline,        &
          Time_k                    = omi_time(iloop),                &
          SpacecraftAltitude_k      = omi_auraalt(iloop),             &
          Latitude_k                = omi_latitude(1:nxtrack,iloop),  &
          Longitude_k               = omi_longitude(1:nxtrack,iloop), &
          SolarZenithAngle_k        = omi_szenith(1:nxtrack,iloop),   &
          SolarAzimuthAngle_k       = omi_sazimuth(1:nxtrack,iloop),  &
          ViewingZenithAngle_k      = omi_vzenith(1:nxtrack,iloop),   &
          ViewingAzimuthAngle_k     = omi_vazimuth(1:nxtrack,iloop),  &
          TerrainHeight_k           = omi_height(1:nxtrack,iloop),    &
          GroundPixelQualityFlags_k = omi_geoflg(1:nxtrack,iloop),    &
          XTrackQualityFlags_k      = omi_xtrflg_l1b(1:nxtrack,iloop)   ) ! gga row anomaly
     CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_READ_L1B_FILE, &
          modulename//f_sep//"L1Br_getGEOline failed.", vb_lev_default, errstat )

     IF ( locerrstat /= omi_s_success ) omi_radiance_errstat(iloop) = pge_errstat_error

     ! --------------------------------------------------------
     ! Expand XTR Quality Flags to something easily parse-able.
     ! Assign Missing Values where necessary. gga
     ! --------------------------------------------------------
     CALL convert_xtqualflag_info ( &
          nxtrack, omi_xtrflg_l1b(1:nxtrack,iloop), omi_xtrflg(1:nxtrack,iloop) )
     
     ! --------------------------------------------------------------
     ! Check for missing data and reinitialize to PGEs local MissVals
     ! --------------------------------------------------------------
     WHERE ( omi_height(1:nxtrack,iloop) <= i2_missval )
        omi_height(1:nxtrack,iloop) = i2_missval
     ENDWHERE
     WHERE ( ABS(omi_latitude(1:nxtrack,iloop)) > max_latitude )
        omi_latitude(1:nxtrack,iloop) = r4_missval
     ENDWHERE
     WHERE ( ABS(omi_longitude(1:nxtrack,iloop)) > max_longitude )
        omi_longitude(1:nxtrack,iloop) = r4_missval
     ENDWHERE

     ! ---------------------------------------------------------------------
     ! For the Zenith Angles we only correct those values < 0, since we want
     ! to maintain the information of the value of SZA even if it is out of
     ! the bounds required for the computation of AMFs
     ! ---------------------------------------------------------------------
      WHERE ( &
           omi_szenith(1:nxtrack,iloop) < min_zenith )
         omi_szenith(1:nxtrack,iloop) = r4_missval
      ENDWHERE
      WHERE ( omi_vzenith(1:nxtrack,iloop) < min_zenith )
         omi_vzenith(1:nxtrack,iloop) = r4_missval
      ENDWHERE
     WHERE ( &
          omi_sazimuth(1:nxtrack,iloop) < min_azimuth .OR. &
          omi_sazimuth(1:nxtrack,iloop) > max_azimuth        )
        omi_sazimuth(1:nxtrack,iloop) = r4_missval
     ENDWHERE
     WHERE ( &
          omi_vazimuth(1:nxtrack,iloop) < min_azimuth .OR. &
          omi_vazimuth(1:nxtrack,iloop) > max_azimuth        )
        omi_vazimuth(1:nxtrack,iloop) = r4_missval
     ENDWHERE


     ! ----------------------------------------------------------------
     ! Relative azimuth angles (requires some checks and adjustments).
     ! It also has become a Geolocation Field rather than a Data Field.
     ! ----------------------------------------------------------------
     ! (1) Map [-180, +180] to [0, 360]
     tmp_sazm(1:nxtrack) = omi_sazimuth(1:nxtrack,iloop)
     tmp_vazm(1:nxtrack) = omi_vazimuth(1:nxtrack,iloop)
     WHERE ( tmp_sazm(1:nxtrack) /= r4_missval .AND. &
          tmp_sazm(1:nxtrack) < 0.0_r4 .AND. tmp_sazm(1:nxtrack) >= -180.0_r4 )
        tmp_sazm(1:nxtrack) = 360.0_r4 + tmp_sazm(1:nxtrack)
     ENDWHERE
     WHERE ( tmp_vazm(1:nxtrack) /= r4_missval .AND. &
          tmp_vazm(1:nxtrack) < 0.0_r4 .AND. tmp_vazm(1:nxtrack) >= -180.0_r4 )
        tmp_vazm(1:nxtrack) = 360.0_r4 + tmp_vazm(1:nxtrack)
     ENDWHERE
     ! (2) Compute relative azimuth angle (RELATIVE means absolute value)
     WHERE ( tmp_sazm(1:nxtrack) >= 0.0_r4 .AND. tmp_vazm(1:nxtrack) >= 0.0_r4 )
        omi_razimuth(1:nxtrack,iloop) = ABS( tmp_sazm(1:nxtrack) - tmp_vazm(1:nxtrack) )
     ENDWHERE
     WHERE ( omi_razimuth(1:nxtrack,iloop) > 180.0_r4 )
        omi_razimuth(1:nxtrack,iloop) = 360.0_r4 - omi_razimuth(1:nxtrack,iloop)
     ENDWHERE
  
     ! -------------------------------------------------------------------------
     ! Compute zenith and relative azimuth angles at PGE TOA (from control file)
     ! -------------------------------------------------------------------------
     ! --------------------------------------------------------------------------------
     ! NOTE: In the GOME data products all angles are given at the Spacecraft, while
     ! in the OMI data product, angles seem to be given at the surface. This difference
     ! requires a change of the adjustment.
     ! --------------------------------------------------------------------------------
     ! OMI
     IF ( zatmos > 0.0_r8 ) & 
          CALL angle_sat2toa ( &
          earth_radius_avg, 0.0_r4, zatmos, nxtrack, &
          omi_szenith(1:nxtrack,iloop), omi_vzenith(1:nxtrack,iloop), omi_razimuth(1:nxtrack,iloop))
     !! GOME
     !CALL angle_sat2toa ( &
     !     earth_radius_avg, eos_aura_avgalt, zatmos, nxtrack, &
     !     omi_szenith(1:nxtrack,iloop), omi_vzenith(1:nxtrack,iloop), omi_razimuth(1:nxtrack,iloop))

     ! ----------------------------------------------
     ! Get radiances associated with wavelength range
     ! ----------------------------------------------
     IF ( omi_radiance_errstat(iloop) /= pge_errstat_error ) THEN
        locerrstat = L1Br_getSIGline ( omi_data_block, blockline,   &
             Signal_k            = tmp_spc(1:nwavel_ccd,1:nxtrack), &
             SignalPrecision_k   = tmp_prc(1:nwavel_ccd,1:nxtrack), &
             PixelQualityFlags_k = tmp_flg(1:nwavel_ccd,1:nxtrack), &
             Wavelength_k        = tmp_wvl(1:nwavel_ccd,1:nxtrack), &
             Nwl_k               = nwl                               )
        CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_READ_L1B_FILE, &
             modulename//f_sep//"L1Br_getSIGline failed.", vb_lev_default, errstat )

        IF ( locerrstat /= omi_s_success ) omi_radiance_errstat(iloop) = pge_errstat_error

        DO ix = 1, nxtrack
           imin = omi_ccdpix_selection(ix,1)
           imax = omi_ccdpix_selection(ix,4)
           icnt = imax - imin + 1
           omi_radiance_wavl(1:icnt,ix,iloop) = REAL ( tmp_wvl(imin:imax,ix), KIND=r8 )
           omi_radiance_spec(1:icnt,ix,iloop) = REAL ( tmp_spc(imin:imax,ix), KIND=r8 )
           omi_radiance_prec(1:icnt,ix,iloop) = REAL ( tmp_prc(imin:imax,ix), KIND=r8 )
           omi_radiance_qflg(1:icnt,ix,iloop) =        tmp_flg(imin:imax,ix)
           omi_nwav_rad     (       ix,iloop) = icnt
        END DO

     END IF

  END DO lineloop

  ! --------------------------
  ! Close data block structure
  ! --------------------------
  locerrstat = L1Br_close ( omi_data_block )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_warning, OMSAO_W_CLOS_L1B_FILE, &
       modulename//f_sep//"L1Br_close failed.", vb_lev_omidebug, errstat )

  ! --------------------------------------------
  ! Return if we don't have and error free lines
  ! --------------------------------------------
  ! (isn't this a bit pointless at this location in the subroutine? tpk note to himself)
  ! ------------------------------------------------------------------------------------
  IF ( ALL ( omi_radiance_errstat(0:nloop-1) /= pge_errstat_ok ) ) RETURN

  RETURN
END SUBROUTINE omi_read_radiance_lines


SUBROUTINE omi_read_glint_ice_flags ( l1bfile, nx, nt, snow_ice_flg, glint_flg, errstat )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: i2_missval
  USE OMSAO_variables_module,  ONLY: verb_thresh_lev
  USE OMSAO_omidata_module,    ONLY: omi_radiance_swathname
  USE OMSAO_errstat_module
  USE hdfeos4_parameters
  USE L1B_Reader_class

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=*), INTENT (IN) :: l1bfile
  INTEGER (KIND=i4), INTENT (IN) :: nx, nt

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i2), DIMENSION (nx,0:nt-1), INTENT (OUT)   :: glint_flg, snow_ice_flg
  INTEGER (KIND=i4),                        INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)               :: locerrstat, iline
  INTEGER (KIND=i2), DIMENSION (nx) :: geoflg, land_water_flg

  TYPE (L1B_block_type)  :: omi_data_block

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_read_radiance_lines'

  locerrstat = pge_errstat_ok

  ! -----------------------------------------------------------
  ! Open data block called 'omi_data_block' with default size of 100 lines
  ! -----------------------------------------------------------
  locerrstat = L1Br_open ( omi_data_block, l1bfile, omi_radiance_swathname, nL=nt )
  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_F_OPEN_L1B_FILE, &
       modulename//f_sep//"L1Br_open failed.", vb_lev_default, errstat )
  IF ( errstat >= pge_errstat_error ) RETURN

  ! ---------------------------------------------
  ! Read geo flag and convert to snow/glint flags
  ! ------------------------------------------------------------------------
  ! NOTE: GLINT_FLG and SNOW_ICE_FLG are defined on the whole swath because
  !       they are being used in the AMF computation routine. Hence we need
  !       to save the full array. 
  ! ------------------------------------------------------------------------
  DO iline = 0, nt-1

     geoflg = i2_missval

     locerrstat = L1Br_getGEOline (                                      &
          omi_data_block, iline, GroundPixelQualityFlags_k = geoflg(1:nx) )
     CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_READ_L1B_FILE, &
          modulename//f_sep//"L1Br_getGEOline failed.", vb_lev_default, errstat )
  
  ! After change suggested by jed, gga
!!$     CALL convert_gpqualflag_info (   &
!!$          nx,                         &
!!$          geoflg        (1:nx),       &
!!$          land_water_flg(1:nx),       &
!!$          glint_flg     (1:nx,iline), &
!!$          snow_ice_flg  (1:nx,iline)    )

  ! Bits 0-3 are land/water -- not used here
     ! land_water_flg = iand (geoflg, 15_i2)

  ! Bit 4 is glint
     glint_flg(1:nx,iline) = iand (ishft(geoflg, -4), 1_i2)

  ! Bits 8-14 are snow/ice
     snow_ice_flg(1:nx,iline)  = iand (ishft(geoflg, -8), 127_i2)    

  END DO

  ! --------------------------
  ! Close data block structure
  ! --------------------------
  locerrstat = L1Br_close ( omi_data_block )

  CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_warning, OMSAO_W_CLOS_L1B_FILE, &
       modulename//f_sep//"L1Br_close failed.", vb_lev_omidebug, errstat )

  RETURN
END SUBROUTINE omi_read_glint_ice_flags


SUBROUTINE convert_gpqualflag_info ( &
     nxtrack, omi_geoflg, land_water_flg, glint_flg, snow_ice_flg )

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack
  INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (IN) :: omi_geoflg

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (OUT) :: land_water_flg, glint_flg, snow_ice_flg

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i2),                PARAMETER      :: nbyte = 16
  INTEGER (KIND=i2), DIMENSION (7), PARAMETER      :: seven_byte = (/ 1, 2, 4, 8, 16, 32, 64 /)
  INTEGER (KIND=i2)                                :: i
  INTEGER (KIND=i2), DIMENSION (nxtrack)           :: tmp_flg
  INTEGER (KIND=i2), DIMENSION (nxtrack,0:nbyte-1) :: tmp_bytes

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  land_water_flg = 0 ; glint_flg = 0 ; snow_ice_flg = 0

  ! -----------------------------------------------
  ! Save input variable in TMP_FLG for modification
  ! -----------------------------------------------
  tmp_flg(1:nxtrack) = omi_geoflg(1:nxtrack)  ;  tmp_bytes = 0

  ! -------------------------------------------------------------------
  ! CAREFUL: Only 15 flags/positions (0:14) can be returned or else the
  !          conversion will result in a numeric overflow.
  ! -------------------------------------------------------------------
  CALL convert_2bytes_to_16bits ( &
       nbyte-1, nxtrack, tmp_flg(1:nxtrack), tmp_bytes(1:nxtrack,0:nbyte-2) )

  ! ------------------------------
  ! The Glint flag is easy: Byte 4
  ! ------------------------------
  glint_flg(1:nxtrack) = tmp_bytes(1:nxtrack,4)

  ! ------------------------------------------------------------------
  ! Land/Water and Ice require a bit more work. The BIT slices must be
  ! multiplied with the corresponding powers of 2. The sum over this
  ! product is the information we seek.
  ! ------------------------------------------------------------------
  DO i = 1, nxtrack
     land_water_flg(i) = SUM(tmp_bytes(i,0:3 )*seven_byte(1:4))
     snow_ice_flg  (i) = SUM(tmp_bytes(i,8:14)*seven_byte(1:7))
  END DO

  RETURN
END SUBROUTINE convert_gpqualflag_info

SUBROUTINE convert_xtqualflag_info ( nxtrack, omi_xtrflg_l1b, omi_xtrflg )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: i1_missval, i2_missval
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i1), DIMENSION (nxtrack), INTENT (INOUT) :: omi_xtrflg_l1b

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (OUT) :: omi_xtrflg

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4),                       PARAMETER    :: nbit = 16
  INTEGER (KIND=i1), DIMENSION (0:2),      PARAMETER    :: three_bit = (/ 1, 2, 4 /)
  INTEGER (KIND=i2), DIMENSION (3:nbit-1), PARAMETER    :: add_value = (/ 10, 30, 100, 1000, 10000, &
                                                                           0,  0,   0,    0,     0, &
                                                                           0,  0,   0 /)
  INTEGER (KIND=i4)                                     :: i, j
  INTEGER (KIND=i2), DIMENSION (nxtrack)                :: tmp_flg
  INTEGER (KIND=i2), DIMENSION (nxtrack,0:nbit-1)       :: tmp_bits

  ! -----------------------------------------------------------------------
  ! Initialize output quantities. We can't initialize to "I2_MISSVAL" since
  ! we will be recursively adding values to OMI_XTRFLG and hence have to
  ! start out from Zero.
  ! -----------------------------------------------------------------------
  omi_xtrflg = 0_i2

  ! --------------------------------------------------------
  ! Save input variable in TMP_FLG for modification; in that
  ! process, perform a INT8 --> UINT8 conversion.
  ! --------------------------------------------------------
  tmp_bits = 0
  DO i = 1, nXtrack     

     IF ( omi_xtrflg_l1b(i) > -127_i1 .AND. omi_xtrflg_l1b(i) < 0_i1 ) THEN
        tmp_flg(i) = INT ( omi_xtrflg_l1b(i), KIND=i2 ) + 256_i2
     ELSE
        tmp_flg(i) = INT ( omi_xtrflg_l1b(i), KIND=i2 )
     END IF

     ! -----------------------------------------------------------
     ! The code below fails on 32 bit platforms due to "255" being
     ! outside the range of 8bit Integers. On 64 bit platforms it
     ! works perfectly fine.
     ! -----------------------------------------------------------
     !IF ( omi_xtrflg_l1b(i) > -127_i1 .AND. omi_xtrflg_l1b(i) < 0_i1 ) THEN
     !   tmp_flg(i) = INT ( IAND(omi_xtrflg_l1b(i),255), KIND=i2 )
     !ELSE
     !   tmp_flg(i) = INT ( omi_xtrflg_l1b(i), KIND=i2 )
     !END IF

  END DO

  ! gga 
  ! -----------------------------------------------
  ! Save input variable in TMP_FLG for modification
  ! -----------------------------------------------
  !tmp_flg(1:nxtrack) = omi_xtrflg_l1b(1:nxtrack)  ;  tmp_bits = 0
  ! gga
  CALL convert_2bytes_to_16bits ( &
       nbit, nxtrack, tmp_flg(1:nxtrack), tmp_bits(1:nxtrack,0:nbit-1) )
  ! ------------------------------------------------------------------
  ! The Row Anomaly Flags, Bits 0-2
  ! ------------------------------------------------------------------
  DO i = 1, nxtrack
     omi_xtrflg(i) = INT ( SUM(tmp_bits(i,0:2 )*three_bit(0:2)), KIND=i2 )
  END DO

  ! ----------------------------------------------
  ! Add the other bit flags:
  !
  !  Bit  Effect                      Added Value
  !   3   Reserved for future use        10
  !   4   wavelength-shift               30
  !   5   blockage                      100
  !   6   stray sunlight               1000
  !   7   stray earth radiance        10000
  ! ----------------------------------------------
  DO i = 1, nxtrack
     IF ( omi_xtrflg(i) < 0_i2 ) THEN
        omi_xtrflg(i)     = i2_missval
        omi_xtrflg_l1b(i) = i1_missval
     ELSE
        omi_xtrflg(i) = omi_xtrflg(i) + SUM(INT(tmp_bits(i,3:nbit-1),KIND=i2) * add_value(3:nbit-1))
     END IF
  END DO

  RETURN
END SUBROUTINE convert_xtqualflag_info

SUBROUTINE omi_xtract_swathname ( l1bfile, l1bchan, omiswath )

  USE OMSAO_precision_module

  USE OMSAO_variables_module,  ONLY : verb_thresh_lev, l1b_channel
  USE OMSAO_parameters_module, ONLY : maxchlen
  USE hdfeos4_parameters

  IMPLICIT NONE

  ! --------------
  ! Input Variable
  ! --------------
  CHARACTER (LEN=*), INTENT (IN) :: l1bfile
  CHARACTER (LEN=3), INTENT (IN) :: l1bchan

  ! ----------------
  ! Output variables
  ! ----------------
  CHARACTER (LEN=*), INTENT (OUT) :: omiswath

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)      :: is, ie
  INTEGER   (KIND=i4)      :: swfid, nswath, strbufsize, xswath
  CHARACTER (LEN=maxchlen) :: swathlist

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=20), PARAMETER :: modulename = 'omi_xtract_swathname'

  ! ------------------
  ! External Functions
  ! ------------------
  INTEGER (KIND=i4) :: swinqswath

  ! --------------------------
  ! Initialize OUTPUT variable
  ! --------------------------
  omiswath = '?'

  ! ---------------------------------------------------------
  ! Inquire about the swaths in the current L1b radiance file
  ! ---------------------------------------------------------
  swfid  = SWOpen     ( l1bfile, DFACC_READ )
  nswath = SWInqswath ( l1bfile, swathlist, strbufsize )
  xswath = SWClose    ( swfid )

  ! --------------------------------------------------------------------
  ! Extract the swath name we need. Either there is one one swath in the
  ! file (VIS) or there are two (UV-1, UV-2)
  ! --------------------------------------------------------------------
  SELECT CASE ( nswath )
  CASE ( 1 )
     is = 1 ; ie = strbufsize
  CASE ( 2 )
     CALL find_endstring ( strbufsize, swathlist, 1, ie )
     SELECT CASE ( l1bchan )
     CASE ( 'UV1' )
        is = 1 ; ie = ie-1
     CASE ( 'UV2' )
        is = ie+1 ; ie = strbufsize
     CASE DEFAULT
        ! Nothing to do here except to fold
     END SELECT
  CASE DEFAULT
     ! Nothing to do here except to fold.
  END SELECT

  omiswath = TRIM(ADJUSTL(swathlist(is:ie)))

  RETURN
END SUBROUTINE omi_xtract_swathname
