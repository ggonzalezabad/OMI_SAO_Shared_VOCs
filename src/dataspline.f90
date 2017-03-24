SUBROUTINE dataspline ( xtrack_pix, n_radwvl, curr_rad_wvl, n_max_rspec, errstat )

  ! ---------------------------------------------------------------------
  ! Explanation of subroutine arguments:
  !
  !    xtrack_pix ........... current cross-track pixel number
  !    n_radwvl ............. number of radiance wavelenghts
  !    curr_rad_wvl ......... radiance wavelengths
  !    n_max_rspec .......... maximum number of reference spectra points
  !    errstat .............. error status returned from the subroutine
  ! ---------------------------------------------------------------------

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: max_rs_idx, mxs_idx, max_calfit_idx, &
       solar_idx, refspec_strings, hwe_idx, asy_idx, comm_idx, us1_idx, &
       us2_idx
  USE OMSAO_parameters_module, ONLY: zerospec_string, solar_i0_scd, &
       yn_i0_spc
  USE OMSAO_variables_module, ONLY: refspecs_original, common_mode_spec, &
       database, fitvar_rad_init, lo_radbnd, up_radbnd, yn_use_labslitfunc, &
       yn_solar_i0, pge_idx
  USE OMSAO_omidata_module, ONLY : omi_solcal_pars
  USE OMSAO_slitfunction_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                       INTENT (IN) :: xtrack_pix, n_radwvl, n_max_rspec
  REAL    (KIND=r8), DIMENSION (n_radwvl), INTENT (IN) :: curr_rad_wvl

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  LOGICAL                                       :: yn_full_range
  INTEGER (KIND=i4)                             :: idx, npts, locerrstat, iii, nsol, ios, k
  REAL    (KIND=r8)                             :: DU_load
  REAL    (KIND=r8), DIMENSION (n_max_rspec)    :: tmp_spec, tmp_wavl
  REAL    (KIND=r8), DIMENSION (n_radwvl)       :: dbase_loc
  REAL    (KIND=r8), DIMENSION (:), ALLOCATABLE :: solar_spc, solar_wvl, solar_conv, xsec_i0_spc

  CHARACTER (LEN=11), PARAMETER :: modulename = 'dataspline'

  locerrstat = pge_errstat_ok
  
  nsol = refspecs_original(solar_idx)%nPoints
  IF ( nsol > 0 ) THEN
     ALLOCATE ( solar_spc  (1:nsol), STAT=ios )
     ALLOCATE ( solar_wvl  (1:nsol), STAT=ios )
     ALLOCATE ( solar_conv (1:nsol), STAT=ios )
     ALLOCATE ( xsec_i0_spc(1:nsol), STAT=ios )
  END IF

  ! --------------------------------------
  ! Compute correction for Solar I0 effect
  ! --------------------------------------
  IF ( yn_solar_i0 ) THEN
     idx  = solar_idx
     iii  = max_calfit_idx + (idx-1)*mxs_idx

     solar_wvl(1:nsol) = refspecs_original(idx)%RefSpecWavs(1:nsol)
     solar_spc(1:nsol) = refspecs_original(idx)%RefSpecData(1:nsol)
     CALL convolve_data (                                                              &
          xtrack_pix, nsol,  solar_wvl(1:nsol), solar_spc(1:nsol), yn_use_labslitfunc, &
          omi_solcal_pars(hwe_idx,xtrack_pix), omi_solcal_pars(asy_idx,xtrack_pix),    &
          solar_conv(1:nsol), errstat )
  END IF
  ! ---------------------------------------------------------------------
  ! Load results into the database array. The order of the spectra is
  ! determined by the molecule indices in OMSAO_indices_module. Only 
  ! those spectra that are requested for read-in are interpolated.
  ! Note that this might be different from the actual fitting parameters
  ! to be varied in the fit - a non-zero but constant fitting parameter
  ! will still require a spectrum to make a contribution. However, if all
  ! fitting parameters, including the upper and lower bounds, are ZERO,
  ! then there is no way that spectrum can contribute, and it therefore
  ! doesn't have to be interpolated.
  !
  ! If the original reference spectrum does not cover the wavelength
  ! range of the current radiance wavelength, we interpolate only the
  ! part that is covered and set the rest to Zero.
  ! ---------------------------------------------------------------------

  ! ------------------------
  ! Spline Reference Spectra
  ! ------------------------
  DO idx = 1, max_rs_idx

     ! --------------------------------------------------------------------------
     ! We don't spline the solar reference spectrum and the undersampling spectra
     ! (computed and assigned to DATABASE in the UNDERSAMPLE subroutine), and the
     ! common mode spectrum (which hasn't been defined yet anyway).
     ! --------------------------------------------------------------------------
     IF ( (idx == solar_idx) .OR. (idx == us1_idx) .OR. (idx == us2_idx) ) CYCLE

     iii = max_calfit_idx + (idx-1)*mxs_idx

     IF ( refspecs_original(idx)%nPoints /= 0                                          .AND. &
          INDEX (TRIM(ADJUSTL(refspecs_original(idx)%FileName)), zerospec_string) == 0 .AND. &
          ( ANY (fitvar_rad_init(iii+1:iii+mxs_idx) /= 0) .OR. &
            ANY (lo_radbnd      (iii+1:iii+mxs_idx) /= 0) .OR. &
            ANY (up_radbnd      (iii+1:iii+mxs_idx) /= 0)        ) ) THEN

        ! -----------------------------------------------
        ! Define short-hand for number of spectral points
        ! -----------------------------------------------
        npts             = refspecs_original(idx)%nPoints
        tmp_wavl(1:npts) = refspecs_original(idx)%RefSpecWavs(1:npts)
        tmp_spec(1:npts) = refspecs_original(idx)%RefSpecData(1:npts)

        ! -------------------------------------------------------------
        ! Solar I0 correction, Yes or No?
        ! Check if it is a corrected species (set in params module) CCM
        ! -------------------------------------------------------------
        IF ( yn_solar_i0 .AND. yn_i0_spc(idx)) THEN
				
           DU_load = solar_i0_scd(idx)
           
           !print*,'idx:',idx
           !print*,'solar_i0_scd(idx):',solar_i0_scd(idx)
           ! --------------------------------------------------------------------
           ! 1: Interpolate cross sections to solar reference spectrum wavelength
           ! --------------------------------------------------------------------
           CALL interpolation ( &
                modulename//'--'//TRIM(ADJUSTL(refspecs_original(idx)%FileName)), &
                npts, tmp_wavl(1:npts), tmp_spec(1:npts),                         &
                nsol, solar_wvl(1:nsol), xsec_i0_spc(1:nsol),                     &
                'fillvalue', 0.0_r8, yn_full_range, locerrstat )
           ! ---------------------
           ! 2: Undo normalization
           ! ---------------------
           xsec_i0_spc(1:nsol) = xsec_i0_spc(1:nsol) * refspecs_original(idx)%NormFactor
           ! ------------------------------------
           ! 3: Compute attenuated solar spectrum
           ! ------------------------------------
           tmp_spec(1:nsol) = solar_spc(1:nsol) * EXP(-xsec_i0_spc(1:nsol)*DU_load)
           ! ------------------------------
           ! 4: Convolve with slit function
           ! ------------------------------
           CALL convolve_data (                                                            &
                xtrack_pix, nsol, solar_wvl(1:nsol), tmp_spec(1:nsol), yn_use_labslitfunc, &
                omi_solcal_pars(hwe_idx,xtrack_pix), omi_solcal_pars(asy_idx,xtrack_pix),  &
                xsec_i0_spc(1:nsol), errstat )
           ! -----------------------------------
           ! 5: Compute corrected cross sections
           ! -----------------------------------
           WHERE ( solar_conv(1:nsol) > 0.0_r8 .AND. xsec_i0_spc(1:nsol) > 0.0_r8 )
              xsec_i0_spc(1:nsol) = -LOG( xsec_i0_spc(1:nsol) / solar_conv(1:nsol) ) / DU_load
           END WHERE
           ! ---------------------
           ! 6: Redo normalization
           ! ---------------------
           xsec_i0_spc(1:nsol) = xsec_i0_spc(1:nsol) / refspecs_original(idx)%NormFactor
           ! ------------------------------------------------------------
           ! Copy values to the variables used in the interpolation below
           ! ------------------------------------------------------------
           npts             = nsol
           tmp_spec(1:npts) = xsec_i0_spc(1:npts)
           tmp_wavl(1:npts) = solar_wvl  (1:npts)
        ELSE
           ! --------------------------------------------------------------------------
           ! All other cross sections just need to be convolved with the slit function,
           ! EXCEPT for the Common Mode spectra.
           ! --------------------------------------------------------------------------
           IF ( (idx == comm_idx) ) THEN
              tmp_spec(1:npts) = common_mode_spec%RefSpecData(xtrack_pix,1:npts)
           ELSE
              CALL convolve_data (                                                           &
                   xtrack_pix, npts, tmp_wavl(1:npts), tmp_spec(1:npts), yn_use_labslitfunc, &
                   omi_solcal_pars(hwe_idx,xtrack_pix), omi_solcal_pars(asy_idx,xtrack_pix), &
                   tmp_spec(1:npts), errstat )
              !DO k = 1, npts
              !   WRITE (idx+10,'(0P1F15.5,1PE20.10)') tmp_wavl(k), tmp_spec(k)*refspecs_original(idx)%NormFactor
              !END DO
           END IF

        END IF

        ! ----------------------------------------------------------------------------
        ! Call interpolation and check returned error status. WARNING status indicates
        ! missing parts of the interpolated spectrum, while ERROR status indicates a
        ! more serious condition that requires termination.
        ! ----------------------------------------------------------------------------
        CALL interpolation ( &
             modulename//'--'//TRIM(ADJUSTL(refspecs_original(idx)%FileName)), &
             npts, tmp_wavl(1:npts), tmp_spec(1:npts),                         &
             n_radwvl, curr_rad_wvl(1:n_radwvl), dbase_loc(1:n_radwvl),        &
             'fillvalue', 0.0_r8, yn_full_range, locerrstat )

        database(idx, 1:n_radwvl) = dbase_loc(1:n_radwvl)
        CALL error_check ( &
             locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL_REFSPEC, &
             modulename//f_sep//TRIM(ADJUSTL(refspec_strings(idx))), vb_lev_default, errstat )
        IF ( .NOT. yn_full_range ) THEN
           CALL error_check ( &
             0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
             modulename//f_sep//TRIM(ADJUSTL(refspec_strings(idx))), vb_lev_develop, errstat )
        END IF
     END IF
  END DO        

  IF ( nsol > 0 ) THEN
     IF ( ALLOCATED ( solar_spc   ) )  DEALLOCATE ( solar_spc   )
     IF ( ALLOCATED ( solar_wvl   ) )  DEALLOCATE ( solar_wvl   )
     IF ( ALLOCATED ( solar_conv  ) )  DEALLOCATE ( solar_conv  )
     IF ( ALLOCATED ( xsec_i0_spc ) )  DEALLOCATE ( xsec_i0_spc )
  END IF

  !STOP
  RETURN
END SUBROUTINE dataspline


SUBROUTINE convolve_data (                                     &
     xtrack_pix, npts, wvl_in, spec_in, yn_labslit, hw1e, asy, &
     spec_conv, errstat )

  USE OMSAO_precision_module
  USE OMSAO_slitfunction_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                   INTENT (IN) :: xtrack_pix, npts
  REAL    (KIND=r8),                   INTENT (IN) :: hw1e, asy
  LOGICAL,                             INTENT (IN) :: yn_labslit
  REAL    (KIND=r8), DIMENSION (npts), INTENT (IN) :: spec_in, wvl_in

  ! ----------------
  ! Output variables
  ! ----------------
  REAL (KIND=r8), DIMENSION (npts), INTENT (OUT) :: spec_conv

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)           :: errstat, i
  CHARACTER (LEN=13), PARAMETER :: modulename = 'convolve_data'


  errstat = pge_errstat_ok

  ! -----------------------------------------------------------------
  ! Either laboratory slit function (tabulated) or Gaussian (fitted).
  ! -----------------------------------------------------------------
  IF ( yn_labslit ) THEN
     CALL omi_slitfunc_convolve (                     &
          xtrack_pix, npts, wvl_in(1:npts), spec_in(1:npts), spec_conv(1:npts), errstat )
     CALL error_check ( &
          errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_INTERPOL, &
          modulename//f_sep//'Convolution', vb_lev_default, errstat )
  ELSE
     ! -----------------------------------------------------------------
     ! Here is the Gaussian branch. We need to make sure that we use the
     ! slit function information for the current pixel. There is also no
     ! need to save the convolved spectrum, because we have to convolve
     ! for each and every pixel due to the varying slit function.
     ! -----------------------------------------------------------------
     CALL asymmetric_gaussian_sf (                                           &
          npts, hw1e, asy, wvl_in(1:npts), spec_in(1:npts), spec_conv(1:npts) )
  END IF

  RETURN
END SUBROUTINE convolve_data
