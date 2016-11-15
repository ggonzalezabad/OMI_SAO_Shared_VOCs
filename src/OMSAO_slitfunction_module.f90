MODULE OMSAO_slitfunction_module

  ! =================================================================
  !
  ! This module defines variables associated with error handling. It
  ! also loads/includes all (SDPTK) files that define error messages
  ! and generally deal with error handling.
  !
  ! =================================================================

  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_indices_module,    ONLY: omi_slitfunc_lun
  USE OMSAO_variables_module,  ONLY: omi_slitfunc_fname, l1b_channel
  USE OMSAO_omidata_module,    ONLY: nxtrack_max
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! --------------------------------------------------------------------
  ! Some maximum dimension definitions, basically used to circumvent the
  ! need of allocating arrays dynamically.
  ! --------------------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_max_sf_col = 23, n_max_sf_hyp = 5, n_max_sf_pol = 7
  INTEGER (KIND=i4), PARAMETER :: n_max_sf_pts = 201, n_sf_segs = 3

  ! -------------------------------------------------------------------
  ! Variables and arrays to hold the tabulated slit function parameters
  ! -------------------------------------------------------------------
  INTEGER (KIND=i4)                           :: n_sf_row, n_sf_firstrow, n_sf_colwvl, n_sf_tabwvl
  INTEGER (KIND=i4)                           :: n_sf_pol, n_sf_hyp, n_sf_pts
  INTEGER (KIND=i4), DIMENSION (2, n_sf_segs) :: sf_idx
  REAL    (KIND=r8), DIMENSION (n_max_sf_col) :: sf_colwvl
  REAL    (KIND=r8), DIMENSION (n_max_sf_pts) :: sf_tabwvl
  REAL    (KIND=r8), DIMENSION (n_max_sf_col, n_max_sf_pol, n_max_sf_hyp, n_sf_segs) :: sf_hyppar


  ! ------------------------------------------------------------
  ! Character strings to look for in the slit function data file
  ! ------------------------------------------------------------
  INTEGER   (KIND=i4),  PARAMETER, PRIVATE :: lstr = 25
  CHARACTER (LEN=lstr), PARAMETER, PRIVATE :: &
       sf_uv2_str = 'UV2 channel slit function', sf_vis_str = 'VIS channel slit function'


  ! --------------------------------------------------
  ! Tabulated slit function dimensions, maxmum values
  ! --------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_sf_cmax = 290, n_sf_rmax = 301

  INTEGER (KIND=i4), DIMENSION (nxtrack_max), PARAMETER :: sf_xtrack_center_rows = (/ &
       53,   61,  69,  77,  85,  93, 101, 109, 117, 125, 133, 141, 149, 157, 165, &
       173, 181, 189, 197, 205, 213, 221, 229, 237, 245, 253, 261, 269, 277, 285, &
       293, 301, 309, 317, 325, 333, 341, 349, 357, 365, 373, 381, 389, 397, 405, &
       413, 421, 429, 437, 445, 453, 461, 469, 477, 485, 493, 501, 509, 517, 525   /)


  ! ----------------------------------------------------------------
  ! Cross-track position mapping for the true Spatial Zoom pixels 
  ! (not rebinned to 30 "global size" cross-track pixels) has to
  ! take into account the different binning factors (8 in global,
  ! 4 in spatial zoom). The 1-60 spatial zoom cross-track positions
  ! are mapped onto positions 16-45 in global mode as follows: 
  !
  ! ipos_g = INT(i/8),                   i = 0, 489
  ! ipos_z = FLOOR((REAL(i)-120.)/4.)+1, i = 0, 489
  ! 
  ! where i = 0, 498 are the CCD detector elements that are being
  ! used.
  ! 
  ! ZOOM2GLOBAL_MAP(j), for j = 1, 60, holds the global cross-track
  ! position (and hence slit function table entry) for zoom
  ! cross-track position j.
  ! ----------------------------------------------------------------
  INTEGER (KIND=i4), DIMENSION (nxtrack_max), PARAMETER :: zoom2global_map = (/ &
       16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24,      &
       24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32,      &
       33, 33, 34, 34, 35, 35, 36, 36, 37, 37, 38, 38, 39, 39, 40, 40, 41,      &
       41, 42, 42, 43, 43, 44, 44, 45, 45 /)



  ! ---------------------------------------------------------------------
  ! The following two quantities are used to determine whether we need to
  ! reconvolve the solar spectrum. Only if either SHIFT or SQUEEZE have 
  ! changed from one iteration to the other is a reconvolution necessary.
  ! ---------------------------------------------------------------------
  REAL (KIND=r8) :: saved_shift = -1.0E+30_r8, saved_squeeze = -1.0E+30_r8


CONTAINS

  SUBROUTINE omi_slitfunc_read ( errstat )

    IMPLICIT NONE

    ! ---------------------------------------------------------------------
    ! Explanation of subroutine arguments:
    !
    !    errstat .............. error status returned from the subroutine
    ! ---------------------------------------------------------------------

    CHARACTER (LEN=17), PARAMETER :: modulename = 'omi_slitfunc_read'

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! --------------
    ! Local variable
    ! --------------
    INTEGER   (KIND=i4)  :: locerrstat, file_read_stat, funit, version
    INTEGER   (KIND=i4)  :: i, j, k
    INTEGER   (KIND=i4)  :: i_sf_ll, i_sf_l, i_sf_r, i_sf_rr
    CHARACTER (LEN=lstr) :: sf_curr_str, lastline

    ! ---------------------------------
    ! External OMI and Toolkit routines
    ! ---------------------------------
    INTEGER (KIND=i4) :: pgs_io_gen_openf, pgs_io_gen_closef

    locerrstat = pge_errstat_ok

    ! --------------------------------------
    ! Open file with tabluated slit function
    ! --------------------------------------
    version = 1
    locerrstat = PGS_IO_GEN_OPENF ( &
         omi_slitfunc_lun, PGSd_IO_Gen_RSeqFrm, 0, funit, version )
    CALL error_check ( &
         locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_OPEN_REFSPEC_FILE, &
         modulename//f_sep//TRIM(ADJUSTL(omi_slitfunc_fname)), vb_lev_default, locerrstat )
    IF (  locerrstat /= pge_errstat_ok ) THEN
       errstat = MAX(errstat, locerrstat); locerrstat = PGS_IO_GEN_CLOSEF ( funit ) ; RETURN
    END IF

    ! ----------------------------------------------------------------
    ! Set landmark string depending on in which channel we are fitting
    ! ----------------------------------------------------------------
    SELECT CASE ( l1b_channel )
    CASE ( 'UV2' )
       sf_curr_str = sf_uv2_str
    CASE ( 'VIS' )
       sf_curr_str = sf_vis_str
    END SELECT

    ! -------------------------------
    ! Skip to channel landmark string
    ! -------------------------------
    CALL skip_to_filemark ( funit, sf_curr_str, lastline, file_read_stat )
    CALL error_check ( &
         file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_E_READ_REFSPEC_FILE, &
         modulename//f_sep//sf_curr_str, vb_lev_default, locerrstat )
    IF (  locerrstat /= pge_errstat_ok ) THEN
       errstat = MAX(errstat, locerrstat) ; locerrstat = PGS_IO_GEN_CLOSEF ( funit ) ; RETURN
    END IF

    ! ---------------------------
    ! Read the slit function data
    ! ---------------------------

    ! ----------------------------
    ! Dimensions and column values
    ! ----------------------------
    READ (UNIT=funit, FMT=*, IOSTAT=file_read_stat) n_sf_row, n_sf_firstrow
    READ (UNIT=funit, FMT=*, IOSTAT=file_read_stat) n_sf_colwvl
    READ (UNIT=funit, FMT=*, IOSTAT=file_read_stat) sf_colwvl(1:n_sf_colwvl)
    READ (UNIT=funit, FMT=*, IOSTAT=file_read_stat) n_sf_tabwvl
    READ (UNIT=funit, FMT=*, IOSTAT=file_read_stat) sf_tabwvl(1:n_sf_tabwvl)
    READ (UNIT=funit, FMT=*, IOSTAT=file_read_stat) i_sf_ll, i_sf_l, i_sf_r, i_sf_rr
    READ (UNIT=funit, FMT=*, IOSTAT=file_read_stat) n_sf_pol,  n_sf_hyp

    ! --------------------------------------------------------
    ! Compute number of slit function points, i.e., the number
    ! of columns for which the slit function will be computed.
    ! --------------------------------------------------------
    n_sf_pts = i_sf_rr - i_sf_ll + 1

    ! ----------------------------------------------------
    ! Save the left and right bounding indices in an array
    ! ----------------------------------------------------
    sf_idx(1:2,1) = (/ i_sf_ll,  i_sf_l-1 /)
    sf_idx(1:2,2) = (/ i_sf_l,   i_sf_r   /)
    sf_idx(1:2,3) = (/ i_sf_r+1, i_sf_rr  /)

    ! --------------------------------------------------------------
    ! Skip the string header and read the block of hyper-parameters.
    ! There are three segments/blocks to read.
    ! --------------------------------------------------------------
    DO k = 1, n_sf_segs
       READ (UNIT=funit, FMT=*, IOSTAT=file_read_stat) lastline
       DO j = 1, n_sf_pol
          DO i = 1, n_sf_colwvl
             READ (UNIT=funit, FMT=*, IOSTAT=file_read_stat) sf_hyppar(i,j,1:n_sf_hyp,k)
          END DO
       END DO
    END DO

    ! --------------------------
    ! Close data file and return
    ! --------------------------
    locerrstat = PGS_IO_GEN_CLOSEF ( funit )

    RETURN
  END SUBROUTINE omi_slitfunc_read


  SUBROUTINE omi_slitfunc_compute ( row, nwvl, wvl, npts, sf_wvals, sf_profiles )

    IMPLICIT NONE

    ! ------------------------------------------------------------------------------
    ! Explanation of subroutine arguments:
    !
    !    row .................. the row value (just one)
    !    nwvl ................. number of center wavelengths for which to compute
    !                           slit function profiles
    !    wvl .................. the center wavelenght values
    !    npts ................. number of points in the slit function profiles
    !    sf_wvals ............. wavelenght values for the slit function profiles
    !    sf_profiles .......... the slit function profiles
    ! ------------------------------------------------------------------------------

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                   INTENT (IN) :: row, nwvl, npts
    REAL    (KIND=r8), DIMENSION (nwvl), INTENT (IN) :: wvl

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r8), DIMENSION (nwvl, npts) :: sf_wvals, sf_profiles

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                        :: i, j, k, k0, k1, k2
    INTEGER (KIND=i4)                                        :: locerrstat
    REAL    (KIND=r8)                                        :: tmprow
    REAL    (KIND=r8), DIMENSION (nwvl)                      :: tmpwvl, tmphyp
    REAL    (KIND=r8), DIMENSION (n_sf_hyp)                  :: xhyp
    REAL    (KIND=r8), DIMENSION (nwvl, n_sf_pol, n_sf_segs) :: wvlpar
    REAL    (KIND=r8), DIMENSION (n_sf_pol, npts)            :: xpol

    ! -------------------------------
    ! Name of the function/subroutine
    ! -------------------------------
    CHARACTER (LEN=20), PARAMETER :: modulename = 'omi_slitfunc_compute'


    ! ----------------------------------------------------------------------
    ! Compose  row and array of column values, making sure that they stay 
    ! within the  bounds of the tabulated slit function parameters. Constant
    ! extension is used beyod the lowest and highest tabulated values.
    ! ----------------------------------------------------------------------
    ! First the row
    ! --------------
    tmprow = MAX ( MIN ( REAL(row,KIND=r8 ), REAL(n_sf_row+n_sf_firstrow-1,KIND=r8) ), &
                   REAL(n_sf_firstrow, KIND=r8) )

    ! --------------------------
    ! Then the wavelength values
    ! --------------------------
    tmpwvl(1:nwvl) = wvl(1:nwvl)
    WHERE ( tmpwvl(1:nwvl) < sf_colwvl(1) )
       tmpwvl(1:nwvl) = sf_colwvl(1)
    ENDWHERE
    WHERE ( tmpwvl(1:nwvl) > sf_colwvl(n_sf_colwvl) )
       tmpwvl(1:nwvl) = sf_colwvl(n_sf_colwvl)
    ENDWHERE

    ! --------------------------------------------------------------------
    ! Compute the powers of the hyper-parameter and polynomial base values
    ! --------------------------------------------------------------------
    DO i = 1, n_sf_hyp
       xhyp(i) = tmprow**(i-1)
    END DO

    DO i = 1, n_sf_pol
       xpol(i, 1:npts) = sf_tabwvl(1:npts)**(i-1)
    END DO

    ! -----------------------------------------------------------------------------------
    ! Compose the polynomial parameters (a function of column)  from the hyper-parameters
    ! (a function of rows): Interpolate the hyper-parameters to the destination columns
    ! and sum to the particular row value the slit function is to be computed for. Save
    ! the results in a temporary parameter array.
    ! -----------------------------------------------------------------------------------
    wvlpar = 0.0_r8
    DO k = 1, n_sf_segs
       DO i = 1, n_sf_pol
          DO j = 1, n_sf_hyp
             locerrstat = pge_errstat_ok
             CALL ezspline_1d_interpolation ( &
                  n_sf_colwvl, sf_colwvl(1:n_sf_colwvl), sf_hyppar(1:n_sf_colwvl,i,j,k), &
                  nwvl, tmpwvl(1:nwvl), tmphyp(1:nwvl), locerrstat )
             wvlpar(1:nwvl,i,k) = wvlpar(1:nwvl,i,k) + tmphyp(1:nwvl)*xhyp(j)
          END DO
       END DO
    END DO

    ! --------------------------------------------------------------------------
    ! Sum over the polynomial coefficients to create the slit function profiles.
    ! --------------------------------------------------------------------------
    k0 = (n_sf_tabwvl-1) / 2
    sf_profiles = 0.0_r8

    DO j = 1, nwvl
       DO k = 1, n_sf_segs
          k1 = sf_idx(1,k)+k0 ;  k2 = sf_idx(2,k)+k0
          DO i = 1, n_sf_pol
             sf_profiles(j,k1:k2) = sf_profiles(j,k1:k2) + wvlpar(j,i,k)*xpol(i,k1:k2)
          END DO
       END DO
       IF ( MAXVAL(sf_profiles(j,1:npts)) > 0.0_r8 ) &
            sf_profiles(j,1:npts) = sf_profiles(j,1:npts) / MAXVAL(sf_profiles(j,1:npts))
       WHERE ( sf_profiles(j,1:npts) < 0.0_r8 )
          sf_profiles(j,1:npts) = 0.0_r8
       ENDWHERE
    END DO

    sf_wvals   = 0.0_r8
    DO i = 1, nwvl
       sf_wvals(i,1:npts) =  sf_tabwvl(1:npts) + wvl(i)
    END DO

    RETURN
  END SUBROUTINE omi_slitfunc_compute


  SUBROUTINE omi_slitfunc_convolve ( xtrack_pix, nwvl, wvl, spec, spec_conv, errstat )

    IMPLICIT NONE

    ! ---------------------------------------------------------------------
    ! Explanation of subroutine arguments:
    !
    !    xtrack_pix ........... current cross-track pixel number
    !    nwvl ................. number of points in current spectrum
    !    wvl .................. wavelengths of current spectrum
    !    spec ................. current spectrum
    !    spec_conv ............ convolved spectrum
    !    errstat .............. error status returned from the subroutine
    ! ---------------------------------------------------------------------

    ! ---------------------------------------
    ! gga, several changes to fix convolution
    ! ---------------------------------------

    CHARACTER (LEN=21), PARAMETER :: modulename = 'omi_slitfunc_convolve'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                  INTENT (IN) :: nwvl, xtrack_pix
    REAL    (KIND=r8), DIMENSION(nwvl), INTENT (IN) :: wvl, spec
    
    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION(nwvl), INTENT (OUT) :: spec_conv
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, row
    INTEGER (KIND=i4) :: l, j, j1, j2, k1, k2, k1ext, k2ext
    REAL    (KIND=r8)                                :: sf_area, sf_err, delvar, swvl, ewvl
    REAL    (KIND=r8), DIMENSION (nwvl, n_sf_tabwvl) :: sf_wvals, sf_profiles
    REAL    (KIND=r8), DIMENSION (nwvl+1000)              :: wvlext, specext, convtmpext
    REAL    (KIND=r8), DIMENSION (nwvl+1000, n_sf_tabwvl) :: sf_wvalsext, sf_profilesext
    LOGICAL                                          :: yn_full_range
    INTEGER (KIND=i4)                                :: ntmp, nswvl, newvl, nextwvl
    ! The MAX(,) dimension prevents array bound problems for small values of NWVL
    REAL    (KIND=r8), DIMENSION (MAX(n_sf_tabwvl,nwvl)) :: sfwvl_tmp, sfpro_tmp

    ! ---------------------------------------------------
    ! Initialize important output and temporary variables
    ! ---------------------------------------------------
    locerrstat   = pge_errstat_ok
    sf_area = 0.0_r8     ; delvar  = 0.0_r8
    swvl    = 0.0_r8     ; ewvl    = 0.0_r8
    sf_wvals = 0.0_r8    ; sf_profiles = 0.0_r8
    wvlext = 0.0_r8      ; specext = 0.0_r8
    sf_wvalsext = 0.0_r8 ; sf_profilesext = 0.0_r8

    ! -------------------------------------------------------------------
    ! First things first: Compute all the slit function profiles required
    ! for the convolution of the spectrum. Since the parameterization of
    ! the slit function is row-based, we first have to retrieve current
    ! row value, a function of the cross track position.
    ! -------------------------------------------------------------------
    row = sf_xtrack_center_rows(xtrack_pix)

    ! --------------------------------------------------------------
    ! Read the slit function for the wavelength we are interested at
    ! --------------------------------------------------------------
    CALL omi_slitfunc_compute ( row, nwvl, wvl(1:nwvl), n_sf_tabwvl, &
         sf_wvals(1:nwvl,1:n_sf_tabwvl), sf_profiles(1:nwvl,1:n_sf_tabwvl) )

    ! ----------------------------------------------------------
    ! No we work how much further the edges go from the wvl grid
    ! ----------------------------------------------------------
    ! End points for sf_wvals(1,1) and sf_wvals(nwvl,n_sf_tabwvl)
    ! Slit function space
    ! -------------------
    swvl = sf_wvals(1,1)
    ewvl = sf_wvals(nwvl,n_sf_tabwvl)

    ! ----------------------------------------------------------------
    ! Now extend wvl to wvlext with a similar spacing to wvl grid from
    ! wvl(1) to swvl and wvl(nwvl) to ewvl
    ! ----------------------------------------------------------------
    delvar = wvl(2)-wvl(1); nswvl = 2*((wvl(1)-swvl) / delvar)
    DO l = 1, nswvl
       wvlext(l)  = wvl(1) - delvar * (nswvl-l+1)
       specext(l) = spec(1)
    END DO
    DO l = nswvl+1, nswvl+1+nwvl-1
       wvlext(l)  = wvl(l-nswvl)
       specext(l) = spec(l-nswvl)
    END DO
    delvar = wvl(nwvl)-wvl(nwvl-1); newvl = 2*((ewvl-wvl(nwvl)) / delvar)
    DO l = nswvl+nwvl+1,nswvl+nwvl+1+newvl
       wvlext(l)  = wvl(nwvl) + delvar * (l-(nswvl+nwvl+1)+1)
       specext(l) = spec(nwvl)
    END DO
   
    nextwvl = nswvl + nwvl + newvl

    ! ----------------------------------------------------------------
    ! Compute again the slit functions but this time with the extended
    ! grid
    ! ----------------------------------------------------------------
    CALL omi_slitfunc_compute ( row, nextwvl, wvlext(1:nextwvl), n_sf_tabwvl, &
         sf_wvalsext(1:nextwvl,1:n_sf_tabwvl), sf_profilesext(1:nextwvl,1:n_sf_tabwvl) )

    ! -----------------------------------------------------------------------------
    ! In the following loop over the wavelengths in the input spectrum, we perform
    ! the following steps:
    ! (1) Interpolation of the slit function to the spectrum wavelength:
    !     (a) Find indices J1 and J2 in the slit function wavelength array SF_WVALS
    !         bounding the non-zero part of the slit function from left and right;
    !     (b) Find indices K1 and K2 in the spectrum wavelength array WVL such that
    !         SF_WVALS(J1) <= WVL(K1)  < WVL(K2) <= SF_WAVLS(J2)
    !     (c) Interpolation.
    ! (2) Normalize the interpolated slit function to AREA=1.
    ! (3) Convolve the spectrum.
    ! -----------------------------------------------------------------------------
    ! gga
    spec_conv = 0.0_r8

    DO l = 1, nextwvl

       convtmpext = 0.0_r8
       ! ----------------------------------------------
       ! Find indices J1 and J2; add 1 at either end to
       ! assure that we include the ZERO values.
       ! -----------------------------------------------------------------
       ! J1 is the minimum position of SF > 0 minus 1, but not less than 1
       ! -----------------------------------------------------------------
       j1  = 1
       get_sf_start: DO j = 2, n_sf_tabwvl
          IF ( sf_profilesext(l,j) > 0.0_r8 ) THEN
             j1 = j-1 ; EXIT get_sf_start
          END IF
       END DO get_sf_start
       ! ------------------------------------------------------------------------
       ! J2 is the maximum position of SF > 0 plus 1, but not more than N_SF_NPTS
       ! ------------------------------------------------------------------------
       j2  = n_sf_tabwvl
       get_sf_end: DO j = n_sf_tabwvl-1, 1, -1
          IF ( sf_profilesext(l,j) > 0.0_r8 ) THEN
             j2 = j+1 ; EXIT get_sf_end
          END IF
       END DO get_sf_end
       IF ( j2 <= 0 ) j2 = n_sf_tabwvl

       ! ----------------------------------------------------------------------------------
       ! Now find the bounding indices K1ext and K2ext in the extended spectrum wavelengths
       ! ----------------------------------------------------------------------------------
       CALL array_locate_r8 ( nextwvl, wvlext(1:nextwvl), sf_wvalsext(l,j1), 'GE', k1ext )
       CALL array_locate_r8 ( nextwvl, wvlext(1:nextwvl), sf_wvalsext(l,j2), 'LE', k2ext )
       IF (k1ext == k2ext) THEN
          CYCLE
       END IF

       ! -----------------------------------------------------------------------------
       ! Now interpolate the slit function to the spectrum wavlengths. We are doing it
       ! this way disregard of the number of spectral points in either array (i.e., we
       ! do not interpolate the lower resolution grid to the higher resolution one)
       ! because we don't want to worry about introducing any additional undersampling
       ! effects. This may or may not be an issue but better safe than sorry.
       ! -----------------------------------------------------------------------------
       ! Copy to temporary arrays to avoid creation of same in call to subroutine
       ! ------------------------------------------------------------------------
       ntmp              = j2-j1+1
       sfwvl_tmp(1:ntmp) = sf_wvalsext   (l,j1:j2)
       sfpro_tmp(1:ntmp) = sf_profilesext(l,j1:j2)

       ! -------------------------------------------------------
       ! Interpolation of the extended SF to obtain the area for
       ! normalization applied below
       ! -------------------------------------------------------
       CALL interpolation ( &
            modulename, ntmp, sfwvl_tmp(1:ntmp), sfpro_tmp(1:ntmp),                           &
            k2ext-k1ext+1, wvlext(k1ext:k2ext), convtmpext(k1ext:k2ext), 'endpoints', 0.0_r8, &
            yn_full_range, locerrstat )

       ! ------------------------------------------------
       ! If the SF doesn't contribute between nwvl(1) and 
       ! nwvl(nwvl) then we can skip this l.
       ! ------------------------------------------------
       IF (wvlext(k2ext) < wvl(1) .OR. wvlext(k1ext) > wvl(nwvl)) CYCLE

       ! --------------------------------------------------------
       ! Renormalize the slit function to AREA = 1, then convolve
       ! --------------------------------------------------------
       sf_area = SUM(convtmpext(k1ext:k2ext))
       IF ( sf_area == 0.0_r8 ) sf_area = 1.0_r8

       ! -----------------------
       ! Save to temporal arrays
       ! -----------------------
       convtmpext(k1ext:k2ext) = convtmpext(k1ext:k2ext) / sf_area

       ! --------------------------------------------------------------
       ! Find indices  that contribute to the actual fitting window wvl
       ! --------------------------------------------------------------
       CALL array_locate_r8 ( nwvl, wvl(1:nwvl), wvlext(k1ext), 'GE', k1 )
       CALL array_locate_r8 ( nwvl, wvl(1:nwvl), wvlext(k2ext), 'LE', k2 )

       ! -----------------------------------------------------------------
       ! Find corresponding indices to wvl(k1) and wvl(k2) in wvlext array
       ! -----------------------------------------------------------------
       k1ext = MINLOC(ABS(wvlext-wvl(k1)),1)
       k2ext = MINLOC(ABS(wvlext-wvl(k2)),1)
      
       ! ----------------------------------------------------------------------------
       ! Done normalization, we perform the convolution of the spectrum with the slit
       ! function.
       ! ----------------------------------------------------------------------------
       spec_conv(k1:k2) = spec_conv(k1:k2) + specext(l) * convtmpext(k1ext:k2ext)

    END DO
    ! gga 
    
    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE omi_slitfunc_convolve

  SUBROUTINE asymmetric_gaussian_sf ( npoints, hw1e, e_asym, wvlarr, specarr, specmod)

    ! =========================================================================
    !
    ! Convolves input spectrum with an asymmetric Gaussian slit function of
    ! specified HW1E (half-width at 1/e intensity) and asymmetry factor E_ASYM.
    !
    ! The asymetric Gaussian g(x) is defined as
    !                   _                                   _
    !                  |               x^2                   |
    !      g(x) =  EXP | - --------------------------------- |
    !                  |_   (hw1e * (1 + SIGN(x)*e_asym))^2 _|
    !
    ! g(x) becomes symmetric for E_ASYM = 0.
    !
    ! =========================================================================


    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: npoints
    REAL    (KIND=r8),                      INTENT (IN) :: hw1e, e_asym
    REAL    (KIND=r8), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr
    
    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION (npoints), INTENT (OUT) :: specmod

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                        :: i, j, nslit, sslit, eslit
    REAL    (KIND=r8)                        :: slitsum, sliterr, cwvl, lwvl, rwvl
    REAL    (KIND=r8), DIMENSION (3*npoints) :: spc_temp, wvl_temp, sf_val, xtmp, ytmp

    REAL (KIND=r8) :: signdp
    EXTERNAL signdp

    ! --------------------------------------------------------
    ! Initialize output variable (default for "no convolution"
    ! --------------------------------------------------------
    specmod(1:npoints) = specarr(1:npoints)
  
    ! -----------------------------------------------
    ! No Gaussian convolution if Halfwidth @ 1/e is 0
    ! -----------------------------------------------
    IF ( hw1e == 0.0_r8 ) RETURN


    ! ------------------------------------------------------------------------
    ! One temporary variable is SPC_TEMP, which is three times the size of
    ! SPEC. For the convolution routine to work (hopefully) in each and
    ! every case, we reflect the spectrum at its end points to always have a
    ! fully filled slit function. But this causes some real index headaches
    ! when the slit function wraps around at the ends. Performing the mirror
    ! imaging before we get to the convolution helps to keep things a little
    ! more simple.
    !
    ! Note that this approach is the same as for the pre-tabulated OMI lab
    ! slit function. It is adopted here because now we not only convolve the
    ! solar spectrum, but also any higher resolution reference cross sections,
    ! and these may not necessarily be equidistant in wavelength.
    ! ------------------------------------------------------------------------
    spc_temp(npoints+1:2*npoints) = specarr(1:npoints)
    wvl_temp(npoints+1:2*npoints) = wvlarr (1:npoints)
    DO i = 1, npoints
       spc_temp(npoints+1-i) = specarr(i)
       wvl_temp(npoints+1-i) = 2.0_r8*wvlarr(1)-wvlarr(i) -0.001_r8
       spc_temp(2*npoints+i) = specarr(npoints+1-i)
       wvl_temp(2*npoints+i) = 2.0_r8*wvlarr(npoints)-wvlarr(npoints+1-i) +0.001_r8
    END DO

    ! ------------------------------------------------------------------------
    ! We now compute the asymmetric Gaussian for every point in the spectrum.
    ! Starting from the center point, we go outwards and stop accumulating
    ! points when both sides are less than 0.001 of the maximum slit function.
    ! Since we are starting at the center wavelength, this can be set to 1.0.
    ! Remember that the original wavelength array is now located at indices
    ! NPOINTS+1:2*NPOINTS
    ! ------------------------------------------------------------------------
    DO i = 1, npoints
       sf_val = 0.0_r8
       cwvl = wvl_temp(npoints+i)

       sf_val(npoints+i) = 1.0_r8
       getslit: DO j = 1, npoints
!!$ gga         sslit = npoints+i-j ; lwvl = wvl_temp(sslit) - cwvl 
!!$ gga         eslit = npoints+i+j ; rwvl = wvl_temp(eslit) - cwvl
          sslit = npoints+i-j ; lwvl = cwvl - wvl_temp(sslit)
          eslit = npoints+i+j ; rwvl = cwvl - wvl_temp(eslit)

 !         sf_val(sslit) = EXP(-lwvl**2 / ( hw1e * (1.0_r8 + signdp(lwvl)*e_asym) )**2)
 !        sf_val(eslit) = EXP(-rwvl**2 / ( hw1e * (1.0_r8 + signdp(rwvl)*e_asym) )**2)
          sf_val(sslit) = EXP(-(ABS(lwvl/  hw1e ))**e_asym)
          sf_val(eslit) = EXP(-(ABS(rwvl/  hw1e ))**e_asym)
          IF ( sf_val(sslit) < 0.0005_r8 .AND. sf_val(sslit) < 0.0005_r8 ) EXIT getslit
       END DO getslit

       ! ----------------------------------
       ! The number of slit function points
       ! ----------------------------------
       nslit = eslit - sslit + 1
       ! ----------------------------------------------------------------
       ! Compute the norm of the slitfunction. It should be close to 1
       ! already, but making sure doesn't hurt.
       ! ----------------------------------------------------------------
       xtmp(1:nslit) = wvl_temp(sslit:eslit)-cwvl
       ytmp(1:nslit) = sf_val  (sslit:eslit)
       CALL cubint ( &
            nslit, xtmp(1:nslit), ytmp(1:nslit), 1, nslit, slitsum, sliterr)

       IF ( slitsum > 0.0_r8 ) sf_val(sslit:eslit) = sf_val(sslit:eslit) / slitsum

       ! ---------------------------------------------------------------------
       ! Prepare array for integration: Multiply slit function values with the
       ! spectrum array to be convolved.
       ! ---------------------------------------------------------------------
       ytmp(1:nslit) = sf_val(sslit:eslit) * spc_temp(sslit:eslit)

       ! ----------------------------------------------------------
       ! Folding (a.k.a. integration) of spectrum and slit function
       ! ----------------------------------------------------------
       CALL cubint ( &
            nslit, xtmp(1:nslit), ytmp(1:nslit), 1, nslit, specmod(i), sliterr)
    END DO

    RETURN
  END SUBROUTINE asymmetric_gaussian_sf

END MODULE OMSAO_slitfunction_module
