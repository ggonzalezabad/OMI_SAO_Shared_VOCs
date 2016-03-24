SUBROUTINE asym_gauss ( npoints, hw1e, e_asym, wvlarr, specarr, specmod)

  ! =========================================================================
  !
  ! Convolves input spectrum with an asymmetric Gaussian slit function of
  ! specified HW1E (half-width at 1/e intensity) and asymmetry factor E_ASYM.
  !
  ! The asymetric Gaussian g(x) is defined as
  !                   _                                  _
  !                  |               x^2                  |
  !      g(x) =  EXP | - -------------------------------- |
  !                  |_   (hw1e * (1 + SIGN(x)*e_aym))^2 _|
  !
  ! g(x) becomes symmetric for E_ASYM = 0.
  !
  ! =========================================================================

  USE OMSAO_precision_module

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER (KIND=i4),                      INTENT (IN) :: npoints
  REAL    (KIND=r8),                      INTENT (IN) :: hw1e, e_asym
  REAL    (KIND=r8), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=r8), DIMENSION (npoints), INTENT (OUT) :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER (KIND=i4)                      :: i, j, j1, j2, nslit, mslit, sslit, eslit, locerrstat
  REAL    (KIND=r8)                      :: delwvl, slitsum,  maxslit
  REAL    (KIND=r8), DIMENSION (npoints) :: slit, locwvl, locsli
  INTEGER (KIND=i4), DIMENSION (npoints) :: idx

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

  ! --------------------------------------------------------------
  ! Find the number of spectral points that fall within a Gaussian
  ! slit function with values >= 0.001. Remember that we have an
  ! asymmetric Gaussian, so we create a wavelength array symmetric
  ! around 0. The spacing is provided by the equidistant WVLARR.
  ! --------------------------------------------------------------
  delwvl = wvlarr(2) - wvlarr(1) ; locsli = 0.0_r8
  DO i = -npoints/2, npoints/2
     j = MIN ( i + npoints/2 + 1, npoints )
     locwvl(j) = delwvl * REAL(i,KIND=r8)
     locsli(j) = &
          EXP(-(locwvl(j))**2 / ( hw1e * (1.0_r8 + signdp(locwvl(j))*e_asym) )**2)
  END DO

  ! ------------------------------------------------------------------------
  ! Find the array entries that mark the region of SLIT >= 0.0001*MAX(SLIT).
  ! We start by initializing the auxilliary IDX as IDX(n) = n. We then set
  ! all those entries to 0 that are outside the accepted minimum value of 
  ! the slit function. The start and end indices simply follow as the 
  ! minimum and maximum values of IDX where IDX /= 0.
  ! ------------------------------------------------------------------------
  maxslit = MAXVAL(locsli(1:npoints))
  sslit = 0 ; eslit = 0 ; idx = (/ (i, i = 1, npoints) /)
  WHERE ( locsli < 0.001_r8*maxslit )
     idx = 0
  END WHERE
  sslit = MINVAL(MINLOC( idx, MASK = (idx /= 0) ))
  eslit = MAXVAL(MAXLOC( idx, MASK = (idx /= 0) ))

  ! ----------------------------------------------------------
  ! Compute number of slit function points, and the final slit
  ! function array; the latter will be normalized to 1.
  ! ----------------------------------------------------------
  nslit = eslit - sslit + 1
  slit(1:nslit) = locsli(sslit:eslit)

  ! ------------------------------------------------------------------
  ! The original summation - now superceeded by the integration scheme
  ! ------------------------------------------------------------------
  !  slitsum = SUM(slit(1:nslit))
  ! ------------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! Folding (a.k.a. integration) of solar spectrum and slit function
  ! ----------------------------------------------------------------
  CALL davint ( &
       locwvl(sslit:eslit), slit(1:nslit), nslit, locwvl(sslit), locwvl(nslit), &
       slitsum, locerrstat)
  IF ( slitsum > 0.0_r8 ) slit(1:nslit) = slit(1:nslit) / slitsum

  ! --------------------------------------------------------
  ! Find the maximum of the slit. Note that since we may be 
  ! with an asymmetric Gaussian, the maximum must not always
  ! be in the center of the Gaussian function.
  ! --------------------------------------------------------
  mslit = MAXVAL( MAXLOC ( slit(1:nslit) ) )

  specmod(1:npoints) = 0.0_r8 !specarr(1:npoints)
  !IF ( slitsum > 0.0_r8 ) specmod(1:npoints) = specmod(1:npoints) / slitsum

  ! ---------------------------------------------------------------
  ! Convolve spectrum. First do the middle part, where we have full
  ! overlap coverage of the slit function. Again, remember the
  ! asymmetry of the Gaussian, which makes impossible a simple
  ! 50-50 division of the summation interval.
  ! ---------------------------------------------------------------

  ! Make a local copy of the NSLIT spectrum points to be convolved
  ! with the slit function. The spectrum points to be convolved are
  ! arranged such that the updated index corresponds to the maximum
  ! of the slit function (MSLIT). For simplicity we reflect the
  ! spectrum at the array end points.

  ! ----------------------------------------------------
  ! Loop over all points of the spectrum to be convolved
  ! ----------------------------------------------------
  specmod(1:npoints) = 0.0_r8 !specarr(1:npoints)
  DO i = 1, npoints
     locsli(1:nslit) = 0.0_r8
     ! First do the right half of the slit function
     DO j = mslit, nslit
        j1 = i+(j-mslit) ; IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
        locsli(j) = specarr(j1)
        idx(j) = j1
     END DO
     ! Now the left half of the slit function
     DO j = mslit-1, 1, -1
        j2 = i-((mslit-1)-j) - 1 ; IF ( j2 < 1 ) j2 = ABS(j2) + 2
        locsli(j) = specarr(j2)
        idx(j) = j2
     END DO
     specmod(i) = DOT_PRODUCT(slit(1:nslit), locsli(1:nslit))
     locsli(1:nslit) =locsli(1:nslit) * slit(1:nslit)

     CALL davint ( &
          locwvl(sslit:eslit), locsli(1:nslit), nslit, locwvl(sslit), locwvl(eslit), &
          specmod(i), locerrstat)
  END DO

  RETURN
END SUBROUTINE asym_gauss


SUBROUTINE gauss (wvlarr, specarr, specmod, npoints, hw1e)

  !     convolves input spectrum with gaussian slit function of specified
  !     hw1e (half-width at 1/e intensity).
  
  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER (KIND=i4),                      INTENT (IN)    :: npoints
  REAL    (KIND=r8),                      INTENT (IN)    :: hw1e
  REAL    (KIND=r8), DIMENSION (npoints), INTENT (IN)    :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=r8), DIMENSION (npoints), INTENT (OUT) :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER (KIND=i4)                      :: nhi, nlo, i, j, nslit
  REAL    (KIND=r8)                      :: emult, delwvl, slitsum, slit0
  REAL    (KIND=r8), DIMENSION (npoints) :: slit

  ! ----------------------------------------------------------------
  ! Initialization of output variable (default for "no convolution")
  ! ----------------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)

  ! --------------------------------------
  ! No convolution if halfwidth @ 1/e is 0
  ! --------------------------------------
  IF ( hw1e == 0.0_r8 ) RETURN

  emult = -1.0_r8 / ( hw1e * hw1e )
  delwvl = wvlarr(2) - wvlarr(1)

  !     Apply slit function convolution

  !     Calculate slit function values out to 0.001 times x0 value,
  !     normalize so that sum = 1.

  slitsum = 1.0_r8  ;  slit0 = 1.0_r8
  i = 1  ;  nslit = 0
  numslit: DO WHILE ( nslit < npoints )
     slit (i) = EXP (emult * (delwvl * REAL(i,KIND=r8))**2)
     slitsum = slitsum + 2.0_r8 * slit (i)
     IF (slit (i) / slit0 <= 0.001_r8 .OR. i == npoints ) THEN
        nslit = i
        EXIT numslit
     END IF
     i = i + 1
  END DO numslit
  
  slit0 = slit0 / slitsum
  slit(1:nslit) = slit(1:nslit) / slitsum

  !     Convolve spectrum.  reflect at endpoints.
  specmod(1:npoints) = slit0 * specarr(1:npoints)
  DO i = 1, npoints
     DO j = 1, nslit
        nlo = MAXVAL( (/ 1, ABS(i - j) /) )
        nhi = i + j
        IF ( nhi > npoints ) nhi = MAXVAL ( (/ 1, 2*npoints-nhi /) )
        specmod(i) = specmod(i) + slit(j) * ( specarr(nlo) + specarr(nhi) )
     END DO
  END DO

  RETURN
END SUBROUTINE gauss

