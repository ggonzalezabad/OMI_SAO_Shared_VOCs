SUBROUTINE subtract_cubic ( locwvl, npts, ll_rad, lu_rad )

  USE OMSAO_precision_module
  USE OMSAO_indices_module,        ONLY : max_rs_idx, solar_idx
  USE OMSAO_parameters_module,     ONLY : max_spec_pts, doas_npol, elsunc_np, elsunc_nw
  USE OMSAO_variables_module,      ONLY : cubic_x, cubic_y, cubic_w, database
  USE OMSAO_elsunc_fitting_module, ONLY : elsunc

  IMPLICIT NONE

  EXTERNAL cubic_specfit

  INTEGER (KIND=i4),                   INTENT (IN) :: npts, ll_rad, lu_rad
  REAL    (KIND=r8), DIMENSION (npts), INTENT (IN) :: locwvl

  INTEGER (KIND=i4)                   :: i, nlower, nupper, nfitted
  REAL    (KIND=r8)                   :: locavg, chisq
  REAL    (KIND=r8), DIMENSION (npts) :: x, ptemp, sig

  ! ================
  ! ELSUNC variables
  ! ================
  INTEGER (KIND=i4)                                     :: exval
  INTEGER (KIND=i4)                                     :: elbnd
  INTEGER (KIND=i4), DIMENSION (elsunc_np)              :: p
  REAL    (KIND=r8), DIMENSION (elsunc_nw)              :: w
  REAL    (KIND=r8), DIMENSION (doas_npol)              :: blow, bupp
  REAL    (KIND=r8), DIMENSION (max_spec_pts)           :: f
  REAL    (KIND=r8), DIMENSION (max_spec_pts,doas_npol) :: dfda
  REAL    (KIND=r8), DIMENSION (doas_npol)              :: par


  ! ======================
  ! Assign fitting weights
  ! ======================
  sig = 1.0_r8

  ! -------------------------------------------------------------------------
  ! Find limits for polynomial fitting, with ~1 nm overlap
  ! -------------------------------------------------------------------------
  ! ARRAY_LOCATE is the preferred way to calculate the bounds; however, this 
  ! particular  section of the code needs still to be tested. tpk 09 Feb 2007
  ! -------------------------------------------------------------------------
  nlower = MINVAL ( MINLOC ( locwvl(1:npts), MASK=(locwvl(1:npts) >= locwvl(ll_rad)-1.0_r8) ) )
  nupper = MAXVAL ( MAXLOC ( locwvl(1:npts), MASK=(locwvl(1:npts) <= locwvl(lu_rad)+1.0_r8) ) )
  !CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl(ll_rad)-1.0_r8, 'GE', nlower )
  !CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl(lu_rad)+1.0_r8, 'LE', nupper )

  nfitted = nupper - nlower + 1

  ! Find average position over fitted region
  locavg = SUM ( locwvl(1+nlower-1:nfitted+nlower-1) ) / REAL ( nfitted, KIND=r8 )

  ! Load temporary position file: re-define positions in order to fit
  ! about mean position
  DO i = 1, nfitted
     ptemp(i) = locwvl(i+nlower-1) - locavg
  END DO

  !     Load and fit database spectra nos. 2-11
  DO i = 1, max_rs_idx
     IF ( i /= solar_idx ) THEN
        ! ===============================================================
        ! ELBND: 0 = unconstrained
        !        1 = all variables have same lower bound
        !        else: lower and upper bounds must be supplied by the use
        ! ===============================================================  
        elbnd = 0  ;  exval = 0
        p   = -1 ; p(1) = 0  ;  p(3) = 5  ; w = -1.0
        blow(1:doas_npol) = 0.0_r8  ;  bupp(1:doas_npol) = 0.0_r8
  
        cubic_x(1:nfitted) = ptemp(1:nfitted)
        cubic_y(1:nfitted) = database(i, 1+nlower-1:nfitted+nlower-1)
        cubic_w(1:nfitted) = sig(1:nfitted)

        par = 0.0_r8 ; f = 0.0_r8 ; dfda = 0.0_r8
        CALL elsunc ( &
             par, doas_npol, nfitted, nfitted, cubic_specfit, elbnd, blow(1:doas_npol), &
             bupp(1:doas_npol), p, w, exval, f(1:nfitted), dfda(1:nfitted,1:doas_npol) )
        chisq = SUM  ( f(1:nfitted)**2 ) ! This gives the same CHI**2 as the NR routines
        x(1:npts) = locwvl(1:npts) - locavg

        cubic_x(1:npts) = x(1:npts)
        exval = 3
        CALL cubic_specfit ( &
             par(1:doas_npol), doas_npol, cubic_y(1:npts), npts, exval, dfda, 0 )
        database(i,1:npts) = database(i,1:npts) - cubic_y(1:npts)
     END IF
  END DO

  RETURN
END SUBROUTINE subtract_cubic

SUBROUTINE subtract_cubic_meas ( locwvl, npoints, locspec, ll_rad, lu_rad )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,     ONLY : max_spec_pts, doas_npol, elsunc_np, elsunc_nw
  USE OMSAO_variables_module,      ONLY : cubic_x, cubic_y, cubic_w
  USE OMSAO_elsunc_fitting_module, ONLY : elsunc

  IMPLICIT NONE

  EXTERNAL cubic_specfit

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                      INTENT (IN) :: npoints, ll_rad, lu_rad
  REAL    (KIND=r8), DIMENSION (npoints), INTENT (IN) :: locwvl

  ! ------------------
  ! Modified variables
  ! ------------------
  REAL (KIND=r8), DIMENSION (npoints), INTENT (INOUT) :: locspec

  ! ---------------
  ! Local variables
  ! ---------------
  REAL    (KIND=r8), DIMENSION (npoints)             :: tmp, ptmp, sig
  REAL    (KIND=r8), DIMENSION (doas_npol)           :: r
  INTEGER (KIND=i4)                                  :: i, nlower, nupper, nfitted
  REAL    (KIND=r8)                                  :: locavg
  REAL    (KIND=r8), DIMENSION (npoints)             :: x

  ! ================
  ! ELSUNC variables
  ! ================
  INTEGER (KIND=i4)                                :: exval
  REAL    (KIND=r8)                                :: chisq2
  INTEGER (KIND=i4)                                :: elbnd
  INTEGER (KIND=i4), DIMENSION (11)                :: p
  REAL    (KIND=r8), DIMENSION (6)                 :: w
  REAL    (KIND=r8), DIMENSION (doas_npol)         :: blow, bupp
  REAL    (KIND=r8), DIMENSION (npoints)           :: f
  REAL    (KIND=r8), DIMENSION (npoints,doas_npol) :: dfda


  ! ======================
  ! Assign fitting weights
  ! ======================
  sig = 1.0_r8

  ! -------------------------------------------------------------------------
  !     Find limits for polynomial fitting, with ~1 nm overlap
  ! -------------------------------------------------------------------------
  ! ARRAY_LOCATE is the preferred way to calculate the bounds; however, this 
  ! particular  section of the code needs still to be tested. tpk 09 Feb 2007
  ! -------------------------------------------------------------------------
  nlower = MINVAL(MINLOC( locwvl(1:npoints), MASK=(locwvl(1:npoints) >= locwvl(ll_rad)-1.0_r8) ))
  nupper = MAXVAL(MAXLOC( locwvl(1:npoints), MASK=(locwvl(1:npoints) <= locwvl(lu_rad)+1.0_r8) ))
  !CALL array_locate_r8 ( npoints, locwvl(1:npoints), locwvl(ll_rad)-1.0_r8, 'GE', nlower )
  !CALL array_locate_r8 ( npoints, locwvl(1:npoints), locwvl(lu_rad)+1.0_r8, 'LE', nupper )
  nfitted = nupper - nlower + 1

  !     Find average position over fitted region
  locavg = SUM ( locwvl(1+nlower-1:nfitted+nlower-1) ) / REAL ( nfitted, KIND=r8 )

  !     Load temporary position file: re-define positions in order to fit
  !     about mean position
  DO i = 1, nfitted
     ptmp(i) = locwvl(i+nlower-1) - locavg
  END DO

  !     Load and fit spectrum
  tmp(1:nfitted) = locspec(1+nlower-1:nfitted+nlower-1)

  ! ===============================================================
  ! ELBND: 0 = unconstrained
  !        1 = all variables have same lower bound
  !        else: lower and upper bounds must be supplied by the use
  ! ===============================================================  
  elbnd = 0  ;  exval = 0
  p   = -1 ; p(1) = 0  ;  p(3) = 5  ; w = -1.0
  blow(1:doas_npol) = 0.0_r8  ;  bupp(1:doas_npol) = 0.0_r8
  
  cubic_x(1:nfitted) = ptmp(1:nfitted)
  cubic_y(1:nfitted) =  tmp(1:nfitted)
  cubic_w(1:nfitted) =  sig(1:nfitted)

  r = 0.0_r8 ; f = 0.0_r8 ; dfda = 0.0_r8

  CALL elsunc ( &
       r, doas_npol, nfitted, nfitted, cubic_specfit, elbnd, blow(1:doas_npol), &
       bupp(1:doas_npol), p, w, exval, f(1:nfitted), dfda(1:nfitted,1:doas_npol) )
  chisq2 = SUM  ( f(1:nfitted)**2 ) ! This gives the same CHI**2 as the NR routines

  ! Re-load spec with high-pass filtered data, over whole  spectral region
  x(1:npoints) = locwvl(1:npoints) - locavg

  cubic_x(1:npoints) = x(1:npoints)
  exval = 3
  CALL cubic_specfit ( &
       r(1:doas_npol), doas_npol, cubic_y(1:npoints), npoints, exval, dfda, 0 )
  locspec(1:npoints) = locspec(1:npoints) - cubic_y(1:npoints)


  RETURN
END SUBROUTINE subtract_cubic_meas
