SUBROUTINE specfit (                                       &
     nfitvar, fitvar, nspecpts, lowbnd, uppbnd, max_itnum, &
     covar, fitspec, fitres, exval, itnum, fitfunc           )

  USE OMSAO_precision_module
  USE OMSAO_indices_module,        ONLY: elsunc_userdef
  USE OMSAO_parameters_module,     ONLY: elsunc_np, elsunc_nw, forever
  USE OMSAO_variables_module,      ONLY: &
       tol, epsrel, epsabs, epsx, num_fitfunc_calls, num_fitfunc_jacobi, max_fitfunc_calls
  USE OMSAO_elsunc_fitting_module, ONLY: elsunc

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER (KIND=i4),                       INTENT (IN) :: nfitvar, nspecpts, max_itnum
  REAL    (KIND=r8), DIMENSION (nfitvar),  INTENT (IN) :: lowbnd, uppbnd

  ! ==================
  ! Modified variables
  ! ==================
  REAL (KIND=r8), DIMENSION (nfitvar),  INTENT (INOUT) :: fitvar

  ! ================
  ! Output variables
  ! ================
  INTEGER (KIND=i4),                              INTENT (OUT) :: exval, itnum
  REAL    (KIND=r8), DIMENSION (nfitvar,nfitvar), INTENT (OUT) :: covar
  REAL    (KIND=r8), DIMENSION (nspecpts),        INTENT (OUT) :: fitres, fitspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER (KIND=i4)                               :: elbnd
  INTEGER (KIND=i4), DIMENSION (elsunc_np)        :: p
  REAL    (KIND=r8), DIMENSION (elsunc_nw)        :: w
  REAL    (KIND=r8), DIMENSION (nfitvar)          :: blow, bupp
  REAL    (KIND=r8), DIMENSION (nspecpts)         :: f
  REAL    (KIND=r8), DIMENSION (nspecpts,nfitvar) :: dfda

  EXTERNAL fitfunc

  ! ===============================================================
  ! ELBND: 0 = unconstrained
  !        1 = all variables have same lower bound
  !        else: lower and upper bounds must be supplied by the use
  ! ===============================================================  
  elbnd = elsunc_userdef

  exval = 0 ; itnum = 0
  
  p   = -1    ;  p(1)   = 0  ;  p(3) = max_itnum
  w   = -1.0  ;  w(1:4) = (/ tol,  epsrel,  epsabs,  epsx /)

  blow(1:nfitvar) = lowbnd(1:nfitvar)
  bupp(1:nfitvar) = uppbnd(1:nfitvar)

  ! ---------------------------------------------------------------------------------
  ! Reset to ZERO the number of calls to the fitting function, NUM_FITFUNC_CALLS.
  ! This variable is incremented with each call to the fitting function and 
  ! checked agains MAX_FITFUNC_CALLS. Exceeded function calls lead to uncomputability
  ! followed by termination of the iteration process.
  ! ---------------------------------------------------------------------------------
  num_fitfunc_calls = 0 ; num_fitfunc_jacobi = 0
  max_fitfunc_calls = max_itnum * nfitvar * nfitvar ! Empiric value
  IF ( max_itnum < 0 ) max_fitfunc_calls = forever

  CALL elsunc ( &
       fitvar(1:nfitvar), nfitvar, nspecpts, nspecpts, fitfunc, elbnd, blow(1:nfitvar), &
       bupp(1:nfitvar), p, w, exval, f(1:nspecpts), dfda(1:nspecpts,1:nfitvar) )

  ! -----------------------------
  ! Save the number of iterations
  ! -----------------------------
  itnum = p(6)

  ! ------------------------------------------------------------------
  ! Call to ELSUNC fitting function to obtain complete fitted spectrum
  ! ------------------------------------------------------------------
  CALL fitfunc ( fitvar(1:nfitvar), nfitvar, fitspec(1:nspecpts), nspecpts, 3, dfda, 0 )

  ! ---------------------------------------------------------------
  ! Compute fitting residual
  ! FITRES is the negative of the returned function F = Model-Data.
  ! ---------------------------------------------------------------
  fitres(1:nspecpts) = -f(1:nspecpts)

  ! Covariance matrix.
  ! ------------------
  covar(1:nfitvar,1:nfitvar) = dfda(1:nfitvar,1:nfitvar)

  RETURN
END SUBROUTINE specfit
