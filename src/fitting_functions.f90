SUBROUTINE specfit_func_sol ( fitvar, nfitvar, ymod, npoints, ctrl, dyda, mdy )

  !
  ! Calculates the Solar spectrum and its derivatives for ELSUNC
  !
  ! NOTE: the variable DYDA as required here is the transpose of that 
  !       rquired for the Numerical Recipes
  !

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : elsunc_parsoob_eval, elsunc_infloop_eval
  USE OMSAO_variables_module,  ONLY : &
       yn_smooth, yn_doas, fitwavs, fitweights, currspec, sol_wav_avg, &
       lobnd, upbnd, num_fitfunc_calls, num_fitfunc_jacobi, max_fitfunc_calls

  IMPLICIT NONE

  INTEGER (KIND=i4),                              INTENT (IN)    :: nfitvar, npoints, mdy
  INTEGER (KIND=i4),                              INTENT (INOUT) :: ctrl
  REAL    (KIND=r8), DIMENSION (nfitvar),         INTENT (IN)    :: fitvar
  REAL    (KIND=r8), DIMENSION (npoints),         INTENT (INOUT) :: ymod
  REAL    (KIND=r8), DIMENSION (npoints,nfitvar), INTENT (INOUT) :: dyda

  REAL    (KIND=r8), DIMENSION (npoints) :: locwvl

  locwvl(1:npoints) = fitwavs(1:npoints)
  SELECT CASE ( ABS ( ctrl ) )
  CASE ( 1 )
     ! -------------------------------------------------
     ! Count the number of calls to the fitting function
     ! and terminate if we exceed the allowed maximum.
     ! -------------------------------------------------
     num_fitfunc_calls = num_fitfunc_calls + 1
     IF ( num_fitfunc_calls > max_fitfunc_calls ) THEN
        ctrl = INT(elsunc_infloop_eval, KIND=i4)
        RETURN
     END IF

     ! ------------------------------------------------
     ! Check for any fitting variables "out of bounds".
     ! ------------------------------------------------
     IF ( ANY(fitvar(1:nfitvar) < lobnd(1:nfitvar)) .OR. &
          ANY(fitvar(1:nfitvar) > upbnd(1:nfitvar)) ) THEN
        ctrl = INT(elsunc_parsoob_eval, KIND=i4)
        RETURN
     END IF
     ! -----------------------------------------------------------------------
     ! Calculate the weighted difference between fitted and measured spectrum.
     ! -----------------------------------------------------------------------
     CALL spectrum_solar ( &
          npoints, nfitvar, sol_wav_avg, locwvl(1:npoints), ymod(1:npoints), &
          fitvar(1:nfitvar) )
     ymod(1:npoints) = ( ymod(1:npoints) - currspec(1:npoints) ) * fitweights(1:npoints)

  CASE ( 2 )
     ! -------------------------------------------------
     ! Count the number of calls to the fitting function
     ! with request for the Jacobian and terminate if we
     ! exceed the allowed maximum (just to be safe!).
     ! -------------------------------------------------
     num_fitfunc_jacobi = num_fitfunc_jacobi + 1
     IF ( num_fitfunc_jacobi > max_fitfunc_calls ) THEN
        ctrl = INT(elsunc_infloop_eval, KIND=i4)
        RETURN
     END IF
     ! ---------------------------------------------------------------------
     ! The following sets up ELSUNC for numerical computation of the fitting
     ! function derivative. It is faster and more flexible than the original
     ! "manual" (AUTODIFF) scheme, and gives better fitting uncertainties.
     ! ---------------------------------------------------------------------
     ctrl = 0; RETURN

  CASE ( 3 )
     ! Calculate the spectrum, without weighting
     CALL spectrum_solar ( &
          npoints, nfitvar, sol_wav_avg, locwvl(1:npoints), ymod(1:npoints), &
          fitvar(1:nfitvar))

  CASE DEFAULT
     !WRITE (*,'(A,I4)') &
     !     "ERROR in function ELSUNC_SPECFIT_FUNC. Don't know how to handle CTRL = ", ctrl
  END SELECT

  RETURN
END SUBROUTINE specfit_func_sol


SUBROUTINE specfit_func ( fitvar, nfitvar, ymod, npoints, ctrl, dyda, mdy )

  !
  !     Calculates the spectrum and its derivatives for ELSUNC
  !
  ! NOTE: the variable DYDA as required here is the transpose of that 
  !       rquired for the Numerical Recipes
  !

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : elsunc_infloop_eval
  USE OMSAO_variables_module, ONLY : database, yn_doas, yn_smooth, &
       rad_wav_avg, fitwavs, fitweights, currspec, lobnd, upbnd, &
       num_fitfunc_calls, num_fitfunc_jacobi, max_fitfunc_calls
  USE OMSAO_radiance_ref_module, ONLY: yn_reference_fit

  IMPLICIT NONE

  INTEGER (KIND=i4),                              INTENT (IN)    :: nfitvar, npoints, mdy
  INTEGER (KIND=i4),                              INTENT (INOUT) :: ctrl
  REAL    (KIND=r8), DIMENSION (nfitvar),         INTENT (INOUT) :: fitvar
  REAL    (KIND=r8), DIMENSION (npoints),         INTENT (INOUT) :: ymod
  REAL    (KIND=r8), DIMENSION (npoints,nfitvar), INTENT (INOUT) :: dyda

  REAL    (KIND=r8), DIMENSION (npoints)       :: locwvl

  locwvl(1:npoints) = fitwavs(1:npoints)

  SELECT CASE ( ABS ( ctrl ) )
  CASE ( 1 )
     ! -------------------------------------------------
     ! Count the number of calls to the fitting function
     ! and terminate if we exceed the allowed maximum.
     ! -------------------------------------------------
     num_fitfunc_calls = num_fitfunc_calls + 1
     IF ( num_fitfunc_calls > max_fitfunc_calls ) THEN
        ctrl = INT(elsunc_infloop_eval, KIND=i4)
        RETURN
     END IF
     ! -------------------------------------------------------------
     ! Check whether any of the fitting variables are out of bounds.
     ! If so, indicate uncomputability and return. "Uncomputability"
     ! is indicated by setting  CTRL to < -10.
     ! (NOTE that this is slightly different from the description in
     ! the ELSUNC fitting module).
     ! -------------------------------------------------------------
     IF ( ANY(fitvar(1:nfitvar) < lobnd(1:nfitvar)) .OR. &
          ANY(fitvar(1:nfitvar) > upbnd(1:nfitvar)) ) THEN
        ! DON'T USE THIS! ctrl = INT(elsunc_parsoob_eval, KIND=i4)
        ctrl = -1
        RETURN
     END IF

     ! -----------------------------------------------------------------------
     ! Calculate the weighted difference between fitted and measured spectrum.
     ! -----------------------------------------------------------------------
     CALL spectrum_earthshine ( &
          npoints, nfitvar, rad_wav_avg, locwvl(1:npoints), ymod(1:npoints), &
          fitvar(1:nfitvar), database, yn_doas )
     !IF ( .NOT. yn_reference_fit ) THEN
     !   WRITE (99,'(3I6)') num_fitfunc_calls, nfitvar, npoints
     !   WRITE (99,'(1P100(E15.5:))') fitvar(1:nfitvar)
     !   DO i = 1, npoints
     !      WRITE (99,'(0PF12.4, 1P3E15.5)') locwvl(i), currspec(i), ymod(i), fitweights(i)
     !   END DO
     !END IF

     ymod(1:npoints) = ( ymod(1:npoints) - currspec(1:npoints) ) * fitweights(1:npoints)

  CASE ( 2 )
     ! -------------------------------------------------
     ! Count the number of calls to the fitting function
     ! with request for the Jacobian and terminate if we
     ! exceed the allowed maximum (just to be safe!).
     ! -------------------------------------------------
     num_fitfunc_jacobi = num_fitfunc_jacobi + 1
     IF ( num_fitfunc_jacobi > max_fitfunc_calls ) THEN
        ctrl = INT(elsunc_infloop_eval, KIND=i4)
        RETURN
     END IF
     ! ---------------------------------------------------------------------
     ! The following sets up ELSUNC for numerical computation of the fitting
     ! function derivative. It is faster and more flexible than the original
     ! "manual" (AUTODIFF) scheme, and gives better fitting uncertainties.
     ! ---------------------------------------------------------------------
     ctrl = 0; RETURN

  CASE ( 3 )
     ! Calculate the spectrum, without weighting
     CALL spectrum_earthshine ( &
          npoints, nfitvar, rad_wav_avg, locwvl(1:npoints), ymod(1:npoints), &
          fitvar(1:nfitvar), database, yn_doas )

  CASE DEFAULT
     !WRITE (*,'(A,I4)') &
     !     "ERROR in function ELSUNC_SPECFIT_FUNC. Don't know how to handle CTRL = ", ctrl
  END SELECT

  RETURN
END SUBROUTINE specfit_func


SUBROUTINE specfit_func_o3exp ( fitvar, nfitvar, ymod, npoints, ctrl, dyda, mdy )

  !
  !     Calculates the spectrum and its derivatives for ELSUNC
  !
  ! NOTE: the variable DYDA as required here is the transpose of that 
  !       rquired for the Numerical Recipes
  !

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : elsunc_infloop_eval
  USE OMSAO_variables_module, ONLY : &
       database, yn_doas, yn_smooth, rad_wav_avg, fitwavs, fitweights, currspec, &
       lobnd, upbnd, num_fitfunc_calls, num_fitfunc_jacobi, max_fitfunc_calls
  USE OMSAO_radiance_ref_module, ONLY: yn_reference_fit

  IMPLICIT NONE

  INTEGER (KIND=i4),                              INTENT (IN)    :: nfitvar, npoints, mdy
  INTEGER (KIND=i4),                              INTENT (INOUT) :: ctrl
  REAL    (KIND=r8), DIMENSION (nfitvar),         INTENT (INOUT) :: fitvar
  REAL    (KIND=r8), DIMENSION (npoints),         INTENT (INOUT) :: ymod
  REAL    (KIND=r8), DIMENSION (npoints,nfitvar), INTENT (INOUT) :: dyda

  REAL    (KIND=r8), DIMENSION (npoints)       :: locwvl

  locwvl(1:npoints) = fitwavs(1:npoints)

  SELECT CASE ( ABS ( ctrl ) )
  CASE ( 1 )
     ! -------------------------------------------------
     ! Count the number of calls to the fitting function
     ! and terminate if we exceed the allowed maximum.
     ! -------------------------------------------------
     num_fitfunc_calls = num_fitfunc_calls + 1
     IF ( num_fitfunc_calls > max_fitfunc_calls ) THEN
        ctrl = INT(elsunc_infloop_eval, KIND=i4)
        RETURN
     END IF
     ! -------------------------------------------------------------
     ! Check whether any of the fitting variables are out of bounds.
     ! If so, indicate uncomputability and return. "Uncomputability"
     ! is indicated by setting  CTRL to < -10.
     ! (NOTE that this is slightly different from the description in
     ! the ELSUNC fitting module).
     ! -------------------------------------------------------------
     IF ( ANY(fitvar(1:nfitvar) < lobnd(1:nfitvar)) .OR. &
          ANY(fitvar(1:nfitvar) > upbnd(1:nfitvar)) ) THEN
        ! DON'T USE THIS! ctrl = INT(elsunc_parsoob_eval, KIND=i4)
        ctrl = -1
        RETURN
     END IF

     ! -----------------------------------------------------------------------
     ! Calculate the weighted difference between fitted and measured spectrum.
     ! -----------------------------------------------------------------------
     CALL spectrum_earthshine_o3exp ( &
          npoints, nfitvar, rad_wav_avg, locwvl(1:npoints), ymod(1:npoints), &
          fitvar(1:nfitvar), database, yn_doas )

     !WRITE (*,'(1P100(E12.4:))') fitvar(1:nfitvar)
     !IF ( .NOT. yn_reference_fit ) THEN
     !   WRITE (99,'(3I6)') num_fitfunc_calls, nfitvar, npoints
     !   WRITE (99,'(1P100(E15.5:))') fitvar(1:nfitvar)
     !   DO i = 1, npoints
     !      WRITE (99,'(0PF12.4, 1P3E15.5)') locwvl(i), currspec(i), ymod(i), fitweights(i)
     !   END DO
     !END IF

     ymod(1:npoints) = ( ymod(1:npoints) - currspec(1:npoints) ) * fitweights(1:npoints)
  CASE ( 2 )
     ! -------------------------------------------------
     ! Count the number of calls to the fitting function
     ! with request for the Jacobian and terminate if we
     ! exceed the allowed maximum (just to be safe!).
     ! -------------------------------------------------
     num_fitfunc_jacobi = num_fitfunc_jacobi + 1
     IF ( num_fitfunc_jacobi > max_fitfunc_calls ) THEN
        ctrl = INT(elsunc_infloop_eval, KIND=i4)
        RETURN
     END IF
     ! ---------------------------------------------------------------------
     ! The following sets up ELSUNC for numerical computation of the fitting
     ! function derivative. It is faster and more flexible than the original
     ! "manual" (AUTODIFF) scheme, and gives better fitting uncertainties.
     ! ---------------------------------------------------------------------
     ctrl = 0; RETURN

  CASE ( 3 )
     ! Calculate the spectrum, without weighting
     CALL spectrum_earthshine_o3exp ( &
          npoints, nfitvar, rad_wav_avg, locwvl(1:npoints), ymod(1:npoints), &
          fitvar(1:nfitvar), database, yn_doas )

  CASE DEFAULT
     !WRITE (*,'(A,I4)') &
     !     "ERROR in function ELSUNC_SPECFIT_FUNC. Don't know how to handle CTRL = ", ctrl
  END SELECT

  RETURN
END SUBROUTINE specfit_func_o3exp


SUBROUTINE cubic_func ( x, afunc, ma )

  ! ***************************************************
  !
  !   Computes a third-order polynomial. Used in LFIT
  !
  ! ***************************************************

  USE OMSAO_precision_module

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER (KIND=i4), INTENT (IN) :: ma
  REAL    (KIND=r8), INTENT (IN) :: x

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=r8), DIMENSION (ma), INTENT (OUT) :: afunc
  
  ! ===============
  ! Local variables
  ! ===============
  INTEGER (KIND=i4) :: i

  afunc(1) = 1.0_r8
  DO i = 2, ma
     afunc(i) = afunc(i-1) * x
  END DO

  RETURN
END SUBROUTINE cubic_func

SUBROUTINE cubic_specfit ( a, na, y, m, ctrl, dyda, mdy )

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY : cubic_x, cubic_y, cubic_w

  IMPLICIT NONE

  ! Input parameters
  ! ================
  INTEGER (KIND=i4),                  INTENT (IN)  :: na, m, mdy
  REAL    (KIND=r8), DIMENSION (na),  INTENT (IN)  :: a

  ! Modified parameters
  ! ===================
  INTEGER (KIND=i4), INTENT (INOUT) :: ctrl

  ! Output parameters
  ! =================
  REAL (KIND=r8), DIMENSION (m),    INTENT (OUT)  :: y
  REAL (KIND=r8), DIMENSION (m,na), INTENT (OUT)  :: dyda

  ! Local variables
  ! ===============
  REAL (KIND=r8), DIMENSION (m) :: x, y0

  x  = cubic_x(1:m)
  y0 = a(1) + a(2)*x + a(3)*x*x + a(4)*x*x*x


  SELECT CASE ( ABS(ctrl) )
  CASE ( 1 )
     y  = ( y0 - cubic_y(1:m) ) / cubic_w(1:m)
  CASE ( 2 )
     dyda = 0.0_r8
     dyda(1:m,1) = 1.0_r8
     dyda(1:m,2) = x(1:m)
     dyda(1:m,3) = x(1:m)*x(1:m)
     dyda(1:m,4) = x(1:m)*x(1:m)*x(1:m)
  CASE ( 3 )
     ! This CASE is included to get the complete fitted spectrum
     y  = y0
  CASE DEFAULT
     !WRITE (*, '(A,I3)') "Don't know how to handle CTRL = ", ctrl
  END SELECT

  RETURN
END SUBROUTINE cubic_specfit

