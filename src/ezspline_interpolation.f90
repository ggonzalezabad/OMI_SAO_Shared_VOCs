SUBROUTINE ezspline_1d_interpolation ( n_in, x_in, y_in, n_out, x_out, y_out, errstat )

  USE EZspline_obj
  USE EZspline  
  USE OMSAO_precision_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                    INTENT (IN) :: n_in, n_out
  REAL    (KIND=r8), DIMENSION (n_in),  INTENT (IN) :: x_in, y_in
  REAL    (KIND=r8), DIMENSION (n_out), INTENT (IN) :: x_out

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4),                    INTENT (INOUT) :: errstat
  REAL    (KIND=r8), DIMENSION (n_out), INTENT (OUT)   :: y_out

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER (KIND=i4), PARAMETER :: ezs_r8 = SELECTED_REAL_KIND(12,100) ! real*8 kind

  REAL    (KIND=ezs_r8), DIMENSION (n_in)  :: f ! independent variable and function
  REAL    (KIND=ezs_r8), DIMENSION (n_out) :: z1, fz
  INTEGER (KIND=i4),     DIMENSION (2)     :: BCS1(2)
  TYPE (EZspline1_r8)                      :: spline_o ! 1-d EZspline object

  INTEGER (KIND=i4) :: locerrstat

  locerrstat     = pge_errstat_ok
  y_out(1:n_out) = 0.0_r8


  BCS1 = (/  0,  0 /) ! not a knot

  ! --------------------------
  ! Initialize/allocate memory
  ! --------------------------
  CALL EZspline_init (spline_o, n_in, BCS1, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     CALL Ezspline_free (spline_o, locerrstat)
     CALL error_check ( locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
          "EZspline_Init", vb_lev_default, errstat )
     RETURN
  END IF

  ! ---------------------------
  ! Set explicitely spline_o%x1
  ! ---------------------------
  spline_o%x1 = REAL(x_in(1:n_in), KIND=ezs_r8)

  ! ---------------------
  ! Assign function value
  ! ---------------------
  f(1:n_in) = REAL ( y_in(1:n_in), KIND=ezs_r8 )

  ! ----------------------------
  ! Setting up interpolation ...
  ! ----------------------------
  CALL EZspline_setup(spline_o, f, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     CALL Ezspline_free (spline_o, locerrstat)
     CALL error_check ( locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
          "EZspline_Setup", vb_lev_default, errstat )
     RETURN
  END IF

  ! -------------------
  ! Array interpolation
  ! -------------------
  z1 = REAL ( x_out, KIND=ezs_r8 )

  CALL EZspline_interp (spline_o, n_out, z1, fz, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     CALL Ezspline_free (spline_o, locerrstat)
     CALL error_check ( locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
          "EZspline_Interp", vb_lev_default, errstat )
     RETURN
  END IF

  ! -----------------------------------------------
  ! Assign interpolated results to output variables
  ! -----------------------------------------------
  y_out(1:n_out) = REAL ( fz(1:n_out), KIND=r8 )

  ! ---------------------------
  ! Clean up and free up memory
  ! ---------------------------
  CALL Ezspline_free (spline_o, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     CALL error_check ( locerrstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_INTERPOL, &
          "EZspline_free", vb_lev_default, errstat )
     RETURN
  END IF

  RETURN
END SUBROUTINE ezspline_1d_interpolation

SUBROUTINE ezspline_1d_setup_only ( n_in, x_in, y_in, spline_o, errstat )

  USE EZspline_obj
  USE EZspline  
  USE OMSAO_precision_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                    INTENT (IN) :: n_in
  REAL    (KIND=r8), DIMENSION (n_in),  INTENT (IN) :: x_in, y_in

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4),   INTENT (INOUT) :: errstat
  TYPE (EZspline1_r8), INTENT (OUT)   :: spline_o ! 1-d EZspline object

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER (KIND=i4), PARAMETER :: ezs_r8 = SELECTED_REAL_KIND(12,100) ! real*8 kind

  REAL    (KIND=ezs_r8), DIMENSION (n_in)  :: f ! independent variable and function
  INTEGER (KIND=i4),     DIMENSION (2)     :: BCS1(2)
  INTEGER (KIND=i4)                        :: locerrstat

  locerrstat     = pge_errstat_ok


  BCS1 = (/  0,  0 /) ! not a knot

  ! --------------------------
  ! Initialize/allocate memory
  ! --------------------------
  CALL EZspline_init (spline_o, n_in, BCS1, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     CALL Ezspline_free (spline_o, locerrstat)
     CALL error_check ( locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
          "EZspline_Init", vb_lev_default, errstat )
     RETURN
  END IF

  ! ---------------------------
  ! Set explicitely spline_o%x1
  ! ---------------------------
  spline_o%x1 = REAL(x_in(1:n_in), KIND=ezs_r8)

  ! ---------------------
  ! Assign function value
  ! ---------------------
  f(1:n_in) = REAL ( y_in(1:n_in), KIND=ezs_r8 )

  ! ----------------------------
  ! Setting up interpolation ...
  ! ----------------------------
  CALL EZspline_setup(spline_o, f, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     CALL Ezspline_free (spline_o, locerrstat)
     CALL error_check ( locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
          "EZspline_Setup", vb_lev_default, errstat )
     RETURN
  END IF


  RETURN
END SUBROUTINE ezspline_1d_setup_only


SUBROUTINE ezspline_1d_ipol_only ( spline_o, n_out, x_out, y_out, errstat )

  USE EZspline_obj
  USE EZspline  
  USE OMSAO_precision_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  TYPE (EZspline1_r8),                  INTENT (IN) :: spline_o ! 1-d EZspline object
  INTEGER (KIND=i4),                    INTENT (IN) :: n_out
  REAL    (KIND=r8), DIMENSION (n_out), INTENT (IN) :: x_out

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4),                    INTENT (INOUT) :: errstat
  REAL    (KIND=r8), DIMENSION (n_out), INTENT (OUT)   :: y_out

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER (KIND=i4), PARAMETER :: ezs_r8 = SELECTED_REAL_KIND(12,100) ! real*8 kind

  REAL    (KIND=ezs_r8), DIMENSION (n_out) :: z1, fz
  INTEGER (KIND=i4),     DIMENSION (2)     :: BCS1(2)
  INTEGER (KIND=i4)                        :: locerrstat

  locerrstat     = pge_errstat_ok
  y_out(1:n_out) = 0.0_r8


  BCS1 = (/  0,  0 /) ! not a knot

  ! -------------------
  ! Array interpolation
  ! -------------------
  z1 = REAL ( x_out, KIND=ezs_r8 )

  CALL EZspline_interp (spline_o, n_out, z1, fz, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     CALL Ezspline_free (spline_o, locerrstat)
     CALL error_check ( locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
          "EZspline_Interp", vb_lev_default, errstat )
     RETURN
  END IF

  ! -----------------------------------------------
  ! Assign interpolated results to output variables
  ! -----------------------------------------------
  y_out(1:n_out) = REAL ( fz(1:n_out), KIND=r8 )


  RETURN
END SUBROUTINE ezspline_1d_ipol_only


SUBROUTINE ezspline_1d_memfree ( spline_o, errstat )

  USE EZspline_obj
  USE EZspline  
  USE OMSAO_precision_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ----------------------
  ! Input/output variables
  ! ----------------------
  TYPE (EZspline1_r8), INTENT (INOUT) :: spline_o ! 1-d EZspline object
  INTEGER (KIND=i4),   INTENT (INOUT) :: errstat

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER (KIND=i4) :: locerrstat

  locerrstat     = pge_errstat_ok
  ! ---------------------------
  ! Clean up and free up memory
  ! ---------------------------
  CALL Ezspline_free (spline_o, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     CALL error_check ( locerrstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_INTERPOL, &
          "EZspline_free", vb_lev_default, errstat )
     RETURN
  END IF

  RETURN
END SUBROUTINE ezspline_1d_memfree


SUBROUTINE ezspline_2d_interpolation ( n1, n2, x1, x2, z_in, m1, m2, y1, y2, z_out, errstat )

  USE EZspline_obj
  USE EZspline  
  USE OMSAO_precision_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                    INTENT (IN) :: n1, n2, m1, m2
  REAL    (KIND=r8), DIMENSION (n1),    INTENT (IN) :: x1
  REAL    (KIND=r8), DIMENSION (n2),    INTENT (IN) :: x2
  REAL    (KIND=r8), DIMENSION (n1,n2), INTENT (IN) :: z_in
  REAL    (KIND=r8), DIMENSION (m1),    INTENT (IN) :: y1
  REAL    (KIND=r8), DIMENSION (m2),    INTENT (IN) :: y2

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4),                    INTENT (OUT) :: errstat
  REAL    (KIND=r8), DIMENSION (m1,m2), INTENT (OUT) :: z_out


  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER (KIND=i4), PARAMETER :: ezs_r8 = SELECTED_REAL_KIND(12,100) ! real*8 kind

  REAL    (KIND=ezs_r8), DIMENSION (n1,n2) :: f   ! =z_in

  REAL    (KIND=ezs_r8), DIMENSION (m1)    :: ey1 ! =z1
  REAL    (KIND=ezs_r8), DIMENSION (m2)    :: ey2 ! =z2
  REAL    (KIND=ezs_r8), DIMENSION (m1,m2) :: fz  ! =z_out

  INTEGER (KIND=i4),  DIMENSION (2) :: BCS1(2), BCS2(2)
  TYPE (EZspline2_r8)               :: spline_o ! 2-d EZspline object

  INTEGER (KIND=i4) :: locerrstat

  locerrstat = pge_errstat_ok

  BCS1 = (/  0,  0 /) ! not a knot
  BCS2 = (/  0,  0 /) ! not a knot

  ! --------------------------
  ! Initialize/allocate memory
  ! --------------------------
  CALL EZspline_init (spline_o, n1, n2, BCS1, BCS2, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     errstat = pge_errstat_error; CALL Ezspline_free (spline_o, locerrstat); RETURN
  END IF

  ! ---------------------------
  ! Set explicitely spline_o%x1
  ! ---------------------------
  spline_o%x1 = REAL(x1(1:n1), KIND=ezs_r8)
  spline_o%x2 = REAL(x2(1:n2), KIND=ezs_r8)

  ! ---------------------
  ! Assign function value
  ! ---------------------
  f(1:n1,1:n2) = REAL ( z_in(1:n1,1:n2), KIND=ezs_r8 )

  ! ----------------------------
  ! Setting up interpolation ...
  ! ----------------------------
  CALL EZspline_setup(spline_o, f, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     errstat = pge_errstat_error; CALL Ezspline_free (spline_o, locerrstat); RETURN
  END IF

  ey1 = REAL ( y1, KIND=ezs_r8 )
  ey2 = REAL ( y2, KIND=ezs_r8 )

  CALL EZspline_interp (spline_o, m1, m2, ey1, ey2, fz, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     errstat = pge_errstat_error; CALL Ezspline_free (spline_o, locerrstat); RETURN
  END IF

  ! -----------------------------------------------
  ! Assign interpolated results to output variables
  ! -----------------------------------------------
  z_out(1:m1,1:m2) = REAL (fz(1:m1,1:m2), KIND=r8)

  ! ---------------------------
  ! Clean up and free up memory
  ! ---------------------------
  CALL Ezspline_free (spline_o, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     errstat = pge_errstat_warning; RETURN
  END IF

  RETURN
END SUBROUTINE ezspline_2d_interpolation

SUBROUTINE ezspline_3d_interpolation ( &
     n1, n2, n3, x1, x2, x3, z_in, m1, m2, m3, y1, y2, y3, z_out, errstat )

  USE EZspline_obj
  USE EZspline  
  USE OMSAO_precision_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                       INTENT (IN) :: n1, n2, n3, m1, m2, m3
  REAL    (KIND=r8), DIMENSION (n1),       INTENT (IN) :: x1
  REAL    (KIND=r8), DIMENSION (n2),       INTENT (IN) :: x2
  REAL    (KIND=r8), DIMENSION (n3),       INTENT (IN) :: x3
  REAL    (KIND=r8), DIMENSION (n1,n2,n3), INTENT (IN) :: z_in
  REAL    (KIND=r8), DIMENSION (m1),       INTENT (IN) :: y1
  REAL    (KIND=r8), DIMENSION (m2),       INTENT (IN) :: y2
  REAL    (KIND=r8), DIMENSION (m3),       INTENT (IN) :: y3

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4),                       INTENT (OUT) :: errstat
  REAL    (KIND=r8), DIMENSION (m1,m2,m3), INTENT (OUT) :: z_out

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER (KIND=i4), PARAMETER :: ezs_r8 = SELECTED_REAL_KIND(12,100) ! real*8 kind

  REAL    (KIND=ezs_r8), DIMENSION (n1,n2,n3) :: f   ! =z_in

  REAL    (KIND=ezs_r8), DIMENSION (m1)       :: ey1 ! =z1
  REAL    (KIND=ezs_r8), DIMENSION (m2)       :: ey2 ! =z2
  REAL    (KIND=ezs_r8), DIMENSION (m3)       :: ey3 ! =z3
  REAL    (KIND=ezs_r8), DIMENSION (m1,m2,m3) :: fz  ! =z_out

  INTEGER (KIND=i4),  DIMENSION (2) :: bcs1(2), bcs2(2), bcs3(2)
  TYPE (EZspline3_r8)               :: spline_o ! 3-d EZspline object

  INTEGER (KIND=i4) :: locerrstat

  locerrstat = pge_errstat_ok

  bcs1 = (/  0,  0 /) ! not a knot
  bcs2 = (/  0,  0 /) ! not a knot
  bcs3 = (/  0,  0 /) ! not a knot

  ! --------------------------
  ! Initialize/allocate memory
  ! --------------------------
  CALL EZspline_init (spline_o, n1, n2, n3, bcs1, bcs2, bcs3, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     errstat = pge_errstat_error; CALL Ezspline_free (spline_o, locerrstat); RETURN
  END IF

  !spline_o%isHermite = 1

  ! ---------------------------
  ! Set explicitely spline_o%x1
  ! ---------------------------
  spline_o%x1 = REAL(x1(1:n1), KIND=ezs_r8)
  spline_o%x2 = REAL(x2(1:n2), KIND=ezs_r8)
  spline_o%x3 = REAL(x3(1:n3), KIND=ezs_r8)

  ! ---------------------
  ! Assign function value
  ! ---------------------
  f(1:n1,1:n2,1:n3) = REAL ( z_in(1:n1,1:n2,1:n3), KIND=ezs_r8 )

  ! ----------------------------
  ! Setting up interpolation ...
  ! ----------------------------
  CALL EZspline_setup(spline_o, f, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     errstat = pge_errstat_error; CALL Ezspline_free (spline_o, locerrstat); RETURN
  END IF

  ey1 = REAL ( y1, KIND=ezs_r8 )
  ey2 = REAL ( y2, KIND=ezs_r8 )
  ey3 = REAL ( y3, KIND=ezs_r8 )

  CALL EZspline_interp (spline_o, m1, m2, m3, ey1, ey2, ey3, fz, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     errstat = pge_errstat_error; CALL Ezspline_free (spline_o, locerrstat); RETURN
  END IF

  ! -----------------------------------------------
  ! Assign interpolated results to output variables
  ! -----------------------------------------------
  z_out(1:m1,1:m2,1:m3) = REAL (fz(1:m1,1:m2,1:m3), KIND=r8)

  ! ---------------------------
  ! Clean up and free up memory
  ! ---------------------------
  CALL Ezspline_free (spline_o, locerrstat)
  IF ( locerrstat /= pge_errstat_ok ) THEN
     errstat = pge_errstat_warning; RETURN
  END IF

  RETURN
END SUBROUTINE ezspline_3d_interpolation
