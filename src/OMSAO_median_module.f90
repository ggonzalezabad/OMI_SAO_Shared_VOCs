MODULE OMSAO_median_module

  USE OMSAO_precision_module, ONLY: i4, r8

  IMPLICIT NONE

CONTAINS


  ! --------------------------------------------------------------------
  ! REAL FUNCTION  Median() :
  !    This function receives an array X of N entries, copies its value
  !    to a local array Temp(), sorts Temp() and computes the median.
  !    The returned value is of REAL type.
  ! --------------------------------------------------------------------

  REAL (KIND=r8) FUNCTION median (n, x)
  
    IMPLICIT  NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                INTENT (IN) :: n
    REAL    (KIND=r8), DIMENSION (n), INTENT (IN) :: x

    ! ---------------
    ! Local variables
    ! ---------------
    REAL    (KIND=r8), DIMENSION (n) :: tmp

    tmp(1:n) = x(1:n)               ! make a copy

    CALL sort(n, tmp)               ! sort the copy
    IF (MOD(n,2) == 0) THEN           ! compute the median
       median = (tmp(n/2) + tmp(n/2+1)) / 2.0_r8
    ELSE
       median = tmp(n/2+1)
    END IF

    RETURN
  END FUNCTION  median



  ! --------------------------------------------------------------------
  ! REAL FUNCTION  WeightMean() :
  !    This function receives an array X of N entries, computes the
  !    standard deviation, and subsequently the mean of those values
  !    that lie within M sigma of the regular mean values
  !    The returned value is of REAL type.
  ! --------------------------------------------------------------------

  REAL (KIND=r8) FUNCTION weightmean (n, x, m)
  
    IMPLICIT  NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                INTENT (IN) :: n, m
    REAL    (KIND=r8), DIMENSION (n), INTENT (IN) :: x

    ! ---------------
    ! Local variables
    ! ---------------
    REAL    (KIND=r8) :: mean, smean, sdev
    INTEGER (KIND=i4) :: i, k

    weightmean = 0.0_r8

    IF ( n < 0 ) RETURN

    mean = SUM(x(1:n)) / REAL(n, KIND=r8)
    sdev = 0.0_r8
    DO i = 1, n
       sdev = sdev + (x(i)-mean)*(x(i)-mean)
    END DO
    sdev = SQRT( sdev ) / REAL(n, KIND=r8)

    smean = 0.0_r8 ; k = 0
    DO i = 1, n
       IF ( ( x(i) >= mean - REAL(m,KIND=r8)*sdev) .AND. &
            ( x(i) <= mean + REAL(m,KIND=r8)*sdev)        ) THEN
          smean = smean + x(i)
          k     = k + 1
       END IF
    END DO
    IF ( k > 0 ) smean = smean / REAL(k, KIND=r8)

    weightmean = smean

    RETURN
  END FUNCTION  weightmean


  ! --------------------------------------------------------------------
  ! INTEGER FUNCTION  FindMinimum():
  !    This function returns the location of the minimum in the section
  ! between Start and End.
  ! --------------------------------------------------------------------

  INTEGER (KIND=i4) FUNCTION FindMinimumPosition(n, x, start, end) RESULT (minpos)

    IMPLICIT  NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                 INTENT (IN) :: n
    REAL    (KIND=r8), DIMENSION(1:n), INTENT (IN) :: x
    INTEGER (KIND=i4),                 INTENT (IN) :: start, end

    ! ---------------
    ! Result variable
    ! ---------------
    INTEGER (KIND=i4) :: location

    ! ---------------
    ! Local variables
    ! ---------------
    REAL    (KIND=r8) :: minimum
    INTEGER (KIND=i4) :: i

    minimum  = x(start)             ! assume the first is the min
    location = start                ! record its position
    DO i = start+1, end             ! start with next elements
       IF ( x(i) < minimum ) THEN   !   if x(i) less than the min?
          minimum  = x(i)           !      Yes, a new minimum found
          location = i              !      record its position
       END IF
    END DO

    minpos = location

  END FUNCTION  FindMinimumPosition


  ! --------------------------------------------------------------------
  ! SUBROUTINE  Swap():
  !    This subroutine swaps the values of its two formal arguments.
  ! --------------------------------------------------------------------

  SUBROUTINE  swap(a, b)

    IMPLICIT  NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), INTENT (INOUT) :: a, b

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: tmp

    tmp = a
    a   = b
    b   = tmp

    RETURN
  END SUBROUTINE  swap


  ! --------------------------------------------------------------------
  ! SUBROUTINE  Sort():
  !    This subroutine receives an array x() and sorts it 
  !    into ascending order.
  ! --------------------------------------------------------------------

  SUBROUTINE  sort(n, x)

    IMPLICIT  NONE
    
    ! --------------
    ! Input variable
    ! --------------
    INTEGER (KIND=i4), INTENT (IN) :: n

    ! -----------------
    ! Modified variable
    ! -----------------
    REAL (KIND=r8), DIMENSION(n), INTENT (INOUT) :: x

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: i, minpos


    DO i = 1, n-1                                   ! except for the last
       minpos = FindMinimumPosition ( n, x, i, n )  ! find min from this to last
       CALL swap(x(i), x(minpos))                   ! swap this and the minimum
    END DO

    RETURN
  END SUBROUTINE  sort

END MODULE  OMSAO_median_module
