!*DECK DAVINT
SUBROUTINE DAVINT (X, Y, N, XLO, XUP, ANS, IERR)
  !***BEGIN PROLOGUE  DAVINT
  !***PURPOSE  Integrate a function tabulated at arbitrarily spaced
  !            abscissas using overlapping parabolas.
  !***LIBRARY   SLATEC
  !***CATEGORY  H2A1B2
  !***TYPE      DOUBLE PRECISION (AVINT-S, DAVINT-D)
  !***KEYWORDS  INTEGRATION, QUADRATURE, TABULATED DATA
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract
  !         DAVINT integrates a function tabulated at arbitrarily spaced
  !         abscissas.  The limits of integration need not coincide
  !         with the tabulated abscissas.
  !
  !         A method of overlapping parabolas fitted to the data is used
  !         provided that there are at least 3 abscissas between the
  !         limits of integration.  DAVINT also handles two special cases.
  !         If the limits of integration are equal, DAVINT returns a
  !         result of zero regardless of the number of tabulated values.
  !         If there are only two function values, DAVINT uses the
  !         trapezoid rule.
  !
  !     Description of Parameters
  !         The user must dimension all arrays appearing in the call list
  !              X(N), Y(N)
  !
  !         Input--
  !      X    - DOUBLE PRECISION array of abscissas, which must be in
  !             increasing order.
  !      Y    - DOUBLE PRECISION array of function values. i.e.,
  !                Y(I)=FUNC(X(I))
  !      N    - The integer number of function values supplied.
  !                N .GE. 2 unless XLO = XUP.
  !      XLO  - DOUBLE PRECISION lower limit of integration
  !      XUP  - DOUBLE PRECISION upper limit of integration.  Must have
  !              XLO.LE.XUP
  !
  !         Output--
  !      ANS  - Double Precision computed approximate value of integral
  !      IERR - A status code
  !           --Normal Code
  !                =1 Means the requested integration was performed.
  !           --Abnormal Codes
  !                =2 Means XUP was less than XLO.
  !                =3 Means the number of X(I) between XLO and XUP
  !                   (inclusive) was less than 3 and neither of the two
  !                   special cases described in the abstract occurred.
  !                   No integration was performed.
  !                =4 Means the restriction X(I+1).GT.X(I) was violated.
  !                =5 Means the number N of function values was .lt. 2.
  !                   ANS is set to zero if IERR=2,3,4,or 5.
  !
  !    DAVINT is documented completely in SC-M-69-335
  !    Original program from *Numerical Integration* by Davis & Rabinowitz
  !    Adaptation and modifications by Rondall E Jones.
  !
  !***REFERENCES  R. E. Jones, Approximate integrator of functions
  !                 tabulated at arbitrarily spaced abscissas,
  !                 Report SC-M-69-335, Sandia Laboratories, 1969.
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   690901  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DAVINT
  !

  USE OMSAO_precision_module
  IMPLICIT NONE

  INTEGER I, IERR, INLFT, INRT, ISTART, ISTOP, N
  DOUBLE PRECISION A, ANS, B, C, CA, CB, CC, FL, FR, R3, RP5,  &
       SLOPE, SUM, SYL, SYL2, SYL3, SYU, SYU2, SYU3, TERM1, TERM2,  &
       TERM3, X1, X12, X13, X2, X23, X3, XLO, XUP

  DOUBLE PRECISION x(n), y(n)
  !DIMENSION X(*),Y(*)

  !     BEGIN BLOCK PERMITTING ...EXITS TO 190
  !        BEGIN BLOCK PERMITTING ...EXITS TO 180
  !***FIRST EXECUTABLE STATEMENT  DAVINT
  IERR = 1
  ANS = 0.0D0
  IF (XLO .GT. XUP) GO TO 160
  IF (XLO .EQ. XUP) GO TO 150
  IF (N .GE. 2) GO TO 10
  IERR = 5
  CALL XERMSG ('SLATEC', 'DAVINT',  &
       'LESS THAN TWO FUNCTION VALUES WERE SUPPLIED.',  &
       4, 1)
  !     ...............EXIT
  GO TO 190
10 CONTINUE
  DO I = 2, N
     !        ............EXIT
     IF (X(I) .LE. X(I-1)) GO TO 180
     !                 ...EXIT
     IF (X(I) .GT. XUP) GO TO 30
  END DO
30 CONTINUE
  IF (N .GE. 3) GO TO 40
  !
  !                    SPECIAL N=2 CASE
  SLOPE = (Y(2) - Y(1))/(X(2) - X(1))
  FL = Y(1) + SLOPE*(XLO - X(1))
  FR = Y(2) + SLOPE*(XUP - X(2))
  ANS = 0.5D0*(FL + FR)*(XUP - XLO)
  !     ...............EXIT
  GO TO 190
40 CONTINUE
  IF (X(N-2) .GE. XLO) GO TO 50
  IERR = 3
  CALL XERMSG ('SLATEC', 'DAVINT',  &
       'THERE WERE LESS THAN THREE FUNCTION VALUES ' //  &
       'BETWEEN THE LIMITS OF INTEGRATION.', 4, 1)
  !     ...............EXIT
  GO TO 190
50 CONTINUE
  IF (X(3) .LE. XUP) GO TO 60
  IERR = 3
  CALL XERMSG ('SLATEC', 'DAVINT',  &
       'THERE WERE LESS THAN THREE FUNCTION VALUES ' //  &
       'BETWEEN THE LIMITS OF INTEGRATION.', 4, 1)
  !     ...............EXIT
  GO TO 190
60 CONTINUE
  I = 1
70 IF (X(I) .GE. XLO) GO TO 80
  I = I + 1
  GO TO 70
80 CONTINUE
  INLFT = I
  I = N
90 IF (X(I) .LE. XUP) GO TO 100
  I = I - 1
  GO TO 90
100 CONTINUE
  INRT = I
  IF ((INRT - INLFT) .GE. 2) GO TO 110
  IERR = 3
  CALL XERMSG ('SLATEC', 'DAVINT',  &
       'THERE WERE LESS THAN THREE FUNCTION VALUES ' //  &
       'BETWEEN THE LIMITS OF INTEGRATION.', 4, 1)
  !     ...............EXIT
  GO TO 190
110 CONTINUE
  ISTART = INLFT
  IF (INLFT .EQ. 1) ISTART = 2
  ISTOP = INRT
  IF (INRT .EQ. N) ISTOP = N - 1
  !
  R3 = 3.0D0
  RP5 = 0.5D0
  SUM = 0.0D0
  SYL = XLO
  SYL2 = SYL*SYL
  SYL3 = SYL2*SYL
  !
  DO I = ISTART, ISTOP
     X1 = X(I-1)
     X2 = X(I)
     X3 = X(I+1)
     X12 = X1 - X2
     X13 = X1 - X3
     X23 = X2 - X3
     TERM1 = Y(I-1)/(X12*X13)
     TERM2 = -Y(I)/(X12*X23)
     TERM3 = Y(I+1)/(X13*X23)
     A = TERM1 + TERM2 + TERM3
     B = -(X2 + X3)*TERM1 - (X1 + X3)*TERM2  &
          - (X1 + X2)*TERM3
     C = X2*X3*TERM1 + X1*X3*TERM2 + X1*X2*TERM3
     IF (I .GT. ISTART) GO TO 120
     CA = A
     CB = B
     CC = C
     GO TO 130
120  CONTINUE
     CA = 0.5D0*(A + CA)
     CB = 0.5D0*(B + CB)
     CC = 0.5D0*(C + CC)
130  CONTINUE
     SYU = X2
     SYU2 = SYU*SYU
     SYU3 = SYU2*SYU
     SUM = SUM + CA*(SYU3 - SYL3)/R3  &
          + CB*RP5*(SYU2 - SYL2) + CC*(SYU - SYL)
     CA = A
     CB = B
     CC = C
     SYL = SYU
     SYL2 = SYU2
     SYL3 = SYU3
  END DO
140 CONTINUE
  SYU = XUP
  ANS = SUM + CA*(SYU**3 - SYL3)/R3  &
       + CB*RP5*(SYU**2 - SYL2) + CC*(SYU - SYL)
150 CONTINUE
  GO TO 170
160 CONTINUE
  IERR = 2
  CALL XERMSG ('SLATEC', 'DAVINT',  &
       'THE UPPER LIMIT OF INTEGRATION WAS NOT GREATER ' //  &
       'THAN THE LOWER LIMIT.', 4, 1)
170 CONTINUE
  !     ......EXIT
  GO TO 190
180 CONTINUE
  IERR = 4
  CALL XERMSG ('SLATEC', 'DAVINT',  &
       'THE ABSCISSAS WERE NOT STRICTLY INCREASING.  MUST HAVE ' //  &
       'X(I-1) .LT. X(I) FOR ALL I.', 4, 1)
190 CONTINUE
  RETURN
END SUBROUTINE DAVINT

!*DECK FDUMP
SUBROUTINE FDUMP
  !***BEGIN PROLOGUE  FDUMP
  !***PURPOSE  Symbolic dump (should be locally written).
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3
  !***TYPE      ALL (FDUMP-A)
  !***KEYWORDS  ERROR, XERMSG
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !        ***Note*** Machine Dependent Routine
  !        FDUMP is intended to be replaced by a locally written
  !        version which produces a symbolic dump.  Failing this,
  !        it should be replaced by a version which prints the
  !        subprogram nesting list.  Note that this dump must be
  !        printed on each of up to five files, as indicated by the
  !        XGETUA routine.  See XSETUA and XGETUA for details.
  !
  !     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  FDUMP
  !***FIRST EXECUTABLE STATEMENT  FDUMP

  USE OMSAO_precision_module
  IMPLICIT NONE

  RETURN
END SUBROUTINE FDUMP

!*DECK I1MACH
INTEGER FUNCTION I1MACH (I)
  !***BEGIN PROLOGUE  I1MACH
  !***PURPOSE  Return integer machine dependent constants.
  !***LIBRARY   SLATEC
  !***CATEGORY  R1
  !***TYPE      INTEGER (I1MACH-I)
  !***KEYWORDS  MACHINE CONSTANTS
  !***AUTHOR  Fox, P. A., (Bell Labs)
  !           Hall, A. D., (Bell Labs)
  !           Schryer, N. L., (Bell Labs)
  !***DESCRIPTION
  !
  !   I1MACH can be used to obtain machine-dependent parameters for the
  !   local machine environment.  It is a function subprogram with one
  !   (input) argument and can be referenced as follows:
  !
  !        K = I1MACH(I)
  !
  !   where I=1,...,16.  The (output) value of K above is determined by
  !   the (input) value of I.  The results for various values of I are
  !   discussed below.
  !
  !   I/O unit numbers:
  !     I1MACH( 1) = the standard input unit.
  !     I1MACH( 2) = the standard output unit.
  !     I1MACH( 3) = the standard punch unit.
  !     I1MACH( 4) = the standard error message unit.
  !
  !   Words:
  !     I1MACH( 5) = the number of bits per integer storage unit.
  !     I1MACH( 6) = the number of characters per integer storage unit.
  !
  !   Integers:
  !     assume integers are represented in the S-digit, base-A form
  !
  !                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
  !
  !                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
  !     I1MACH( 7) = A, the base.
  !     I1MACH( 8) = S, the number of base-A digits.
  !     I1MACH( 9) = A**S - 1, the largest magnitude.
  !
  !   Floating-Point Numbers:
  !     Assume floating-point numbers are represented in the T-digit,
  !     base-B form
  !                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
  !
  !                where 0 .LE. X(I) .LT. B for I=1,...,T,
  !                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
  !     I1MACH(10) = B, the base.
  !
  !   Single-Precision:
  !     I1MACH(11) = T, the number of base-B digits.
  !     I1MACH(12) = EMIN, the smallest exponent E.
  !     I1MACH(13) = EMAX, the largest exponent E.
  !
  !   Double-Precision:
  !     I1MACH(14) = T, the number of base-B digits.
  !     I1MACH(15) = EMIN, the smallest exponent E.
  !     I1MACH(16) = EMAX, the largest exponent E.
  !
  !   To alter this function for a particular environment, the desired
  !   set of DATA statements should be activated by removing the C from
  !   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
  !   checked for consistency with the local operating system.
  !
  !***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
  !                 a portable library, ACM Transactions on Mathematical
  !                 Software 4, 2 (June 1978), pp. 177-188.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   891012  Added VAX G-floating constants.  (WRB)
  !   891012  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900618  Added DEC RISC constants.  (WRB)
  !   900723  Added IBM RS 6000 constants.  (WRB)
  !   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
  !           (RWC)
  !   910710  Added HP 730 constants.  (SMR)
  !   911114  Added Convex IEEE constants.  (WRB)
  !   920121  Added SUN -r8 compiler option constants.  (WRB)
  !   920229  Added Touchstone Delta i860 constants.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !   920625  Added Convex -p8 and -pd8 compiler option constants.
  !           (BKS, WRB)
  !   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
  !   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
  !           options.  (DWL, RWC and WRB).
  !***END PROLOGUE  I1MACH
  !

  USE OMSAO_precision_module
  IMPLICIT NONE

  INTEGER :: i

  INTEGER IMACH(16),OUTPUT
  SAVE IMACH
  EQUIVALENCE (IMACH(4),OUTPUT)

  !     ==================================================================
  !     MACHINE CONSTANTS, some  modeled on those for for the DEC ALPHA,
  !     others computed using F90 intrinsic routines. Note that the latter
  !     make this subroutine incompatible with non-F90/F95 compilers.
  !     ==================================================================
  !     ( T. Kurosu, 17 June, 2001)
  !     ==================================================================

  IMACH( 1) = 5  !?
  IMACH( 2) = 6  !?
  IMACH( 3) = 6  !?
  IMACH( 4) = 6  !?
  IMACH( 5) = BIT_SIZE(1)
  IMACH( 6) = 4  !?
  IMACH( 7) = RADIX(1)
  IMACH( 8) = DIGITS(1)
  IMACH( 9) = HUGE(1)
  IMACH(10) = RADIX(1.0)
  IMACH(11) = DIGITS(1.0)
  IMACH(12) = MINEXPONENT(1.0)
  IMACH(13) = MAXEXPONENT(1.0)
  IMACH(14) = DIGITS(1.0D0)
  IMACH(15) = MINEXPONENT(1.0D0)
  IMACH(16) = MAXEXPONENT(1.0D0)

  !***FIRST EXECUTABLE STATEMENT  I1MACH
  IF (I < 1  .OR.  I > 16) GO TO 10
  !
  I1MACH = IMACH(I)
  RETURN
  !
10 CONTINUE
  WRITE (UNIT = OUTPUT, FMT = 9000)
9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
  !
  !     CALL FDUMP
  !
  STOP
END FUNCTION I1MACH

!*DECK J4SAVE
INTEGER FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
  !***BEGIN PROLOGUE  J4SAVE
  !***SUBSIDIARY
  !***PURPOSE  Save or recall global variables needed by error
  !            handling routines.
  !***LIBRARY   SLATEC (XERROR)
  !***TYPE      INTEGER (J4SAVE-I)
  !***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract
  !        J4SAVE saves and recalls several global variables needed
  !        by the library error handling routines.
  !
  !     Description of Parameters
  !      --Input--
  !        IWHICH - Index of item desired.
  !                = 1 Refers to current error number.
  !                = 2 Refers to current error control flag.
  !                = 3 Refers to current unit number to which error
  !                    messages are to be sent.  (0 means use standard.)
  !                = 4 Refers to the maximum number of times any
  !                     message is to be printed (as set by XERMAX).
  !                = 5 Refers to the total number of units to which
  !                     each error message is to be written.
  !                = 6 Refers to the 2nd unit for error messages
  !                = 7 Refers to the 3rd unit for error messages
  !                = 8 Refers to the 4th unit for error messages
  !                = 9 Refers to the 5th unit for error messages
  !        IVALUE - The value to be set for the IWHICH-th parameter,
  !                 if ISET is .TRUE. .
  !        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
  !                 given the value, IVALUE.  If ISET=.FALSE., the
  !                 IWHICH-th parameter will be unchanged, and IVALUE
  !                 is a dummy parameter.
  !      --Output--
  !        The (old) value of the IWHICH-th parameter will be returned
  !        in the function value, J4SAVE.
  !
  !***SEE ALSO  XERMSG
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900205  Minor modifications to prologue.  (WRB)
  !   900402  Added TYPE section.  (WRB)
  !   910411  Added KEYWORDS section.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  J4SAVE

  USE OMSAO_precision_module
  IMPLICIT NONE

  LOGICAL ISET
  INTEGER               :: iwhich, ivalue
  INTEGER, DIMENSION(9) :: IPARAM
  SAVE IPARAM
  DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
  DATA IPARAM(5)/1/
  DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
  !***FIRST EXECUTABLE STATEMENT  J4SAVE
  J4SAVE = IPARAM(IWHICH)
  IF (ISET) IPARAM(IWHICH) = IVALUE
  RETURN
END FUNCTION J4SAVE

!*DECK XERCNT
SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
  !***BEGIN PROLOGUE  XERCNT
  !***SUBSIDIARY
  !***PURPOSE  Allow user control over handling of errors.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XERCNT-A)
  !***KEYWORDS  ERROR, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract
  !        Allows user control over handling of individual errors.
  !        Just after each message is recorded, but before it is
  !        processed any further (i.e., before it is printed or
  !        a decision to abort is made), a call is made to XERCNT.
  !        If the user has provided his own version of XERCNT, he
  !        can then override the value of KONTROL used in processing
  !        this message by redefining its value.
  !        KONTRL may be set to any value from -2 to 2.
  !        The meanings for KONTRL are the same as in XSETF, except
  !        that the value of KONTRL changes only for this message.
  !        If KONTRL is set to a value outside the range from -2 to 2,
  !        it will be moved back into that range.
  !
  !     Description of Parameters
  !
  !      --Input--
  !        LIBRAR - the library that the routine is in.
  !        SUBROU - the subroutine that XERMSG is being called from
  !        MESSG  - the first 20 characters of the error message.
  !        NERR   - same as in the call to XERMSG.
  !        LEVEL  - same as in the call to XERMSG.
  !        KONTRL - the current value of the control flag as set
  !                 by a call to XSETF.
  !
  !      --Output--
  !        KONTRL - the new value of KONTRL.  If KONTRL is not
  !                 defined, it will remain at its original value.
  !                 This changed value of control affects only
  !                 the current occurrence of the current message.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900206  Routine changed from user-callable to subsidiary.  (WRB)
  !   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
  !           names, changed routine name from XERCTL to XERCNT.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERCNT

  USE OMSAO_precision_module
  IMPLICIT NONE

  INTEGER :: kontrl, level, nerr
  CHARACTER*(*) LIBRAR, SUBROU, MESSG
  !***FIRST EXECUTABLE STATEMENT  XERCNT
  RETURN
END SUBROUTINE XERCNT

!*DECK XERHLT
SUBROUTINE XERHLT (MESSG)
  !***BEGIN PROLOGUE  XERHLT
  !***SUBSIDIARY
  !***PURPOSE  Abort program execution and print error message.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XERHLT-A)
  !***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract
  !        ***Note*** machine dependent routine
  !        XERHLT aborts the execution of the program.
  !        The error message causing the abort is given in the calling
  !        sequence, in case one needs it for printing on a dayfile,
  !        for example.
  !
  !     Description of Parameters
  !        MESSG is as in XERMSG.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900206  Routine changed from user-callable to subsidiary.  (WRB)
  !   900510  Changed calling sequence to delete length of character
  !           and changed routine name from XERABT to XERHLT.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERHLT

  USE OMSAO_precision_module
  IMPLICIT NONE

  CHARACTER*(*) MESSG
  !***FIRST EXECUTABLE STATEMENT  XERHLT
  STOP
END SUBROUTINE XERHLT

!*DECK XERMSG
SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
  !***BEGIN PROLOGUE  XERMSG
  !***PURPOSE  Process error messages for SLATEC and other libraries.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XERMSG-A)
  !***KEYWORDS  ERROR MESSAGE, XERROR
  !***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
  !***DESCRIPTION
  !
  !   XERMSG processes a diagnostic message in a manner determined by the
  !   value of LEVEL and the current value of the library error control
  !   flag, KONTRL.  See subroutine XSETF for details.
  !
  !    LIBRAR   A character constant (or character variable) with the name
  !             of the library.  This will be 'SLATEC' for the SLATEC
  !             Common Math Library.  The error handling package is
  !             general enough to be used by many libraries
  !             simultaneously, so it is desirable for the routine that
  !             detects and reports an error to identify the library name
  !             as well as the routine name.
  !
  !    SUBROU   A character constant (or character variable) with the name
  !             of the routine that detected the error.  Usually it is the
  !             name of the routine that is calling XERMSG.  There are
  !             some instances where a user callable library routine calls
  !             lower level subsidiary routines where the error is
  !             detected.  In such cases it may be more informative to
  !             supply the name of the routine the user called rather than
  !             the name of the subsidiary routine that detected the
  !             error.
  !
  !    MESSG    A character constant (or character variable) with the text
  !             of the error or warning message.  In the example below,
  !             the message is a character constant that contains a
  !             generic message.
  !
  !                   CALL XERMSG ('SLATEC', 'MMPY',
  !                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
  !                  *3, 1)
  !
  !             It is possible (and is sometimes desirable) to generate a
  !             specific message--e.g., one that contains actual numeric
  !             values.  Specific numeric values can be converted into
  !             character strings using formatted WRITE statements into
  !             character variables.  This is called standard Fortran
  !             internal file I/O and is exemplified in the first three
  !             lines of the following example.  You can also catenate
  !             substrings of characters to construct the error message.
  !             Here is an example showing the use of both writing to
  !             an internal file and catenating character strings.
  !
  !                   CHARACTER*5 CHARN, CHARL
  !                   WRITE (CHARN,10) N
  !                   WRITE (CHARL,10) LDA
  !                10 FORMAT(I5)
  !                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
  !                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
  !                  *   CHARL, 3, 1)
  !
  !             There are two subtleties worth mentioning.  One is that
  !             the // for character catenation is used to construct the
  !             error message so that no single character constant is
  !             continued to the next line.  This avoids confusion as to
  !             whether there are trailing blanks at the end of the line.
  !             The second is that by catenating the parts of the message
  !             as an actual argument rather than encoding the entire
  !             message into one large character variable, we avoid
  !             having to know how long the message will be in order to
  !             declare an adequate length for that large character
  !             variable.  XERMSG calls XERPRN to print the message using
  !             multiple lines if necessary.  If the message is very long,
  !             XERPRN will break it into pieces of 72 characters (as
  !             requested by XERMSG) for printing on multiple lines.
  !             Also, XERMSG asks XERPRN to prefix each line with ' *  '
  !             so that the total line length could be 76 characters.
  !             Note also that XERPRN scans the error message backwards
  !             to ignore trailing blanks.  Another feature is that
  !             the substring '$$' is treated as a new line sentinel
  !             by XERPRN.  If you want to construct a multiline
  !             message without having to count out multiples of 72
  !             characters, just use '$$' as a separator.  '$$'
  !             obviously must occur within 72 characters of the
  !             start of each line to have its intended effect since
  !             XERPRN is asked to wrap around at 72 characters in
  !             addition to looking for '$$'.
  !
  !    NERR     An integer value that is chosen by the library routine's
  !             author.  It must be in the range -99 to 999 (three
  !             printable digits).  Each distinct error should have its
  !             own error number.  These error numbers should be described
  !             in the machine readable documentation for the routine.
  !             The error numbers need be unique only within each routine,
  !             so it is reasonable for each routine to start enumerating
  !             errors from 1 and proceeding to the next integer.
  !
  !    LEVEL    An integer value in the range 0 to 2 that indicates the
  !             level (severity) of the error.  Their meanings are
  !
  !            -1  A warning message.  This is used if it is not clear
  !                that there really is an error, but the user's attention
  !                may be needed.  An attempt is made to only print this
  !                message once.
  !
  !             0  A warning message.  This is used if it is not clear
  !                that there really is an error, but the user's attention
  !                may be needed.
  !
  !             1  A recoverable error.  This is used even if the error is
  !                so serious that the routine cannot return any useful
  !                answer.  If the user has told the error package to
  !                return after recoverable errors, then XERMSG will
  !                return to the Library routine which can then return to
  !                the user's routine.  The user may also permit the error
  !                package to terminate the program upon encountering a
  !                recoverable error.
  !
  !             2  A fatal error.  XERMSG will not return to its caller
  !                after it receives a fatal error.  This level should
  !                hardly ever be used; it is much better to allow the
  !                user a chance to recover.  An example of one of the few
  !                cases in which it is permissible to declare a level 2
  !                error is a reverse communication Library routine that
  !                is likely to be called repeatedly until it integrates
  !                across some interval.  If there is a serious error in
  !                the input such that another step cannot be taken and
  !                the Library routine is called again without the input
  !                error having been corrected by the caller, the Library
  !                routine will probably be called forever with improper
  !                input.  In this case, it is reasonable to declare the
  !                error to be fatal.
  !
  !    Each of the arguments to XERMSG is input; none will be modified by
  !    XERMSG.  A routine may make multiple calls to XERMSG with warning
  !    level messages; however, after a call to XERMSG with a recoverable
  !    error, the routine should return to the user.  Do not try to call
  !    XERMSG with a second recoverable error after the first recoverable
  !    error because the error package saves the error number.  The user
  !    can retrieve this error number by calling another entry point in
  !    the error handling package and then clear the error number when
  !    recovering from the error.  Calling XERMSG in succession causes the
  !    old error number to be overwritten by the latest error number.
  !    This is considered harmless for error numbers associated with
  !    warning messages but must not be done for error numbers of serious
  !    errors.  After a call to XERMSG with a recoverable error, the user
  !    must be given a chance to call NUMXER or XERCLR to retrieve or
  !    clear the error number.
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
  !***REVISION HISTORY  (YYMMDD)
  !   880101  DATE WRITTEN
  !   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
  !           THERE ARE TWO BASIC CHANGES.
  !           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
  !               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
  !               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
  !               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
  !               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
  !               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
  !               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
  !               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
  !           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
  !               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
  !               OF LOWER CASE.
  !   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
  !           THE PRINCIPAL CHANGES ARE
  !           1.  CLARIFY COMMENTS IN THE PROLOGUES
  !           2.  RENAME XRPRNT TO XERPRN
  !           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
  !               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
  !               CHARACTER FOR NEW RECORDS.
  !   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
  !           CLEAN UP THE CODING.
  !   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
  !           PREFIX.
  !   891013  REVISED TO CORRECT COMMENTS.
  !   891214  Prologue converted to Version 4.0 format.  (WRB)
  !   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
  !           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
  !           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
  !           XERCTL to XERCNT.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERMSG

  USE OMSAO_precision_module
  IMPLICIT NONE

  INTEGER, EXTERNAL :: j4save
  INTEGER :: maxmes, i, ltemp, mkntrl, lkntrl, llevel, lerr, kount, kdummy, level, nerr
  CHARACTER*(*) LIBRAR, SUBROU, MESSG
  CHARACTER*8 XLIBR, XSUBR
  CHARACTER*72  TEMP
  CHARACTER*20  LFIRST
  !***FIRST EXECUTABLE STATEMENT  XERMSG
  LKNTRL = J4SAVE (2, 0, .FALSE.)
  MAXMES = J4SAVE (4, 0, .FALSE.)
  !
  !       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
  !       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
  !          SHOULD BE PRINTED.
  !
  !       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
  !          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
  !          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
  !
  IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.  &
       LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
     CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //  &
          'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//  &
          'JOB ABORT DUE TO FATAL ERROR.', 72)
     CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
     CALL XERHLT (' ***XERMSG -- INVALID INPUT')
     RETURN
  ENDIF
  !
  !       RECORD THE MESSAGE.
  !
  I = J4SAVE (1, NERR, .TRUE.)
  CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
  !
  !       HANDLE PRINT-ONCE WARNING MESSAGES.
  !
  IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
  !
  !       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
  !
  XLIBR  = LIBRAR
  XSUBR  = SUBROU
  LFIRST = MESSG
  LERR   = NERR
  LLEVEL = LEVEL
  CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
  !
  LKNTRL = MAX(-2, MIN(2,LKNTRL))
  MKNTRL = ABS(LKNTRL)
  !
  !       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
  !       ZERO AND THE ERROR IS NOT FATAL.
  !
  IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
  IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
  IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
  IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
  !
  !       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
  !       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
  !       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
  !       IS NOT ZERO.
  !
  IF (LKNTRL .NE. 0) THEN
     TEMP(1:21) = 'MESSAGE FROM ROUTINE '
     I = MIN(LEN(SUBROU), 16)
     TEMP(22:21+I) = SUBROU(1:I)
     TEMP(22+I:33+I) = ' IN LIBRARY '
     LTEMP = 33 + I
     I = MIN(LEN(LIBRAR), 16)
     TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
     TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
     LTEMP = LTEMP + I + 1
     CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
  ENDIF
  !
  !       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
  !       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
  !       FROM EACH OF THE FOLLOWING THREE OPTIONS.
  !       1.  LEVEL OF THE MESSAGE
  !              'INFORMATIVE MESSAGE'
  !              'POTENTIALLY RECOVERABLE ERROR'
  !              'FATAL ERROR'
  !       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
  !              'PROG CONTINUES'
  !              'PROG ABORTED'
  !       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
  !           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
  !           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
  !              'TRACEBACK REQUESTED'
  !              'TRACEBACK NOT REQUESTED'
  !       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
  !       EXCEED 74 CHARACTERS.
  !       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
  !
  IF (LKNTRL .GT. 0) THEN
     !
     !       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
     !
     IF (LEVEL .LE. 0) THEN
        TEMP(1:20) = 'INFORMATIVE MESSAGE,'
        LTEMP = 20
     ELSEIF (LEVEL .EQ. 1) THEN
        TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
        LTEMP = 30
     ELSE
        TEMP(1:12) = 'FATAL ERROR,'
        LTEMP = 12
     ENDIF
     !
     !       THEN WHETHER THE PROGRAM WILL CONTINUE.
     !
     IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.  &
          (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
        TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
        LTEMP = LTEMP + 14
     ELSE
        TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
        LTEMP = LTEMP + 16
     ENDIF
     !
     !       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
     !
     IF (LKNTRL .GT. 0) THEN
        TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
        LTEMP = LTEMP + 20
     ELSE
        TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
        LTEMP = LTEMP + 24
     ENDIF
     CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
  ENDIF
  !
  !       NOW SEND OUT THE MESSAGE.
  !
  CALL XERPRN (' *  ', -1, MESSG, 72)
  !
  !       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
  !          TRACEBACK.
  !
  IF (LKNTRL .GT. 0) THEN
     WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
     DO I=16,22
        IF (TEMP(I:I) .NE. ' ') GO TO 20
     END DO
10   CONTINUE
     !
20   CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
     CALL FDUMP
  ENDIF
  !
  !       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
  !
  IF (LKNTRL .NE. 0) THEN
     CALL XERPRN (' *  ', -1, ' ', 72)
     CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
     CALL XERPRN ('    ',  0, ' ', 72)
  ENDIF
  !
  !       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
  !       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
  !
30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
  !
  !       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
  !       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
  !       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
  !
  IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
     IF (LEVEL .EQ. 1) THEN
        CALL XERPRN  &
             (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
     ELSE
        CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
     ENDIF
     CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
     CALL XERHLT (' ')
  ELSE
     CALL XERHLT (MESSG)
  ENDIF
  RETURN
END SUBROUTINE XERMSG

!*DECK XERPRN
SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
  !***BEGIN PROLOGUE  XERPRN
  !***SUBSIDIARY
  !***PURPOSE  Print error messages processed by XERMSG.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XERPRN-A)
  !***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
  !***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
  !***DESCRIPTION
  !
  ! This routine sends one or more lines to each of the (up to five)
  ! logical units to which error messages are to be sent.  This routine
  ! is called several times by XERMSG, sometimes with a single line to
  ! print and sometimes with a (potentially very long) message that may
  ! wrap around into multiple lines.
  !
  ! PREFIX  Input argument of type CHARACTER.  This argument contains
  !         characters to be put at the beginning of each line before
  !         the body of the message.  No more than 16 characters of
  !         PREFIX will be used.
  !
  ! NPREF   Input argument of type INTEGER.  This argument is the number
  !         of characters to use from PREFIX.  If it is negative, the
  !         intrinsic function LEN is used to determine its length.  If
  !         it is zero, PREFIX is not used.  If it exceeds 16 or if
  !         LEN(PREFIX) exceeds 16, only the first 16 characters will be
  !         used.  If NPREF is positive and the length of PREFIX is less
  !         than NPREF, a copy of PREFIX extended with blanks to length
  !         NPREF will be used.
  !
  ! MESSG   Input argument of type CHARACTER.  This is the text of a
  !         message to be printed.  If it is a long message, it will be
  !         broken into pieces for printing on multiple lines.  Each line
  !         will start with the appropriate prefix and be followed by a
  !         piece of the message.  NWRAP is the number of characters per
  !         piece; that is, after each NWRAP characters, we break and
  !         start a new line.  In addition the characters '$$' embedded
  !         in MESSG are a sentinel for a new line.  The counting of
  !         characters up to NWRAP starts over for each new line.  The
  !         value of NWRAP typically used by XERMSG is 72 since many
  !         older error messages in the SLATEC Library are laid out to
  !         rely on wrap-around every 72 characters.
  !
  ! NWRAP   Input argument of type INTEGER.  This gives the maximum size
  !         piece into which to break MESSG for printing on multiple
  !         lines.  An embedded '$$' ends a line, and the count restarts
  !         at the following character.  If a line break does not occur
  !         on a blank (it would split a word) that word is moved to the
  !         next line.  Values of NWRAP less than 16 will be treated as
  !         16.  Values of NWRAP greater than 132 will be treated as 132.
  !         The actual line length will be NPREF + NWRAP after NPREF has
  !         been adjusted to fall between 0 and 16 and NWRAP has been
  !         adjusted to fall between 16 and 132.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  I1MACH, XGETUA
  !***REVISION HISTORY  (YYMMDD)
  !   880621  DATE WRITTEN
  !   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
  !           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
  !           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
  !           SLASH CHARACTER IN FORMAT STATEMENTS.
  !   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
  !           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
  !           LINES TO BE PRINTED.
  !   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
  !           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
  !   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
  !   891214  Prologue converted to Version 4.0 format.  (WRB)
  !   900510  Added code to break messages between words.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERPRN

  USE OMSAO_precision_module
  IMPLICIT NONE

  INTEGER, EXTERNAL :: i1mach
  CHARACTER*(*) PREFIX, MESSG
  INTEGER NPREF, NWRAP, lwrap, lpref, i, idelta, nextc, lpiece, lenmsg, n
  CHARACTER*148 CBUFF
  INTEGER IU(5), NUNIT
  CHARACTER*2 NEWLIN
  PARAMETER (NEWLIN = '$$')
  !***FIRST EXECUTABLE STATEMENT  XERPRN
  CALL XGETUA(IU,NUNIT)
  !
  !       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
  !       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
  !       ERROR MESSAGE UNIT.
  !
  N = I1MACH(4)
  DO I=1,NUNIT
     IF (IU(I) .EQ. 0) IU(I) = N
  END DO
10 CONTINUE
  !
  !       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
  !       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
  !       THE REST OF THIS ROUTINE.
  !
  IF ( NPREF /= 0 ) THEN
     LPREF = LEN(PREFIX)
  ELSE
     LPREF = NPREF
  ENDIF
  LPREF = MIN(16, LPREF)
  IF (LPREF /= 0) CBUFF(1:LPREF) = PREFIX
  !
  !       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
  !       TIME FROM MESSG TO PRINT ON ONE LINE.
  !
  LWRAP = MAX(16, MIN(132, NWRAP))
  !
  !       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
  !
  LENMSG = LEN(MESSG)
  N = LENMSG
  DO I=1,N
     IF (MESSG(LENMSG:LENMSG) /= ' ') GO TO 30
     LENMSG = LENMSG - 1
  END DO
20 CONTINUE
30 CONTINUE
  !
  !       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
  !
  IF (LENMSG .EQ. 0) THEN
     CBUFF(LPREF+1:LPREF+1) = ' '
     DO I=1,NUNIT
        WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
     END DO
40   CONTINUE
     RETURN
  ENDIF
  !
  !       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
  !       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
  !       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
  !       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
  !
  !       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
  !       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
  !       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
  !       OF THE SECOND ARGUMENT.
  !
  !       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
  !       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
  !       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
  !       POSITION NEXTC.
  !
  !       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
  !                       REMAINDER OF THE CHARACTER STRING.  LPIECE
  !                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
  !                       WHICHEVER IS LESS.
  !
  !       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
  !                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
  !                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
  !                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
  !                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
  !                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
  !                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
  !                       SHOULD BE INCREMENTED BY 2.
  !
  !       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
  !
  !       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
  !                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
  !                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
  !                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
  !                       AT THE END OF A LINE.
  !
  NEXTC = 1
50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
  IF (LPIECE .EQ. 0) THEN
     !
     !       THERE WAS NO NEW LINE SENTINEL FOUND.
     !
     IDELTA = 0
     LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
     IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
        DO I=LPIECE+1,2,-1
           IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
              LPIECE = I-1
              IDELTA = 1
              GOTO 54
           ENDIF
        END DO
52      CONTINUE
     ENDIF
54   CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC = NEXTC + LPIECE + IDELTA
  ELSEIF (LPIECE .EQ. 1) THEN
     !
     !       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
     !       DON'T PRINT A BLANK LINE.
     !
     NEXTC = NEXTC + 2
     GO TO 50
  ELSEIF (LPIECE .GT. LWRAP+1) THEN
     !
     !       LPIECE SHOULD BE SET DOWN TO LWRAP.
     !
     IDELTA = 0
     LPIECE = LWRAP
     DO I=LPIECE+1,2,-1
        IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
           LPIECE = I-1
           IDELTA = 1
           GOTO 58
        ENDIF
     END DO
56   CONTINUE
58   CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC = NEXTC + LPIECE + IDELTA
  ELSE
     !
     !       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
     !       WE SHOULD DECREMENT LPIECE BY ONE.
     !
     LPIECE = LPIECE - 1
     CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC  = NEXTC + LPIECE + 2
  ENDIF
  !
  !       PRINT
  !
  DO I=1,NUNIT
     WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
  END DO
60 CONTINUE
  !
  IF (NEXTC .LE. LENMSG) GO TO 50
  RETURN
END SUBROUTINE XERPRN

!*DECK XERSVE
SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,  &
     ICOUNT)
  !***BEGIN PROLOGUE  XERSVE
  !***SUBSIDIARY
  !***PURPOSE  Record that an error has occurred.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3
  !***TYPE      ALL (XERSVE-A)
  !***KEYWORDS  ERROR, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  ! *Usage:
  !
  !        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
  !        CHARACTER * (len) LIBRAR, SUBROU, MESSG
  !
  !        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
  !
  ! *Arguments:
  !
  !        LIBRAR :IN    is the library that the message is from.
  !        SUBROU :IN    is the subroutine that the message is from.
  !        MESSG  :IN    is the message to be saved.
  !        KFLAG  :IN    indicates the action to be performed.
  !                      when KFLAG > 0, the message in MESSG is saved.
  !                      when KFLAG=0 the tables will be dumped and
  !                      cleared.
  !                      when KFLAG < 0, the tables will be dumped and
  !                      not cleared.
  !        NERR   :IN    is the error number.
  !        LEVEL  :IN    is the error severity.
  !        ICOUNT :OUT   the number of times this message has been seen,
  !                      or zero if the table has overflowed and does not
  !                      contain this message specifically.  When KFLAG=0,
  !                      ICOUNT will not be altered.
  !
  ! *Description:
  !
  !   Record that this error occurred and possibly dump and clear the
  !   tables.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  I1MACH, XGETUA
  !***REVISION HISTORY  (YYMMDD)
  !   800319  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900413  Routine modified to remove reference to KFLAG.  (WRB)
  !   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
  !           sequence, use IF-THEN-ELSE, make number of saved entries
  !           easily changeable, changed routine name from XERSAV to
  !           XERSVE.  (RWC)
  !   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERSVE

  USE OMSAO_precision_module
  IMPLICIT NONE

  INTEGER, EXTERNAL :: i1mach
  INTEGER, PARAMETER :: LENTAB=10
  INTEGER               :: &
       nunit, nmsg, kflag, kountx, lenmsg, iunit, kunit, i, nerr, level, icount
  INTEGER, DIMENSION(5) :: LUN
  CHARACTER*(*) LIBRAR, SUBROU, MESSG
  CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
  CHARACTER*20 MESTAB(LENTAB), MES

  INTEGER, DIMENSION(lentab) :: NERTAB, LEVTAB, KOUNT

  SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
  DATA KOUNTX/0/, NMSG/0/
  !***FIRST EXECUTABLE STATEMENT  XERSVE
  !
  IF (KFLAG.LE.0) THEN
     !
     !        Dump the table.
     !
     IF (NMSG.EQ.0) RETURN
     !
     !        Print to each unit.
     !
     CALL XGETUA (LUN, NUNIT)
     DO KUNIT = 1,NUNIT
        IUNIT = LUN(KUNIT)
        IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
        !
        !           Print the table header.
        !
        WRITE (IUNIT,9000)
        !
        !           Print body of table.
        !
        DO I = 1,NMSG
           WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),  &
                NERTAB(I),LEVTAB(I),KOUNT(I)
        END DO
10      CONTINUE
        !
        !           Print number of other errors.
        !
        IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
        WRITE (IUNIT,9030)
     END DO
20   CONTINUE
     !
     !        Clear the error tables.
     !
     IF (KFLAG.EQ.0) THEN
        NMSG = 0
        KOUNTX = 0
     ENDIF
  ELSE
     !
     !        PROCESS A MESSAGE...
     !        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
     !        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
     !
     LIB = LIBRAR
     SUB = SUBROU
     MES = MESSG
     DO I = 1,NMSG
        IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.  &
             MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.  &
             LEVEL.EQ.LEVTAB(I)) THEN
           KOUNT(I) = KOUNT(I) + 1
           ICOUNT = KOUNT(I)
           RETURN
        ENDIF
     END DO
30   CONTINUE
     !
     IF (NMSG.LT.LENTAB) THEN
        !
        !           Empty slot found for new message.
        !
        NMSG = NMSG + 1
        LIBTAB(I) = LIB
        SUBTAB(I) = SUB
        MESTAB(I) = MES
        NERTAB(I) = NERR
        LEVTAB(I) = LEVEL
        KOUNT (I) = 1
        ICOUNT    = 1
     ELSE
        !
        !           Table is full.
        !
        KOUNTX = KOUNTX+1
        ICOUNT = 0
     ENDIF
  ENDIF
  RETURN
  !
  !     Formats.
  !
9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /  &
          ' LIBRARY    SUBROUTINE MESSAGE START             NERR',  &
          '     LEVEL     COUNT')
9010 FORMAT (1X,A,3X,A,3X,A,3I10)
9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
9030 FORMAT (1X)
END SUBROUTINE XERSVE

!*DECK XGETUA
SUBROUTINE XGETUA (IUNITA, N)
  !***BEGIN PROLOGUE  XGETUA
  !***PURPOSE  Return unit number(s) to which error messages are being
  !            sent.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XGETUA-A)
  !***KEYWORDS  ERROR, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract
  !        XGETUA may be called to determine the unit number or numbers
  !        to which error messages are being sent.
  !        These unit numbers may have been set by a call to XSETUN,
  !        or a call to XSETUA, or may be a default value.
  !
  !     Description of Parameters
  !      --Output--
  !        IUNIT - an array of one to five unit numbers, depending
  !                on the value of N.  A value of zero refers to the
  !                default unit, as defined by the I1MACH machine
  !                constant routine.  Only IUNIT(1),...,IUNIT(N) are
  !                defined by XGETUA.  The values of IUNIT(N+1),...,
  !                IUNIT(5) are not defined (for N .LT. 5) or altered
  !                in any way by XGETUA.
  !        N     - the number of units to which copies of the
  !                error messages are being sent.  N will be in the
  !                range from 1 to 5.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  J4SAVE
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XGETUA

  USE OMSAO_precision_module
  IMPLICIT NONE

  INTEGER, EXTERNAL :: j4save
  INTEGER :: n, i, index
  INTEGER, DIMENSION (5) :: IUNITA
  !***FIRST EXECUTABLE STATEMENT  XGETUA
  N = J4SAVE(5,0,.FALSE.)
  DO I=1,N
     INDEX = I+4
     IF (I.EQ.1) INDEX = 3
     IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
  END DO
30 CONTINUE
  RETURN
END SUBROUTINE XGETUA

!*DECK DPOLFT
SUBROUTINE DPOLFT (N, X, Y, W, MAXDEG, NDEG, EPS, R, IERR, A)
  !***BEGIN PROLOGUE  DPOLFT
  !***PURPOSE  Fit discrete data in a least squares sense by polynomials
  !            in one variable.
  !***LIBRARY   SLATEC
  !***CATEGORY  K1A1A2
  !***TYPE      DOUBLE PRECISION (POLFIT-S, DPOLFT-D)
  !***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
  !***AUTHOR  Shampine, L. F., (SNLA)
  !           Davenport, S. M., (SNLA)
  !           Huddleston, R. E., (SNLL)
  !***DESCRIPTION
  !
  !     Abstract
  !
  !     Given a collection of points X(I) and a set of values Y(I) which
  !     correspond to some function or measurement at each of the X(I),
  !     subroutine  DPOLFT  computes the weighted least-squares polynomial
  !     fits of all degrees up to some degree either specified by the user
  !     or determined by the routine.  The fits thus obtained are in
  !     orthogonal polynomial form.  Subroutine  DP1VLU  may then be
  !     called to evaluate the fitted polynomials and any of their
  !     derivatives at any point.  The subroutine  DPCOEF  may be used to
  !     express the polynomial fits as powers of (X-C) for any specified
  !     point C.
  !
  !     The parameters for  DPOLFT  are
  !
  !     Input -- All TYPE REAL variables are DOUBLE PRECISION
  !         N -      the number of data points.  The arrays X, Y and W
  !                  must be dimensioned at least  N  (N .GE. 1).
  !         X -      array of values of the independent variable.  These
  !                  values may appear in any order and need not all be
  !                  distinct.
  !         Y -      array of corresponding function values.
  !         W -      array of positive values to be used as weights.  If
  !                  W(1) is negative,  DPOLFT  will set all the weights
  !                  to 1.0, which means unweighted least squares error
  !                  will be minimized.  To minimize relative error, the
  !                  user should set the weights to:  W(I) = 1.0/Y(I)**2,
  !                  I = 1,...,N .
  !         MAXDEG - maximum degree to be allowed for polynomial fit.
  !                  MAXDEG  may be any non-negative integer less than  N.
  !                  Note -- MAXDEG  cannot be equal to  N-1  when a
  !                  statistical test is to be used for degree selection,
  !                  i.e., when input value of  EPS  is negative.
  !         EPS -    specifies the criterion to be used in determining
  !                  the degree of fit to be computed.
  !                  (1)  If  EPS  is input negative,  DPOLFT  chooses the
  !                       degree based on a statistical F test of
  !                       significance.  One of three possible
  !                       significance levels will be used:  .01, .05 or
  !                       .10.  If  EPS=-1.0 , the routine will
  !                       automatically select one of these levels based
  !                       on the number of data points and the maximum
  !                       degree to be considered.  If  EPS  is input as
  !                       -.01, -.05, or -.10, a significance level of
  !                       .01, .05, or .10, respectively, will be used.
  !                  (2)  If  EPS  is set to 0.,  DPOLFT  computes the
  !                       polynomials of degrees 0 through  MAXDEG .
  !                  (3)  If  EPS  is input positive,  EPS  is the RMS
  !                       error tolerance which must be satisfied by the
  !                       fitted polynomial.  DPOLFT  will increase the
  !                       degree of fit until this criterion is met or
  !                       until the maximum degree is reached.
  !
  !     Output -- All TYPE REAL variables are DOUBLE PRECISION
  !         NDEG -   degree of the highest degree fit computed.
  !         EPS -    RMS error of the polynomial of degree  NDEG .
  !         R -      vector of dimension at least NDEG containing values
  !                  of the fit of degree  NDEG  at each of the  X(I) .
  !                  Except when the statistical test is used, these
  !                  values are more accurate than results from subroutine
  !                  DP1VLU  normally are.
  !         IERR -   error flag with the following possible values.
  !             1 -- indicates normal execution, i.e., either
  !                  (1)  the input value of  EPS  was negative, and the
  !                       computed polynomial fit of degree  NDEG
  !                       satisfies the specified F test, or
  !                  (2)  the input value of  EPS  was 0., and the fits of
  !                       all degrees up to  MAXDEG  are complete, or
  !                  (3)  the input value of  EPS  was positive, and the
  !                       polynomial of degree  NDEG  satisfies the RMS
  !                       error requirement.
  !             2 -- invalid input parameter.  At least one of the input
  !                  parameters has an illegal value and must be corrected
  !                  before  DPOLFT  can proceed.  Valid input results
  !                  when the following restrictions are observed
  !                       N .GE. 1
  !                       0 .LE. MAXDEG .LE. N-1  for  EPS .GE. 0.
  !                       0 .LE. MAXDEG .LE. N-2  for  EPS .LT. 0.
  !                       W(1)=-1.0  or  W(I) .GT. 0., I=1,...,N .
  !             3 -- cannot satisfy the RMS error requirement with a
  !                  polynomial of degree no greater than  MAXDEG .  Best
  !                  fit found is of degree  MAXDEG .
  !             4 -- cannot satisfy the test for significance using
  !                  current value of  MAXDEG .  Statistically, the
  !                  best fit found is of order  NORD .  (In this case,
  !                  NDEG will have one of the values:  MAXDEG-2,
  !                  MAXDEG-1, or MAXDEG).  Using a higher value of
  !                  MAXDEG  may result in passing the test.
  !         A -      work and output array having at least 3N+3MAXDEG+3
  !                  locations
  !
  !     Note - DPOLFT  calculates all fits of degrees up to and including
  !            NDEG .  Any or all of these fits can be evaluated or
  !            expressed as powers of (X-C) using  DP1VLU  and  DPCOEF
  !            after just one call to  DPOLFT .
  !
  !***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
  !                 Curve fitting by polynomials in one variable, Report
  !                 SLA-74-0270, Sandia Laboratories, June 1974.
  !***ROUTINES CALLED  DP1VLU, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   740601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900911  Added variable YP to DOUBLE PRECISION declaration.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !   920527  Corrected erroneous statements in DESCRIPTION.  (WRB)
  !***END PROLOGUE  DPOLFT
  INTEGER N, I,IDEGF,IERR,J,JP1,JPAS,K1,K1PJ,K2,K2PJ,K3,K3PI,K4,  &
       K4PI,K5,K5PI,KSIG,M,MAXDEG,MOP1,NDEG,NDER,NFAIL
  DOUBLE PRECISION TEMD1,TEMD2
  DOUBLE PRECISION A(*),DEGF,DEN,EPS,ETST,F,FCRIT,R(*),SIG,SIGJ,  &
       SIGJM1,SIGPAS,TEMP,X(*),XM,Y(*),YP,W(*),W1,W11
  DOUBLE PRECISION CO(4,3)
  SAVE CO
  DATA  CO(1,1), CO(2,1), CO(3,1), CO(4,1), CO(1,2), CO(2,2),    & 
       CO(3,2), CO(4,2), CO(1,3), CO(2,3), CO(3,3),             &
       CO(4,3)/-13.086850D0,-2.4648165D0,-3.3846535D0,-1.2973162D0, &
       -3.3381146D0,-1.7812271D0,-3.2578406D0,-1.6589279D0, &
       -1.6282703D0,-1.3152745D0,-3.2640179D0,-1.9829776D0/
  !***FIRST EXECUTABLE STATEMENT  DPOLFT
  M = ABS(N)
  IF (M .EQ. 0) GO TO 30
  IF (MAXDEG .LT. 0) GO TO 30
  A(1) = MAXDEG
  MOP1 = MAXDEG + 1
  IF (M .LT. MOP1) GO TO 30
  IF (EPS .LT. 0.0D0 .AND.  M .EQ. MOP1) GO TO 30
  XM = M
  ETST = EPS*EPS*XM
  IF (W(1) .LT. 0.0D0) GO TO 2
  DO 1 I = 1,M
     IF (W(I) .LE. 0.0D0) GO TO 30
1 END DO
  GO TO 4
2 DO 3 I = 1,M
     W(I) = 1.0D0
3 END DO
4 IF (EPS .GE. 0.0D0) GO TO 8
  !
  ! DETERMINE SIGNIFICANCE LEVEL INDEX TO BE USED IN STATISTICAL TEST FOR
  ! CHOOSING DEGREE OF POLYNOMIAL FIT
  !
  IF (EPS .GT. (-.55D0)) GO TO 5
  IDEGF = M - MAXDEG - 1
  KSIG = 1
  IF (IDEGF .LT. 10) KSIG = 2
  IF (IDEGF .LT. 5) KSIG = 3
  GO TO 8
5 KSIG = 1
  IF (EPS .LT. (-.03D0)) KSIG = 2
  IF (EPS .LT. (-.07D0)) KSIG = 3
  !
  ! INITIALIZE INDEXES AND COEFFICIENTS FOR FITTING
  !
8 K1 = MAXDEG + 1
  K2 = K1 + MAXDEG
  K3 = K2 + MAXDEG + 2
  K4 = K3 + M
  K5 = K4 + M
  DO 9 I = 2,K4
     A(I) = 0.0D0
9 END DO
  W11 = 0.0D0
  IF (N .LT. 0) GO TO 11
  !
  ! UNCONSTRAINED CASE
  !
  DO 10 I = 1,M
     K4PI = K4 + I
     A(K4PI) = 1.0D0
     W11 = W11 + W(I)
10 END DO
  GO TO 13
  !
  ! CONSTRAINED CASE
  !
11 DO 12 I = 1,M
     K4PI = K4 + I
     W11 = W11 + W(I)*A(K4PI)**2
12 END DO
  !
  ! COMPUTE FIT OF DEGREE ZERO
  !
13 TEMD1 = 0.0D0
  DO 14 I = 1,M
     K4PI = K4 + I
     TEMD1 = TEMD1 + W(I)*Y(I)*A(K4PI)
14 END DO
  TEMD1 = TEMD1/W11
  A(K2+1) = TEMD1
  SIGJ = 0.0D0
  DO 15 I = 1,M
     K4PI = K4 + I
     K5PI = K5 + I
     TEMD2 = TEMD1*A(K4PI)
     R(I) = TEMD2
     A(K5PI) = TEMD2 - R(I)
     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
15 END DO
  J = 0
  !
  ! SEE IF POLYNOMIAL OF DEGREE 0 SATISFIES THE DEGREE SELECTION CRITERION
  !
  IF (EPS) 24,26,27
  !
  ! INCREMENT DEGREE
  !
16 J = J + 1
  JP1 = J + 1
  K1PJ = K1 + J
  K2PJ = K2 + J
  SIGJM1 = SIGJ
  !
  ! COMPUTE NEW B COEFFICIENT EXCEPT WHEN J = 1
  !
  IF (J .GT. 1) A(K1PJ) = W11/W1
  !
  ! COMPUTE NEW A COEFFICIENT
  !
  TEMD1 = 0.0D0
  DO 18 I = 1,M
     K4PI = K4 + I
     TEMD2 = A(K4PI)
     TEMD1 = TEMD1 + X(I)*W(I)*TEMD2*TEMD2
18 END DO
  A(JP1) = TEMD1/W11
  !
  ! EVALUATE ORTHOGONAL POLYNOMIAL AT DATA POINTS
  !
  W1 = W11
  W11 = 0.0D0
  DO 19 I = 1,M
     K3PI = K3 + I
     K4PI = K4 + I
     TEMP = A(K3PI)
     A(K3PI) = A(K4PI)
     A(K4PI) = (X(I)-A(JP1))*A(K3PI) - A(K1PJ)*TEMP
     W11 = W11 + W(I)*A(K4PI)**2
19 END DO
  !
  ! GET NEW ORTHOGONAL POLYNOMIAL COEFFICIENT USING PARTIAL DOUBLE
  ! PRECISION
  !
  TEMD1 = 0.0D0
  DO 20 I = 1,M
     K4PI = K4 + I
     K5PI = K5 + I
     TEMD2 = W(I)*((Y(I)-R(I))-A(K5PI))*A(K4PI)
     TEMD1 = TEMD1 + TEMD2
20 END DO
  TEMD1 = TEMD1/W11
  A(K2PJ+1) = TEMD1
  !
  ! UPDATE POLYNOMIAL EVALUATIONS AT EACH OF THE DATA POINTS, AND
  ! ACCUMULATE SUM OF SQUARES OF ERRORS.  THE POLYNOMIAL EVALUATIONS ARE
  ! COMPUTED AND STORED IN EXTENDED PRECISION.  FOR THE I-TH DATA POINT,
  ! THE MOST SIGNIFICANT BITS ARE STORED IN  R(I) , AND THE LEAST
  ! SIGNIFICANT BITS ARE IN  A(K5PI) .
  !
  SIGJ = 0.0D0
  DO I = 1,M !21
     K4PI = K4 + I
     K5PI = K5 + I
     TEMD2 = R(I) + A(K5PI) + TEMD1*A(K4PI)
     R(I) = TEMD2
     A(K5PI) = TEMD2 - R(I)
     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
  END DO !21 CONTINUE
  !
  ! SEE IF DEGREE SELECTION CRITERION HAS BEEN SATISFIED OR IF DEGREE
  ! MAXDEG  HAS BEEN REACHED
  !
  IF (EPS) 23,26,27
  !
  ! COMPUTE F STATISTICS  (INPUT EPS .LT. 0.)
  !
23 IF (SIGJ .EQ. 0.0D0) GO TO 29
  DEGF = M - J - 1
  DEN = (CO(4,KSIG)*DEGF + 1.0D0)*DEGF
  FCRIT = (((CO(3,KSIG)*DEGF) + CO(2,KSIG))*DEGF + CO(1,KSIG))/DEN
  FCRIT = FCRIT*FCRIT
  F = (SIGJM1 - SIGJ)*DEGF/SIGJ
  IF (F .LT. FCRIT) GO TO 25
  !
  ! POLYNOMIAL OF DEGREE J SATISFIES F TEST
  !
24 SIGPAS = SIGJ
  JPAS = J
  NFAIL = 0
  IF (MAXDEG .EQ. J) GO TO 32
  GO TO 16
  !
  ! POLYNOMIAL OF DEGREE J FAILS F TEST.  IF THERE HAVE BEEN THREE
  ! SUCCESSIVE FAILURES, A STATISTICALLY BEST DEGREE HAS BEEN FOUND.
  !
25 NFAIL = NFAIL + 1
  IF (NFAIL .GE. 3) GO TO 29
  IF (MAXDEG .EQ. J) GO TO 32
  GO TO 16
  !
  ! RAISE THE DEGREE IF DEGREE  MAXDEG  HAS NOT YET BEEN REACHED  (INPUT
  ! EPS = 0.)
  !
26 IF (MAXDEG .EQ. J) GO TO 28
  GO TO 16
  !
  ! SEE IF RMS ERROR CRITERION IS SATISFIED  (INPUT EPS .GT. 0.)
  !
27 IF (SIGJ .LE. ETST) GO TO 28
  IF (MAXDEG .EQ. J) GO TO 31
  GO TO 16
  !
  ! RETURNS
  !
28 IERR = 1
  NDEG = J
  SIG = SIGJ
  GO TO 33
29 IERR = 1
  NDEG = JPAS
  SIG = SIGPAS
  GO TO 33
30 IERR = 2
  CALL XERMSG ('SLATEC', 'DPOLFT', 'INVALID INPUT PARAMETER.', 2,  &
       1)
  GO TO 37
31 IERR = 3
  NDEG = MAXDEG
  SIG = SIGJ
  GO TO 33
32 IERR = 4
  NDEG = JPAS
  SIG = SIGPAS
  !
33 A(K3) = NDEG
  !
  ! WHEN STATISTICAL TEST HAS BEEN USED, EVALUATE THE BEST POLYNOMIAL AT
  ! ALL THE DATA POINTS IF  R  DOES NOT ALREADY CONTAIN THESE VALUES
  !
  IF(EPS .GE. 0.0  .OR.  NDEG .EQ. MAXDEG) GO TO 36
  NDER = 0
  DO I = 1,M
     CALL DP1VLU (NDEG,NDER,X(I),R(I),(/YP/),A)
  END DO
  ! 35     CONTINUE
36 EPS = SQRT(SIG/XM)
37 RETURN
END SUBROUTINE DPOLFT

!*DECK DP1VLU
SUBROUTINE DP1VLU (L, NDER, X, YFIT, YP, A)
  !***BEGIN PROLOGUE  DP1VLU
  !***PURPOSE  Use the coefficients generated by DPOLFT to evaluate the
  !            polynomial fit of degree L, along with the first NDER of
  !            its derivatives, at a specified point.
  !***LIBRARY   SLATEC
  !***CATEGORY  K6
  !***TYPE      DOUBLE PRECISION (PVALUE-S, DP1VLU-D)
  !***KEYWORDS  CURVE FITTING, LEAST SQUARES, POLYNOMIAL APPROXIMATION
  !***AUTHOR  Shampine, L. F., (SNLA)
  !           Davenport, S. M., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract
  !
  !     The subroutine  DP1VLU  uses the coefficients generated by  DPOLFT
  !     to evaluate the polynomial fit of degree  L , along with the first
  !     NDER  of its derivatives, at a specified point.  Computationally
  !     stable recurrence relations are used to perform this task.
  !
  !     The parameters for  DP1VLU  are
  !
  !     Input -- ALL TYPE REAL variables are DOUBLE PRECISION
  !         L -      the degree of polynomial to be evaluated.  L  may be
  !                  any non-negative integer which is less than or equal
  !                  to  NDEG , the highest degree polynomial provided
  !                  by  DPOLFT .
  !         NDER -   the number of derivatives to be evaluated.  NDER
  !                  may be 0 or any positive value.  If NDER is less
  !                  than 0, it will be treated as 0.
  !         X -      the argument at which the polynomial and its
  !                  derivatives are to be evaluated.
  !         A -      work and output array containing values from last
  !                  call to  DPOLFT .
  !
  !     Output -- ALL TYPE REAL variables are DOUBLE PRECISION
  !         YFIT -   value of the fitting polynomial of degree  L  at  X
  !         YP -     array containing the first through  NDER  derivatives
  !                  of the polynomial of degree  L .  YP  must be
  !                  dimensioned at least  NDER  in the calling program.
  !
  !***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
  !                 Curve fitting by polynomials in one variable, Report
  !                 SLA-74-0270, Sandia Laboratories, June 1974.
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   740601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DP1VLU
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  INTEGER I,IC,ILO,IN,INP1,IUP,K1,K1I,K2,K3,K3P1,K3PN,K4,K4P1,K4PN,  &
       KC,L,LM1,LP1,MAXORD,N,NDER,NDO,NDP1,NORD
  DOUBLE PRECISION A(*),CC,DIF,VAL,X,YFIT,YP(*)
  CHARACTER*8 XERN1, XERN2
  !***FIRST EXECUTABLE STATEMENT  DP1VLU
  IF (L .LT. 0) GO TO 12
  NDO = MAX(NDER,0)
  NDO = MIN(NDO,L)
  MAXORD = A(1) + 0.5D0
  K1 = MAXORD + 1
  K2 = K1 + MAXORD
  K3 = K2 + MAXORD + 2
  NORD = A(K3) + 0.5D0
  IF (L .GT. NORD) GO TO 11
  K4 = K3 + L + 1
  IF (NDER .LT. 1) GO TO 2
  DO I = 1,NDER
     YP(I) = 0.0D0
  END DO
2 IF (L .GE. 2) GO TO 4
 IF (L .EQ. 1) GO TO 3
 !
 ! L IS 0
 !
 VAL = A(K2+1)
 GO TO 10
 !
 ! L IS 1
 !
3 CC = A(K2+2)
 VAL = A(K2+1) + (X-A(2))*CC
 IF (NDER .GE. 1) YP(1) = CC
 GO TO 10
 !
 ! L IS GREATER THAN 1
 !
4 NDP1 = NDO + 1
 K3P1 = K3 + 1
 K4P1 = K4 + 1
 LP1 = L + 1
 LM1 = L - 1
 ILO = K3 + 3
 IUP = K4 + NDP1
 DO I = ILO,IUP
    A(I) = 0.0D0
 END DO
 DIF = X - A(LP1)
 KC = K2 + LP1
 A(K4P1) = A(KC)
 A(K3P1) = A(KC-1) + DIF*A(K4P1)
 A(K3+2) = A(K4P1)
 !
 ! EVALUATE RECURRENCE RELATIONS FOR FUNCTION VALUE AND DERIVATIVES
 !
 DO I = 1,LM1
    IN = L - I
    INP1 = IN + 1
    K1I = K1 + INP1
    IC = K2 + IN
    DIF = X - A(INP1)
    VAL = A(IC) + DIF*A(K3P1) - A(K1I)*A(K4P1)
    IF (NDO .LE. 0) GO TO 8
    DO N = 1,NDO
       K3PN = K3P1 + N
       K4PN = K4P1 + N
       YP(N) = DIF*A(K3PN) + N*A(K3PN-1) - A(K1I)*A(K4PN)
    END DO
    !
    ! SAVE VALUES NEEDED FOR NEXT EVALUATION OF RECURRENCE RELATIONS
    !
    DO N = 1,NDO
       K3PN = K3P1 + N
       K4PN = K4P1 + N
       A(K4PN) = A(K3PN)
       A(K3PN) = YP(N)
    END DO
8   A(K4P1) = A(K3P1)
    A(K3P1) = VAL
 END DO
 !
 ! NORMAL RETURN OR ABORT DUE TO ERROR
 !
10 YFIT = VAL
 RETURN
 !
11 WRITE (XERN1, '(I8)') L
 WRITE (XERN2, '(I8)') NORD
 CALL XERMSG ('SLATEC', 'DP1VLU',  &
      'THE ORDER OF POLYNOMIAL EVALUATION, L = ' // XERN1 //  &
      ' REQUESTED EXCEEDS THE HIGHEST ORDER FIT, NORD = ' // XERN2 //  &
      ', COMPUTED BY DPOLFT -- EXECUTION TERMINATED.', 8, 2)
 RETURN
 !
12 CALL XERMSG ('SLATEC', 'DP1VLU',  &
        'INVALID INPUT PARAMETER.  ORDER OF POLYNOMIAL EVALUATION ' //  &
        'REQUESTED IS NEGATIVE.', 2, 2)
 RETURN
END SUBROUTINE DP1VLU
