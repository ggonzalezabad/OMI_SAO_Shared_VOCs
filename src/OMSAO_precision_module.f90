MODULE OMSAO_precision_module

  USE ISO_C_BINDING, ONLY: C_LONG
  IMPLICIT NONE

  ! =====================================================
  ! Define KIND variables for single and double precision
  ! =====================================================

  INTEGER, PARAMETER :: i1 = SELECTED_INT_KIND(2**1)
  INTEGER, PARAMETER :: i2 = SELECTED_INT_KIND(2**2)
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(2**3)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(2**4)

  ! ---------------------------------------------------------------------
  ! HDF-EOS5 requires long integers (i.e., KIND=8) for some of the
  ! variables in the Fortran wrappers for the C interface when compiling
  ! with 64bit compilers. 32bit compilers require KIND=4 for those
  ! variables and will fail if they are declared as long integers. The
  ! easiest way to satisfy those orthogonal requirements is to defive a
  ! precision variable "i48" that is either 4 or 8 dpending on the 
  ! compiler that is being used, and to place that in a preprocessor
  ! clause. We then only have to make sure that
  ! (a) All relevant variables are declared as KIND=i48, which are mainly
  !     dimensions in calls to HDF-EOS5 routines; and
  ! (b) That the compiler directives include "-Dtrue64"
  ! ---------------------------------------------------------------------
!#ifdef true64
!  INTEGER, PARAMETER :: i48 = SELECTED_INT_KIND(2**4)
!#else
!  INTEGER, PARAMETER :: i48 = i4
!#endif

  INTEGER, PARAMETER :: r4 = KIND(1.0)
  INTEGER, PARAMETER :: r8 = KIND(1.0D0)

END MODULE OMSAO_precision_module

