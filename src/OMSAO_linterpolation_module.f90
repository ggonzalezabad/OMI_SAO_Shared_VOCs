module OMSAO_lininterpolation_module

  !F90---------------------------------------------------------------------------!
  !
  !Description:
  ! This module consists of linear interpolation functions for 1 to 7 dimensions.
  ! Lininterpol is the only interface that can be called from ouside the module.
  !
  !Input:
  ! linInterpol:
  ! dim1, .., dimn            : number of nodes in dimension i
  ! dimscale1, .., dimscalen  : dimension scale of the nodes
  ! datand                    : data to be interpolated
  ! x1, .., xn                : value of each dimension to interpolate to
  ! node1_in, ..,noden_in     ! nodes for all dimensions (optional input)
  ! status                    : return status: 0 if succes; -1 if fail
  !
  !Output:
  !
  !Input/Output:
  ! linInterpol:
  ! returns the interpolated value.
  !
  !Revision History:
  ! nodes are optional input,   Roeland van Oss, 26/02/2003
  ! Changed real into REAL(SP) and introduced USE constants, ONLY : SP
  ! Roeland van Oss,,17/05/2002
  ! Version 1.0 07/01/2002
  ! Developed by J.P. Veefkind, KNMI (veefkind@knmi.nl)
  !
  !OMI-KNMI header:
  !
  ! CfA Revision history:
  ! Modified to suit OMSAO OMI applications (gga June 2013)
  !  Removed USE constants, ONLE : SP
  !END---------------------------------------------------------------------------!

  implicit none

  INTEGER, PARAMETER :: SP = 8

  interface linInterpol
    module procedure interpol1D
    module procedure interpol2D
    module procedure interpol3D
    module procedure interpol4D
    module procedure interpol5D
    module procedure interpol6D
    module procedure interpol7D
  end interface

  PUBLIC :: GetNode
  private :: monotonDecrease, findNode
  private :: interpol1D, interpol2D, interpol3D, interpol4D, interpol5D, interpol6D, interpol7D

contains

  !==============================================================================!
  logical function monotonIncrease(dimScale,ndim)
    !
    ! checks if array is monotonically increasing
    !
    integer, intent(in) :: ndim
    REAL(SP), dimension(ndim), intent(in) :: dimscale

    if ( minval( dimScale(2:ndim) - dimScale(1:ndim-1)) <= 0. ) then
      monotonIncrease = .FALSE.
    else
      monotonIncrease = .TRUE.
    end if

  end function  monotonIncrease

  !==============================================================================!
  logical function monotonDecrease(dimScale,ndim)
    !
    ! checks if array is monotonically increasing
    !
    integer, intent(in) :: ndim
    REAL(SP), dimension(ndim), intent(in) :: dimscale

    if ( minval( dimScale(2:ndim) - dimScale(1:ndim-1)) >= 0. ) then
      monotonDecrease = .FALSE.
    else
      monotonDecrease = .TRUE.
    end if

  end function  monotonDecrease

  integer function findNode(dimScale, ndim, x)

    integer, intent(in) :: ndim
    REAL(SP), dimension(ndim), intent(in) :: dimscale
    REAL(SP), intent(in) :: x
    integer, dimension(1) :: nearestLocation

    ! check if the dimension scale is monotonically increasing or decreasing
    if (.NOT.monotonIncrease(dimScale,ndim) .and. &
      .NOT.monotonDecrease(dimScale,ndim)) then
      findnode = -2
      return
    end if

    ! deal with two exceptions:
    ! (1) x = dimScale(1)
    ! (2) x = dimScale(ndim)
    if ( x == dimScale(1) ) then
      findNode = 1
      return
    end if

    if ( x == dimScale(ndim) ) then
      findNode = ndim -1
      return
    end if

    ! find the nearest grid point near x
    nearestLocation=minloc(abs(dimscale-x))

    if (dimscale(nearestLocation(1)) > x ) then
      findNode = nearestLocation(1) - 1
    else
      findNode = nearestLocation(1)
    end if

    if (findNode == 0) findnode = -1  ! extrapolation

    if (findNode == ndim) findnode = -1 ! extrapolation

    return

  end function findNode

  !==============================================================================!
  REAL(SP) function interpol1D(dim1, dimScale1, data1D, x1, node1_in, status)

    integer, parameter :: ndim=1

    integer, intent(in) :: dim1
    REAL(SP), dimension(dim1), intent(in) :: dimScale1
    REAL(SP), dimension(dim1), intent(in) :: data1D
    REAL(SP), intent(in) :: x1
    integer, intent(in), OPTIONAL :: node1_in
    integer, intent(out), OPTIONAL :: status

    integer :: node1
    REAL(SP), dimension(ndim) :: distance
    REAL(SP), dimension(2**ndim) :: weights, values
    integer :: i
    integer :: npoint
    status = 0 ! set to success

    !adapted by RF van Oss:
    IF ( PRESENT(node1_in) ) THEN
      !nodes are given on input
      node1 = node1_in
    ELSE
      node1 = findNode(dimScale1, dim1, x1)
    ENDIF
    if ( node1 <= 0 ) then ! error getting nodes, exit with error
      status = -1
      interpol1D = -huge(interpol1D)
      return
    end if

    npoint = 0
    do i = node1, node1 + 1
      npoint = npoint + 1

      distance(1) = (x1 - dimscale1(i)) / (dimscale1(node1+1) - dimscale1(node1))

      distance = 1. - abs(distance)

      weights(npoint) = product(distance)

      values(npoint) = data1D(i)

    end do

    ! normalize weights
    weights(1:npoint) =  weights(1:npoint) / sum(weights(1:npoint))

    ! determine output
    interpol1D = 0.
    do i=1, npoint
      interpol1D = interpol1D + weights(i) * values(i)
    end do

    return

  end function interpol1D

  !==============================================================================!
  REAL(SP) function interpol2D(dim1, dim2 , &
      dimScale1, dimScale2, &
      data2D, x1, x2, node1_in, node2_in, status)

    integer, parameter :: ndim=2

    integer, intent(in) :: dim1, dim2
    REAL(SP), dimension(dim1), intent(in) :: dimScale1
    REAL(SP), dimension(dim2), intent(in) :: dimScale2
    REAL(SP), dimension(dim1,dim2), intent(in) :: data2D
    REAL(SP), intent(in) :: x1,x2
    integer, intent(in), OPTIONAL :: node1_in, node2_in
    integer, intent(out), OPTIONAL :: status

    integer :: node1, node2
    REAL(SP), dimension(ndim) :: distance
    REAL(SP), dimension(2**ndim) :: weights, values
    integer :: i, j
    integer :: npoint

    status = 0 ! set to success

    !adapted by RF van Oss:
    IF ( PRESENT(node1_in) ) THEN
      !nodes are given on input
      node1 = node1_in
      node2 = node2_in
    ELSE
      node1 = findNode(dimScale1, dim1, x1)
      node2 = findNode(dimScale2, dim2, x2)
    ENDIF

    if ( min(node1,node2) <= 0 ) then ! error getting nodes, exit with error
      status = -1
      interpol2D = -huge(interpol2D)
      return
    end if

    npoint = 0
    do i = node1, node1 + 1
      do j=node2, node2 + 1

        npoint = npoint + 1

        distance(1) = (x1 - dimscale1(i)) / (dimscale1(node1+1) - dimscale1(node1))

        distance(2) = (x2 - dimscale2(j)) / (dimscale2(node2+1) - dimscale2(node2))

        distance = 1. - abs(distance)

        weights(npoint) = product(distance)

        values(npoint) = data2D(i,j)

      end do
    end do

    ! normalize weights
    weights(1:npoint) =  weights(1:npoint) / sum(weights(1:npoint))

    ! determine output
    interpol2D = 0.
    do i=1, npoint
      interpol2D = interpol2D + weights(i) * values(i)
    end do

    return

  end function interpol2D

  !==============================================================================!
  REAL(SP) function interpol3D(dim1, dim2, dim3 , &
      dimScale1, dimScale2, dimScale3, &
      data3D, x1, x2, x3, node1_in, node2_in, node3_in, status)

    integer, parameter :: ndim=3

    integer, intent(in) :: dim1, dim2, dim3
    REAL(SP), dimension(dim1), intent(in) :: dimScale1
    REAL(SP), dimension(dim2), intent(in) :: dimScale2
    REAL(SP), dimension(dim3), intent(in) :: dimScale3
    REAL(SP), dimension(dim1,dim2,dim3), intent(in) :: data3D
    REAL(SP), intent(in) :: x1,x2,x3
    integer, intent(in), OPTIONAL :: node1_in, node2_in, node3_in
    integer, intent(out), OPTIONAL :: status

    integer :: node1, node2,node3
    REAL(SP), dimension(ndim) :: distance
    REAL(SP), dimension(2**ndim) :: weights, values
    integer :: i, j, k
    integer :: npoint

    status = 0 ! set to success

    !adapted by RF van Oss:
    IF ( PRESENT(node1_in) ) THEN
      !nodes are given on input
      node1 = node1_in
      node2 = node2_in
      node3 = node3_in
    ELSE
      node1 = findNode(dimScale1, dim1, x1)
      node2 = findNode(dimScale2, dim2, x2)
      node3 = findNode(dimScale3, dim3, x3)
    ENDIF

    if ( min(node1,node2,node3) <= 0 ) then ! error getting nodes, exit with error
      status = -1
      interpol3D = -huge(interpol3D)
      return
    end if

    npoint = 0
    do i = node1, node1 + 1
      do j=node2, node2 + 1
        do k=node3, node3 + 1

          npoint = npoint + 1

          distance(1) = (x1 - dimscale1(i)) / (dimscale1(node1+1) - dimscale1(node1))

          distance(2) = (x2 - dimscale2(j)) / (dimscale2(node2+1) - dimscale2(node2))

          distance(3) = (x3 - dimscale3(k)) / (dimscale3(node3+1) - dimscale3(node3))

          distance = 1. - abs(distance)

          weights(npoint) = product(distance)

          values(npoint) = data3D(i,j,k)

        end do
      end do
    end do

    ! normalize weights
    weights(1:npoint) =  weights(1:npoint) / sum(weights(1:npoint))

    ! determine output
    interpol3D = 0.
    do i=1, npoint
      interpol3D = interpol3D + weights(i) * values(i)
    end do

    return

  end function interpol3D

  !==============================================================================!
  REAL(SP) function interpol4D(dim1, dim2, dim3 , dim4, &
      dimScale1, dimScale2, dimScale3, dimScale4, &
      data4D, x1, x2, x3, x4,&
      node1_in, node2_in, node3_in, node4_in,  &
      status)

    integer, parameter :: ndim=4

    integer, intent(in) :: dim1, dim2, dim3, dim4
    REAL(SP), dimension(dim1), intent(in) :: dimScale1
    REAL(SP), dimension(dim2), intent(in) :: dimScale2
    REAL(SP), dimension(dim3), intent(in) :: dimScale3
    REAL(SP), dimension(dim4), intent(in) :: dimScale4
    REAL(SP), dimension(dim1,dim2,dim3,dim4), intent(in) :: data4D
    REAL(SP), intent(in) :: x1,x2,x3,x4
    integer, intent(in), OPTIONAL :: node1_in, node2_in, node3_in, node4_in
    integer, intent(out), OPTIONAL :: status

    integer :: node1, node2,node3,node4
    REAL(SP), dimension(ndim) :: distance
    REAL(SP), dimension(2**ndim) :: weights, values
    integer :: i, j, k,l
    integer :: npoint

    status = 0 ! set to success

    !adapted by RF van Oss:
    IF ( PRESENT(node1_in) ) THEN
      !nodes are given on input
      node1 = node1_in
      node2 = node2_in
      node3 = node3_in
      node4 = node4_in
    ELSE
      node1 = findNode(dimScale1, dim1, x1)
      node2 = findNode(dimScale2, dim2, x2)
      node3 = findNode(dimScale3, dim3, x3)
      node4 = findNode(dimScale4, dim4, x4)
    ENDIF

    if ( min(node1,node2,node3,node4) <= 0 ) then ! error getting nodes, exit with error
      status = -1
      interpol4D = -huge(interpol4D)
      return
    end if

    npoint = 0
    do i = node1, node1 + 1
      do j=node2, node2 + 1
        do k=node3, node3 + 1
          do l=node4, node4 + 1

            npoint = npoint + 1

            distance(1) = (x1 - dimscale1(i)) / (dimscale1(node1+1) - dimscale1(node1))

            distance(2) = (x2 - dimscale2(j)) / (dimscale2(node2+1) - dimscale2(node2))

            distance(3) = (x3 - dimscale3(k)) / (dimscale3(node3+1) - dimscale3(node3))

            distance(4) = (x4 - dimscale4(l)) / (dimscale4(node4+1) - dimscale4(node4))

            distance = 1. - abs(distance)

            weights(npoint) = product(distance)

            values(npoint) = data4D(i,j,k,l)

          end do
        end do
      end do
    end do

    ! normalize weights
    weights(1:npoint) =  weights(1:npoint) / sum(weights(1:npoint))

    ! determine output
    interpol4D = 0.
    do i=1, npoint
      interpol4D = interpol4D + weights(i) * values(i)
    end do

    return

  end function interpol4D

  !==============================================================================!
  REAL(SP) function interpol5D( &
      dim1, dim2, dim3 , dim4, dim5,&
      dimScale1, dimScale2, dimScale3, dimScale4, dimscale5, &
      data5D, &
      x1, x2, x3, x4, x5, &
      node1_in, node2_in, node3_in, node4_in, node5_in, &
      status)

    integer, parameter :: ndim=5

    !INPUT
    integer, intent(in) :: dim1, dim2, dim3, dim4, dim5
    REAL(SP), dimension(dim1), intent(in) :: dimScale1
    REAL(SP), dimension(dim2), intent(in) :: dimScale2
    REAL(SP), dimension(dim3), intent(in) :: dimScale3
    REAL(SP), dimension(dim4), intent(in) :: dimScale4
    REAL(SP), dimension(dim5), intent(in) :: dimScale5
    REAL(SP), dimension(dim1,dim2,dim3,dim4,dim5), intent(in) :: data5D
    REAL(SP), intent(in) :: x1,x2,x3,x4,x5
    integer, intent(in), OPTIONAL :: node1_in, node2_in, node3_in, node4_in, node5_in
    integer, intent(out), OPTIONAL  :: status

    !LOCAL
    integer :: node1, node2, node3, node4, node5
    REAL(SP), dimension(ndim) :: distance
    REAL(SP), dimension(2**ndim) :: weights, values
    integer :: i, j, k, l, m
    integer :: npoint

    status = 0 ! set to success

    !adapted by RF van Oss:
    IF ( PRESENT(node1_in) ) THEN
      !nodes are given on input
      node1 = node1_in
      node2 = node2_in
      node3 = node3_in
      node4 = node4_in
      node5 = node5_in
    ELSE
      node1 = findNode(dimScale1, dim1, x1)
      node2 = findNode(dimScale2, dim2, x2)
      node3 = findNode(dimScale3, dim3, x3)
      node4 = findNode(dimScale4, dim4, x4)
      node5 = findNode(dimScale5, dim5, x5)
    ENDIF

    if ( min(node1,node2,node3,node4,node5) <= 0 ) then
      !invalid nodes, exit with error
      status = -1
      interpol5D = -huge(interpol5D)
      return
    end if

    if ( min(dim1-node1,dim2-node2,dim3-node3,dim4-node4,dim5-node5) <= 0 ) then
      !invalid nodes, exit with error
      status = -2
      interpol5D = -huge(interpol5D)
      return
    end if

    npoint = 0
    do i = node1, node1 + 1
      do j=node2, node2 + 1
        do k=node3, node3 + 1
          do l=node4, node4 + 1
            do m=node5, node5 + 1

              npoint = npoint + 1

              distance(1) = (x1 - dimscale1(i)) / (dimscale1(node1+1) - dimscale1(node1))

              distance(2) = (x2 - dimscale2(j)) / (dimscale2(node2+1) - dimscale2(node2))

              distance(3) = (x3 - dimscale3(k)) / (dimscale3(node3+1) - dimscale3(node3))

              distance(4) = (x4 - dimscale4(l)) / (dimscale4(node4+1) - dimscale4(node4))

              distance(5) = (x5 - dimscale5(m)) / (dimscale5(node5+1) - dimscale5(node5))

              distance = 1. - abs(distance)

              weights(npoint) = product(distance)

              values(npoint) = data5D(i,j,k,l,m)

            end do
          end do
        end do
      end do
    end do

    ! normalize weights
    weights(1:npoint) =  weights(1:npoint) / sum(weights(1:npoint))

    ! determine output
    interpol5D = 0.
    do i=1, npoint
      interpol5D = interpol5D + weights(i) * values(i)
    end do

    return

  end function interpol5D

  !==============================================================================!
  REAL(SP) function interpol6D(dim1, dim2, dim3 , dim4, dim5, dim6, &
      dimScale1, dimScale2, dimScale3, dimScale4, dimscale5, dimscale6, &
      data6D, &
      x1, x2, x3, x4, x5, x6, &
      node1_in, node2_in, node3_in, node4_in, node5_in, node6_in,&
      status)

    integer, parameter :: ndim=6

    integer, intent(in) :: dim1, dim2, dim3, dim4, dim5, dim6
    REAL(SP), dimension(dim1), intent(in) :: dimScale1
    REAL(SP), dimension(dim2), intent(in) :: dimScale2
    REAL(SP), dimension(dim3), intent(in) :: dimScale3
    REAL(SP), dimension(dim4), intent(in) :: dimScale4
    REAL(SP), dimension(dim5), intent(in) :: dimScale5
    REAL(SP), dimension(dim6), intent(in) :: dimScale6
    REAL(SP), dimension(dim1,dim2,dim3,dim4,dim5,dim6), intent(in) :: data6D
    REAL(SP), intent(in) :: x1,x2,x3,x4,x5,x6
    integer, intent(in), OPTIONAL :: node1_in, node2_in, node3_in, node4_in, node5_in, node6_in
    integer, intent(out), OPTIONAL  :: status

    integer :: node1, node2,node3,node4,node5,node6
    REAL(SP), dimension(ndim) :: distance
    REAL(SP), dimension(2**ndim) :: weights, values
    integer :: i, j, k, l, m, n
    integer :: npoint

    status = 0 ! set to success

    !adapted by RF van Oss:
    IF ( PRESENT(node1_in) ) THEN
      !nodes are given on input
      node1 = node1_in
      node2 = node2_in
      node3 = node3_in
      node4 = node4_in
      node5 = node5_in
      node6 = node6_in
    ELSE
      node1 = findNode(dimScale1, dim1, x1)
      node2 = findNode(dimScale2, dim2, x2)
      node3 = findNode(dimScale3, dim3, x3)
      node4 = findNode(dimScale4, dim4, x4)
      node5 = findNode(dimScale5, dim5, x5)
      node6 = findNode(dimScale6, dim6, x6)
    ENDIF

    if ( min(node1,node2,node3,node4,node5,node6) <= 0 ) then ! error getting nodes, exit with error
      status = -1
      interpol6D = -huge(interpol6D)
      return
    end if

    npoint = 0
    do i = node1, node1 + 1
      do j = node2, node2 + 1
        do k = node3, node3 + 1
          do l = node4, node4 + 1
            do m = node5, node5 + 1
              do n = node6, node6 + 1

                npoint = npoint + 1

                distance(1) = (x1 - dimscale1(i)) / (dimscale1(node1+1) - dimscale1(node1))

                distance(2) = (x2 - dimscale2(j)) / (dimscale2(node2+1) - dimscale2(node2))

                distance(3) = (x3 - dimscale3(k)) / (dimscale3(node3+1) - dimscale3(node3))

                distance(4) = (x4 - dimscale4(l)) / (dimscale4(node4+1) - dimscale4(node4))

                distance(5) = (x5 - dimscale5(m)) / (dimscale5(node5+1) - dimscale5(node5))

                distance(6) = (x6 - dimscale6(n)) / (dimscale6(node6+1) - dimscale6(node6))

                distance = 1. - abs(distance)

                weights(npoint) = product(distance)

                values(npoint) = data6D(i,j,k,l,m,n)

              end do
            end do
          end do
        end do
      end do
    end do

    ! normalize weights
    weights(1:npoint) =  weights(1:npoint) / sum(weights(1:npoint))

    ! determine output
    interpol6D = 0.
    do i=1, npoint
      interpol6D = interpol6D + weights(i) * values(i)
    end do

    return

  end function interpol6D

  !==============================================================================!
  REAL(SP) function interpol7D(dim1, dim2, dim3 , dim4, dim5, dim6, dim7,&
      dimScale1, dimScale2, dimScale3, dimScale4, dimscale5, dimscale6, dimscale7,&
      data7D, &
      x1, x2, x3, x4, x5, x6, x7, &
      node1_in, node2_in, node3_in, node4_in, node5_in, node6_in, node7_in, &
      status)

    integer, parameter :: ndim=7

    integer, intent(in) :: dim1, dim2, dim3, dim4, dim5, dim6, dim7
    REAL(SP), dimension(dim1), intent(in) :: dimScale1
    REAL(SP), dimension(dim2), intent(in) :: dimScale2
    REAL(SP), dimension(dim3), intent(in) :: dimScale3
    REAL(SP), dimension(dim4), intent(in) :: dimScale4
    REAL(SP), dimension(dim5), intent(in) :: dimScale5
    REAL(SP), dimension(dim6), intent(in) :: dimScale6
    REAL(SP), dimension(dim6), intent(in) :: dimScale7
    REAL(SP), dimension(dim1,dim2,dim3,dim4,dim5,dim6,dim7), intent(in) :: data7D
    REAL(SP), intent(in) :: x1,x2,x3,x4,x5,x6,x7
    integer, intent(in), OPTIONAL :: node1_in, node2_in, node3_in, node4_in, node5_in, node6_in, node7_in
    integer, intent(out), OPTIONAL :: status

    integer :: node1, node2,node3,node4,node5,node6,node7
    REAL(SP), dimension(ndim) :: distance
    REAL(SP), dimension(2**ndim) :: weights, values
    integer :: i, j, k, l, m, n, o
    integer :: npoint

    status = 0 ! set to success

    !adapted by RF van Oss:
    IF ( PRESENT(node1_in) ) THEN
      !nodes are given on input
      node1 = node1_in
      node2 = node2_in
      node3 = node3_in
      node4 = node4_in
      node5 = node5_in
      node6 = node6_in
      node7 = node7_in
    ELSE
      node1 = findNode(dimScale1, dim1, x1)
      node2 = findNode(dimScale2, dim2, x2)
      node3 = findNode(dimScale3, dim3, x3)
      node4 = findNode(dimScale4, dim4, x4)
      node5 = findNode(dimScale5, dim5, x5)
      node6 = findNode(dimScale6, dim6, x6)
      node7 = findNode(dimScale7, dim7, x7)
    ENDIF

    if ( min(node1,node2,node3,node4,node5,node6,node7) <= 0 ) then ! error getting nodes, exit with error
      status = -1
      interpol7D = -huge(interpol7D)
      return
    end if

    npoint = 0
    do i = node1, node1 + 1
      do j = node2, node2 + 1
        do k = node3, node3 + 1
          do l = node4, node4 + 1
            do m = node5, node5 + 1
              do n = node6, node6 + 1
                do o = node7, node7 + 1

                  npoint = npoint + 1

                  distance(1) = (x1 - dimscale1(i)) / (dimscale1(node1+1) - dimscale1(node1))

                  distance(2) = (x2 - dimscale2(j)) / (dimscale2(node2+1) - dimscale2(node2))

                  distance(3) = (x3 - dimscale3(k)) / (dimscale3(node3+1) - dimscale3(node3))

                  distance(4) = (x4 - dimscale4(l)) / (dimscale4(node4+1) - dimscale4(node4))

                  distance(5) = (x5 - dimscale5(m)) / (dimscale5(node5+1) - dimscale5(node5))

                  distance(6) = (x6 - dimscale6(n)) / (dimscale6(node6+1) - dimscale6(node6))

                  distance(7) = (x7 - dimscale7(o)) / (dimscale7(node7+1) - dimscale7(node7))

                  distance = 1. - abs(distance)

                  weights(npoint) = product(distance)

                  values(npoint) = data7D(i,j,k,l,m,n,o)

                end do
              end do
            end do
          end do
        end do
      end do
    end do

    ! normalize weights
    weights(1:npoint) =  weights(1:npoint) / sum(weights(1:npoint))

    ! determine output
    interpol7D = 0.
    do i=1, npoint
      interpol7D = interpol7D + weights(i) * values(i)
    end do

    return

  end function interpol7D

  !==============================================================================!

  SUBROUTINE GetNode(dimScale, x, Node, LowerOrUpper)

    ! This subroutine is partly based on (private) function findNode in
    ! this Module .
    ! doesn't use minloc(abs()), but instead starts from the first index and
    !  exits when the node has been found.
    !In case LowerOrUpper = 'Lower',:
    ! GetNode is the index of the largest value in the increasing array dimscale
    ! that is still smaller than x.
    !In case LowerOrUpper = 'Upper':
    ! GetNode is the index of the smallest value in the increasing array dimscale
    ! that is still larger than x.
    !
    !  Roeland van Oss, May 2002

    !INPUT
    REAL(SP), DIMENSION( : ), INTENT(IN) :: dimscale
    REAL(SP), INTENT(IN) :: x
    CHARACTER(*) :: LowerOrUpper

    !OUTPUT
    INTEGER, INTENT(OUT) :: Node

    !LOCAL
    INTEGER :: is, xdim, sdim

    !Find dimension of dimscale
    sdim = SIZE( dimscale)

    ! check if the dimension scale is monotonically increasing as it should
    IF (.NOT. monotonIncrease(dimScale,sdim) ) THEN
      Node = -1
      RETURN
    END IF

    IF ( LowerOrUpper == 'Lower' ) THEN

      is = 0
      DO
        is = is + 1
        IF ( (dimScale(is) > x) .OR. (is > sdim) ) EXIT
      END DO
      Node = is - 1

      !ERROR: x is smaller than the smallest value in dimscale
      IF ( Node == 0    ) Node = -2

      !ERROR: x is larger than the largest value in dimscale
      IF ( Node == sdim ) Node = -3

    ELSE IF ( LowerOrUpper == 'Upper' ) THEN

      is = sdim + 1
      DO
        is = is - 1
        IF ( (dimScale(is) < x) .OR. (is < 1) ) EXIT
      END DO
      Node = is + 1

      !ERROR: x is smaller than the smallest value in dimscale
      IF ( Node == 1        ) Node = -2

      !ERROR: x is larger than the largest value in dimscale
      IF ( Node == sdim + 1 ) Node = -3

    ELSE

    ENDIF

    RETURN

  END SUBROUTINE GetNode

end module OMSAO_lininterpolation_module

