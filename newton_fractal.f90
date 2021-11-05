! Newton fractal
!##############################################################################!

module newton

  implicit none

  contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horner( params, x, res )

    ! variable declarations
    !--------------------------------------------------------------------------!

    ! input arguments
    !..........................................................................!

    real, dimension ( : ), intent ( in ) :: params

    complex( 8 ), intent ( in ) :: x

    ! output arguments
    !..........................................................................!

    complex( 8 ), intent( out ) :: res

    integer :: i

    res = 0.0
    do i = size ( params ), 1, -1
      res = res * x + params ( i )
    end do

  end subroutine horner

  !----------------------------------------------------------------------------!

  subroutine polyder( params, params_der )

    real, dimension ( : ), intent ( in ) :: params

    real, dimension ( : ), intent ( out ) :: params_der

    real, dimension ( : ), allocatable :: arange, params_reverse

    integer :: m, i

    m = size( params )

    allocate( arange( m ), params_reverse( m ) )

    arange = [ ( i, i = m, 1, -1 ) ]
    params_reverse = params( m : 1 : -1 )

    params_der( : ) = 0.
    params_der( 2 : ) = params_reverse( : m -1 )
    params_der = params_der * arange
    params_der = params_der( m : 1 : -1 )

  end subroutine polyder

  !----------------------------------------------------------------------------!

  subroutine main( params, roots, max_iter, out, convergence )

    implicit none

    ! variable declarations
    !--------------------------------------------------------------------------!

    ! input arguments
    !..........................................................................!

    real( 4 ), dimension ( : ), intent ( in ) :: params

    integer( 4 ), intent( in ) :: max_iter

    complex( 8 ), dimension( : ), intent( in ) :: roots

    ! output arguments
    !..........................................................................!

    integer( 4 ), intent( inout ) :: out( :, : ), convergence( :, : )

    ! internal variables
    !..........................................................................!

    complex( 8 ) :: z, numerator, denominator

    integer( 4 ) :: N, m, i, j, s

    real( 4 ) :: convergence_threshold

    real( 4 ), dimension ( : ), allocatable :: params_der, params_reverse

    ! variable initializations
    !--------------------------------------------------------------------------!

    N = size( out, 1 )
    m = size( params )

    convergence_threshold = 1e-14

    allocate( params_der( m ), params_reverse( m ) )

    params_reverse = params( m : 1 : -1 )

    call polyder( params_reverse, params_der )

    out = 0.
    convergence = max_iter

    ! main loop
    !--------------------------------------------------------------------------!

    !$omp parallel do default(none) &
    !$omp& shared(max_iter, N, out, params_reverse, params_der, roots, convergence, convergence_threshold) &
    !$omp& private(s, i, j, z, numerator, denominator)

    do i = 1, N
      do j = 1, N

        z = cmplx( 4. * float( i ) / ( N - 1. ) - 2, 4. * float( j ) / ( N - 1. ) - 2 )

        do s = 1, max_iter

          call horner( params_reverse, z, numerator )
          call horner( params_der, z, denominator )

          z = z - numerator / denominator

          if ( minval( abs( z - roots ) ) .LE. convergence_threshold ) then
            convergence( i, j ) = s
            exit
          end if

        end do

        out( i, j ) = minloc( abs( z - roots ), dim = 1 )

      end do
    end do

  end subroutine main

  !----------------------------------------------------------------------------!

end module newton

!##############################################################################!