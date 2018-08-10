module numFort
  use kinds
  use quadpack
  use LAPACK95
  implicit none

  real(DP), parameter :: Infty = huge(1.0_DP)
  integer  :: neval, ifail
  real(DP) :: errEstimate

  interface linspace
     module procedure linspace,linspace_real
  end interface linspace
  interface integral
     module procedure integral, integralToInfty, integralOf, integralBreakPts
  end interface integral
  interface rk4
     module procedure rk4,rk4_2,rk4_step
  end interface rk4
  interface splinefit
     module procedure splinefit,splinefit_coeff
  end interface splinefit

contains

  !---------------------------------------------------------------------!
  !                                                                     !
  !                               MeshGrid                              !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine meshgrid(x,y,XX,YY)
    use kinds
    implicit none

    real(DP),dimension(:),intent(in)    :: x
    real(DP),dimension(:),intent(in)    :: y
    real(DP),dimension(:,:),intent(out) :: XX,YY

    integer :: i,j,rows,cols

    rows = size(x)
    cols = size(y)    

    do j = 1,cols
       XX(j,1:rows) = x(1:rows)
    end do
    do j = 1,rows
       YY(1:cols,j) = y(1:cols)
    end do

  end subroutine meshgrid

  !---------------------------------------------------------------------!
  !                                                                     !
  !                           Cubic Spline Fit                          !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine splinefit(x,y,intpts,intvals)
    use kinds
    use LAPACK95
    implicit none 

    real(DP),dimension(:),intent(in)    :: x,y
    real(DP),dimension(:),intent(in)    :: intpts
    real(DP),dimension(:),intent(out)   :: intvals

    integer                             :: N
    real(DP),dimension(:,:),allocatable :: A,B
    integer ,dimension(:)  ,allocatable :: ipiv
    integer                             :: i,j

    N = size(x)
    allocate(A(N,N),B(N,1),ipiv(N))

    B(:,1) = y(:)
    A(:,1) = 1.0_DP
    A(:,2) = x(:)
    do i = 1,N
       do j = 3,N
          A(i,j) = abs(x(i)-x(j-1))
       end do
    end do
    call getrf(A,ipiv)
    call getrs(A,ipiv,B)    
    do i = 1,size(intpts)
       intvals(i) = B(1,1) + B(2,1)*intpts(i)+&
            & sum( B(3:N,1)*abs(intpts(i)-x(2:N-1)) )
    end do

    deallocate(A,B,ipiv)

  end subroutine splinefit

  subroutine splinefit_coeff(x,y,c)
    use kinds
    use LAPACK95
    implicit none 

    real(DP),dimension(:),intent(in)    :: x,y
    real(DP),dimension(:),intent(out)   :: c

    integer                             :: N
    real(DP),dimension(:,:),allocatable :: A,B
    integer ,dimension(:)  ,allocatable :: ipiv
    integer                             :: i,j

    N = size(x)
    allocate(A(N,N),B(N,1),ipiv(N))

    B(:,1) = y(:)
    A(:,1) = 1.0_DP
    A(:,2) = x(:)
    do i = 1,N
       do j = 3,N
          A(i,j) = abs(x(i)-x(j-1))
       end do
    end do
    call getrf(A,ipiv)
    call getrs(A,ipiv,B)       
    c(:) = B(:,1)

    deallocate(A,B,ipiv)

  end subroutine splinefit_coeff


  !---------------------------------------------------------------------!
  !                                                                     !
  !                  Calculate cubic splines value at x                 !
  !                                                                     !
  !---------------------------------------------------------------------!

  function splinevals(c,xj,x)
    use kinds
    implicit none

    real(DP),dimension(:),intent(in) :: c,xj
    real(DP),intent(in) :: x
    real(DP) :: splinevals

    integer :: N

    N = size(c)
    splinevals = c(1)+c(2)*x+sum(c(3:N)*abs(x-xj(2:N-1)))

  end function splinevals

  !---------------------------------------------------------------------!
  !                                                                     !
  !            Polynomial fitting and calculation Algorithim            !
  !                                                                     !
  !---------------------------------------------------------------------!

  function polycal(N,c,x)
    implicit none 

    integer ,intent(in)                :: N
    real(DP),dimension(N+1),intent(in) :: c
    real(DP)               ,intent(in) :: x

    real(DP)                           :: polycal
    real(DP),dimension(N+1)            :: xx
    integer ,dimension(N+1)            :: indexx
    integer                            :: i

    do i = 0,N
       indexx(i+1) = N-i
    end do
    xx = x**indexx
    polycal = sum(c*xx)

  end function polycal

  subroutine polyfit(x,y,N,c)    
    use kinds
    use LAPACK95
    implicit none

    integer ,intent(in)              :: N
    real(DP),dimension(:),intent(in) :: x,y    
    real(DP),dimension(N+1),intent(out) :: c

    real(DP),dimension(N+1,N+1) :: A
    real(DP),dimension(N+1,1) :: B
    integer ,dimension(N+1)   :: ipiv

    integer :: i,j,L

    if ( N > size(x) ) then
       write(*,'(a)') "Error, order of polynomial must be less then the number of entered points"
    else

       L = N+1
       do i = 1,L
          do j = 1,L          
             A(i,j) = x(i)**(L-j)
          end do
          B(i,1) = y(i)
       end do

       call getrf(A,ipiv)
       call getrs(A,ipiv,B)

       c(:) = B(:,1)

    end if

  end subroutine polyfit

  !---------------------------------------------------------------------!
  !                                                                     !
  !                  Runge Kutta 4th order (1 Equation)                 !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine rk4(t0,tN,y0,N,f,t,y)
    use kinds
    implicit none

    real(DP),intent(in) :: t0,tN,y0
    integer ,intent(in) :: N
    interface
       function f(t,y)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: t,y
       end function f
    end interface
    real(DP),dimension(N),intent(out) :: y,t

    real(DP) :: h,k1,k2,k3,k4
    integer  :: i

    h = (tN-t0)/N
    t(1) = t0
    t(N) = tN
    do i = 2,N-1
       t(i) = t(0) + h*(i-1)
    end do
    y(:) = 0
    y(1) = y0
    do i = 1,N-1
       k1 = f(t(i),y(i))
       k2 = f(t(i)+0.5_dp*h,y(i)+0.5_dp*h*k1)
       k3 = f(t(i)+0.5_dp*h,y(i)+0.5_dp*h*k2)
       k4 = f(t(i)+h,y(i)+k3*h)       
       y(i+1) = y(i) + (1/6.0_dp)*(k1+2*k2+2*k3+k4)*h

    end do

  end subroutine rk4

  subroutine rk4_step(h,t0,y0,f,y)
    use kinds
    implicit none

    real(DP),intent(in) :: h,t0,y0
    interface
       function f(t,y)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: t,y
       end function f
    end interface
    real(DP),intent(out) :: y

    real(DP) :: k1,k2,k3,k4

    k1 = f(t0,y0)
    k2 = f(t0+0.5_dp*h,y0+0.5_dp*h*k1)
    k3 = f(t0+0.5_dp*h,y0+0.5_dp*h*k2)
    k4 = f(t0+h,y0+k3*h)       
    y = y0 + (1/6.0_dp)*(k1+2*k2+2*k3+k4)*h

  end subroutine rk4_step


  !---------------------------------------------------------------------!
  !                                                                     !
  !                 Runge Kutta 4th Order (2 Equations)                 !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine rk4_2(h,t0,y0,f,g,y1,y2)
    use kinds
    implicit none

    real(DP),dimension(2) ,intent(in) :: y0
    real(DP),intent(in)               :: h,t0
    interface
       function f(t,y1,y2)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: t,y1,y2
       end function f
    end interface
    interface
       function g(t,y1,y2)
         use kinds
         implicit none
         real(DP)             :: g
         real(DP) ,intent(in) :: t,y1,y2
       end function g
    end interface
    real(DP),intent(out) :: y1,y2

    real(DP) :: k1,k2,k3,k4,l1,l2,l3,l4

    k1 = f(t0,y0(1),y0(2))
    l1 = g(t0,y0(1),y0(2))
    k2 = f(t0+0.5_dp*h,y0(1)+0.5_dp*h*k1,y0(2)+0.5_dp*h*l1)
    l2 = g(t0+0.5_dp*h,y0(1)+0.5_dp*h*k1,y0(2)+0.5_dp*h*l1)
    k3 = f(t0+0.5_dp*h,y0(1)+0.5_dp*h*k2,y0(2)+0.5_dp*h*l2)
    l3 = g(t0+0.5_dp*h,y0(1)+0.5_dp*h*k2,y0(2)+0.5_dp*h*l2)
    k4 = f(t0+h,y0(1)+k3*h,y0(2)+l3*h)
    l4 = g(t0+h,y0(1)+k3*h,y0(2)+l3*h)       
    y1 = y0(1) + (1/6.0_dp)*(k1+2*k2+2*k3+k4)*h
    y2 = y0(2) + (1/6.0_dp)*(l1+2*l2+2*l3+l4)*h

  end subroutine rk4_2

  !---------------------------------------------------------------------!
  !                                                                     !
  !                        Find first zero guess                        !
  !                                                                     !
  !---------------------------------------------------------------------!

  function GuessZero(fvals)
    use Kinds
    implicit none
    real(DP),dimension(:),intent(in) :: fvals
    integer                          :: GuessZero

    integer                          :: jj

    do jj=3,size(fvals)
       if ( sign( fvals(jj-1)/abs(fvals(jj-1)) , fvals(jj)/abs(fvals(jj)) ) .ne. &
            sign( fvals(jj-2)/abs(fvals(jj-2)) , fvals(jj-1)/abs(fvals(jj-1)) ) ) exit
    end do

    GuessZero = jj

  end function GuessZero

  !---------------------------------------------------------------------!
  !                                                                     !
  !                           Newtons Method                            !
  !                                                                     !
  !---------------------------------------------------------------------!

  function newton1D(fn,guess)
    use Kinds
    implicit none
    real(DP) :: newton1D
    interface
       function fn(x)
         use kinds
         implicit none
         real(DP)             :: fn
         real(DP), intent(in) :: x
       end function fn
    end interface
    real(DP) :: guess
    ! Local Variables
    real(DP) :: newt_tol = 1e-6
    real(DP) :: h ! = 1e-6_DP
    real(DP) :: fx, fxfh,fxbh, dfdx 
    integer :: attempt
    integer :: attempt_limit = 20
    logical :: verbose = .false.

    if (verbose) write(*,*) 'Newton-Raphson Method'
    h = guess * 1.0e-6_DP
    newton1D = guess
    attempt = 1

    do
       fx = fn(newton1D)
       if (abs(fx) <= newt_tol) then
          if (verbose) then
             write(*,'(a,i3,a)') 'Success after', attempt, ' attempts.'
             write(*,*) 'Zero  = ', newton1D
          end if
          exit
       else if (attempt >= attempt_limit) then
          write(*,*) 'Failed, change initial guess or increase attempt limit.'
          exit
       end if

       fxfh = fn(newton1D + h/2)
       fxbh = fn(newton1D - h/2)
       dfdx = (fxfh-fxbh) / h
       newton1D = newton1D - fx / dfdx

       attempt = attempt + 1
    end do
  end function newton1D


  !---------------------------------------------------------------------!
  !                                                                     !
  !                          Linspace Function                          !
  !                                                                     !
  !---------------------------------------------------------------------!

  function linspace(start,finish,N)
    integer                :: N
    real(DP), intent(in)   :: start, finish
    real(DP), dimension(N) :: linspace

    real(DP)               :: int
    integer                :: i

    linspace(1) = start
    linspace(N) = finish
    int  = (real(finish)-start)/(N-1)

    do i=2,N-1
       linspace(i)=linspace(i-1)+int
    end do

  end function linspace

  function linspace_real(start,finish,N)
    integer                :: N
    real, intent(in)       :: start, finish
    real, dimension(N)     :: linspace_real

    real(DP)               :: int
    integer                :: i

    linspace_real(1) = start
    linspace_real(N) = finish
    int  = (real(finish)-start)/(N-1)

    do i=2,N-1
       linspace_real(i)=linspace_real(i-1)+int
    end do

  end function linspace_real

  function linspace_int(start,finish,N)
    integer                :: N
    integer, intent(in)       :: start, finish
    integer, dimension(N)     :: linspace_int

    real(DP)               :: int
    integer                :: i

    linspace_int(1) = start
    linspace_int(N) = finish
    int  = (real(finish)-start)/(N-1)

    do i=2,N-1
       linspace_int(i)=linspace_int(i-1)+int
    end do

  end function linspace_int

  !---------------------------------------------------------------------!
  !                                                                     !
  !                          Finite Difference                          !
  !                                                                     !
  !---------------------------------------------------------------------!

  function deriv(f,x0)
    use Kinds
    implicit none
    real(DP)                 :: deriv
    interface
       function f(x)
         use kinds
         implicit none
         real(DP)            :: f
         real(DP),intent(in) :: x
       end function f
    end interface
    real(DP),intent(in)      :: x0
    real(DP)                 :: h,fxfh,fxbh

    h     = x0*1e-6_DP
    fxfh  = f(x0+h/2)
    fxbh  = f(x0-h/2)    
    deriv = (fxfh-fxbh)/h

  end function deriv

  !---------------------------------------------------------------------!
  !                                                                     !
  !                      QuadPack Integral wrapper                      !
  !                                                                     !
  !---------------------------------------------------------------------!

  function integral(f, a, b, absErr, relErr)
    implicit none
    real(DP)                  :: integral
    interface                                   ! Interfaces: Sections 5.11, 5.18
       function f(x)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: x
       end function f
    end interface
    real(DP), intent(in)      :: a, b, absErr, relErr

    !  Local variables
    !
    real(DP)                  :: integResult, bound=0.0_DP
    integer                   :: inf

    !  Determine if the limits include infinity and call qagi if nessary
    !
    if (a == -Infty) then
       if (b == Infty) then
          inf=2
          call qagi(f, bound, inf, absErr, relErr, integResult, errEstimate, neval, ifail)
          if ( ifail /= 0 ) then
             write(*,*) 'Warning from qagi: the error code is ', ifail
          end if
       else
          inf = -1
          bound = b
          call qagi(f, bound, inf, absErr, relErr, integResult, errEstimate, neval, ifail)
          if ( ifail /= 0 ) then
             write(*,*) 'Warning from qagi: the error code is ', ifail
          end if
       end if
    else 
       if (b == Infty) then
          inf = 1
          bound = a
          call qagi(f, bound, inf, absErr, relErr, integResult, errEstimate, neval, ifail)
          if ( ifail /= 0 ) then
             write(*,*) 'Warning from qagi: the error code is ', ifail
          end if
       else
          call qags(f, a, b, absErr, relErr, integResult, errEstimate, neval, ifail)
          if ( ifail /= 0 ) then
             write(*,*) 'Warning from qags: the error code is ', ifail
          end if
       end if
    end if
    integral = integResult

  end function integral


  function integralToInfty(f, bound, absErr, relErr)
    implicit none
    real(DP)                  :: integralToInfty
    interface
       function f(x)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: x
       end function f
    end interface
    real(DP), intent(in)      :: bound, absErr, relErr

    !  Local variables
    !
    real(DP)                  :: integResult
    integer, parameter        :: inf=1

    call qagi(f, bound, inf, absErr, relErr, integResult, errEstimate, neval, ifail)
    if ( ifail /= 0 ) then
       write(*,*) 'Warning from qagi: the error code is ', ifail
    end if
    integralToInfty = integResult

  end function integralToInfty


  function integralOf(f, absErr, relErr)
    implicit none
    real(DP)                  :: integralOf
    interface
       function f(x)
         implicit none
         integer,  parameter  :: DP = kind(1.0d0)
         real(DP)             :: f
         real(DP), intent(in) :: x
       end function f
    end interface
    real(DP), intent(in)      :: absErr, relErr

    !  Local variables
    !
    real(DP)                  :: integResult, bound=0.0d0
    integer, parameter        :: inf=2

    call qagi(f, bound, inf, absErr, relErr, integResult, errEstimate, neval, ifail)
    if ( ifail /= 0 ) then
       write(*,*) 'Warning from qagi: the error code is ', ifail
    end if
    integralOf = integResult

  end function integralOf


  function integralBreakPts(f, a, b, absErr, relErr, nBreakPts, BreakPts)
    implicit none
    real(DP)                  :: integralBreakPts
    interface
       function f(x)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: x
       end function f
    end interface
    real(DP), intent(in)      :: a, b, absErr, relErr
    integer,  intent(in)      :: nBreakPts
    real(DP), intent(in), dimension(nBreakPts) :: BreakPts

    !  Local variables
    !
    real(DP)                         :: integResult
    real(DP), dimension(nBreakPts+2) :: BreakPtsP2        ! Automatic array.  Similar to allocatable arrays.

    BreakPtsP2(1:nBreakPts) = BreakPts(1:nBreakPts)       ! Array section limits are required here.

    call qagp(f, a, b, nBreakPts+2, BreakPtsP2, absErr, relErr, integResult, errEstimate, neval, ifail)
    if ( ifail /= 0 ) then
       write(*,*) 'Warning from qagp: the error code is ', ifail
    end if
    integralBreakPts = integResult

  end function integralBreakPts

end module numFort
