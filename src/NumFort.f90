!*****************************************************************************!
!                               NUMFORT                                       !
!              Authors: Anthony Kalaitzis, Curtis Abell,                      !
!                 Derek Leinweber and John Burkardt                           !
!                                                                             !
! Numerical library for fortran. Use the following module titles              !
! to Move Around easily.                                                      !
!                                                                             !
! 1. QuadPack                                                                 !
! 2. NumFort                                                                  !
!                                                                             !
! Below contains a list of the function names inside numFort for easier       !
! navigation purposes. Jump to the bottom for numFort.                        !
!                                                                             !
! - factorial                                                                 !
! - bessel                                                                    !
! - deriv                                                                     !
! - Trace                                                                     !
! - inv                                                                       !
! - meshgrid                                                                  !
! - FFT                                                                       !
! - splinefit                                                                 !
! - splineval                                                                 !
! - polyfit                                                                   !
! - polyval                                                                   !
! - polyInt                                                                   !
! - minimize                                                                  !
! - linspace                                                                  !
! - PrintTime                                                                 !
! - guessZero                                                                 !
! - Newton1D                                                                  !
! - Euler                                                                     !
! - rk4                                                                       !
! - trapZ                                                                     !
! - integral                                                                  !
! - integralPV                                                                !
! - writeData                                                                 !
! - LoopTimer                                                                 !
! - pyplot                                                                    !
!                                                                             !
!*****************************************************************************!

!---------------------------------------------------------------------!
!                                                                     !
!                      NumFort numerical library                      !
!                                                                     !
!---------------------------------------------------------------------!

module numFort
  use kinds
  use quadpack
  use LAPACK95
  use minFun
  implicit none

  interface linspace
     module procedure linspace,linspaceReal,linspaceInt
  end interface linspace
  interface integral
     module procedure integral, integralToInfty, integralOf, integralBreakPts
  end interface integral
  interface rk4
     module procedure rk4,rk4N
  end interface rk4
  interface GuessZero
     module procedure GuessZero,GuessZeroF
  end interface GuessZero
  interface Trace
     module procedure TraceDP,TraceSP,TraceComplexDP,TraceComplexSP
  end interface Trace
  interface inv
     module procedure invMatSP, invMatDP, invMatComplexSP, invMatComplexDP
  end interface inv
  interface pyplot
     module procedure pyplotN,pyplotXY,pyplotXYZW
  end interface pyplot
  interface Euler
     module procedure euler,eulerND
  end interface Euler
  interface writeData
     module procedure writeDataXY,writeDataXYZW,writeDataN,writeDataPXY,writeDataPXYZW,writeDataPN
  end interface writeData
  interface trapZ
     module procedure trapZ,trapZZ
  end interface trapZ

contains

  !---------------------------------------------------------------------!
  !                                                                     !
  !                            FFT Algorithim                           !
  !                                                                     !
  !---------------------------------------------------------------------!

  ! In place Cooley-Tukey FFT
  recursive subroutine fft(x)
    complex(DP), dimension(:), intent(inout)  :: x
    complex(DP)                               :: t

    complex(DP), dimension(:), allocatable    :: even, odd
    integer :: N
    integer :: i

    N=size(x)

    if(N .le. 1) return

    allocate(odd((N+1)/2))
    allocate(even(N/2))

    ! divide
    odd  = x(1:N:2)
    even = x(2:N:2)

    ! conquer
    call fft(odd)
    call fft(even)

    ! combine
    do i=1,N/2
       t=exp(cmplx(0.0_DP,-2.0_DP*pi*real(i-1,DP)/real(N,DP),DP))*even(i)
       x(i)     = odd(i) + t
       x(i+N/2) = odd(i) - t
    end do

    deallocate(odd)
    deallocate(even)

  end subroutine fft

  !---------------------------------------------------------------------!
  !                                                                     !
  !                         Write data to a file                        !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine writeDataXY(x,y,title)
    use Kinds
    implicit none
    real(DP),dimension(:),intent(in)     :: x,y
    character(len=*),intent(in),optional :: title
    integer :: ii,N

    N = size(x)
    if (present(title)) then
       open(110,file=trim(title),action="write", &
            & status="replace",form="formatted")
    else
       open(110,file="data.dat",action="write", &
            & status="replace",form="formatted")
    end if
    do ii = 1,N
       write(110,'(2e30.15e4)') x(ii), y(ii)
    end do

    close(110)

  end subroutine writeDataXY

  subroutine writeDataXYZW(x,y,z,w,title)
    use Kinds
    implicit none
    real(DP),dimension(:),intent(in)     :: x,y,z,w
    character(len=*),intent(in),optional :: title
    integer :: ii,N

    N = size(x)
    if (present(title)) then
       open(110,file=trim(title),action="write", &
            & status="replace",form="formatted")
    else
       open(110,file="data.dat",action="write", &
            & status="replace",form="formatted")
    end if
    do ii = 1,N
       write(110,'(4e30.15e4)') x(ii), y(ii), z(ii), w(ii)
    end do

    close(110)

  end subroutine writeDataXYZW

  subroutine writeDataN(x,title)
    use Kinds
    implicit none
    real(DP),dimension(:,:),intent(in)     :: x
    character(len=*),intent(in),optional   :: title
    character(len=14)                      :: fmt
    integer :: ii,N

    N = size(x,dim=2)
    if ( N < 10) then
       write(fmt,'(a1,i1,a9)') '(', N, 'e30.15e4)'
    else
       write(fmt,'(a1,i2,a9)') '(', N, 'e30.15e4)'
    end if

    if (present(title)) then
       open(110,file=trim(title),action="write", &
            & status="replace",form="formatted")
    else
       open(110,file="data.dat",action="write", &
            & status="replace",form="formatted")
    end if
    do ii = 1,size(x,dim=1)
       write(110,fmt) x(ii,:)
    end do

    close(110)

  end subroutine writeDataN

  !---------------------------------------------------------------------!
  !                                                                     !
  !                        Write data with a path                       !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine writeDataPXY(x,y,pwd,title)
    use Kinds
    implicit none
    real(DP),dimension(:),intent(in)     :: x,y
    character(len=*),intent(in)          :: title,pwd
    integer :: ii,N

    N = size(x)
    open(110,file=trim(pwd)//trim(title),action="write", &
         & status="replace",form="formatted")
    do ii = 1,N
       write(110,'(2e30.15e4)') x(ii), y(ii)
    end do

    close(110)

  end subroutine writeDataPXY

  subroutine writeDataPXYZW(x,y,z,w,pwd,title)
    use Kinds
    implicit none
    real(DP),dimension(:),intent(in)     :: x,y,z,w
    character(len=*),intent(in)          :: title,pwd
    integer :: ii,N

    N = size(x)
    open(110,file=trim(pwd)//trim(title),action="write", &
         & status="replace",form="formatted")
    do ii = 1,N
       write(110,'(4e30.15e4)') x(ii), y(ii), z(ii), w(ii)
    end do

    close(110)

  end subroutine writeDataPXYZW

  subroutine writeDataPN(x,pwd,title)
    use Kinds
    implicit none
    real(DP),dimension(:,:),intent(in)     :: x
    character(len=*),intent(in)            :: title,pwd
    character(len=14)                      :: fmt
    integer :: ii,N

    N = size(x,dim=2)
    if ( N < 10) then
       write(fmt,'(a1,i1,a9)') '(', N, 'e30.15e4)'
    else
       write(fmt,'(a1,i2,a9)') '(', N, 'e30.15e4)'
    end if

    open(110,file=trim(pwd)//trim(title),action="write", &
         & status="replace",form="formatted")
    do ii = 1,size(x,dim=1)
       write(110,fmt) x(ii,:)
    end do

    close(110)

  end subroutine writeDataPN

  !---------------------------------------------------------------------!
  !                                                                     !
  !             evaluate Bessel function at x of order 0<n<3            !
  !                                                                     !
  !---------------------------------------------------------------------!

  function bessel(x,n) result(val)
    use kinds
    implicit none

    real(DP),intent(in) :: x
    integer ,intent(in) :: n
    real(DP)            :: val

    if( n > 5 .or. n < 0 ) then
       write(*,'(a)') "Error: n must be an integer between 0 and 5"
    else

       select case(n)
       case(0)
          val = sin(x)/x
       case(1)
          val = sin(x)/x**2-cos(x)/x
       case(2)
          val = (3/x**2-1)*sin(x)/x -3*cos(x)/x**2
       case(3)
          val = (15/x**3-6/x)*sin(x)/x-(15/x**2-1)*sin(x)/x
       case(4)
          val = (sin(x)+4*cos(x)/x-12*sin(x)/x**2-24*cos(x)/x**3+24*sin(x)/x**4)/x
       case(5)
          val = (cos(x)-5*sin(x)/x-20*cos(x)/x**2+60*sin(x)/x**3+120*cos(x)/x**4-120*sin(x)/x**5)/x
       end select
    end if

  end function bessel

  !---------------------------------------------------------------------!
  !                                                                     !
  !                   Calculate the trace of a Matrix                   !
  !                                                                     !
  !---------------------------------------------------------------------!

  function TraceDP(M)
    use kinds
    implicit none

    real(DP),dimension(:,:),intent(in) :: M
    real(DP)                           :: TraceDP

    integer :: rows,cols,i

    rows = size(M(:,1))
    cols = size(M(1,:))
    TraceDP = 0
    if ( rows .ne. cols ) then
       write(*,'(a)') "Error, not a square matrix"
    else
       do i = 1,rows
          TraceDP = TraceDP + M(i,i)
       end do
    end if

  end function TraceDP

  function TraceSP(M)
    use kinds
    implicit none

    real(SP),dimension(:,:),intent(in) :: M
    real(SP)                           :: TraceSP

    integer :: rows,cols,i

    rows = size(M(:,1))
    cols = size(M(1,:))
    TraceSP = 0
    if ( rows .ne. cols ) then
       write(*,'(a)') "Error, not a square matrix"
    else
       do i = 1,rows
          TraceSP = TraceSP + M(i,i)
       end do
    end if

  end function TraceSP

  function TraceComplexDP(M)
    use kinds
    implicit none

    complex(DP),dimension(:,:),intent(in) :: M
    complex(DP)                           :: TraceComplexDP

    integer :: rows,cols,i

    rows = size(M(:,1))
    cols = size(M(1,:))
    TraceComplexDP = 0.0_DP
    if ( rows .ne. cols ) then
       write(*,'(a)') "Error, not a square matrix"
    else
       do i = 1,rows
          TraceComplexDP = TraceComplexDP + M(i,i)
       end do
    end if

  end function TraceComplexDP

  function TraceComplexSP(M)
    use kinds
    implicit none

    complex(SP),dimension(:,:),intent(in) :: M
    complex(SP)                           :: TraceComplexSP

    integer :: rows,cols,i

    rows = size(M(:,1))
    cols = size(M(1,:))
    TraceComplexSP = 0.0_DP
    if ( rows .ne. cols ) then
       write(*,'(a)') "Error, not a square matrix"
    else
       do i = 1,rows
          TraceComplexSP = TraceComplexSP + M(i,i)
       end do
    end if

  end function TraceComplexSP

  !---------------------------------------------------------------------!
  !                                                                     !
  !                  Calculate the inverse of a matrix                  !
  !                                                                     !
  !---------------------------------------------------------------------!

  function invMatSP(A) result(B)
    real(SP), intent(in) :: A(:,:)
    real(SP)             :: B(size(A,1),size(A,2))
    ! Local variables
    integer              :: n

    n = size(A,1)
    select case (n)
    case (1)
       if (A(1,1).eq.0.0) then
          write(*,*) 'Matrix is numerically singular!'
       end if
       B(:,:) = 1.0_SP / A(:,:)
    case (2)
       B(:,:) = invMatSP2(A)
    case (3)
       B(:,:) = invMatSP3(A)
    case (4)
       B(:,:) = invMatSP4(A)
    case (5:)
       B(:,:) = invMatSPN(A)
    end select

  end function invMatSP

  function invMatDP(A) result(B)
    real(DP), intent(in) :: A(:,:)
    real(DP)             :: B(size(A,1),size(A,2))
    ! Local variables
    integer              :: n

    n = size(A,1)
    select case (n)
    case (1)
       if (A(1,1).eq.0.0_DP) then
          write(*,*) 'Matrix is numerically singular!'
       end if
       B(:,:) = 1.0_DP / A(:,:)
    case (2)
       B(:,:) = invMatDP2(A)
    case (3)
       B(:,:) = invMatDP3(A)
    case (4)
       B(:,:) = invMatDP4(A)
    case (5:)
       B(:,:) = invMatDPN(A)
    end select

  end function invMatDP

  function invMatComplexSP(A) result(B)
    complex(SP), intent(in) :: A(:,:)
    complex(SP)             :: B(size(A,1),size(A,2))
    ! Local variables
    integer              :: n

    n = size(A,1)
    select case (n)
    case (1)
       if (A(1,1).eq.cmplx(0.0,0.0)) then
          write(*,*) 'Matrix is numerically singular!'
       end if
       B(:,:) = 1.0 / A(:,:)
    case (2)
       B(:,:) = invMatCmplxSP2(A)
    case (3)
       B(:,:) = invMatCmplxSP3(A)
    case (4)
       B(:,:) = invMatCmplxSP4(A)
    case (5:)
       B(:,:) = invMatCmplxSPN(A)
    end select

  end function invMatComplexSP

  function invMatComplexDP(A) result(B)
    complex(DP), intent(in) :: A(:,:)
    complex(DP)             :: B(size(A,1),size(A,2))
    ! Local variables
    integer              :: n

    n = size(A,1)
    select case (n)
    case (1)
       if (A(1,1).eq.cmplx(0.0_DP,0.0_DP)) then
          write(*,*) 'Matrix is numerically singular!'
       end if
       B(:,:) = 1.0_DP / A(:,:)
    case (2)
       B(:,:) = invMatCmplxDP2(A)
    case (3)
       B(:,:) = invMatCmplxDP3(A)
    case (4)
       B(:,:) = invMatCmplxDP4(A)
    case (5:)
       B(:,:) = invMatCmplxDPN(A)
    end select

  end function invMatComplexDP

  function invMatDP2(A) result(B)
    ! Performs a direct calculation of the inverse of a 2×2 matrix.
    real(dp), intent(in) :: A(2,2)   ! Matrix
    real(dp)             :: B(2,2)   ! Inverse matrix
    real(dp)             :: det

    ! Calculate the inverse determinant of the matrix
    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    if (det.eq.0.0_DP) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  A(2,2)
    B(2,1) = -A(2,1)
    B(1,2) = -A(1,2)
    B(2,2) =  A(1,1)
    B(:,:) =  B(:,:)/det
  end function invMatDP2

  function invMatDP3(A) result(B)
    ! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(dp), intent(in) :: A(3,3)   ! Matrix
    real(dp)             :: B(3,3)   ! Inverse matrix
    real(dp)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    if (det.eq.0.0_DP) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    B(:,:) =  B(:,:)/det
  end function invMatDP3

  function invMatDP4(A) result(B)
    ! Performs a direct calculation of the inverse of a 4×4 matrix.
    real(dp), intent(in) :: A(4,4)   ! Matrix
    real(dp)             :: B(4,4)   ! Inverse matrix
    real(dp)             :: det

    ! Calculate the inverse determinant of the matrix
    det = &
         1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    if (det.eq.0.0_DP) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) = (A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = (A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = (A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = (A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = (A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = (A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = (A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = (A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = (A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = (A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = (A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = (A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = (A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = (A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = (A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = (A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    B(:,:) = B(:,:)*det
  end function invMatDP4

  function invMatDPN(A) result(Ainv)

    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(size(A,1),size(A,2)) :: Ainv

    real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dgetrf(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call dgetri(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function invMatDPN

  function invMatSP2(A) result(B)
    ! Performs a direct calculation of the inverse of a 2×2 matrix.
    real(SP), intent(in) :: A(2,2)   ! Matrix
    real(SP)             :: B(2,2)   ! Inverse matrix
    real(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    if (det.eq.0.0) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  A(2,2)
    B(2,1) = -A(2,1)
    B(1,2) = -A(1,2)
    B(2,2) =  A(1,1)
    B(:,:) =  B(:,:)/det
  end function invMatSP2

  function invMatSP3(A) result(B)
    ! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(SP), intent(in) :: A(3,3)   ! Matrix
    real(SP)             :: B(3,3)   ! Inverse matrix
    real(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    if (det.eq.0.0) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    B(:,:) =  B(:,:)/det
  end function invMatSP3

  function invMatSP4(A) result(B)
    ! Performs a direct calculation of the inverse of a 4×4 matrix.
    real(SP), intent(in) :: A(4,4)   ! Matrix
    real(SP)             :: B(4,4)   ! Inverse matrix
    real(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = &
         1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    if (det.eq.0.0) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) = (A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = (A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = (A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = (A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = (A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = (A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = (A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = (A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = (A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = (A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = (A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = (A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = (A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = (A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = (A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = (A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    B(:,:) = B(:,:)*det
  end function invMatSP4

  function invMatSPN(A) result(Ainv)

    real(SP), dimension(:,:), intent(in) :: A
    real(SP), dimension(size(A,1),size(A,2)) :: Ainv

    real(SP), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dgetrf(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call dgetri(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function invMatSPN

  function invMatCmplxSP2(A) result(B)
    ! Performs a direct calculation of the inverse of a 2×2 matrix.
    complex(SP), intent(in) :: A(2,2)   ! Matrix
    complex(SP)             :: B(2,2)   ! Inverse matrix
    complex(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2) - A(1,2)*A(2,1))

    if (det.eq.cmplx(0.0,0.0)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  A(2,2)
    B(2,1) = -A(2,1)
    B(1,2) = -A(1,2)
    B(2,2) =  A(1,1)
    B(:,:) =  B(:,:)/det
  end function invMatCmplxSP2

  function invMatCmplxSP3(A) result(B)
    ! Performs a direct calculation of the inverse of a 3×3 matrix.
    complex(SP), intent(in) :: A(3,3)   ! Matrix
    complex(SP)             :: B(3,3)   ! Inverse matrix
    complex(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    if (det.eq.cmplx(0.0,0.0)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    B(:,:) =  B(:,:)/det
  end function invMatCmplxSP3

  function invMatCmplxSP4(A) result(B)
    ! Performs a direct calculation of the inverse of a 4×4 matrix.
    complex(SP), intent(in) :: A(4,4)   ! Matrix
    complex(SP)             :: B(4,4)   ! Inverse matrix
    complex(SP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = &
         (A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    if (det.eq.cmplx(0.0,0.0)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) = (A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = (A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = (A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = (A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = (A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = (A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = (A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = (A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = (A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = (A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = (A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = (A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = (A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = (A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = (A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = (A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    B(:,:) = B(:,:)/det
  end function invMatCmplxSP4

  function invMatCmplxSPN(A) result(Ainv)
    complex(SP), dimension(:,:), intent(in) :: A
    complex(SP), dimension(size(A,1),size(A,2)) :: Ainv

    complex(SP), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call zgetrf(n, n, Ainv, n, ipiv, info)

    if (info .ne. 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call zgetri(n, Ainv, n, ipiv, work, n, info)

    if (info .ne. 0) then
       stop 'Matrix inversion failed!'
    end if
  end function invMatCmplxSPN

  function invMatCmplxDP2(A) result(B)
    ! Performs a direct calculation of the inverse of a 2×2 matrix.
    complex(DP), intent(in) :: A(2,2)   ! Matrix
    complex(DP)             :: B(2,2)   ! Inverse matrix
    complex(DP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2) - A(1,2)*A(2,1))

    if (det.eq.cmplx(0.0_DP,0.0_DP)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  A(2,2)
    B(2,1) = -A(2,1)
    B(1,2) = -A(1,2)
    B(2,2) =  A(1,1)
    B(:,:) =  B(:,:)/det
  end function invMatCmplxDP2

  function invMatCmplxDP3(A) result(B)
    ! Performs a direct calculation of the inverse of a 3×3 matrix.
    complex(DP), intent(in) :: A(3,3)   ! Matrix
    complex(DP)             :: B(3,3)   ! Inverse matrix
    complex(DP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    if (det.eq.cmplx(0.0_DP,0.0_DP)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    B(:,:) =  B(:,:)/det
  end function invMatCmplxDP3

  function invMatCmplxDP4(A) result(B)
    ! Performs a direct calculation of the inverse of a 4×4 matrix.
    complex(DP), intent(in) :: A(4,4)   ! Matrix
    complex(DP)             :: B(4,4)   ! Inverse matrix
    complex(DP)             :: det

    ! Calculate the inverse determinant of the matrix
    det = &
         (A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    if (det.eq.cmplx(0.0_DP,0.0_DP)) then
       write(*,*) 'Matrix is numerically singular!'
    end if

    ! Calculate the inverse of the matrix
    B(1,1) = (A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = (A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = (A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = (A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = (A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = (A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = (A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = (A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = (A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = (A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = (A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = (A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = (A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = (A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = (A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = (A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    B(:,:) = B(:,:)/det
  end function invMatCmplxDP4

  function invMatCmplxDPN(A) result(Ainv)
    complex(DP), dimension(:,:), intent(in) :: A
    complex(DP), dimension(size(A,1),size(A,2)) :: Ainv

    complex(DP), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call zgetrf(n, n, Ainv, n, ipiv, info)

    if (info .ne. 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call zgetri(n, Ainv, n, ipiv, work, n, info)

    if (info .ne. 0) then
       stop 'Matrix inversion failed!'
    end if
  end function invMatCmplxDPN


  !---------------------------------------------------------------------!
  !                                                                     !
  !                      Useful Mathamatical functions                  !
  !                                                                     !
  !---------------------------------------------------------------------!

  function identityMat(n) result(id)
    implicit none

    integer, intent(in) :: n
    integer, dimension(n,n) :: id
    integer :: i

    id = 0
    do i = 1,n
       id(i,i) = 1
    end do

  end function identityMat

  function factorial(n)
    implicit none

    integer,intent(in) :: n
    integer            :: factorial

    integer            :: i

    factorial = 1

    do i = 1,n
       factorial = factorial*i
    end do

  end function factorial

  function sgnn(x)
    use kinds
    implicit none

    real(DP), intent(in) :: x
    real(DP)             :: sgnn

    if( x .eq. 0.0_DP ) then
       sgnn = 1.0_DP
    else
       sgnn = abs(x)/x
    end if

  end function sgnn


  !---------------------------------------------------------------------!
  !                                                                     !
  !                               MeshGrid                              !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! creates a unique lattice of points given 2 vectors x,y of length    !
  ! N and M. Can be used for 3d plots. See Matlab meshgrid              !
  ! documentation for more details                                      !
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
  ! Fits a cubic spline to input data. First version fits and outputs   !
  ! values at specified interpolated points. Next version returns       !
  ! the coefficients of the cubic spline fit which can then be passed   !
  ! to the next function to evaluate the spline at a certain x.         !
  !                                                                     !
  ! Uses Lapack for necessary linear algebra                            !
  !---------------------------------------------------------------------!

  function splinefit(x,y,order) result(c)
    use kinds
    use LAPACK95
    implicit none

    real(DP),dimension(:),intent(in)    :: x,y
    integer              ,intent(in)    :: order
    real(DP),dimension(size(x))         :: c

    integer                             :: N
    real(DP),dimension(:,:),allocatable :: A,Ainv
    integer                             :: i,j

    N = size(x)
    allocate(A(N,N),Ainv(N,N))

    A(:,1) = 1.0_DP
    A(:,2) = x(:)
    do i = 1,N
       do j = 3,N
          A(i,j) = abs(x(i)-x(j-1))**order
       end do
    end do

    Ainv = inv(A)
    c = matmul(Ainv,y)

    deallocate(A,Ainv)

  end function splinefit

  function splineVal(c,xj,x,order)
    use kinds
    implicit none

    real(DP),dimension(:),intent(in) :: c,xj
    real(DP),intent(in)              :: x
    integer              ,intent(in) :: order
    real(DP)                         :: splineval

    real(DP), dimension(size(c)-2)   :: dist,diff
    integer                          :: Lc,Lxj

    Lc  = size(c)
    Lxj = Lc-1
    diff = (abs( x-xj(2:Lxj) ))**order
    dist  = c(3:Lc)*diff
    splineval = c(1)+c(2)*x+sum(dist)

  end function splineVal

  function splineDeriv(c,xj,x,order) result(val)
    use kinds
    implicit none

    real(DP),dimension(:),intent(in) :: c,xj
    real(DP),intent(in)              :: x
    integer              ,intent(in) :: order
    real(DP)                         :: val

    integer, dimension(size(c)-2)    :: signn
    real(DP), dimension(size(c)-2)   :: vals
    integer :: ii,Lc,Lxj

    Lc  = size(c)
    Lxj = Lc-1

    do ii = 1,Lc-2
       if( x-xj(ii+1) .ge. 0) then
          signn(ii) = 1
       else
          signn(ii) = -1
       end if
    end do

    vals = (abs(x-xj(2:Lxj)))**(order-1)
    val = c(2)+order*sum(c(3:Lc)*signn*vals)

  end function splineDeriv

  !---------------------------------------------------------------------!
  !                                                                     !
  !                               PolyFit                               !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! N dimensional polynomial fitting algorithim which outputs the       !
  ! coefficitents of the fit. Can then be passed into polycal to        !
  ! evaluate the fit at a certain point.                                !
  !                                                                     !
  ! Uses Lapack for necessary linear algebra                            !
  !---------------------------------------------------------------------!

  function polyVal(c,x) result(val)
    use kinds
    implicit none

    real(DP),dimension(:),intent(in) :: c
    real(DP)             ,intent(in) :: x

    real(DP)                           :: val
    integer                            :: i,powers(size(c))

    powers = linspace(0,size(c)-1,size(c))
    val = sum(c*x**powers)

  end function polyVal


  function polyfit(vx, vy, d)
    use kinds
    use Lapack95
    implicit none
    integer, intent(in)                   :: d
    real(dp), dimension(d+1)              :: polyfit
    real(dp), dimension(:), intent(in)    :: vx, vy

    real(dp), dimension(:,:), allocatable :: X
    real(dp), dimension(:,:), allocatable :: XT
    real(dp), dimension(:,:), allocatable :: XTX

    integer :: i, j

    integer     :: n, lda, lwork
    integer :: info
    integer, dimension(:), allocatable :: ipiv
    real(dp), dimension(:), allocatable :: work

    n = d+1
    lda = n
    lwork = n

    allocate(ipiv(n))
    allocate(work(lwork))
    allocate(XT(n, size(vx)))
    allocate(X(size(vx), n))
    allocate(XTX(n, n))

    ! prepare the matrix
    do i = 0, d
       do j = 1, size(vx)
          X(j, i+1) = vx(j)**i
       end do
    end do

    XT  = transpose(X)
    XTX = matmul(XT, X)

    ! calls to LAPACK subs DGETRF and DGETRI
    call DGETRF(n, n, XTX, lda, ipiv, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if
    call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if

    polyfit = matmul( matmul(XTX, XT), vy)

    deallocate(ipiv)
    deallocate(work)
    deallocate(X)
    deallocate(XT)
    deallocate(XTX)

  end function polyfit

  !---------------------------------------------------------------------!
  !                                                                     !
  !                         Polynomial Integral                         !
  !                                                                     !
  !---------------------------------------------------------------------!

  function polyInt(c,a,b) result(val)
    use kinds
    implicit none

    real(DP), dimension(:),intent(in) :: c
    real(DP)              ,intent(in) :: a,b
    real(DP)                          :: val

    real(DP) :: aval,bval
    integer  :: ii

    bval = 0
    aval = 0

    do ii = 1,size(c)
       bval = bval+c(ii)*(b**ii)/ii
       aval = aval+c(ii)*(a**ii)/ii
    end do

    val = bval - aval

  end function polyInt

  !---------------------------------------------------------------------!
  !                                                                     !
  !                            Eulers method                            !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Eulers method for solving 1 or N coupled DE's                       !
  !---------------------------------------------------------------------!

  function Euler(f,h,t0,y0)
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
    real(DP) :: Euler

    Euler = h*f(t0,y0)+y0

  end function Euler

  function eulerND(f,h,t0,y0)
    use kinds
    implicit none

    real(DP),intent(in)              :: t0,h
    real(DP),dimension(:),intent(in) :: y0
    interface
       function f(t0,y0,nEq)
         use kinds
         implicit none
         integer,intent(in)    :: nEq
         real(DP),intent(in)   :: t0,y0(nEq)
         real(DP),dimension(nEq) :: f
       end function f
    end interface
    real(DP),dimension(size(y0)) :: eulerND
    integer :: nEq

    nEq = size(y0)
    eulerND = h*f(t0,y0,nEq)+y0

  end function eulerND

  !---------------------------------------------------------------------!
  !                                                                     !
  !                       Runge Kutta 4th Order                         !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Runge Kutta calculation for input function(s). For a single step    !
  ! or a range of specified values. Algorithims for a single DE or      !
  ! N coupled DE's.                                                     !
  !---------------------------------------------------------------------!

  function rk4(f,h,t0,y0)
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
    real(DP) :: rk4

    real(DP) :: k1,k2,k3,k4

    k1 = f(t0,y0)
    k2 = f(t0+0.5_dp*h,y0+0.5_dp*h*k1)
    k3 = f(t0+0.5_dp*h,y0+0.5_dp*h*k2)
    k4 = f(t0+h,y0+k3*h)
    rk4 = y0 + (1/6.0_dp)*(k1+2*k2+2*k3+k4)*h

  end function rk4

  function rk4N(f,h,t0,y0)
    use kinds
    implicit none

    real(DP),intent(in)              :: t0,h
    real(DP),dimension(:),intent(in) :: y0
    interface
       function f(t0,y0,nEq)
         use kinds
         implicit none
         integer,intent(in)    :: nEq
         real(DP),intent(in)   :: t0,y0(nEq)
         real(DP),dimension(nEq) :: f
       end function f
    end interface
    real(DP),dimension(size(y0)) :: rk4N

    real(DP),dimension(size(y0)) :: k1,k2,k3,k4
    integer :: nEq

    nEq = size(y0)
    k1 = f(t0,y0,nEq)
    k2 = f(t0 + 0.5_DP*h,y0 + 0.5_DP*h*k1,nEq)
    k3 = f(t0 + 0.5_DP*h,y0 + 0.5_DP*h*k2,nEq)
    k4 = f(t0 + h,y0 + h*k3,nEq)

    rk4N = y0 + (1/6.0_DP)*(k1 + 2*k2 + 2*k3 + k4)*h

  end function rk4N

  !---------------------------------------------------------------------!
  !                                                                     !
  !                               Guess Zero                            !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Returns the points where a function changes sign. Can be an input   !
  ! function of an input set of values.                                 !
  !---------------------------------------------------------------------!

  function guessZero(fvals)
    use Kinds
    implicit none

    real(DP),dimension(:),intent(in) :: fvals
    integer                          :: GuessZero

    integer  :: j,sign,newsign

    sign = int(fvals(1)/abs(fvals(1)))
    do j = 2,size(fvals)
       newsign = int(fvals(j)/abs(fvals(j)))
       if ( newsign .ne. sign) exit
    end do
    if ( j .eq. size(fvals) ) then
       write(*,'(a)') "ERROR did not see any sign change/ zero"
       write(*,'(a)') "maybe try again with more points"
    else
       GuessZero = j
    end if

  end function GuessZero

  function guessZeroF(f,a,b)
    use kinds
    implicit none

    real(DP),intent(in) :: a,b
    real(DP)            :: GuessZeroF
    interface
       function f(x)
         use kinds
         implicit none
         real(DP),intent(in) :: x
         real(DP)            :: f
       end function f
    end interface

    real(DP) :: h,fval
    integer :: i,sign,newsign,N

    h = a*1e-1
    N = int((b-a)/h)
    fval = f(a)
    sign = int(fval/abs(fval))

    do i = 1,N
       fval = f(a+i*h)
       newsign = int(fval/abs(fval))
       if ( newsign .ne. sign ) exit
    end do
    if ( i .eq. N ) then
       write(*,*) "Error Max itterations reached"
    else
       GuessZeroF = a+i*h
    end if

  end function GuessZeroF

  !---------------------------------------------------------------------!
  !                                                                     !
  !                         Newtons Method                              !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Newtons method for an input function of 1 variable                  !
  !                                                                     !
  !                 CODE WRITTEN BY CURTIS ABELL                        !
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
  !                               Linspace                              !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Creates a linear space of points. See matlab documentation for      !
  ! further details.                                                    !
  !---------------------------------------------------------------------!

  function linspace(start,finish,N)
    use kinds
    implicit none

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

  function linspaceReal(start,finish,N)
    use kinds
    implicit none

    integer                :: N
    real, intent(in)       :: start, finish
    real, dimension(N)     :: linspaceReal

    real(DP)               :: int
    integer                :: i

    linspaceReal(1) = start
    linspaceReal(N) = finish
    int  = (real(finish)-start)/(N-1)

    do i=2,N-1
       linspaceReal(i)=linspaceReal(i-1)+int
    end do

  end function linspaceReal

  function linspaceInt(start,finish,N)
    use kinds
    implicit none

    integer                :: N
    integer, intent(in)    :: start, finish
    integer, dimension(N)  :: linspaceInt

    real(DP)               :: int
    integer                :: i

    linspaceInt(1) = start
    linspaceInt(N) = finish
    int  = (real(finish)-start)/(N-1)

    do i=2,N-1
       linspaceInt(i)=linspaceInt(i-1)+int
    end do

  end function linspaceInt

  function logspace(x0,Lpow,Rpow,npts) result(vec)
    use kinds
    implicit none

    real(DP), intent(in) :: x0,Lpow,Rpow
    integer , intent(in) :: npts
    real(DP) , dimension(npts) :: vec

    real(DP) :: powers(npts)
    integer  :: ii

    powers = linspace(Lpow,Rpow,npts)

    do ii = 1,npts
       vec(ii) = x0*10**powers(ii)
    end do

  end function logspace

  !---------------------------------------------------------------------!
  !                                                                     !
  !                         Iteration loop timer                        !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine LoopTimer(i,L,t0,tF)
    use kinds
    implicit none

    real(DP), intent(in)  :: t0
    real(DP), intent(out) :: tF
    integer , intent(in)  :: i,L

    character(len=50)  :: bar,barNew
    real(DP):: delta
    integer :: ii,timeTotal

    call cpu_time(tF)

    bar = '='
    barNew = bar
    delta = tF-t0

    timeTotal = int(delta*L-delta*i)

    do ii = 1,((real(i)/L)*100)/2-1
       barNew = trim(barNew)//trim(bar)
    end do

    barNew = trim(barNew)//'>'

    write(*,'(a13,i3,a1,i3,a2,a50,a2,f5.1,a11,f5.1,a3)',advance='no') &
         & '| Iteration: ',i,'/',L,' |',barNew,'| ',(real(i)/L)*100,'% | It Len ',delta,&
         & 's |'

    if (timeTotal/60 .le. 9 .and. mod(timeTotal,60) .le. 9)   write(*,'(i2,a2,i1,a2,i1)') (timeTotal/60)/60,':0',(timeTotal/60),':0',mod(timeTotal,60)
    if (timeTotal/60 .le. 9 .and. mod(timeTotal,60) .ge. 10)  write(*,'(i2,a2,i1,a1,i2)') (timeTotal/60)/60,':0',(timeTotal/60),':', mod(timeTotal,60)
    if (timeTotal/60 .ge. 10 .and. mod(timeTotal,60) .le. 9)  write(*,'(i2,a1,i2,a2,i1)') (timeTotal/60)/60,':', (timeTotal/60),':0',mod(timeTotal,60)
    if (timeTotal/60 .ge. 10 .and. mod(timeTotal,60) .ge. 10) write(*,'(i2,a1,i2,a1,i2)') (timeTotal/60)/60,':', (timeTotal/60),':', mod(timeTotal,60)

  end subroutine LoopTimer

  !---------------------------------------------------------------------!
  !                                                                     !
  !                          Finite difference                          !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Numerically calculates the derivative with finite difference        !
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
  !                       Trapezodial Integration                       !
  !                                                                     !
  !---------------------------------------------------------------------!

  function trapZ(f,a,b,relErr) result(val)
    use kinds
    implicit none

    real(DP), intent(in) :: a,b,relErr
    interface
       function f(x)
         use kinds
         implicit none
         real(DP),intent(in) :: x
         real(DP)            :: f
       end function f
    end interface
    real(DP) :: val

    real(DP) :: width,f1,f2
    integer  :: N,i

    val = 0
    width = relErr
    N = int((b-a)/width)

    do i = 2,N+1
       f1 = f(a+width*(i-1))
       f2 = f(a+width*(i-2))
       val = val+(f1+f2)*width/2
    end do

  end function trapZ


  function trapZZ(x,y) result(val)
    use kinds
    implicit none

    real(DP), dimension(:), intent(in) :: x,y
    real(DP) :: val

    real(DP) :: f1,f2,width
    integer  :: N,i

    val = 0
    N = size(x)
    width = abs(x(1)-x(2))

    do i = 2,N
       f1 = y(i)
       f2 = y(i-1)
       val = val+(f1+f2)*width/2
    end do

  end function trapZZ

  !---------------------------------------------------------------------!
  !                                                                     !
  !                     General Minimization Wrapper                    !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine Minimize(TestFun, n, x, e, maxStepErrorScaleFactor &
       &, minimum, doPrinting, maxIters)
    use kinds
    implicit none
    interface
       subroutine TestFun(n, x, f)
         USE kinds
         implicit none
         integer               , intent(in)  :: n
         real(DP), dimension(n), intent(in)  :: x
         real(DP)              , intent(out) :: f
       end subroutine TestFun
    end interface
    integer,                intent(in)    :: n
    real(DP), dimension(n), intent(inout) :: x, e
    real(DP),               intent(in)    :: maxStepErrorScaleFactor
    real(DP),               intent(out)   :: minimum
    logical, optional                     :: doPrinting
    integer, optional                     :: maxIters

    ! Local variables to interface with minf
    ! Written by Derek Leinweber, optional arguements added by Curtis Abell
    !
    ! These variables control how minf works
    !
    ! icon=2 helps to avoid local minima.  Highly recommended.
    ! See minf for details on the other parameters
    !
    integer, parameter           :: icon=2
    integer                      :: iprint, maxit, iterCount
    integer                      :: nwork
    real(DP), dimension(n*(n+3)) :: w

    ! If the optional doPrinting arguement is present and false, minf
    !    will not print any output
    if (present(doPrinting)) then
       if (doPrinting) then
          iprint = 2
       else
          iprint = 0
       end if
    else
       iprint = 2
    end if

    ! Set the max iterations to 100,000 if they are not specified
    if (present(maxIters)) then
       maxit = maxIters
    else
       maxit = 100000
    end if


    nwork = n*(n+3)

    call minf( TestFun, x, e, n, minimum, maxStepErrorScaleFactor, &
         &     iprint, icon, maxit, iterCount, w, nwork )

  end Subroutine Minimize

  !---------------------------------------------------------------------!
  !                                                                     !
  !                      QuadPack Integral wrapper                      !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Wrappers for the Quadpack integral routine shown above.             !
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
    integer  :: neval, ifail
    real(DP) :: errEstimate

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
    integer  :: neval, ifail
    real(DP) :: errEstimate

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
    integer  :: neval, ifail
    real(DP) :: errEstimate

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
    integer  :: neval, ifail
    real(DP) :: errEstimate

    BreakPtsP2(1:nBreakPts) = BreakPts(1:nBreakPts)       ! Array section limits are required here.

    call qagp(f, a, b, nBreakPts+2, BreakPtsP2, absErr, relErr, integResult, errEstimate, neval, ifail)
    if ( ifail /= 0 ) then
       write(*,*) 'Warning from qagp: the error code is ', ifail
    end if
    integralBreakPts = integResult

  end function integralBreakPts

  !---------------------------------------------------------------------!
  !                                                                     !
  !                      Principle Value Integrator                     !
  !                                                                     !
  !---------------------------------------------------------------------!
  ! Evaluates the Cauchy-Principle value integral for f(x)/(x-c)        !
  !---------------------------------------------------------------------!

  function integralPV(f, c, a, b, absErr, relErr)
    use kinds
    implicit none
    real(DP)                  :: integralPV
    interface
       function f(x)
         use kinds
         implicit none
         real(DP)             :: f
         real(DP), intent(in) :: x
       end function f
    end interface
    real(DP), intent(in)      :: c, a, b, absErr, relErr

    ! ---------------------------Local Variables--------------------------
    real(DP)                  :: res, errRequest
    integer                   :: inf ! -1 for -infty, +1 for infty
    ! Used to split the integral if a or b are at plus/minux infty
    !    into semi-infinite integrals and a PV integral. This is done
    !    by moving the upper and/or lower bounds of the integral
    !    away from the pole and splitting the integral into 2 (3).
    real(DP)                  :: intBoundaryHigh, intBoundaryLow
    integer  :: neval, ifail
    real(DP) :: errEstimate


    integralPV = 0.0_DP
    intBoundaryLow = a
    intBoundaryHigh = b

    ! Integral from -Infty to lower bound
    if (a.eq.-Infty) then
       inf = -1
       if (abs(c).lt.2.0_DP) then
          intBoundaryLow = -4.0_DP
       else
          intBoundaryLow = c / 2.0_DP
       end if
       call qagi(f, intBoundaryLow, inf, absErr, relErr, res &
            & , errEstimate, neval, ifail)
       integralPV = integralPV + res
       if ( ifail.ne.0 ) then
          write(*,*) ' Warning from qagi from -Infty: the error code is ', ifail
       end if
    end if

    ! Integral from upper bound to Infty
    if (b.eq.Infty) then
       inf = 1
       if (abs(c).lt.2.0_DP) then
          intBoundaryHigh = 4.0_DP
       else
          intBoundaryHigh = c * 2.0_DP
       end if
       call qagi(fNonSingular, intBoundaryHigh, inf, absErr  &
            & , relErr, res, errEstimate, neval, ifail)
       integralPV = integralPV + res
       if ( ifail.ne.0 ) then
          write(*,*) ' Warning from qagi to Infty: the error code is ', ifail
       end if
    end if

    ! Principle-value integral
    call qawc(f, intBoundaryLow, intBoundaryHigh, c, absErr, relErr &
         & , res, errEstimate, neval, ifail)
    if ( ifail.ne.0 ) then
       write(*,*) ' Warning from qawc: the error code is ', ifail
    end if
    integralPV = integralPV + res

  contains
    ! Function to pass to the semi-infinite integrals
    !    - now includes the (x-c) on the denominator
    function fNonSingular(x)
      use kinds
      implicit none
      real(DP) :: fNonSingular
      real(DP), intent(in) :: x
      fNonSingular = f(x)/(x - c)
    end function fNonSingular

  end function integralPV

  !---------------------------------------------------------------------!
  !                                                                     !
  !                       Multi-Dimensional Pyplots                     !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine pyplotN(x,title,xaxis,yaxis,legend)
    use Kinds
    implicit none
    character(len=*),intent(in),optional              :: xaxis,yaxis,title
    character(len=*),dimension(:),intent(in),optional :: legend
    real(DP),dimension(:,:),intent(in)                :: x

    character(len=10),dimension(size(x(1,:))/2)       :: ld
    character(len=14)                                 :: fmt,xlabel,ylabel,name
    integer                                           :: ii,j,N,L

    N = size(x,dim=2)
    L = size(ld)

    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    name   = ""
    xlabel = ""
    ylabel = ""
    ld     = ""

    if(present(title)) name   = title
    if(present(xaxis)) xlabel = xaxis
    if(present(yaxis)) ylabel = yaxis
    if(present(legend))ld(1:L)= legend(1:L)

    write(100,'(a20)') name
    write(100,'(a15)') xlabel
    write(100,'(a15)') ylabel
    write(fmt,'(a1,i1,a7)') '(', N, 'es20.9)'

    do ii=1,L
       write(100,'(a12)') ld(ii)
    end do

    do ii=1,size(x,dim=1)
       write(101,fmt) x(ii,:)
    end do

    close(100)
    close(101)

  end subroutine pyplotN

  !---------------------------------------------------------------------!
  !                                                                     !
  !                          Standard Pyplots                           !
  !                                                                     !
  !---------------------------------------------------------------------!

  subroutine pyplotXY(x,y,title,xaxis,yaxis)
    use Kinds
    implicit none
    character(len=*),intent(in),optional :: xaxis,yaxis,title
    real(DP),dimension(:),intent(in)     :: x,y

    character(len=14)                    :: xlabel,ylabel,name,ld
    integer                              :: ii

    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    name   = ""
    xlabel = ""
    ylabel = ""
    ld     = "empty"

    if(present(title)) name   = title
    if(present(xaxis)) xlabel = xaxis
    if(present(yaxis)) ylabel = yaxis

    write(100,'(a20)') name
    write(100,'(a15)') xlabel
    write(100,'(a15)') ylabel
    write(100,'(a7)') ld

    do ii=1,size(x)
       write(101,'(2es20.9)') x(ii), y(ii)
    end do

    close(100)
    close(101)

  end subroutine pyplotXY

  subroutine pyplotXYZW(x,y,z,w,title,xaxis,yaxis,legend)
    use Kinds
    implicit none
    character(len=*),dimension(2),intent(in),optional :: legend
    character(len=*),intent(in),optional :: xaxis,yaxis,title
    real(DP),dimension(:),intent(in)     :: x,y,z,w

    character(len=14)                    :: xlabel,ylabel,name
    character(len=14),dimension(2)       :: ld
    integer                              :: ii

    open(100,file="titles.dat",action="write", &
         & status="replace",form="formatted")
    open(101,file="output.dat",action="write", &
         & status="replace",form="formatted")

    name   = ""
    xlabel = ""
    ylabel = ""
    ld     = "empty"

    if(present(title)) name   = title
    if(present(xaxis)) xlabel = xaxis
    if(present(yaxis)) ylabel = yaxis
    if(present(yaxis)) ld     = legend

    write(100,'(a20)') name
    write(100,'(a15)') xlabel
    write(100,'(a15)') ylabel
    write(100,'(a7)') ld
    write(100,'(a7)') ld

    do ii=1,size(x)
       write(101,'(4es20.9)') x(ii), y(ii), z(ii), w(ii)
    end do

    close(100)
    close(101)

  end subroutine pyplotXYZW

end module numFort
