module kinds
  implicit none

  integer, parameter :: I4B = SELECTED_INT_kind(9)
  integer, parameter :: I2B = SELECTED_INT_kind(4)
  integer, parameter :: I1B = SELECTED_INT_kind(2)
  integer, parameter :: SP = kind(1.0)
  integer, parameter :: DP = kind(1.0D0)
  integer, parameter :: SPC = kind((1.0,1.0))
  integer, parameter :: DPC = kind((1.0D0,1.0D0))
  integer, parameter :: LGT = kind(.true.)
  real(DP), parameter :: PI=3.141592653589793238462643383279502884197_sp
  real(DP), parameter :: PIO2=1.57079632679489661923132169163975144209858_sp
  real(DP), parameter :: TWOPI=6.283185307179586476925286766559005768394_sp
  real(DP), parameter :: SQRT2=1.41421356237309504880168872420969807856967_sp
  real(DP), parameter :: EULER=0.5772156649015328606065120900824024310422_sp
  real(DP), parameter :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  real(DP), parameter :: TWOPI_D=6.283185307179586476925286766559005768394_dp
  type sprs2_sp
     integer(I4B) :: n,len
     real(DP), dimension(:), pointer :: val
     integer(I4B), dimension(:), pointer :: irow
     integer(I4B), dimension(:), pointer :: jcol
  end type sprs2_sp
  type sprs2_dp
     integer(I4B) :: n,len
     real(DP), dimension(:), pointer :: val
     integer(I4B), dimension(:), pointer :: irow
     integer(I4B), dimension(:), pointer :: jcol
  end type sprs2_dp

end module kinds
