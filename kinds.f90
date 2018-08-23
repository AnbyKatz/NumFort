module kinds
  implicit none

  integer,  parameter :: SP = kind(1.0)
  integer,  parameter :: DP = kind(1.0D0)
  real(DP), parameter :: PI = 3.141592653589793238462643383279502884197_dp
  real(DP), parameter :: Infty = huge(1.0_DP)

end module kinds
