module hompack_kinds
!! Real kinds and common numeric constants used by `hompack`.

   use iso_fortran_env, only: real64
   implicit none
   private

   public :: dp
   public :: zero, one, pi

   integer, parameter :: dp = real64

   real(dp), parameter :: zero = 0.0_dp
   real(dp), parameter :: one = 1.0_dp
   real(dp), parameter :: pi = 4*atan(one)

end module hompack_kinds
